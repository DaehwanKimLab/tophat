/*
 *  tophat_reports.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/20/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif

#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <getopt.h>
#include <inttypes.h>
#include <boost/thread.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "common.h"
#include "utils.h"
#include "bwt_map.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "fusions.h"
#include "align_status.h"
#include "fragments.h"
#include "wiggles.h"
#include "tokenize.h"
#include "reads.h"
#include "coverage.h"
#include "inserts.h"
#include "bam_merge.h"

using namespace std;
using namespace seqan;
using std::set;

static const JunctionSet empty_junctions;
static const InsertionSet empty_insertions;
static const DeletionSet empty_deletions;
static const FusionSet empty_fusions;
static const Coverage empty_coverage;

// daehwan - this is redundancy, which should be removed.
void get_seqs(istream& ref_stream,
	      RefSequenceTable& rt,
	      bool keep_seqs = true)
{    
  while(ref_stream.good() && !ref_stream.eof())
    {
      RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence();
      string name;
      readMeta(ref_stream, name, Fasta());
      string::size_type space_pos = name.find_first_of(" \t\r");
      if (space_pos != string::npos)
	{
	  name.resize(space_pos);
	}
      fprintf(stderr, "\tLoading %s...", name.c_str());
      seqan::read(ref_stream, *ref_str, Fasta());
      fprintf(stderr, "done\n");
      
      rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
    }
}

struct cmp_read_alignment
{
  bool operator() (const BowtieHit& l, const BowtieHit& r) const
  {
    return l.alignment_score() > r.alignment_score();
  }
};

struct cmp_read_equal
{
  bool operator() (const BowtieHit& l, const BowtieHit& r) const
  {
    if (l.insert_id() != r.insert_id() ||
	l.ref_id() != r.ref_id() || 
	l.ref_id2() != r.ref_id2() ||
	l.left() != r.left() ||
	l.right() != r.right() ||
	l.antisense_align() != r.antisense_align() ||
	l.mismatches() != r.mismatches() ||
	l.edit_dist() != r.edit_dist() ||
	l.cigar().size() != r.cigar().size())
      return false;

    return true;
  }
};

char read_best_alignments(const HitsForRead& hits_for_read,
		HitsForRead& best_hits,
		const JunctionSet& gtf_junctions,
		const JunctionSet& junctions = empty_junctions,
		const InsertionSet& insertions = empty_insertions,
		const DeletionSet& deletions = empty_deletions,
		const FusionSet& fusions = empty_fusions,
		const Coverage& coverage = empty_coverage,
		bool final_report = false,
		boost::mt19937* rng = NULL)
{
	char ret_code=0;
	const vector<BowtieHit>& hits = hits_for_read.hits;
	/* if (hits.size() >= max_multihits * 5) {
		ret_code |= 16;
		return ret_code; //too many hits
	}*/
	unsigned int nhits=0;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		if (hits[i].mismatches() > read_mismatches ||
				hits[i].gap_length() > read_gap_length ||
				hits[i].edit_dist() > read_edit_dist)
			continue;
		BowtieHit hit = hits[i];
		AlignStatus align_status(hit, gtf_junctions,
				junctions, insertions, deletions, fusions, coverage);
		hit.alignment_score(align_status._alignment_score);
		if (report_secondary_alignments || !final_report)
		{
			best_hits.hits.push_back(hit);
		}
		else
		{
			// Is the new status better than the current best one?
			if (best_hits.hits.size() == 0 || cmp_read_alignment()(hit, best_hits.hits[0]))
			{
				best_hits.hits.clear();
				best_hits.hits.push_back(hit);
			}
			else if (!cmp_read_alignment()(best_hits.hits[0], hit)) // is it just as good?
			{
				best_hits.hits.push_back(hit);
			}
		}
	}

	// due to indel alignments, there may be alignments with the same location
	std::sort(best_hits.hits.begin(), best_hits.hits.end());
	vector<BowtieHit>::iterator new_end = std::unique(best_hits.hits.begin(), best_hits.hits.end(), cmp_read_equal());
	best_hits.hits.erase(new_end, best_hits.hits.end());

	if ((report_secondary_alignments || !final_report) && best_hits.hits.size() > 0)
	{
		sort(best_hits.hits.begin(), best_hits.hits.end(), cmp_read_alignment());
	}

	if (final_report) {
		if (best_hits.hits.size() > max_multihits) {
			ret_code |= 16; //read has too many valid mappings
			if (suppress_hits)
				best_hits.hits.clear();
			else {
				//if (best_hits.hits.size() > max_multihits)
				// there may be several alignments with the same alignment scores,
				// all of which we can not include because of this limit.
				// let's pick up some of them randomly up to this many max multihits.
				vector<size_t> tie_indexes;
				int tie_alignment_score = best_hits.hits[max_multihits - 1].alignment_score();
				int count_better_alignments = 0;
				for (size_t i = 0; i < best_hits.hits.size(); ++i)
				{
					int temp_alignment_score = best_hits.hits[i].alignment_score();
					if (temp_alignment_score == tie_alignment_score)
						tie_indexes.push_back(i);
					else if (temp_alignment_score < tie_alignment_score)
						break;
					else
						++count_better_alignments;
				}

				while (count_better_alignments + tie_indexes.size() > max_multihits)
				{
					int random_index = (*rng)() % tie_indexes.size();
					tie_indexes.erase(tie_indexes.begin() + random_index);
				}

				for (size_t i = 0; i < tie_indexes.size(); ++i)
				{
					if (count_better_alignments + i != tie_indexes[i])
						best_hits.hits[count_better_alignments + i] = best_hits.hits[tie_indexes[i]];
				}

				best_hits.hits.erase(best_hits.hits.begin() + max_multihits, best_hits.hits.end());
			}
		}
		if ((nhits=best_hits.hits.size())>0) {
			ret_code |= 1; //read has a valid mapping
			if (nhits>1) ret_code |= 4; //read has multiple mappings
		}
	}
	return ret_code;
}

bool is_fusion_insert_alignment(const BowtieHit& lh, const BowtieHit& rh)
{
  bool left_fusion = lh.fusion_opcode() != FUSION_NOTHING;
  bool right_fusion = rh.fusion_opcode() != FUSION_NOTHING;
  if (left_fusion || right_fusion)
    return true;
  
  if (lh.ref_id() != rh.ref_id())
    return true;
  
  if (lh.ref_id() == rh.ref_id())
    {
      if (lh.antisense_align() == rh.antisense_align())
	return true;
      else
	{
	  int inner_dist = 0;
	  if (lh.antisense_align())
	    // used rh.left() instead of rh.right() for the cases,
	    // where reads overlap with each other or reads span introns
	    inner_dist = lh.left() - rh.left(); 
	  else
	    inner_dist = rh.left() - lh.left();
	  
	  if (inner_dist < 0 || inner_dist > (int)fusion_min_dist)
	    return true;
	}
    }

  return false;
}

bool set_insert_alignment_grade(const BowtieHit& lh, const BowtieHit& rh, const JunctionSet& junctions, InsertAlignmentGrade& grade)
{
  bool left_fusion = lh.fusion_opcode() != FUSION_NOTHING;
  bool right_fusion = rh.fusion_opcode() != FUSION_NOTHING;
  if (left_fusion && right_fusion)
    return false;
  
  bool fusion = is_fusion_insert_alignment(lh, rh);

  grade = InsertAlignmentGrade(lh, rh, junctions, fusion);
 
  return true;
}

struct cmp_pair_alignment
{
  cmp_pair_alignment(const JunctionSet& junctions) :
    _junctions(&junctions) {}
    
  bool operator() (const pair<BowtieHit, BowtieHit>& l, const pair<BowtieHit, BowtieHit>& r) const
  {
    InsertAlignmentGrade gl; set_insert_alignment_grade(l.first, l.second, *_junctions, gl);
    InsertAlignmentGrade gr; set_insert_alignment_grade(r.first, r.second, *_junctions, gr);

    bool better = gr < gl;
    bool worse = gl < gr;

    if (better && !worse)
      return true;
    else	
      return false;
  }

  const JunctionSet* _junctions;
};

struct cmp_pair_equal
{
  bool operator() (const pair<BowtieHit, BowtieHit>& l, const pair<BowtieHit, BowtieHit>& r) const
  {
    if (!cmp_read_equal()(l.first, r.first) || !cmp_read_equal()(l.second, r.second))
	return false;
	
    return true;
  }
};

struct cmp_pair_less
{
  bool operator() (const pair<BowtieHit, BowtieHit>& l, const pair<BowtieHit, BowtieHit>& r) const
  {
    if (l.first < r.first)
      return true;
    else if (r.first < l.first)
      return false;

    if (l.second < r.second)
      return true;
    else if (r.second < l.second)
      return false;
	
    return false;
  }
};

struct SAlignStats {
	int64_t num_aligned_left; //total left reads aligned
	int64_t num_unmapped_left; //total left reads unmapped, or mapped in too many places (!)
	int64_t num_aligned_left_multi; //total left reads mapped in more than 1 place
	int64_t num_aligned_left_xmulti; //total left reads mapped in too many places (> max_multihits)

	int64_t num_aligned_right; //total right reads aligned
	int64_t num_unmapped_right; //total right reads unmapped, or mapped in too many places (!)
	int64_t num_aligned_right_multi; //total right reads in more than 1 place
	int64_t num_aligned_right_xmulti; //total right reads mapped in too many places (> max_multihits)

	int64_t num_aligned_pairs; //total pairs aligned
	int64_t num_aligned_pairs_multi; //total pairs aligned in more than 1 place
	int64_t num_aligned_pairs_disc; //total pairs aligned discordantly only

	int64_t num_aligned_unpaired; //total unpaired reads aligned, when mixed with PE
	int64_t num_unmapped_unpaired; //total right reads unmapped, or mapped in too many places (!)
	int64_t num_aligned_unpaired_multi;
	int64_t num_aligned_unpaired_xmulti; //total right reads mapped in too many places (> max_multihits)

	SAlignStats():num_aligned_left(0), num_unmapped_left(0), num_aligned_left_multi(0), num_aligned_left_xmulti(0), num_aligned_right(0),
			num_unmapped_right(0), num_aligned_right_multi(0), num_aligned_right_xmulti(0), num_aligned_pairs(0), num_aligned_pairs_multi(0),
			num_aligned_pairs_disc(0),  num_aligned_unpaired(0), num_unmapped_unpaired(0),
			num_aligned_unpaired_multi(0), num_aligned_unpaired_xmulti(0) { }
	void add(SAlignStats& a) {
		num_aligned_left+=a.num_aligned_left;
		num_unmapped_left+=a.num_unmapped_left;
		num_aligned_left_multi+=a.num_aligned_left_multi;
		num_aligned_left_xmulti+=a.num_aligned_left_xmulti;
		num_aligned_right+=a.num_aligned_right;
		num_unmapped_right+=a.num_unmapped_right;
		num_aligned_right_multi+=a.num_aligned_right_multi;
		num_aligned_right_xmulti+=a.num_aligned_right_xmulti;
		num_aligned_pairs+=a.num_aligned_pairs;
		num_aligned_pairs_multi+=a.num_aligned_pairs_multi;
		num_aligned_pairs_disc+=a.num_aligned_pairs_disc;
		num_aligned_unpaired+=a.num_aligned_unpaired;
		num_unmapped_unpaired+=a.num_unmapped_unpaired;
		num_aligned_unpaired_multi+=a.num_aligned_unpaired_multi;
		num_aligned_unpaired_xmulti+=a.num_aligned_unpaired_xmulti;
	}

};

char pair_best_alignments(const HitsForRead& left_hits,
			  const HitsForRead& right_hits,
			  InsertAlignmentGrade& best_grade,
			  vector<pair<BowtieHit, BowtieHit> >& best_hits,
			  const JunctionSet& gtf_junctions,
			  const JunctionSet& junctions = empty_junctions,
			  const InsertionSet& insertions = empty_insertions,
			  const DeletionSet& deletions = empty_deletions,
			  const FusionSet& fusions = empty_fusions,
			  const Coverage& coverage = empty_coverage,
			  bool final_report = false,
			  boost::mt19937* rng = NULL
			  )
{
	const vector<BowtieHit>& left = left_hits.hits;
	const vector<BowtieHit>& right = right_hits.hits;
	char ret_code=0;
	/*
	if (left.size() >= max_multihits * 5) {
		ret_code |= 1;
	}
	if (right.size() >= max_multihits * 5) {
		ret_code |= 2;
	}
	if (ret_code)
		return ret_code;
	*/
	unsigned int nhits=0;

	vector<BowtieHit> rhits;

	for (size_t j = 0; j < right.size(); ++j)
	{
		if (right[j].mismatches() > read_mismatches ||
				right[j].gap_length() > read_gap_length ||
				right[j].edit_dist() > read_edit_dist)
			continue;
		nhits++;
		if ((nhits>>1) > max_multihits)
			break;
		BowtieHit rh = right[j];
		AlignStatus align_status(rh, gtf_junctions,
				junctions, insertions, deletions, fusions, coverage);
		rh.alignment_score(align_status._alignment_score);
		rhits.push_back(rh);
	}
    nhits=0;
	for (size_t i = 0; i < left.size(); ++i)
	{
		if (left[i].mismatches() > read_mismatches ||
				left[i].gap_length() > read_gap_length ||
				left[i].edit_dist() > read_edit_dist)
			continue;
		nhits++;
		if ((nhits>>1)>max_multihits) {
		    break;
		}
		BowtieHit lh = left[i];
		AlignStatus align_status(lh, gtf_junctions,
				junctions, insertions, deletions, fusions, coverage);
		lh.alignment_score(align_status._alignment_score);

		for (size_t j = 0; j < rhits.size(); ++j)
		{
			BowtieHit rh = rhits[j];
			InsertAlignmentGrade g;
			bool allowed;
			allowed = set_insert_alignment_grade(lh, rh, final_report ? junctions : gtf_junctions, g);

			// daehwan - for debugging purposes
#if 0
if (lh.insert_id() == 10790262)
{
	fprintf(stderr, "lh %d:%d %s score: %d (from %d) NM: %d\n",
			lh.ref_id(), lh.left(), print_cigar(lh.cigar()).c_str(),
			lh.alignment_score(), left[i].alignment_score(), lh.edit_dist());
	fprintf(stderr, "rh %d:%d %s score: %d (from %d) NM: %d\n",
			rh.ref_id(), rh.left(), print_cigar(rh.cigar()).c_str(),
			rh.alignment_score(), rhits[j].alignment_score(), rh.edit_dist());
	fprintf(stderr, "combined score: %d is_fusion(%d)\n", g.align_score(), g.is_fusion());
}
#endif

			if (!allowed) continue;
			bool new_best_grade=false;
			if (best_grade < g)
			{
				best_grade = g;
				new_best_grade=true;
			}

			if (g.fusion && !fusion_search && !report_discordant_pair_alignments) continue;

			if (report_secondary_alignments || !final_report)
			{
				best_hits.push_back(make_pair(lh, rh));
			}
			else
			{
				// Is the new status better than the current best one?
				// if (best_grade < g)
				 if (new_best_grade)
				{
					best_hits.clear();
					best_hits.push_back(make_pair(lh, rh));
				}
				else if (!(g < best_grade))
				{
					best_hits.push_back(make_pair(lh, rh));
				}
			}
		}//for j in right mate hits
	} //for i in left mate hits

	std::sort(best_hits.begin(), best_hits.end(), cmp_pair_less());

				// daehwan - for debugging purposes
			#if 0
				if (best_hits.size() > 0 && best_hits[0].first.insert_id() == 10790262)
				{
					for (size_t i = 0; i < best_hits.size(); ++i)
					{
						const BowtieHit& lh = best_hits[i].first;
						const BowtieHit& rh = best_hits[i].second;

						fprintf(stderr, "%d %d:%d %s %d:%d %s\n",
								i,
								lh.ref_id(), lh.left(), print_cigar(lh.cigar()).c_str(),
								rh.ref_id(), rh.left(), print_cigar(rh.cigar()).c_str());
					}

					fprintf(stderr, "\n\n\n");
				}
			#endif

				vector<pair<BowtieHit, BowtieHit> >::iterator new_end = std::unique(best_hits.begin(), best_hits.end(), cmp_pair_equal());
				best_hits.erase(new_end, best_hits.end());

				// daehwan - for debugging purposes
			#if 0
				if (best_hits.size() > 0 && best_hits[0].first.insert_id() == 10790262)
				{
					for (size_t i = 0; i < best_hits.size(); ++i)
					{
						const BowtieHit& lh = best_hits[i].first;
						const BowtieHit& rh = best_hits[i].second;

						fprintf(stderr, "%d %d:%d %s %d:%d %s\n",
								i,
								lh.ref_id(), lh.left(), print_cigar(lh.cigar()).c_str(),
								rh.ref_id(), rh.left(), print_cigar(rh.cigar()).c_str());
					}

					fprintf(stderr, "\n\n\n");
				}
			#endif


				if ((report_secondary_alignments || !final_report) && best_hits.size() > 0)
				{
					cmp_pair_alignment cmp(final_report ? junctions : gtf_junctions);
					sort(best_hits.begin(), best_hits.end(), cmp);
					set_insert_alignment_grade(best_hits[0].first, best_hits[0].second, final_report ? junctions : gtf_junctions, best_grade);
				}

				if (final_report)
				{
					nhits = best_hits.size();
					if (nhits>0) {
					  ret_code|=3; //both reads have acceptable mappings
					}
					if (nhits>1) ret_code |= 12; //both reads have multiple valid mappings
					if (nhits>max_multihits) {
						ret_code |= 48; //both reads have too many mappings
						if (suppress_hits) best_hits.clear();
					}

					if (best_hits.size() > max_multihits)
					{
						vector<size_t> tie_indexes;
						InsertAlignmentGrade temp_grade;
						set_insert_alignment_grade(best_hits[max_multihits - 1].first, best_hits[max_multihits - 1].second, junctions, temp_grade);
						int tie_alignment_score = temp_grade.align_score();
						int count_better_alignments = 0;

						for (size_t i = 0; i < best_hits.size(); ++i)
						{
							set_insert_alignment_grade(best_hits[i].first, best_hits[i].second, junctions, temp_grade);
							int temp_alignment_score = temp_grade.align_score();
							if (temp_alignment_score == tie_alignment_score)
								tie_indexes.push_back(i);
							else if (temp_alignment_score < tie_alignment_score)
								break;
							else
								++count_better_alignments;
						}

						while (count_better_alignments + tie_indexes.size() > max_multihits)
						{
							int random_index = (*rng)() % tie_indexes.size();
							tie_indexes.erase(tie_indexes.begin() + random_index);
						}

						for (size_t i = 0; i < tie_indexes.size(); ++i)
						{
							if (count_better_alignments + i != tie_indexes[i])
								best_hits[count_better_alignments + i] = best_hits[tie_indexes[i]];
						}

						best_hits.erase(best_hits.begin() + max_multihits, best_hits.end());
					}
				} //final report

	best_grade.num_alignments = best_hits.size();
	return ret_code;
}

enum FragmentType {FRAG_UNPAIRED, FRAG_LEFT, FRAG_RIGHT};

void add_auxData(vector<string>& auxdata, vector<string>& sam_toks,
                 const RefSequenceTable& rt,const BowtieHit& bh, FragmentType insert_side,
                 int num_hits, const BowtieHit* next_hit, int hitIndex) {
  bool XS_found = false;
  if (sam_toks.size()>11) {
    
    for (size_t i=11;i<sam_toks.size();++i) {
      if (sam_toks[i].find("NH:i:")==string::npos &&
	  sam_toks[i].find("XS:i:")==string::npos)
	auxdata.push_back(sam_toks[i]);

      if (sam_toks[i].find("XS:A:")!=string::npos)
	XS_found = true;
    }
  }
  string aux("NH:i:");
  str_appendInt(aux, num_hits);
  auxdata.push_back(aux);
  if (next_hit) {
    const char* nh_ref_name = "=";
    nh_ref_name = rt.get_name(next_hit->ref_id());
    assert (nh_ref_name != NULL);
    bool same_contig=(next_hit->ref_id()==bh.ref_id());
    aux="CC:Z:";
    aux+= (same_contig ? "=" : nh_ref_name);
    auxdata.push_back(aux);
    aux="CP:i:";
    int nh_gpos=next_hit->left() + 1;
    str_appendInt(aux, nh_gpos);
    auxdata.push_back(aux);
  } //has next_hit
    // FIXME: this code is still a bit brittle, because it contains no
    // consistency check that the mates are on opposite strands - a current protocol
    // requirement, and that the strand indicated by the alignment is consistent
    // with the orientation of the splices (though that should be handled upstream).
  if (!XS_found) {
	  const string xs_f("XS:A:+");
	  const string xs_r("XS:A:-");
	  if (library_type == FR_FIRSTSTRAND) {
		  if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED) {
			  if (bh.antisense_align())
				  auxdata.push_back(xs_f);
			  else
				  auxdata.push_back(xs_r);
		  }
		  else {
			  if (bh.antisense_align())
				  auxdata.push_back(xs_r);
			  else
				  auxdata.push_back(xs_f);
		  }
	  }
	  else if (library_type == FR_SECONDSTRAND)   {
		  if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED){
			  if (bh.antisense_align())
				  auxdata.push_back(xs_r);
			  else
				  auxdata.push_back(xs_f);
		  }
		  else
		  {
			  if (bh.antisense_align())
				  auxdata.push_back(xs_f);
			  else
				  auxdata.push_back(xs_r);
		  }
	  }
  }
  if (hitIndex >= 0)
    {
      string aux("HI:i:");
      str_appendInt(aux, hitIndex);
      auxdata.push_back(aux);
    }
}

bool rewrite_sam_record(GBamWriter& bam_writer,
			const RefSequenceTable& rt,
                        const BowtieHit& bh,
                        const char* bwt_buf,
                        const char* read_alt_name,
                        FragmentType insert_side,
                        int num_hits,
                        const BowtieHit* next_hit,
                        bool primary,
                        int hitIndex)
{
  // Rewrite this hit, filling in the alt name, mate mapping
  // and setting the pair flag
  vector<string> sam_toks;
  tokenize(bwt_buf, "\t", sam_toks);

  string ref_name = sam_toks[2], ref_name2 = "";
  char cigar1[1024] = {0}, cigar2[1024] = {0};
  string left_seq, right_seq, left_qual, right_qual;
  int left1 = -1, left2 = -1;
  bool fusion_alignment = false;
  size_t XF_index = 0;
  for (size_t i = 11; i < sam_toks.size(); ++i)
    {
      string& tok = sam_toks[i];
      if (strncmp(tok.c_str(), "XF", 2) == 0)
	{
	  fusion_alignment = true;
	  XF_index = i;
	  
	  vector<string> tuple_fields;
	  tokenize(tok.c_str(), " ", tuple_fields);
	  vector<string> contigs;
	  tokenize(tuple_fields[1].c_str(), "-", contigs);
	  if (contigs.size() >= 2)
	    {
	      ref_name = contigs[0];
	      ref_name2 = contigs[1];
	    }	
	  
	  extract_partial_hits(bh, tuple_fields[4].c_str(), tuple_fields[5].c_str(),
			       cigar1, cigar2, left_seq, right_seq,
			       left_qual, right_qual, left1, left2);
	  
	  break;
	}

      else if (strncmp(tok.c_str(), "AS", 2) == 0)
	{
	  char AS_score[128] = {0};
	  sprintf(AS_score, "AS:i:%d", min(0, bh.alignment_score()));
	  tok = AS_score;
	}
    }
  
  string qname(read_alt_name);
  size_t slash_pos=qname.rfind('/');
  if (slash_pos!=string::npos)
    qname.resize(slash_pos);
  //read_alt_name as QNAME
  int flag=atoi(sam_toks[1].c_str()); //FLAG
  if (insert_side != FRAG_UNPAIRED) {
    //flag = atoi(sam_toks[1].c_str());
    // mark this as a singleton mate
    flag |= 0x0001;
    if (insert_side == FRAG_LEFT)
      flag |= 0x0040;
    else if (insert_side == FRAG_RIGHT)
      flag |= 0x0080;
    flag |= 0x0008;
    //char flag_buf[64];
    //sprintf(flag_buf, "%d", flag);
    //sam_toks[t] = flag_buf;
  }
  if (!primary)
    flag |= 0x100;
  
  int gpos=isdigit(sam_toks[3][0]) ? atoi(sam_toks[3].c_str()) : 0;
  int mapQ = 50;
  if (num_hits > 1)  {
    double err_prob = 1 - (1.0 / num_hits);
    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
  }
  int tlen =atoi(sam_toks[8].c_str()); //TLEN
  int mate_pos=atoi(sam_toks[7].c_str());

  string rg_aux = "";
  if (!sam_readgroup_id.empty())
      rg_aux = string("RG:Z:") + sam_readgroup_id;

  GBamRecord* bamrec=NULL;
  if (fusion_alignment) {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name.c_str(), left1 + 1, mapQ,
				 cigar1, sam_toks[6].c_str(), mate_pos,
				 tlen, left_seq.c_str(), left_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;

    auxdata.clear();
    sam_toks[XF_index][5] = '2';
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name2.c_str(), left2 + 1, mapQ,
				 cigar2, sam_toks[6].c_str(), mate_pos,
				 tlen, right_seq.c_str(), right_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  } else {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, sam_toks[2].c_str(), gpos, mapQ,
				 sam_toks[5].c_str(), sam_toks[6].c_str(), mate_pos,
				 tlen, sam_toks[9].c_str(), sam_toks[10].c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  }
    
  return true;
}

bool rewrite_sam_record(GBamWriter& bam_writer,
			const RefSequenceTable& rt,
                        const BowtieHit& bh,
                        const char* bwt_buf,
                        const char* read_alt_name,
                        const InsertAlignmentGrade& grade,
                        FragmentType insert_side,
                        const BowtieHit* partner,
                        int num_hits,
                        const BowtieHit* next_hit,
                        bool primary,
                        int hitIndex)
{
  // Rewrite this hit, filling in the alt name, mate mapping
  // and setting the pair flag
  vector<string> sam_toks;
  tokenize(bwt_buf, "\t", sam_toks);
  string qname(read_alt_name);
  size_t slash_pos=qname.rfind('/');
  if (slash_pos!=string::npos)
    qname.resize(slash_pos);
  //read_alt_name as QNAME
  int flag = atoi(sam_toks[1].c_str());
  // 0x0010 (strand of query) is assumed to be set correctly
  // to begin with
  flag |= 0x0001; //it must be paired
  if (insert_side == FRAG_LEFT)
    flag |= 0x0040;
  else if (insert_side == FRAG_RIGHT)
    flag |= 0x0080;
  if (!primary)
    flag |= 0x100;

  string ref_name = sam_toks[2], ref_name2 = "";
  char cigar1[1024] = {0}, cigar2[1024] = {0};
  string left_seq, right_seq, left_qual, right_qual;
  int left1 = -1, left2 = -1;
  bool fusion_alignment = false;
  size_t XF_tok_idx = 11;
  for (; XF_tok_idx < sam_toks.size(); ++XF_tok_idx)
    {
      string& tok = sam_toks[XF_tok_idx];
      if (strncmp(tok.c_str(), "XF", 2) == 0)
	{
	  fusion_alignment = true;
	  
	  vector<string> tuple_fields;
	  tokenize(tok.c_str(), " ", tuple_fields);
	  vector<string> contigs;
	  tokenize(tuple_fields[1].c_str(), "-", contigs);
	  if (contigs.size() >= 2)
	    {
	      ref_name = contigs[0];
	      ref_name2 = contigs[1];
	    }	
	  
	  extract_partial_hits(bh, tuple_fields[4].c_str(), tuple_fields[5].c_str(),
			       cigar1, cigar2, left_seq, right_seq,
			       left_qual, right_qual, left1, left2);
	  
	  break;
	}

      else if (strncmp(tok.c_str(), "AS", 2) == 0)
	{
	  char AS_score[128] = {0};
	  sprintf(AS_score, "AS:i:%d", min(0, bh.alignment_score()));
	  tok = AS_score;
	}
    }
  
  int gpos=isdigit(sam_toks[3][0]) ? atoi(sam_toks[3].c_str()) : 0;
  int mapQ = 50;
  if (grade.num_alignments > 1) {
    double err_prob = 1 - (1.0 / grade.num_alignments);
    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
  }
  int tlen=0; //TLEN
  int mate_pos=atoi(sam_toks[7].c_str());
  string mate_contig = "*", mate_contig2 = "*";
  if (partner) {
    if (partner->ref_id() == bh.ref_id()) {
      mate_contig = "="; //same chromosome
      //TLEN:

      // a read contains its partner
      if (bh.left() <= partner->left() && bh.right() >= partner->right() && bh.right() - bh.left() > partner->right() - partner->left())
	tlen = bh.right() - bh.left();
      else if (partner->left() <= bh.left() && partner->right() >= bh.right() && partner->right() - partner->left() > bh.right() - bh.left())
	tlen = partner->left() - partner->right();
      else
	tlen = bh.left() < partner->left() ? partner->right() - bh.left() : partner->left() - bh.right();
    }
    
    else { //partner on different chromosome/contig
      mate_contig = rt.get_name(partner->ref_id());
    }
    mate_pos = partner->left() + 1;
    if (grade.happy())
      flag |= 0x0002;
    if (partner->antisense_align())
      flag |=  0x0020;

    if (fusion_alignment)
      {
	if (partner->ref_id() == bh.ref_id2())
	  mate_contig2 = "=";
	else
	  mate_contig2 = rt.get_name(partner->ref_id());
      }

    if (fusion_search)
      {
	string cigar_str = print_cigar(partner->cigar());
	char partner_pos[4096];
	if (partner->fusion_opcode() != FUSION_NOTHING)
	  {
	    sprintf(partner_pos, "XP:Z:%s-%s %d %s",
		    rt.get_name(partner->ref_id()),
		    rt.get_name(partner->ref_id2()),
		    partner->left() + 1,
		    cigar_str.c_str());
	  }
	else
	  {
	    sprintf(partner_pos, "XP:Z:%s %d %s",
		    rt.get_name(partner->ref_id()),
		    partner->left() + 1,
		    cigar_str.c_str());
	  }
	sam_toks.push_back(partner_pos);
      }
  }
  else {
    mate_pos = 0;
    flag |= 0x0008;
  }

  string rg_aux = "";
  if (!sam_readgroup_id.empty())
    rg_aux = string("RG:Z:") + sam_readgroup_id;

  GBamRecord* bamrec=NULL;
  if (fusion_alignment) {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name.c_str(), left1 + 1, mapQ,
				 cigar1, mate_contig.c_str(), mate_pos,
				 tlen, left_seq.c_str(), left_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;

    auxdata.clear();
    sam_toks[XF_tok_idx][5] = '2';
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name2.c_str(), left2 + 1, mapQ,
				 cigar2, mate_contig2.c_str(), mate_pos,
				 tlen, right_seq.c_str(), right_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  } else {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    if (rg_aux != "")
      auxdata.push_back(rg_aux);
    bamrec=bam_writer.new_record(qname.c_str(), flag, sam_toks[2].c_str(), gpos, mapQ,
				 sam_toks[5].c_str(), mate_contig.c_str(), mate_pos,
				 tlen, sam_toks[9].c_str(), sam_toks[10].c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  }

  return true;
}

struct lex_hit_sort
{
    lex_hit_sort(const RefSequenceTable& rt, const HitsForRead& hits)
	: _rt(rt), _hits(hits)
    {}
    
    bool operator()(const uint32_t& l, const uint32_t& r) const
    {
        const BowtieHit& lhs = _hits.hits[l];
        const BowtieHit& rhs = _hits.hits[r];

        uint32_t l_id = lhs.ref_id();
        uint32_t r_id = rhs.ref_id();
        if (l_id != r_id)
	  return l_id < r_id;

        return lhs.left() < rhs.left();
    }
    
    const RefSequenceTable& _rt;
    const HitsForRead& _hits;
};

void print_sam_for_single(const RefSequenceTable& rt,
			  const HitsForRead& hits,
			  FragmentType frag_type,
			  const string& alt_name,
			  GBamWriter& bam_writer,
			  boost::mt19937& rng)
{
    //assert(!read.alt_name.empty());
    if (hits.hits.empty())
      return;
    
    lex_hit_sort s(rt, hits);
    vector<uint32_t> index_vector;

    size_t i;
    for (i = 0; i < hits.hits.size(); ++i)
        index_vector.push_back(i);
    
    sort(index_vector.begin(), index_vector.end(), s);

    size_t primaryHit = 0;
    if (!report_secondary_alignments)
      primaryHit = rng() % hits.hits.size();
    
    bool multipleHits = (hits.hits.size() > 1);
    for (i = 0; i < hits.hits.size(); ++i)
      {
	size_t index = index_vector[i];
	bool primary = (index == primaryHit);
	const BowtieHit& bh = hits.hits[index];
	rewrite_sam_record(bam_writer, rt,
			   bh,
			   bh.hitfile_rec().c_str(),
			   alt_name.c_str(),
			   frag_type,
			   hits.hits.size(),
			   (i < hits.hits.size()-1) ? &(hits.hits[index_vector[i+1]]) : NULL,
			   primary,
			   (multipleHits? i: -1));
      }
}

void print_sam_for_pair(const RefSequenceTable& rt,
			const vector<pair<BowtieHit, BowtieHit> >& best_hits,
                        const InsertAlignmentGrade& grade,
                        GBamWriter& bam_writer,
                        const string& left_alt_name,
                        const string& right_alt_name,
			boost::mt19937& rng,
			uint64_t begin_id = 0,
                        uint64_t end_id = std::numeric_limits<uint64_t>::max())
{
    Read left_read;
    Read right_read;
    if (best_hits.empty())
      return;

    size_t i;
    HitsForRead right_hits;
    for (i = 0; i < best_hits.size(); ++i)
      right_hits.hits.push_back(best_hits[i].second);
    
    size_t primaryHit = 0;
    vector<uint32_t> index_vector;
    lex_hit_sort s(rt, right_hits);
    for (i = 0; i < right_hits.hits.size(); ++i)
      index_vector.push_back(i);
    
    sort(index_vector.begin(), index_vector.end(), s);

    if (!report_secondary_alignments)
      primaryHit = rng() % right_hits.hits.size();
    
    bool multipleHits = (best_hits.size() > 1);
    for (i = 0; i < best_hits.size(); ++i)
      {
	size_t index = index_vector[i];
	bool primary = (index == primaryHit);
	const BowtieHit& right_bh = best_hits[index].second;
	const BowtieHit& left_bh = best_hits[index].first;
        
	rewrite_sam_record(bam_writer, rt,
			   right_bh,
			   right_bh.hitfile_rec().c_str(),
			   right_alt_name.c_str(),
			   grade,
			   FRAG_RIGHT,
			   &left_bh,
			   best_hits.size(),
			   (i < best_hits.size() - 1) ? &(best_hits[index_vector[i+1]].second) : NULL,
			   primary, 
			   (multipleHits? i: -1));
	rewrite_sam_record(bam_writer, rt,
			   left_bh,
			   left_bh.hitfile_rec().c_str(),
			   left_alt_name.c_str(),
			   grade,
			   FRAG_LEFT,
			   &right_bh,
			   best_hits.size(),
			   (i < best_hits.size() - 1) ? &(best_hits[index_vector[i+1]].first) : NULL,
			   primary,
			   (multipleHits? i: -1));
      }
}

/**
 * Given all of the hits for a particular read, update the read counts for insertions and deletions.
 * @param hits hits The alignments for a particular read
 * @param insertions Maps from an insertion to the number of supporting reads for that insertion
 * @param deletions Maps from a deletion to the number of supporting reads for that deletion
 */
void update_insertions_and_deletions(const HitsForRead& hits,
                                     InsertionSet& insertions,
                                     DeletionSet& deletions)
{
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      const BowtieHit& bh = hits.hits[i];
      insertions_from_alignment(bh, insertions);
      deletions_from_alignment(bh, deletions);
    }
}

void update_coverage(const HitsForRead& hits,
		     Coverage& coverage)
{
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      const BowtieHit& hit = hits.hits[i];
      const vector<CigarOp>& cigar = hit.cigar();
      unsigned int positionInGenome = hit.left();
      RefID ref_id = hit.ref_id();
      
      for(size_t c = 0; c < cigar.size(); ++c)
	{
	  int opcode = cigar[c].opcode;
	  int length = cigar[c].length;
	  switch(opcode)
	    {
	    case REF_SKIP:
	    case MATCH:
	    case DEL:
	      if (opcode == MATCH)
		coverage.add_coverage(ref_id, positionInGenome, length);
      
	      positionInGenome += length;
	      break;
	    case rEF_SKIP:
	    case mATCH:
	    case dEL:
	      positionInGenome -= length;

	      if (opcode == mATCH)
		coverage.add_coverage(ref_id, positionInGenome + 1, length);
	      break;
	    case FUSION_FF:
	    case FUSION_FR:
	    case FUSION_RF:
	    case FUSION_RR:
		positionInGenome = length;
		ref_id = hit.ref_id2();
	      break;
	    default:
	      break;
	    }	
	}
    }
}


void update_fusions(const HitsForRead& hits,
		    RefSequenceTable& rt,
		    FusionSet& fusions,
		    const FusionSet& fusions_ref = empty_fusions)
{
  if (!fusion_search)
    return;
  
  if (hits.hits.size() > fusion_multireads)
    return;

  bool update_stat = fusions_ref.size() > 0;
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      const BowtieHit& bh = hits.hits[i];

      if (bh.edit_dist() > fusion_read_mismatches)
	continue;

      fusions_from_alignment(bh, fusions, rt, update_stat);

      if (update_stat)
	unsupport_fusions(bh, fusions, fusions_ref);
    }
}

void update_junctions(const HitsForRead& hits,
		      JunctionSet& junctions)
{
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      const BowtieHit& bh = hits.hits[i];
      junctions_from_alignment(bh, junctions);
    }
}

// Extracts junctions from all the SAM hits (based on REF_SKIPs) in the hit file
// resets the stream when finished.
void exclude_hits_on_filtered_junctions(const JunctionSet& junctions, HitsForRead& hits)
{
  HitsForRead remaining;
  remaining.insert_id = hits.insert_id;
  
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      BowtieHit& bh = hits.hits[i];
      if (bh.mismatches() > read_mismatches ||
	  bh.gap_length() > read_gap_length ||
	  bh.edit_dist() > read_edit_dist)
	continue;
      
      bool filter_hit = false;
      if (!bh.contiguous())
	{
	  JunctionSet bh_juncs;
	  junctions_from_alignment(bh, bh_juncs);
	  for (JunctionSet::iterator itr = bh_juncs.begin();
	       itr != bh_juncs.end();
	       itr++)
	    {
	      const Junction& j = itr->first;
	      JunctionSet::const_iterator target = junctions.find(j);
	      if (target == junctions.end() || !target->second.accepted)
		{
		  filter_hit = true;
		  break;
		}
	    }
	}
      if (!filter_hit)
	remaining.hits.push_back(bh);
    }
  hits = remaining;
}

void realign_reads(HitsForRead& hits,
		   const RefSequenceTable& rt,
		   const JunctionSet& junctions,
		   const JunctionSet& rev_junctions,
		   const InsertionSet& insertions,
		   const DeletionSet& deletions,
		   const DeletionSet& rev_deletions,
		   const FusionSet& fusions)
{
  if (color)
    return;

  vector<BowtieHit> additional_hits;
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      BowtieHit& bh = hits.hits[i];
      if (fusion_search && bh.fusion_opcode() != FUSION_NOTHING)
	return;

      const vector<CigarOp>& cigars = bh.cigar();
      int pos = bh.left();
      int refid = bh.ref_id();
      for (size_t j = 0; j < cigars.size(); ++j)
	{
	  const CigarOp& op = cigars[j];
	  if (j == 0 || j == cigars.size() - 1)
	    {
	      // let's do this for MATCH case only,
	      if (op.opcode == MATCH || op.opcode == mATCH)
		{
		  int left1, left2;
		  if (op.opcode == MATCH)
		    {
		      left1 = pos;
		      left2 = pos + op.length - 1;
		    }
		  else
		    {
		      left1 = pos - op.length + 1;
		      left2 = pos;
		    }

		  {
		    static const size_t max_temp_juncs = 5;
		    
		    JunctionSet::const_iterator lb, ub;
		    JunctionSet temp_junctions;
		    if (j == 0)
		      {
			lb = rev_junctions.upper_bound(Junction(refid, left1, 0, true));
			ub = rev_junctions.lower_bound(Junction(refid, left2, left2, false));
			while (lb != ub && lb != rev_junctions.end())
			  {
			    Junction temp_junction = lb->first;
			    temp_junction.left = lb->first.right;
			    temp_junction.right = lb->first.left;
			    temp_junctions[temp_junction] = lb->second;
			    ++lb;

			    if (temp_junctions.size() > max_temp_juncs)
			      break;
			  }
		      }
		    
		    if (j == cigars.size() - 1)
		      {
			int common_right = left2 + max_report_intron_length;
			lb = junctions.upper_bound(Junction(refid, left1, common_right, true));
			ub = junctions.lower_bound(Junction(refid, left2, common_right, false));
			while (lb != ub && lb != junctions.end())
			  {
			    temp_junctions[lb->first] = lb->second;
			    ++lb;

			    if (temp_junctions.size() > max_temp_juncs)
			      break;
			  }
		      }

		    // daehwan - for debugging purposes
		    /*
		    if (bh.insert_id() == 15461 && cigars.size() == 1)
		      {
			printf("%d: %s\n",
			       bh.insert_id(), print_cigar(bh.cigar()).c_str());

			printf("candidate junctions: %d - max junctions: %d\n",
			       temp_junctions.size(),
			       max_temp_juncs);

			JunctionSet::const_iterator junc_iter = temp_junctions.begin();
			for (; junc_iter != temp_junctions.end(); ++junc_iter)
			  {
			    Junction junc = junc_iter->first;
			    fprintf(stderr, "%d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
				    bh.insert_id(), bh.left(), bh.right(),
				    print_cigar(bh.cigar()).c_str(),
				    bh.alignment_score(), bh.edit_dist(),
				    junc.left, junc.right);
			  }
		      }
		    */

		    if (temp_junctions.size() > max_temp_juncs)
		      continue;

		    JunctionSet::const_iterator junc_iter = temp_junctions.begin();
		    for (; junc_iter != temp_junctions.end(); ++junc_iter)
		      {
			Junction junc = junc_iter->first;
			
#if 0
			fprintf(stderr, "%d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
				bh.insert_id(), bh.left(), bh.right(),
				print_cigar(bh.cigar()).c_str(),
				bh.alignment_score(), bh.edit_dist(),
				junc.left, junc.right);
#endif
			
			int new_left = bh.left();
			int intron_length = junc.right - junc.left - 1;
			vector<CigarOp> new_cigars;
			bool anchored = false;
			if (j == 0 && bh.left() > (int)junc.left)
			  {
			    new_left -= intron_length;
			    int before_match_length = junc.left - new_left + 1;;
			    int after_match_length = op.length - before_match_length;
			    
			    if (before_match_length > 0 && after_match_length > 0)
			      {
				anchored = true;
				new_cigars.push_back(CigarOp(MATCH, before_match_length));
				new_cigars.push_back(CigarOp(REF_SKIP, intron_length));
				new_cigars.push_back(CigarOp(MATCH, after_match_length));
				
				new_cigars.insert(new_cigars.end(), cigars.begin() + 1, cigars.end());
			      }
			  }
			else if (j == cigars.size() - 1 && pos < (int)junc.left)
			  {
			    new_cigars.insert(new_cigars.end(), cigars.begin(), cigars.end() - 1);
			    
			    int before_match_length = junc.left - pos + 1;
			    int after_match_length = op.length - before_match_length;
			    
			    if (before_match_length > 0 && after_match_length > 0)
			      {
				anchored = true;
				new_cigars.push_back(CigarOp(MATCH, before_match_length));
				new_cigars.push_back(CigarOp(REF_SKIP, intron_length));
				new_cigars.push_back(CigarOp(MATCH, after_match_length));
			      }
			  }

			if (!anchored)
			  continue;

			BowtieHit new_bh(bh.ref_id(),
					 bh.ref_id2(),
					 bh.insert_id(), 
					 new_left,  
					 new_cigars,
					 bh.antisense_align(),
					 junc.antisense,
					 0, /* mismatches - needs to be recalculated */
					 0, /* edit_dist - needs to be recalculated */
					 0, /* splice_mms - needs to be recalculated */
					 false);

			new_bh.seq(bh.seq());
			new_bh.qual(bh.qual());

			const RefSequenceTable::Sequence* ref_str = rt.get_seq(bh.ref_id());

			if (new_left >= 0 && new_bh.right() <= (int)length(*ref_str))
			  {
			    vector<string> aux_fields;
			    bowtie_sam_extra(new_bh, rt, aux_fields);
			    
			    vector<string>::const_iterator aux_iter = aux_fields.begin();
			    for (; aux_iter != aux_fields.end(); ++aux_iter)
			      {
				const string& aux_field = *aux_iter;
				if (strncmp(aux_field.c_str(), "AS", 2) == 0)
				  {
				    int alignment_score = atoi(aux_field.c_str() + 5);
				    new_bh.alignment_score(alignment_score);
				  }
				else if (strncmp(aux_field.c_str(), "XM", 2) == 0)
				  {
				    int XM_value =  atoi(aux_field.c_str() + 5);
				    new_bh.mismatches(XM_value);
				    new_bh.edit_dist(XM_value + gap_length(new_cigars));
				  }
			      }

			    string NM = "NM:i:";
			    str_appendInt(NM, new_bh.edit_dist());
			    aux_fields.push_back(NM);

			    // replace the previous sam auxiliary fields with the new ones
			    vector<string> sam_toks;
			    tokenize(bh.hitfile_rec().c_str(), "\t", sam_toks);
			    
			    char coord[20] = {0,};
			    sprintf(coord, "%d", new_bh.left() + 1);
			    sam_toks[3] = coord;
			    sam_toks[5] = print_cigar(new_bh.cigar());
			    for (size_t a = 11; a < sam_toks.size(); ++a)
			      {
				string& sam_tok = sam_toks[a];
				for (size_t b = 0; b < aux_fields.size(); ++b)
				  {
				    const string& aux_tok = aux_fields[b];
				    if (strncmp(sam_tok.c_str(), aux_tok.c_str(), 5) == 0)
				      {
					sam_tok = aux_tok;
					break;
				      }
				  }
			      }

			    if (!bh.is_spliced())
			      {
				if (junc.antisense)
				  sam_toks.push_back("XS:A:-");
				else
				  sam_toks.push_back("XS:A:+");
			      }
			    
			    string new_rec = "";
			    for (size_t d = 0; d < sam_toks.size(); ++d)
			      {
				new_rec += sam_toks[d];
				if (d < sam_toks.size() - 1)
				  new_rec += "\t";
			      }
			    
			    new_bh.hitfile_rec(new_rec);
			    
			    if (new_bh.edit_dist() <= bh.edit_dist())
			      additional_hits.push_back(new_bh);

#if 0
			    fprintf(stderr, "\t%d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
				    new_bh.insert_id(), new_bh.left(), new_bh.right(),
				    print_cigar(new_bh.cigar()).c_str(),
				    new_bh.alignment_score(), new_bh.edit_dist(),
				    junc.left, junc.right);
#endif
			  }
		      }
		  }

#if 0
		  {
		    DeletionSet::const_iterator lb, ub;
		    bool use_rev_deletions = (j == 0);
		    const DeletionSet& curr_deletions = (use_rev_deletions ? rev_deletions : deletions);
		    if (use_rev_deletions)
		      {
			lb = curr_deletions.upper_bound(Deletion(refid, left1, 0, true));
			ub = curr_deletions.lower_bound(Deletion(refid, left2, left2, false));
		      }
		    else
		      {
			int common_right = left2 + 100;
			lb = curr_deletions.upper_bound(Deletion(refid, left1, common_right, true));
			ub = curr_deletions.lower_bound(Deletion(refid, left2, common_right, false));
		      }
		  
		    while (lb != curr_deletions.end() && lb != ub)
		      {
			Deletion del = lb->first;
			if (use_rev_deletions)
			  {
			    int temp = del.left;
			    del.left = del.right;
			    del.right = temp;
			  }		      
			
			// daehwan - for debuggin purposes
			/*
			  fprintf(stderr, "(type%d) %d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
			  !use_rev_junctions,
			  bh.insert_id(), bh.left(), bh.right(),
			  print_cigar(bh.cigar()).c_str(),
			  bh.alignment_score(), bh.edit_dist(),
			  junc.left, junc.right);
			*/
			
			int del_length = del.right - del.left - 1;
			int new_left = bh.left();
			if (j == 0)
			  new_left -= del_length;
			
			vector<CigarOp> new_cigars;
			if (j == 0)
			  {
			    int before_match_length = del.left - new_left + 1;;
			    int after_match_length = op.length - before_match_length;
			    
			    if (before_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, before_match_length));
			    new_cigars.push_back(CigarOp(DEL, del_length));
			    if (after_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, after_match_length));
			    
			    new_cigars.insert(new_cigars.end(), cigars.begin() + 1, cigars.end());
			  }
			else
			  {
			    new_cigars.insert(new_cigars.end(), cigars.begin(), cigars.end() - 1);
			    
			    int before_match_length = del.left - pos + 1;
			    int after_match_length = op.length - before_match_length;
			    
			    if (before_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, before_match_length));
			    new_cigars.push_back(CigarOp(DEL, del_length));
			    if (after_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, after_match_length));
			  }
			
			BowtieHit new_bh(bh.ref_id(),
					 bh.ref_id2(),
					 bh.insert_id(), 
					 new_left,  
					 new_cigars,
					 bh.antisense_align(),
					 bh.antisense_splice(),
					 0, /* edit_dist - needs to be recalculated */
					 0, /* splice_mms - needs to be recalculated */
					 false);
			
			new_bh.seq(bh.seq());
			new_bh.qual(bh.qual());
			
			vector<string> aux_fields;
			bowtie_sam_extra(new_bh, rt, aux_fields);
		      
			vector<string>::const_iterator aux_iter = aux_fields.begin();
			for (; aux_iter != aux_fields.end(); ++aux_iter)
			  {
			    const string& aux_field = *aux_iter;
			    if (strncmp(aux_field.c_str(), "AS", 2) == 0)
			      {
				int alignment_score = atoi(aux_field.c_str() + 5);
				new_bh.alignment_score(alignment_score);
			      }
			    else if (strncmp(aux_field.c_str(), "XM", 2) == 0)
			      {
				int XM_value =  atoi(aux_field.c_str() + 5);
				new_bh.edit_dist(XM_value);
			      }
			  }

			vector<string> sam_toks;
			tokenize(bh.hitfile_rec().c_str(), "\t", sam_toks);
			
			char coord[20] = {0,};
			sprintf(coord, "%d", new_bh.left() + 1);
			sam_toks[3] = coord;
			sam_toks[5] = print_cigar(new_bh.cigar());
			for (size_t a = 11; a < sam_toks.size(); ++a)
			  {
			    string& sam_tok = sam_toks[a];
			    for (size_t b = 0; b < aux_fields.size(); ++b)
			      {
				const string& aux_tok = aux_fields[b];
				if (strncmp(sam_tok.c_str(), aux_tok.c_str(), 5) == 0)
				  {
				    sam_tok = aux_tok;
				    break;
				  }
			      }
			  }

			string new_rec = "";
			for (size_t d = 0; d < sam_toks.size(); ++d)
			  {
			    new_rec += sam_toks[d];
			    if (d < sam_toks.size() - 1)
			      new_rec += "\t";
			  }
			
			new_bh.hitfile_rec(new_rec);
			
			if (new_bh.edit_dist() <= bh.edit_dist())
			  additional_hits.push_back(new_bh);
			
			/*
			  fprintf(stderr, "\t%d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
			  new_bh.insert_id(), new_bh.left(), new_bh.right(),
			  print_cigar(new_bh.cigar()).c_str(),
			  new_bh.alignment_score(), new_bh.edit_dist(),
			  junc.left, junc.right);
			*/
			
			++lb;
		      }
		  }

		  {
		    InsertionSet::const_iterator lb, ub;
		    lb = insertions.upper_bound(Insertion(refid, left1, ""));
		    ub = insertions.lower_bound(Insertion(refid, left2, ""));
		  
		    while (lb != insertions.end() && lb != ub)
		      {
			Insertion ins = lb->first;
			
			// daehwan - for debugging purposes
			/*
			  fprintf(stderr, "(type%d) %d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
			  !use_rev_junctions,
			  bh.insert_id(), bh.left(), bh.right(),
			  print_cigar(bh.cigar()).c_str(),
			  bh.alignment_score(), bh.edit_dist(),
			  junc.left, junc.right);
			*/

			int ins_length = ins.sequence.length();
			int new_left = bh.left();
			if (j == 0)
			  new_left -= ins_length;
			
			vector<CigarOp> new_cigars;
			if (j == 0)
			  {
			    int before_match_length = ins.left - new_left + 1;;
			    int after_match_length = op.length - before_match_length - ins_length;
			    
			    if (before_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, before_match_length));
			    new_cigars.push_back(CigarOp(INS, ins_length));
			    if (after_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, after_match_length));
			    
			    new_cigars.insert(new_cigars.end(), cigars.begin() + 1, cigars.end());
			  }
			else
			  {
			    new_cigars.insert(new_cigars.end(), cigars.begin(), cigars.end() - 1);
			    
			    int before_match_length = ins.left - pos + 1;
			    int after_match_length = op.length - before_match_length - ins_length;
			    
			    if (before_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, before_match_length));
			    new_cigars.push_back(CigarOp(INS, ins_length));
			    if (after_match_length > 0)
			      new_cigars.push_back(CigarOp(MATCH, after_match_length));
			  }
			
			BowtieHit new_bh(bh.ref_id(),
					 bh.ref_id2(),
					 bh.insert_id(), 
					 new_left,  
					 new_cigars,
					 bh.antisense_align(),
					 bh.antisense_splice(),
					 0, /* edit_dist - needs to be recalculated */
					 0, /* splice_mms - needs to be recalculated */
					 false);
			
			new_bh.seq(bh.seq());
			new_bh.qual(bh.qual());
			
			vector<string> aux_fields;
			bowtie_sam_extra(new_bh, rt, aux_fields);
		      
			vector<string>::const_iterator aux_iter = aux_fields.begin();
			for (; aux_iter != aux_fields.end(); ++aux_iter)
			  {
			    const string& aux_field = *aux_iter;
			    if (strncmp(aux_field.c_str(), "AS", 2) == 0)
			      {
				int alignment_score = atoi(aux_field.c_str() + 5);
				new_bh.alignment_score(alignment_score);
			      }
			    else if (strncmp(aux_field.c_str(), "XM", 2) == 0)
			      {
				int XM_value =  atoi(aux_field.c_str() + 5);
				new_bh.edit_dist(XM_value);
			      }
			  }
			
			/*
			  fprintf(stderr, "\t%d %d-%d %s (AS:%d XM:%d) with junc %u-%u\n",
			  new_bh.insert_id(), new_bh.left(), new_bh.right(),
			  print_cigar(new_bh.cigar()).c_str(),
			  new_bh.alignment_score(), new_bh.edit_dist(),
			  junc.left, junc.right);
			*/
			
			++lb;
		      }
		  }
#endif
		}
	    }
	  
	  switch(op.opcode)
	    {
	    case REF_SKIP:
	      pos += op.length;
	      break;
	    case rEF_SKIP:
	      pos -= op.length;
	      break;
	    case MATCH:
	    case DEL:
	      pos += op.length;
	      break;
	    case mATCH:
	    case dEL:
	      pos -= op.length;
	      break;
	    case FUSION_FF:
	    case FUSION_FR:
	    case FUSION_RF:
	    case FUSION_RR:
	      pos = op.length;
	      refid = bh.ref_id2();
	      break;
	    default:
	      break;
	    }
	}
    }

  hits.hits.insert(hits.hits.end(), additional_hits.begin(), additional_hits.end());

  std::sort(hits.hits.begin(), hits.hits.end());
  vector<BowtieHit>::iterator new_end = std::unique(hits.hits.begin(), hits.hits.end());
  hits.hits.erase(new_end, hits.hits.end());  
}

class MultipleBAMReader
{
public:
	MultipleBAMReader(ReadTable& it,
			RefSequenceTable& rt,
			const vector<string>& fnames,
			long begin_id,
			long end_id) :
				_it(it),
				_rt(rt),
				_begin_id(begin_id),
				_end_id(end_id),
				_bam_hit_factory(it, rt),
				_bam_merge(NULL)
	{
			// calculate file offsets
			vector<int64_t> offsets;
			for (size_t i = 0; i < fnames.size(); ++i)
			{
				const string& fname = fnames[i];

				vector<uint64_t> next_file_read_ids;

				vector<string> temp_fnames;
				if (fname.substr(fname.length() - 4) == ".bam")
				{
					temp_fnames.push_back(fname);
					next_file_read_ids.push_back(0);
				}
				else
				{
					size_t j = 0;
					while (true)
					{
						char suffix[128];
						sprintf(suffix, "%lu.bam", j);
						string temp_fname = fname + suffix;
						string temp_index_fname = temp_fname + ".index";

						ifstream index_file(temp_index_fname.c_str());
						if (!index_file.is_open())
						{
							next_file_read_ids.push_back(0);
							break;
						}

						temp_fnames.push_back(temp_fname);

						if (j > 0)
						{
							string line;
							int64_t offset = 0;
							uint64_t read_id = 0;
							if (getline(index_file, line))
							{
								istringstream istream(line);
								istream >> read_id >> offset;
								next_file_read_ids.push_back(read_id);
							}
							else
							{
								next_file_read_ids.push_back(0);
							}
						}

						++j;
					}
				}

				for (size_t j = 0; j < temp_fnames.size(); ++j)
				{
					ifstream reads_index_file((temp_fnames[j] + ".index").c_str());
					if (!reads_index_file.is_open())
						continue;

					bool pushed = false;

					int64_t offset = 0, last_offset = 0;
					uint64_t read_id = 0, last_read_id = 0;
					string line;
					while (getline(reads_index_file, line))
					{
						istringstream istream(line);
						istream >> read_id >> offset;
						if (read_id > _begin_id && last_read_id <= _begin_id)
						{
							pushed = true;
							_fnames.push_back(temp_fnames[j]);
							offsets.push_back(last_offset);

	#if 0
							fprintf(stderr, "bet %lu and %lu - %s %lu %ld\n",
									_begin_id, _end_id, temp_fnames[j].c_str(), last_offset, last_read_id);
	#endif

							break;
						}

						last_offset = offset;
						last_read_id = read_id;
					}

					if (!pushed)
					{
						if(next_file_read_ids[j] > _begin_id && last_read_id <= _begin_id)
						{
							pushed = true;
							_fnames.push_back(temp_fnames[j]);
							offsets.push_back(last_offset);

	#if 0
							fprintf(stderr, "2 bet %lu and %lu - %s %lu %ld\n",
									_begin_id, _end_id, temp_fnames[j].c_str(), last_offset, last_read_id);
	#endif
						}
					}

					if (read_id >= _end_id)
						break;

					if (read_id == 0)
					{
						_fnames.push_back(temp_fnames[j]);
						offsets.push_back(0);
					}
				}
			}

			_bam_merge = new BamMerge(_fnames, offsets);
			_bam_hit_factory.set_sam_header(_bam_merge->get_sam_header());
	}

	~MultipleBAMReader()
	{
		if (_bam_merge)
			delete _bam_merge;
	}

	bool next_read_hits(HitsForRead& hits)
	{
		hits.insert_id = 0;

		if (!_bam_merge)
			return false;

		vector<CBamLine> bam_lines;
		while (true)
		{
			if (!_bam_merge->next_bam_lines(bam_lines))
				return false;

			if (bam_lines.size() <= 0)
				return false;

			uint64_t read_id = bam_lines[0].read_id;
			if (read_id >= _begin_id && read_id < _end_id)
			{
				hits.hits.clear();
				for (size_t i = 0; i < bam_lines.size(); ++i)
				{
					CBamLine& bam_line = bam_lines[i];
					BowtieHit bh;

					char seq[MAX_READ_LEN + 1] = {0};
					char qual[MAX_READ_LEN + 1] = {0};

					bool success = _bam_hit_factory.get_hit_from_buf((const char*)bam_line.b, bh, true, NULL, NULL, seq, qual);
					if (success)
					{
						bh.seq(seq);
						bh.qual(qual);

						char* sam_line = bam_format1(_bam_merge->get_sam_header(), bam_line.b);
						bh.hitfile_rec(sam_line);
						free(sam_line);

						hits.insert_id = bh.insert_id();
						hits.hits.push_back(bh);
					}

					bam_line.b_free();
				}
				bam_lines.clear();

				if (hits.hits.size() > 0)
					return true;
			}

			for (size_t i = 0; i < bam_lines.size(); ++i)
				bam_lines[i].b_free();

			bam_lines.clear();

			if (read_id >= _end_id)
				break;
		}

		return false;
	}

private:
	ReadTable& _it;
	RefSequenceTable& _rt;
	vector<string> _fnames;
	uint64_t _begin_id;
	uint64_t _end_id;

	BAMHitFactory _bam_hit_factory;
	BamMerge* _bam_merge;
};


// events include splice junction, indels, and fusion points.
struct ConsensusEventsWorker
{
	void operator()()
	{
		ReadTable it;
		MultipleBAMReader l_hs(it, *rt, left_map_fnames, begin_id, end_id);
		MultipleBAMReader r_hs(it, *rt, right_map_fnames, begin_id, end_id);

		HitsForRead curr_left_hit_group;
		HitsForRead curr_right_hit_group;

		l_hs.next_read_hits(curr_left_hit_group);
		r_hs.next_read_hits(curr_right_hit_group);

		uint32_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		uint32_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);

		// While we still have unreported hits...
		while((curr_left_obs_order != VMAXINT32 || curr_right_obs_order != VMAXINT32) &&
				(curr_left_obs_order < end_id || curr_right_obs_order < end_id))
		{
			// Chew up left singletons
			while (curr_left_obs_order < curr_right_obs_order &&
					curr_left_obs_order < end_id &&
					curr_left_obs_order != VMAXINT32)
			{
				HitsForRead best_hits;
				best_hits.insert_id = curr_left_obs_order;

				// Process hits for left singleton, select best alignments
				read_best_alignments(curr_left_hit_group, best_hits, *gtf_junctions);
				// update_coverage(best_hits, *coverage);
				update_junctions(best_hits, *junctions);
				update_insertions_and_deletions(best_hits, *insertions, *deletions);
				update_fusions(best_hits, *rt, *fusions);

				// Get next hit group
				l_hs.next_read_hits(curr_left_hit_group);
				curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			}

			// Chew up right singletons
			while (curr_left_obs_order > curr_right_obs_order &&
					curr_right_obs_order < end_id &&
					curr_right_obs_order != VMAXINT32)
			{
				HitsForRead best_hits;
				best_hits.insert_id = curr_right_obs_order;

				if (curr_right_obs_order >= begin_id)
				{
					// Process hit for right singleton, select best alignments
					read_best_alignments(curr_right_hit_group, best_hits, *gtf_junctions);
					// update_coverage(best_hits, *coverage);

					update_junctions(best_hits, *junctions);
					update_insertions_and_deletions(best_hits, *insertions, *deletions);
					update_fusions(best_hits, *rt, *fusions);
				}

				// Get next hit group
				r_hs.next_read_hits(curr_right_hit_group);
				curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
			}

			// Since we have both left hits and right hits for this insert,
			// Find the best pairing and print both
			while (curr_left_obs_order == curr_right_obs_order &&
					curr_left_obs_order < end_id &&
					curr_left_obs_order != VMAXINT32)
			{
				vector<pair<BowtieHit, BowtieHit> > best_hits;

				InsertAlignmentGrade grade;
				pair_best_alignments(curr_left_hit_group,
						curr_right_hit_group,
						grade,
						best_hits,
						*gtf_junctions);

				HitsForRead best_left_hit_group; best_left_hit_group.insert_id = curr_left_obs_order;
				HitsForRead best_right_hit_group; best_right_hit_group.insert_id = curr_left_obs_order;

				if (best_hits.size() > 0)
				{
					for (size_t i = 0; i < best_hits.size(); ++i)
					{
						best_left_hit_group.hits.push_back(best_hits[i].first);
						best_right_hit_group.hits.push_back(best_hits[i].second);
					}
				}
				else
				{
					best_left_hit_group.hits = curr_left_hit_group.hits;
					best_right_hit_group.hits = curr_right_hit_group.hits;
				}

				// update_coverage(best_left_hit_group, *coverage);
				update_junctions(best_left_hit_group, *junctions);
				update_insertions_and_deletions(best_left_hit_group, *insertions, *deletions);
				update_fusions(best_left_hit_group, *rt, *fusions);

				// update_coverage(best_right_hit_group, *coverage);
				update_junctions(best_right_hit_group, *junctions);
				update_insertions_and_deletions(best_right_hit_group, *insertions, *deletions);
				update_fusions(best_right_hit_group, *rt, *fusions);

				l_hs.next_read_hits(curr_left_hit_group);
				curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);

				r_hs.next_read_hits(curr_right_hit_group);
				curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
			}
		}
	}

	vector<string> left_map_fnames;
	vector<string> right_map_fnames;

	RefSequenceTable* rt;

	JunctionSet* gtf_junctions;

	uint64_t begin_id;
	uint64_t end_id;
	int64_t left_map_offset;
	int64_t right_map_offset;

	JunctionSet* junctions;
	InsertionSet* insertions;
	DeletionSet* deletions;
	FusionSet* fusions;

	Coverage* coverage;
};
void print_alnStats(SAlignStats& alnStats) {
	string fname(output_dir);
	fname+="/align_summary.txt";
	FILE* f = fopen(fname.c_str(), "w");
	int64_t total_left=alnStats.num_aligned_left+alnStats.num_unmapped_left;
	int64_t total_right=alnStats.num_aligned_right+alnStats.num_unmapped_right;
	int64_t total_unpaired=alnStats.num_aligned_unpaired+alnStats.num_unmapped_unpaired;
	//int64_t accepted_left =alnStats.num_aligned_left -alnStats.num_aligned_left_xmulti; //accepted mappings, < max_multihits
	//int64_t accepted_right=alnStats.num_aligned_right-alnStats.num_aligned_right_xmulti; //accepted right mappings
	string rdn("Left reads");
	if (total_right==0) rdn="Reads";
	fprintf(f, "%s:\n", rdn.c_str());
	fprintf(f, "               Input: %9ld\n", total_left);
	double perc=(100.0*alnStats.num_aligned_left)/total_left;
	fprintf(f, "              Mapped: %9ld (%4.1f%% of input)\n",
			alnStats.num_aligned_left, perc);
	if (alnStats.num_aligned_left) {
		perc=(100.0*alnStats.num_aligned_left_multi)/alnStats.num_aligned_left;
	 fprintf(f,"            of these: %9ld (%4.1f%%) have multiple alignments (%ld have >%d)\n",
			alnStats.num_aligned_left_multi, perc, alnStats.num_aligned_left_xmulti, max_multihits);
	}
  /*perc=(100.0*accepted_left)/total_left;
	fprintf(f, "   Mapped acceptably: %9ld (%4.1f%% of input)\n",
			accepted_left, perc);
  */
	int64_t total_mapped=alnStats.num_aligned_left;
	int64_t total_input=total_left;
	int64_t total_pairs=0;
	if (total_right) {
		fprintf(f, "Right reads:\n");
		fprintf(f, "               Input: %9ld\n", total_right);
		total_input+=total_right;
		perc=(100.0*alnStats.num_aligned_right)/total_right;
		fprintf(f, "              Mapped: %9ld (%4.1f%% of input)\n",
				alnStats.num_aligned_right, perc);
		if (alnStats.num_aligned_right) {
			perc=(100.0* alnStats.num_aligned_right_multi)/alnStats.num_aligned_right;
			fprintf(f,"            of these: %9ld (%4.1f%%) have multiple alignments (%ld have >%d)\n",
					alnStats.num_aligned_right_multi, perc, alnStats.num_aligned_right_xmulti, max_multihits);
		}
		total_mapped+=alnStats.num_aligned_right;
		total_pairs=(total_right<total_left)? total_right : total_left;

		if (total_unpaired) {
				fprintf(f, "Unpaired reads:\n");
				fprintf(f, "               Input: %9ld\n", total_unpaired);
				total_input+=total_unpaired;
				perc=(100.0*alnStats.num_aligned_unpaired)/total_unpaired;
				fprintf(f, "              Mapped: %9ld (%4.1f%% of input)\n",
						alnStats.num_aligned_unpaired, perc);
				if (alnStats.num_aligned_unpaired) {
					perc=(100.0* alnStats.num_aligned_unpaired_multi)/alnStats.num_aligned_unpaired;
					fprintf(f,"            of these: %9ld (%4.1f%%) have multiple alignments (%ld have >%d)\n",
							alnStats.num_aligned_unpaired_multi, perc, alnStats.num_aligned_unpaired_xmulti, max_multihits);
				}
		total_mapped+=alnStats.num_aligned_unpaired;
		}
	}
	perc=(100.0*total_mapped)/total_input;
	fprintf(f, "%4.1f%% overall read alignment rate.\n", perc);
	if (alnStats.num_aligned_pairs) {
		fprintf(f, "\nAligned pairs: %9ld\n", alnStats.num_aligned_pairs);
		perc=(100.0*alnStats.num_aligned_pairs_multi)/alnStats.num_aligned_pairs;
		fprintf(f, "     of these: %9ld (%4.1f%%) have multiple alignments\n",
				alnStats.num_aligned_pairs_multi, perc);
		perc=(100.0*alnStats.num_aligned_pairs_disc)/alnStats.num_aligned_pairs;
		fprintf(f, "          and: %9ld (%4.1f%%) are discordant alignments\n",
				alnStats.num_aligned_pairs_disc, perc);
		perc=(100.0*(alnStats.num_aligned_pairs-alnStats.num_aligned_pairs_disc))/total_pairs;
		fprintf(f, "%4.1f%% concordant pair alignment rate.\n", perc);
	}
	fclose(f);
}

struct CReadProc: public GetReadProc {
	CReadProc(GBamWriter* bamw=NULL, int64_t* um_counter=NULL, int64_t* mm_counter=NULL,
			int64_t* u_um_counter=NULL, int64_t* u_mm_counter=NULL):
		GetReadProc(bamw, um_counter, mm_counter, u_um_counter, u_mm_counter) {}
	virtual bool process(QReadData& rdata, bool& found) { //, bool is_unmapped) {
		//if (is_unmapped || (rdata.trashCode!='M' && !found)) {
		if (rdata.trashCode!='M' && !found) {
			this->writeUnmapped(rdata);
		}
		//
		if (unmapped_counter && !found) {
			if (rdata.trashCode!='M') {
				if (rdata.matenum) (*unmapped_counter)++;
				else if (u_unmapped_counter) (*u_unmapped_counter)++;
			}
			else if (multimapped_counter) {
				if (rdata.matenum) (*multimapped_counter)++;
				else if (u_multimapped_counter) (*u_multimapped_counter)++;
			}
		}
		return true;
	}
};

struct ReportWorker {

	ReportWorker(RefSequenceTable* r=NULL, SAlignStats* s=NULL):  gtf_junctions(NULL), junctions(NULL), rev_junctions(NULL),
			insertions(NULL), deletions(NULL), rev_deletions(NULL), fusions(NULL), coverage(NULL),
			final_junctions(NULL), final_insertions(NULL), final_deletions(NULL), final_fusions(NULL),
			rt(r), begin_id(0), end_id(0),  left_reads_offset(0), left_map_offset(0), right_reads_offset(0),
			right_map_offset(0), alnStats(s), PE_reads(false) {  }

	void write_singleton_alignments(uint64_t curr_obs_order,
	    HitsForRead& curr_hit_group, ReadStream& reads_file,
	    GBamWriter& bam_writer, FragmentType fragment_type,
	    GetReadProc* readProc=NULL,
	    QReadData* gotRead=NULL) {
		int64_t* unmapped_counter = & alnStats->num_unmapped_left;
		int64_t* aligned_counter = & alnStats->num_aligned_left;
		int64_t* aligned_counter_multi = & alnStats->num_aligned_left_multi;
		int64_t* aligned_counter_xmulti = & alnStats->num_aligned_left_xmulti;

		QReadData read;
		if (gotRead==NULL) {
			if (reads_file.getQRead(curr_obs_order, read,
			                   begin_id, end_id, readProc))
                  gotRead=&read;
		}
		if (gotRead==NULL) {
			//Should never happen!
			warn_msg("Error: singleton getQRead() failed for id# %ld.\n", curr_obs_order);
			return;
		}

		if (gotRead->matenum==0) {
			 fragment_type = FRAG_UNPAIRED;
			 if (PE_reads) {
				 unmapped_counter = & alnStats->num_unmapped_unpaired;
				 aligned_counter = & alnStats->num_aligned_unpaired;
				 aligned_counter_multi = & alnStats->num_aligned_unpaired_multi;
				 aligned_counter_xmulti = & alnStats->num_aligned_unpaired_xmulti;
			 }
		} else if (fragment_type == FRAG_RIGHT) {
			unmapped_counter = & alnStats->num_unmapped_right;
			aligned_counter = & alnStats->num_aligned_right;
			aligned_counter_multi = & alnStats->num_aligned_right_multi;
			aligned_counter_xmulti = & alnStats->num_aligned_right_xmulti;
		}

		if (PE_reads && !report_mixed_alignments) {
			//the user only wants paired alignments reported
			/* if (!gotRead) {
				if (curr_hit_group.hits.size()>1) {
					(*aligned_counter_multi)++;
				}
				else {
					(*aligned_counter)++;
				}
			}
			*/
			if (readProc) readProc->writeUnmapped(*gotRead);
			(*unmapped_counter)++;
			return;
		}
		HitsForRead best_hits;
		best_hits.insert_id = curr_obs_order;
		realign_reads(curr_hit_group, *rt, *junctions, *rev_junctions, *insertions,
				*deletions, *rev_deletions, *fusions);

		exclude_hits_on_filtered_junctions(*junctions, curr_hit_group);

		// Process hits for singleton, select best alignments
		const bool final_report = true;
		char map_flags=read_best_alignments(curr_hit_group, best_hits, *gtf_junctions, *junctions,
				*insertions, *deletions, *fusions, *coverage, final_report, &rng);

		char map_code=0;
		if (map_flags & 4) {
			(*aligned_counter_multi)++;
		}
		if (map_flags & 16) {
			map_code='M';
			(*aligned_counter_xmulti)++;
		}
		if (map_flags & 1) (*aligned_counter)++; //acceptable mappings found
		else {
			if (!gotRead->um_written)
				(*unmapped_counter)++;
			if (readProc)
				readProc->writeUnmapped(*gotRead);
		}

		if (best_hits.hits.size() > 0) {
				update_junctions(best_hits, *final_junctions);
				update_insertions_and_deletions(best_hits, *final_insertions,
						*final_deletions);
				update_fusions(best_hits, *rt, *final_fusions, *fusions);
				if (!gotRead->um_written) {
					print_sam_for_single(*rt, best_hits, fragment_type, gotRead->read.alt_name,
							bam_writer, rng);
				}
		}
		else {
			readProc->writeUnmapped(*gotRead);
		}
	}

	void operator()() {
		rng.seed(1);

		ReadTable it;
		GBamWriter bam_writer(bam_output_fname.c_str(), sam_header.c_str());

		PE_reads = right_map_fnames.size() > 0;

		ReadStream left_reads_file(left_reads_fname);
		if (left_reads_file.file() == NULL)
			err_die("Error: cannot open %s for reading\n", left_reads_fname.c_str());
		if (left_reads_file.isBam()) {
			left_reads_file.use_alt_name();
			left_reads_file.ignoreQC();
		}
		if (left_reads_offset > 0)
			left_reads_file.seek(left_reads_offset);

		GBamWriter* left_um_out = new GBamWriter(left_um_fname.c_str(),
		    sam_header.c_str());
		GBamWriter* right_um_out = NULL;

		ReadStream right_reads_file(right_reads_fname);
		if (right_reads_offset > 0)
			right_reads_file.seek(right_reads_offset);

		if (!right_reads_fname.empty()) {
			if (right_reads_file.isBam()) {
				right_reads_file.use_alt_name();
				right_reads_file.ignoreQC();
				right_um_out = new GBamWriter(right_um_fname.c_str(),
				    sam_header.c_str());
			}
		}

		MultipleBAMReader left_hs(it, *rt, left_map_fnames, begin_id, end_id);
		MultipleBAMReader right_hs(it, *rt, right_map_fnames, begin_id, end_id);

		HitsForRead curr_left_hit_group;
		HitsForRead curr_right_hit_group;

		left_hs.next_read_hits(curr_left_hit_group);
		right_hs.next_read_hits(curr_right_hit_group);

		uint64_t curr_left_obs_order = it.observation_order(
		    curr_left_hit_group.insert_id);
		uint64_t curr_right_obs_order = it.observation_order(
		    curr_right_hit_group.insert_id);

		const bool final_report = true;

		// While we still have unreported hits...
		QReadData l_read;
		QReadData r_read;
		int64_t* num_unmapped_unpaired = (PE_reads) ? &(alnStats->num_unmapped_unpaired) :
				&(alnStats->num_unmapped_left);
		int64_t* num_aligned_unpaired_xmulti= (PE_reads) ?  &(alnStats->num_aligned_unpaired_xmulti) :
				&(alnStats->num_aligned_left_xmulti);
		CReadProc l_readProc(left_um_out, &(alnStats->num_unmapped_left),
		    &(alnStats->num_aligned_left_xmulti), num_unmapped_unpaired, num_aligned_unpaired_xmulti);
		CReadProc r_readProc(right_um_out, &(alnStats->num_unmapped_right),
		    &(alnStats->num_aligned_right_xmulti), &(alnStats->num_unmapped_unpaired),
		    &(alnStats->num_aligned_unpaired_xmulti));

		while ((curr_left_obs_order != VMAXINT32
		    || curr_right_obs_order != VMAXINT32)
		    && (curr_left_obs_order < end_id || curr_right_obs_order < end_id)) {

			/*if (curr_left_obs_order >= 3463 || curr_right_obs_order >= 3463) {
				fprintf(stderr, "Debug target reached!\n");
			}*/
			// Chew up left singletons (pairs with right reads unmapped)
			while (curr_left_obs_order < curr_right_obs_order
			    && curr_left_obs_order < end_id && curr_left_obs_order != VMAXINT32) {
				write_singleton_alignments(curr_left_obs_order, curr_left_hit_group,
				    left_reads_file, bam_writer, //*left_um_out,
				    PE_reads ? FRAG_LEFT : FRAG_UNPAIRED, &l_readProc);

				// Get next hit group
				left_hs.next_read_hits(curr_left_hit_group);
				curr_left_obs_order = it.observation_order(
				    curr_left_hit_group.insert_id);
			} //left singletons

			// Chew up right singletons
			while (curr_left_obs_order > curr_right_obs_order
			    && curr_right_obs_order < end_id && curr_right_obs_order != VMAXINT32) {
				write_singleton_alignments(curr_right_obs_order, curr_right_hit_group,
				    right_reads_file, bam_writer, FRAG_RIGHT, &r_readProc);

				// Get next hit group
				right_hs.next_read_hits(curr_right_hit_group);
				curr_right_obs_order = it.observation_order(
				    curr_right_hit_group.insert_id);
			}

			// Since we have both left hits and right hits for this insert,
			//  find the best pairing and print both alignments
			while (curr_left_obs_order == curr_right_obs_order
			    && curr_left_obs_order < end_id && curr_left_obs_order != VMAXINT32) {
				realign_reads(curr_left_hit_group, *rt, *junctions, *rev_junctions,
				    *insertions, *deletions, *rev_deletions, *fusions);
				exclude_hits_on_filtered_junctions(*junctions, curr_left_hit_group);

				realign_reads(curr_right_hit_group, *rt, *junctions, *rev_junctions,
				    *insertions, *deletions, *rev_deletions, *fusions);
				exclude_hits_on_filtered_junctions(*junctions, curr_right_hit_group);

				vector<pair<BowtieHit, BowtieHit> > best_hits;

				bool paired_alignments = curr_left_hit_group.hits.size() > 0
				    && curr_right_hit_group.hits.size() > 0;
				InsertAlignmentGrade grade;
				bool got_left_read = left_reads_file.getQRead(curr_left_obs_order,
						    	l_read, begin_id, end_id, &l_readProc);
				if (!got_left_read) {
					warn_msg("Error: failed to retrieve left read for pair # %ld !\n", curr_left_obs_order);
				}
				bool got_right_read = right_reads_file.getQRead(curr_right_obs_order, r_read,
					    begin_id, end_id, &r_readProc);
				if (!got_right_read) {
					warn_msg("Error: failed to retrieve right read for pair # %ld !\n", curr_right_obs_order);
				}
				char pair_map_flags=0;
				if (paired_alignments) { //there are *some* paired alignments
					pair_map_flags=pair_best_alignments(curr_left_hit_group, curr_right_hit_group, grade,
					    best_hits, *gtf_junctions, *junctions, *insertions, *deletions,
					    *fusions, *coverage, final_report, &rng);
					/*
					if (pair_map_flags & 1) alnStats->num_aligned_left++;
					else alnStats->num_unmapped_left++;
					if (pair_map_flags & 2) alnStats->num_aligned_right++;
					else alnStats->num_unmapped_right++;
					if ((pair_map_flags & 3)==3) //at least one acceptable alignment was found for each read
						alnStats->num_aligned_pairs++;
					char left_map_code=0;
					*/
					/*
					//if (pair_map_flags & 4) alnStats->num_aligned_left_multi++;
					//if (pair_map_flags & 8) alnStats->num_aligned_right_multi++;
					if ((pair_map_flags & 12) == 12) alnStats->num_aligned_pairs_multi++;

					if (pair_map_flags & 16) {
						//left_map_code='M';
						alnStats->num_aligned_left_xmulti++;
					}
					char right_map_code=0;
					if (pair_map_flags & 32) {
						//right_map_code='M';
						alnStats->num_aligned_right_xmulti++;
					}

					//l_read should be written as unmapped if (pair_map_flags & 1)==0
					if (got_left_read && (pair_map_flags & 1)==0) {
						l_readProc.writeUnmapped(l_read);
					}
					// r_read should be written as unmapped if (pair_map_flags & 2)==0
					if (got_right_read && (pair_map_flags & 2)==0) {
						r_readProc.writeUnmapped(r_read);
					}
					*/
					if (((pair_map_flags & 3)==3) && !grade.concordant()) {
						alnStats->num_aligned_pairs_disc++;
					}
					if (report_mixed_alignments) {
						if (best_hits.size() <= 0
						    || (grade.fusion && !fusion_search
						        && !report_discordant_pair_alignments))
							paired_alignments = false;
					}
				}
				if (paired_alignments) { //after some filtering, these paired reads still
					                     // have some alignments together
					HitsForRead best_left_hit_group;
					best_left_hit_group.insert_id = curr_left_obs_order;
					HitsForRead best_right_hit_group;
					best_right_hit_group.insert_id = curr_left_obs_order;

					for (size_t i = 0; i < best_hits.size(); ++i) {
						best_left_hit_group.hits.push_back(best_hits[i].first);
						best_right_hit_group.hits.push_back(best_hits[i].second);
					}

					if (best_hits.size() > 0) {
						if (got_left_read && got_right_read) { //should always be true
							update_junctions(best_left_hit_group, *final_junctions);
							update_insertions_and_deletions(best_left_hit_group,
							    *final_insertions, *final_deletions);
							update_fusions(best_left_hit_group, *rt, *final_fusions,
							    *fusions);

							update_junctions(best_right_hit_group, *final_junctions);
							update_insertions_and_deletions(best_right_hit_group,
							    *final_insertions, *final_deletions);
							update_fusions(best_right_hit_group, *rt, *final_fusions,
							    *fusions);

							pair_support(best_hits, *final_fusions, *fusions);

							print_sam_for_pair(*rt, best_hits, grade, bam_writer,
							    l_read.read.alt_name, r_read.read.alt_name, rng, begin_id, end_id);
							l_read.um_written=true;
							r_read.um_written=true;
						} /*
						else {
							fprintf(stderr, "Error: couldn't get reads for pair #%ld (%d, %d)\n",
									curr_left_obs_order, int(got_left_read), int(got_right_read));
						}*/
					}

					if (best_hits.size()>0) {
						alnStats->num_aligned_left++;
						alnStats->num_aligned_right++;
						alnStats->num_aligned_pairs++;
					}
					else {
						alnStats->num_unmapped_left++;
						alnStats->num_unmapped_right++;
					}
					if (pair_map_flags & 4) alnStats->num_aligned_left_multi++;
					if (pair_map_flags & 8) alnStats->num_aligned_right_multi++;
					if ((pair_map_flags & 12) == 12) alnStats->num_aligned_pairs_multi++;

					if (pair_map_flags & 16) {
						//left_map_code='M';
						alnStats->num_aligned_left_xmulti++;
					}
					if (pair_map_flags & 32) {
						//right_map_code='M';
						alnStats->num_aligned_right_xmulti++;
					}

					//l_read should be written as unmapped if (pair_map_flags & 1)==0
					if (got_left_read && (pair_map_flags & 1)==0) {
						l_readProc.writeUnmapped(l_read);
					}
					// r_read should be written as unmapped if (pair_map_flags & 2)==0
					if (got_right_read && (pair_map_flags & 2)==0) {
						r_readProc.writeUnmapped(r_read);
					}
				}
				else { //paired reads with alignments not paired properly
					if (curr_left_hit_group.hits.size() > 0) {
						write_singleton_alignments(curr_left_obs_order, curr_left_hit_group,
						    left_reads_file, bam_writer, //*left_um_out,
									   PE_reads ? FRAG_LEFT : FRAG_UNPAIRED, &l_readProc,
											   got_left_read ? &l_read : NULL);
					}
					else {
						l_readProc.writeUnmapped(l_read);
						alnStats->num_unmapped_left++;
					}

					if (curr_right_hit_group.hits.size() > 0) {   //only right read mapped
						write_singleton_alignments(curr_right_obs_order,
						    curr_right_hit_group, right_reads_file, bam_writer, //*right_um_out,
									   FRAG_RIGHT, &r_readProc, got_right_read ? &r_read : NULL);
					} else {
						r_readProc.writeUnmapped(r_read);
						alnStats->num_unmapped_right++;
					}
				}

				left_hs.next_read_hits(curr_left_hit_group);
				curr_left_obs_order = it.observation_order(
				    curr_left_hit_group.insert_id);

				right_hs.next_read_hits(curr_right_hit_group);
				curr_right_obs_order = it.observation_order(
				    curr_right_hit_group.insert_id);
			} //both mates have alignments
		} //while we still have unreported hits..

		//print the remaining unmapped reads at the end of each reads' stream

		left_reads_file.getQRead(VMAXINT32, l_read, begin_id,
		    end_id, &l_readProc);
		    //left_um_out, 0, &(alnStats->num_unmapped_left), &(alnStats->num_aligned_left_xmulti));
		if (right_reads_file.file())
			right_reads_file.getQRead(VMAXINT32, r_read, begin_id,
			    end_id, &r_readProc);
			    //right_um_out, 0, &(alnStats->num_unmapped_right), &(alnStats->num_aligned_right_xmulti));

		// pclose (pipe close), which waits for a process to end, seems to conflict with boost::thread::join somehow,
		// resulting in deadlock like behavior.
		delete left_um_out;
		delete right_um_out;
	}

	string bam_output_fname;
	string sam_header_fname;

	string left_reads_fname;
	vector<string> left_map_fnames;
	string right_reads_fname;
	vector<string> right_map_fnames;

	string left_um_fname;
	string right_um_fname;

	JunctionSet* gtf_junctions;

	JunctionSet* junctions;
	JunctionSet* rev_junctions;
	InsertionSet* insertions;
	DeletionSet* deletions;
	DeletionSet* rev_deletions;
	FusionSet* fusions;
	Coverage* coverage;

	JunctionSet* final_junctions;
	InsertionSet* final_insertions;
	DeletionSet* final_deletions;
	FusionSet* final_fusions;

	RefSequenceTable* rt;

	uint64_t begin_id;
	uint64_t end_id;
	int64_t left_reads_offset;
	int64_t left_map_offset;
	int64_t right_reads_offset;
	int64_t right_map_offset;
	//read alignment accounting:
	SAlignStats* alnStats;

	bool PE_reads;

	boost::mt19937 rng;
};

void driver(const string& bam_output_fname,
		istream& ref_stream,
		const vector<string>& left_map_fnames,
		const string& left_reads_fname,
		const vector<string>& right_map_fnames,
		const string& right_reads_fname,
		FILE* junctions_out,
		FILE* insertions_out,
		FILE* deletions_out,
		FILE* fusions_out)
{
	if (!parallel)
		num_threads = 1;

	RefSequenceTable rt(sam_header, true);
	get_seqs(ref_stream, rt, true);

	srandom(1);

	JunctionSet gtf_junctions;
	if (!gtf_juncs.empty())
	{
		char splice_buf[4096];
		FILE* splice_coords = fopen(gtf_juncs.c_str(), "r");
		if (splice_coords)
		{
			while (fgets(splice_buf, sizeof(splice_buf), splice_coords))
			{
				char* nl = strrchr(splice_buf, '\n');
				char* buf = splice_buf;
				if (nl) *nl = 0;

				char* ref_name = get_token((char**)&buf, "\t");
				char* scan_left_coord = get_token((char**)&buf, "\t");
				char* scan_right_coord = get_token((char**)&buf, "\t");
				char* orientation = get_token((char**)&buf, "\t");

				if (!scan_left_coord || !scan_right_coord || !orientation)
				{
					fprintf(stderr,"Error: malformed splice coordinate record in %s\n:%s\n",
							gtf_juncs.c_str(), buf);
					exit(1);
				}

				uint32_t ref_id = rt.get_id(ref_name, NULL, 0);
				uint32_t left_coord = atoi(scan_left_coord);
				uint32_t right_coord = atoi(scan_right_coord);
				bool antisense = *orientation == '-';

				JunctionStats junction_stat;
				junction_stat.gtf_match = true;
				junction_stat.accepted = true;

				gtf_junctions.insert(make_pair<Junction, JunctionStats>(Junction(ref_id, left_coord, right_coord, antisense), junction_stat));
			}
		}
		fprintf(stderr, "Loaded %d GFF junctions from %s.\n", (int)(gtf_junctions.size()), gtf_juncs.c_str());
	}

	vector<uint64_t> read_ids;
	vector<vector<int64_t> > offsets;
	if (num_threads > 1)
	{
		vector<string> fnames;
		if (right_map_fnames.size() > 0)
		{
			fnames.push_back(right_reads_fname);
			fnames.push_back(right_map_fnames.back());
		}
		fnames.push_back(left_reads_fname);
		fnames.push_back(left_map_fnames.back());
		bool enough_data = calculate_offsets(fnames, read_ids, offsets);
		if (!enough_data)
			num_threads = 1;
	}

	vector<JunctionSet> vjunctions(num_threads);
	vector<InsertionSet> vinsertions(num_threads);
	vector<DeletionSet> vdeletions(num_threads);
	vector<FusionSet> vfusions(num_threads);
	vector<Coverage> vcoverages(num_threads);

	vector<boost::thread*> threads;
	for (int i = 0; i < num_threads; ++i)
	{
		ConsensusEventsWorker worker;

		worker.left_map_fnames = left_map_fnames;
		worker.right_map_fnames = right_map_fnames;
		worker.rt = &rt;
		worker.gtf_junctions = &gtf_junctions;

		worker.junctions = &vjunctions[i];
		worker.insertions = &vinsertions[i];
		worker.deletions = &vdeletions[i];
		worker.fusions = &vfusions[i];
		worker.coverage = &vcoverages[i];

		worker.right_map_offset = 0;

		if (i == 0)
		{
			worker.begin_id = 0;
			worker.left_map_offset = 0;
		}
		else
		{
			size_t offsets_size = offsets[i-1].size();

			worker.begin_id = read_ids[i-1];
			worker.left_map_offset = offsets[i-1].back();
			if (offsets_size == 4)
				worker.right_map_offset = offsets[i-1][1];
		}

		worker.end_id = (i+1 < num_threads) ? read_ids[i] : std::numeric_limits<uint64_t>::max();

		if (num_threads > 1 && i + 1 < num_threads)
			threads.push_back(new boost::thread(worker));
		else
			worker();
	}

	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i]->join();
		delete threads[i];
		threads[i] = NULL;
	}
	threads.clear();

	JunctionSet& junctions = vjunctions[0];
	InsertionSet& insertions = vinsertions[0];
	DeletionSet& deletions = vdeletions[0];
	FusionSet& fusions = vfusions[0];
	Coverage& coverage = vcoverages[0];
	for (int i = 1; i < num_threads; ++i)
	{
		merge_with(junctions, vjunctions[i]);
		vjunctions[i].clear();

		merge_with(insertions, vinsertions[i]);
		vinsertions[i].clear();

		merge_with(deletions, vdeletions[i]);
		vdeletions[i].clear();

		merge_with(fusions, vfusions[i]);
		vfusions[i].clear();

		coverage.merge_with(vcoverages[i]);
		vcoverages[i].clear();
	}

	merge_with(junctions, gtf_junctions);
	coverage.calculate_coverage();

	JunctionSet rev_junctions;
	JunctionSet::const_iterator junction_iter = junctions.begin();
	for (; junction_iter != junctions.end(); ++junction_iter)
	{
		const Junction& junction = junction_iter->first;
		Junction rev_junction = junction;
		rev_junction.left = junction.right;
		rev_junction.right = junction.left;
		rev_junctions[rev_junction] = junction_iter->second;
	}

	DeletionSet rev_deletions;
#if 0
DeletionSet::const_iterator deletion_iter = deletions.begin();
for (; deletion_iter != deletions.end(); ++deletion_iter)
{
	const Deletion& deletion = deletion_iter->first;
	Deletion rev_deletion = deletion;
	rev_deletion.left = deletion.right;
	rev_deletion.right = deletion.left;
	rev_deletions[rev_deletion] = deletion_iter->second;
}
#endif

	size_t num_unfiltered_juncs = junctions.size();
	fprintf(stderr, "Loaded %lu junctions\n", (long unsigned int) num_unfiltered_juncs);

	// Read hits, extract junctions, and toss the ones that arent strongly enough supported.
	filter_junctions(junctions, gtf_junctions);

	//size_t num_juncs_after_filter = junctions.size();
	//fprintf(stderr, "Filtered %lu junctions\n",
	//     num_unfiltered_juncs - num_juncs_after_filter);

	/*
			size_t small_overhangs = 0;
			for (JunctionSet::iterator i = junctions.begin(); i != junctions.end(); ++i)
			{
			if (i->second.accepted &&
			(i->second.left_extent < min_anchor_len || i->second.right_extent < min_anchor_len))
			{
			small_overhangs++;
			}
			}

			if (small_overhangs >0)
			fprintf(stderr, "Warning: %lu small overhang junctions!\n", (long unsigned int)small_overhangs);
	 */

	JunctionSet vfinal_junctions[num_threads];
	InsertionSet vfinal_insertions[num_threads];
	DeletionSet vfinal_deletions[num_threads];
	FusionSet vfinal_fusions[num_threads];

	vector<SAlignStats> alnStats(num_threads);

	for (int i = 0; i < num_threads; ++i)
	{
		ReportWorker worker(&rt, &alnStats[i]);

		worker.sam_header_fname = sam_header;
		char filename[1024] = {0};
		sprintf(filename, "%s%d.bam", bam_output_fname.c_str(), i);
		worker.bam_output_fname = filename;
		string tmpoutdir = getFdir(worker.bam_output_fname);
		worker.left_um_fname = tmpoutdir;
		sprintf(filename, "unmapped_left_%d.bam", i);
		worker.left_um_fname+=filename;

		if (right_reads_fname != "")
		{
			sprintf(filename, "unmapped_right_%d.bam", i);
			worker.right_um_fname = tmpoutdir;
			worker.right_um_fname += filename;
		}

		worker.left_reads_fname = left_reads_fname;
		worker.left_map_fnames = left_map_fnames;
		worker.right_reads_fname = right_reads_fname;
		worker.right_map_fnames = right_map_fnames;

		worker.gtf_junctions = &gtf_junctions;
		worker.junctions = &junctions;
		worker.rev_junctions = &rev_junctions;
		worker.insertions = &insertions;
		worker.deletions = &deletions;
		worker.rev_deletions = &rev_deletions;
		worker.fusions = &fusions;
		worker.coverage = &coverage;

		worker.final_junctions = &vfinal_junctions[i];
		worker.final_insertions = &vfinal_insertions[i];
		worker.final_deletions = &vfinal_deletions[i];
		worker.final_fusions = &vfinal_fusions[i];

		//worker.rt = &rt;
		//worker.right_reads_offset = 0;
		//worker.right_map_offset = 0;

		/*if (i == 0)
		{
			worker.begin_id = 0;
			worker.left_reads_offset = 0;
			worker.left_map_offset = 0;
		}
		else */
		if (i != 0)
		{
			size_t offsets_size = offsets[i-1].size();

			worker.begin_id = read_ids[i-1];
			worker.left_reads_offset = offsets[i-1][offsets_size - 2];
			worker.left_map_offset = offsets[i-1].back();
			if (offsets_size == 4)
			{
				worker.right_reads_offset = offsets[i-1][0];
				worker.right_map_offset = offsets[i-1][1];
			}
		}

		worker.end_id = (i+1 < num_threads) ? read_ids[i] : std::numeric_limits<uint64_t>::max();

		if (num_threads > 1 && i + 1 < num_threads)
			threads.push_back(new boost::thread(worker));
		else
			worker();
	} //for each thread

	for (size_t i = 0; i < threads.size(); ++i)
	{
		threads[i]->join();
		delete threads[i];
		threads[i] = NULL;
	}
	threads.clear();

	JunctionSet& final_junctions = vfinal_junctions[0];
	InsertionSet& final_insertions = vfinal_insertions[0];
	DeletionSet& final_deletions = vfinal_deletions[0];
	FusionSet& final_fusions = vfinal_fusions[0];
	for (int i = 1; i < num_threads; ++i)
	{
		alnStats[0].add(alnStats[i]); //merge alignment stats

		merge_with(final_junctions, vfinal_junctions[i]);
		vfinal_junctions[i].clear();

		merge_with(final_insertions, vfinal_insertions[i]);
		vfinal_insertions[i].clear();

		merge_with(final_deletions, vfinal_deletions[i]);
		vfinal_deletions[i].clear();

		merge_with(final_fusions, vfinal_fusions[i]);
		vfinal_fusions[i].clear();
	}

	//small_overhangs = 0;
	for (JunctionSet::iterator i = final_junctions.begin(); i != final_junctions.end();)
	{
		if (i->second.supporting_hits == 0 || i->second.left_extent < 8 || i->second.right_extent < 8)
		{
			final_junctions.erase(i++);
		}
		else
		{
			++i;
		}
	}

	//	if (small_overhangs > 0)
	//		fprintf(stderr, "Warning: %lu small overhang junctions!\n", small_overhangs);
	print_alnStats(alnStats[0]);
	fprintf (stderr, "Printing junction BED track...");
	print_junctions(junctions_out, final_junctions, rt);
	fprintf (stderr, "done\n");

	fprintf (stderr, "Printing insertions...");
	print_insertions(insertions_out, final_insertions,rt);
	fclose(insertions_out);
	fprintf (stderr, "done\n");

	fprintf (stderr, "Printing deletions...");
	print_deletions(deletions_out, final_deletions, rt);
	fclose(deletions_out);
	fprintf (stderr, "done\n");

	if (fusion_search)
	{
		fprintf (stderr, "Printing fusions...");
		print_fusions(fusions_out, final_fusions, rt);
		fclose(fusions_out);
		fprintf (stderr, "done\n");
	}

	fprintf(stderr, "Found %lu junctions from happy spliced reads\n", (long unsigned int)final_junctions.size());
}

void print_usage()
{
	fprintf(stderr, "Usage:   tophat_reports <junctions.bed> <insertions.vcf> <deletions.vcf> <accepted_hits.sam> <left_map1,...,left_mapN> <left_reads.fq>  [right_map1,...,right_mapN] [right_reads.fq]\n");
    
	//	fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

int main(int argc, char** argv)
{
  fprintf(stderr, "tophat_reports v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION);
  fprintf(stderr, "---------------------------------------\n");
  
  //reads_format = FASTQ;
  
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string ref_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string junctions_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string insertions_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string deletions_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string fusions_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string accepted_hits_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_map_filename_list = argv[optind++];
  vector<string> left_map_filenames;
  tokenize(left_map_filename_list, ",", left_map_filenames);

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_reads_filename = argv[optind++];
  string unzcmd=getUnpackCmd(left_reads_filename, false);
  string right_map_filename_list;
  vector<string> right_map_filenames;
  string right_reads_filename;
  
  if (optind < argc)
    {
      right_map_filename_list = argv[optind++];
      tokenize(right_map_filename_list, ",", right_map_filenames);
      if(optind >= argc) {
	print_usage();
	return 1;
      }
      right_reads_filename=argv[optind++];
    }

  ifstream ref_stream(ref_file_name.c_str(), ifstream::in);
  if (!ref_stream.good())
    {
      fprintf(stderr, "Error: cannot open %s for reading\n",
	      ref_file_name.c_str());
      exit(1);
    }
  
  FILE* junctions_file = fopen(junctions_file_name.c_str(), "w");
  if (junctions_file == NULL)
    {
      fprintf(stderr, "Error: cannot open BED file %s for writing\n",
	      junctions_file_name.c_str());
      exit(1);
    }
  
  FILE* insertions_file = fopen(insertions_file_name.c_str(), "w");
  if (insertions_file == NULL)
    {
      fprintf(stderr, "Error: cannot open VCF file %s for writing\n",
	      insertions_file_name.c_str());
      exit(1);
    }
  
  FILE* deletions_file = fopen(deletions_file_name.c_str(), "w");
  if (deletions_file == NULL)
    {
      fprintf(stderr, "Error: cannot open VCF file %s for writing\n",
	      deletions_file_name.c_str());
      exit(1);
    }
  
  FILE* fusions_file = NULL;
  if (fusion_search)
    {
      fusions_file = fopen(fusions_file_name.c_str(), "w");
      if (fusions_file == NULL)
	{
	  fprintf(stderr, "Error: cannot open VCF file %s for writing\n",
		  fusions_file_name.c_str());
	  exit(1);
	}
    }
  
  driver(accepted_hits_file_name,
	 ref_stream,
	 left_map_filenames,
	 left_reads_filename,
	 right_map_filenames,
	 right_reads_filename,
	 junctions_file,
	 insertions_file,
	 deletions_file,
	 fusions_file);
  
  return 0;
}
