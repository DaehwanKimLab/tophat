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

#include <boost/thread.hpp>

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


#include "inserts.h"

using namespace std;
using namespace seqan;
using std::set;

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
      seqan::read(ref_stream, *ref_str, Fasta());
      
      rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
    }
}

void read_best_alignments(const HitsForRead& hits_for_read,
			      FragmentAlignmentGrade& best_grade,
			      HitsForRead& best_hits,
			      const JunctionSet& gtf_junctions)
{
  const vector<BowtieHit>& hits = hits_for_read.hits;
  for (size_t i = 0; i < hits.size(); ++i)
    {
      if (hits[i].edit_dist()>max_read_mismatches) continue;
      
      FragmentAlignmentGrade g(hits[i], gtf_junctions);

      // Is the new status better than the current best one?
      if (best_grade < g)
      {
	if (!report_secondary_alignments)
	  best_hits.hits.clear();
	
        best_grade = g;
        best_hits.hits.push_back(hits[i]);
      }
      else if (!(g < best_grade)) // is it just as good?
      {
        best_grade.num_alignments++;
        best_hits.hits.push_back(hits[i]);
      }
    }
}

void pair_best_alignments(const HitsForRead& left_hits,
			  const HitsForRead& right_hits,
			  InsertAlignmentGrade& best_grade,
			  HitsForRead& left_best_hits,
			  HitsForRead& right_best_hits)
{
  // max mate inner distance (genomic)
  int min_mate_inner_dist = inner_dist_mean - inner_dist_std_dev;
  if (max_mate_inner_dist == -1)
    max_mate_inner_dist = inner_dist_mean + inner_dist_std_dev;
  
  const vector<BowtieHit>& left = left_hits.hits;
  const vector<BowtieHit>& right = right_hits.hits;
  
  for (size_t i = 0; i < left.size(); ++i)
    {
      const BowtieHit& lh = left[i];
      if (lh.edit_dist() > max_read_mismatches) continue;
      for (size_t j = 0; j < right.size(); ++j)
	{
	  const BowtieHit& rh = right[j];
	  if (rh.edit_dist() > max_read_mismatches) continue;

	  bool fusion = false;
	  if (fusion_search)
	    {
	      bool left_fusion = lh.fusion_opcode() != FUSION_NOTHING;
	      bool right_fusion = rh.fusion_opcode() != FUSION_NOTHING;
	      if (left_fusion && right_fusion)
		continue;
	      
	      fusion = left_fusion || right_fusion;
	      if (!fusion && lh.ref_id() != rh.ref_id())
		fusion = true;
	      
	      if (!fusion && lh.ref_id() == rh.ref_id())
		{
		  if (lh.antisense_align() == rh.antisense_align())
		    fusion = true;
		  else
		    {
		      int inter_dist = 0;
		      if (lh.antisense_align())
			inter_dist = lh.left() - rh.right();
		      else
			inter_dist = rh.left() - lh.right();
		      
		      if (inter_dist < -(int)max_insertion_length || inter_dist > (int)fusion_min_dist)
			fusion = true;
		    }
		}
	    }

	  InsertAlignmentGrade g(lh, rh, min_mate_inner_dist, max_mate_inner_dist, fusion);

	  // Is the new status better than the current best one?
	  if (best_grade < g)
	    {
	      if (!report_secondary_alignments)
		{
		  left_best_hits.hits.clear();
		  right_best_hits.hits.clear();
		}
	      
	      best_grade = g;
	      left_best_hits.hits.push_back(lh);
	      right_best_hits.hits.push_back(rh);
	    }
	  else if (!(g < best_grade))
	    {
	      left_best_hits.hits.push_back(lh);
	      right_best_hits.hits.push_back(rh);
	      best_grade.num_alignments++;
	    }
	}
    }
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
                        const FragmentAlignmentGrade& grade,
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
  char rebuf2[2048] = {0};
  char cigar1[256] = {0}, cigar2[256] = {0};
  string left_seq, right_seq, left_qual, right_qual;
  int left = -1, left1 = -1, left2 = -1;
  bool fusion_alignment = false;
  size_t XF_index = 0;
  for (size_t i = 11; i < sam_toks.size(); ++i)
    {
      const string& tok = sam_toks[i];
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
  int mapQ=255;
  if (grade.num_alignments > 1)  {
    double err_prob = 1 - (1.0 / grade.num_alignments);
    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
  }
  int tlen =atoi(sam_toks[8].c_str()); //TLEN
  int mate_pos=atoi(sam_toks[7].c_str());

  GBamRecord* bamrec=NULL;
  if (fusion_alignment) {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name.c_str(), left1 + 1, mapQ,
				 cigar1, sam_toks[6].c_str(), mate_pos,
				 tlen, left_seq.c_str(), left_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;

    auxdata.clear();
    sam_toks[XF_index][5] = '2';
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name2.c_str(), left2 + 1, mapQ,
				 cigar2, sam_toks[6].c_str(), mate_pos,
				 tlen, right_seq.c_str(), right_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  } else {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
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
  char rebuf2[2048] = {0};
  char cigar1[256] = {0}, cigar2[256] = {0};
  string left_seq, right_seq, left_qual, right_qual;
  int left = -1, left1 = -1, left2 = -1;
  bool fusion_alignment = false;
  size_t XF_tok_idx = 11;
  for (; XF_tok_idx < sam_toks.size(); ++XF_tok_idx)
    {
      const string& tok = sam_toks[XF_tok_idx];
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
    }
  
  int gpos=isdigit(sam_toks[3][0]) ? atoi(sam_toks[3].c_str()) : 0;
  int mapQ=255;
  if (grade.num_alignments > 1) {
    double err_prob = 1 - (1.0 / grade.num_alignments);
    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
  }
  int tlen=0; //TLEN
  int mate_pos=atoi(sam_toks[7].c_str());
  if (partner) {
    if (partner->ref_id()==bh.ref_id()) {
      sam_toks[6] = "="; //same chromosome
      //TLEN:
      tlen = bh.left() < partner->left() ? partner->right() - bh.left() :
	partner->left() - bh.right();
    }
    
    else { //partner on different chromosome/contig
      sam_toks[6] = rt.get_name(partner->ref_id());
    }
    mate_pos = partner->left() + 1;
    if (grade.happy())
      flag |= 0x0002;
    if (partner->antisense_align())
      flag |=  0x0020;

    if (partner->fusion_opcode() != FUSION_NOTHING)
      {
	char partner_pos[1024];
	sprintf(partner_pos, "XP:Z:%s-%s %d", rt.get_name(partner->ref_id()), rt.get_name(partner->ref_id2()), partner->left() + 1);
	sam_toks.push_back(partner_pos);
      }
  }
  else {
    sam_toks[6] = "*";
    mate_pos = 0;
    flag |= 0x0008;
  }

  GBamRecord* bamrec=NULL;
  if (fusion_alignment) {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name.c_str(), left1 + 1, mapQ,
				 cigar1, sam_toks[6].c_str(), mate_pos,
				 tlen, left_seq.c_str(), left_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;

    auxdata.clear();
    sam_toks[XF_tok_idx][5] = '2';
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    bamrec=bam_writer.new_record(qname.c_str(), flag, ref_name2.c_str(), left2 + 1, mapQ,
				 cigar2, sam_toks[6].c_str(), mate_pos,
				 tlen, right_seq.c_str(), right_qual.c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
  } else {
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit, hitIndex);
    bamrec=bam_writer.new_record(qname.c_str(), flag, sam_toks[2].c_str(), gpos, mapQ,
				 sam_toks[5].c_str(), sam_toks[6].c_str(), mate_pos,
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
			  const FragmentAlignmentGrade& grade,
			  FragmentType frag_type,
			  Read& read,
			  GBamWriter& bam_writer,
			  FILE* um_out //write unmapped reads here
			  )
{
    assert(!read.alt_name.empty());
    lex_hit_sort s(rt, hits);
    vector<uint32_t> index_vector;

    size_t i;
    for (i = 0; i < hits.hits.size(); ++i)
        index_vector.push_back(i);
    
    sort(index_vector.begin(), index_vector.end(), s);

    size_t primaryHit = (hits.hits.empty()? 0: random() % hits.hits.size());
    bool multipleHits = (hits.hits.size() > 1);
    for (i = 0; i < hits.hits.size(); ++i)
      {
	bool primary = (i == primaryHit);
	size_t index = index_vector[i];
	const BowtieHit& bh = hits.hits[index];
	rewrite_sam_record(bam_writer, rt,
			   bh,
			   bh.hitfile_rec().c_str(),
			   read.alt_name.c_str(),
			   grade,
			   frag_type,
			   hits.hits.size(),
			   (i < hits.hits.size()-1) ? &(hits.hits[index_vector[i+1]]) : NULL,
			   primary,
			   (multipleHits? i: -1));
      }
}

void print_sam_for_pair(const RefSequenceTable& rt,
			const HitsForRead& left_hits,
                        const HitsForRead& right_hits,
                        const InsertAlignmentGrade& grade,
                        ReadStream& left_reads_file,
                        ReadStream& right_reads_file,
                        GBamWriter& bam_writer,
                        FILE* left_um_out,
                        FILE* right_um_out,
			uint64_t begin_id = 0,
                        uint64_t end_id = std::numeric_limits<uint64_t>::max())
{
    assert (left_hits.insert_id == right_hits.insert_id);
    
    Read left_read;
    Read right_read;
    assert (left_hits.hits.size() == right_hits.hits.size() ||
            (left_hits.hits.empty() || right_hits.hits.empty()));
    
    size_t primaryHit = 0;
    vector<uint32_t> index_vector;
    if(right_hits.hits.size() > 0)
    {
      lex_hit_sort s(rt, right_hits);
      for (size_t i = 0; i < right_hits.hits.size(); ++i)
	index_vector.push_back(i);

          sort(index_vector.begin(), index_vector.end(), s);
          primaryHit = random() % right_hits.hits.size();
      }
    else if (left_hits.hits.size() > 0)
    {
          lex_hit_sort s(rt, left_hits);
          for (size_t i = 0; i < left_hits.hits.size(); ++i)
              index_vector.push_back(i);

          sort(index_vector.begin(), index_vector.end(), s);
          primaryHit = random() % left_hits.hits.size();
      }
    
    bool got_left_read = left_reads_file.getRead(left_hits.insert_id, left_read,
						 reads_format, false, left_um_out,
						 false, begin_id, end_id);
    
    bool got_right_read = right_reads_file.getRead(right_hits.insert_id, right_read,
						   reads_format, false, right_um_out,
						   false, begin_id, end_id);
    
    if (left_hits.hits.size() == right_hits.hits.size())
    {
        assert (got_left_read && got_right_read);
        bool multipleHits = (left_hits.hits.size() > 1);
        for (size_t i = 0; i < right_hits.hits.size(); ++i)
        {
            bool primary = (i == primaryHit);
            size_t index = index_vector[i];
            const BowtieHit& right_bh = right_hits.hits[index];
            const BowtieHit& left_bh = left_hits.hits[index];
            
            rewrite_sam_record(bam_writer, rt,
                               right_bh,
                               right_bh.hitfile_rec().c_str(),
                               right_read.alt_name.c_str(),
                               grade,
                               FRAG_RIGHT,
                               &left_bh,
                               right_hits.hits.size(),
                               (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL,
                               primary, 
                               (multipleHits? i: -1));
            rewrite_sam_record(bam_writer, rt,
                               left_bh,
                               left_bh.hitfile_rec().c_str(),
                               left_read.alt_name.c_str(),
                               grade,
                               FRAG_LEFT,
                               &right_bh,
                               left_hits.hits.size(),
                               (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL,
                               primary,
                               (multipleHits? i: -1));
        }
    }
    else if (left_hits.hits.empty())
    { //only right read was mapped properly
        if (right_um_out) {
          //write it in the mapped file with the #MAPPED# flag
          fprintf(right_um_out, "@%s #MAPPED#\n%s\n+\n%s\n", right_read.alt_name.c_str(),
                              right_read.seq.c_str(), right_read.qual.c_str());
          }
        for (size_t i = 0; i < right_hits.hits.size(); ++i)
        {
            bool primary = (i == primaryHit);
            size_t index = index_vector[i];
            const BowtieHit& bh = right_hits.hits[index];
            
            rewrite_sam_record(bam_writer, rt,
                               bh,
                               bh.hitfile_rec().c_str(),
                               right_read.alt_name.c_str(),
                               grade,
                               FRAG_RIGHT,
                               NULL,
                               right_hits.hits.size(),
                               (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL,
                               primary,
                               -1);
            
        }
    }
    else if (right_hits.hits.empty())
    { //only left read was mapped properly
      if (left_um_out) {
        //write it in the mapped file with the #MAPPED# flag
        fprintf(left_um_out, "@%s #MAPPED#\n%s\n+\n%s\n", left_read.alt_name.c_str(), left_read.seq.c_str(),
                            left_read.qual.c_str());
        }

        for (size_t i = 0; i < left_hits.hits.size(); ++i)
        {
            bool primary = (i == primaryHit);
            size_t index = index_vector[i];
            const BowtieHit& bh = left_hits.hits[index];
            rewrite_sam_record(bam_writer, rt,
                               bh,
                               bh.hitfile_rec().c_str(),
                               left_read.alt_name.c_str(),
                               grade,
                               FRAG_LEFT,
                               NULL,
                               left_hits.hits.size(),
                               (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL,
                               primary,			  
                               -1);
            
        }
    }
    else
    {
        assert (false);
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


FusionSet empty_fusions;
void update_fusions(const HitsForRead& hits,
		    RefSequenceTable& rt,
		    FusionSet& fusions,
		    FusionSet& fusions_ref = empty_fusions)
{
  if (hits.hits.size() > fusion_multireads)
    return;

  bool update_stat = fusions_ref.size() > 0;
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      const BowtieHit& bh = hits.hits[i];

      if (bh.edit_dist() > fusion_read_mismatches)
	continue;

      // daehwan
#if 0
      vector<Fusion> new_fusions;
      fusions_from_spliced_hit(bh, new_fusions);
      
      for(size_t i = 0; i < new_fusions.size(); ++i)
	{
	  Fusion fusion = new_fusions[i];
	  if (fusion.left == 110343414 && fusion.right == 80211173 && hits.hits.size() > 1)
	    {
	      for (size_t t = 0; t < hits.hits.size(); ++t)
		{
		  const BowtieHit& lh = hits.hits[t];
		  cout << "insert id: " << lh.insert_id() << endl;
		  cout << (lh.antisense_align() ? "-" : "+") <<  endl;
		  cout << lh.ref_id() << ": " << lh.left() << "-" << lh.right() << endl;
		  cout << lh.ref_id2() << ": " << print_cigar(lh.cigar()) << endl;
		}
	    }
	}
#endif

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

// events include splice junction, indels, and fusion points.
struct ConsensusEventsWorker
{
  void operator()()
  {
    ReadTable it;
    vector<BAMHitFactory*> hit_factories;
    hit_factories.push_back(new BAMHitFactory(it, *rt));
    HitStream l_hs(left_map_fname, hit_factories.back(), false, true, true, true);
    if (left_map_offset > 0)
      l_hs.seek(left_map_offset);

    hit_factories.push_back(new BAMHitFactory(it, *rt));
    HitStream r_hs(right_map_fname, hit_factories.back(), false, true, true, true);
    if (right_map_offset > 0)
      r_hs.seek(right_map_offset);

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
	    FragmentAlignmentGrade grade;
            
	    // Process hits for left singleton, select best alignments
	    read_best_alignments(curr_left_hit_group, grade, best_hits, *gtf_junctions);
	    update_junctions(best_hits, *junctions);
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
		FragmentAlignmentGrade grade;
		// Process hit for right singleton, select best alignments
		read_best_alignments(curr_right_hit_group,grade, best_hits, *gtf_junctions);
		update_junctions(best_hits, *junctions);
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
	    if (curr_left_hit_group.hits.empty())
	      {
		HitsForRead right_best_hits;
		right_best_hits.insert_id = curr_right_obs_order;
                
		FragmentAlignmentGrade grade;
		read_best_alignments(curr_right_hit_group, grade, right_best_hits, *gtf_junctions);
		update_junctions(right_best_hits, *junctions);
		update_fusions(right_best_hits, *rt, *fusions);
	      }
	    else if (curr_right_hit_group.hits.empty())
	      {
		HitsForRead left_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
                
		FragmentAlignmentGrade grade;
		// Process hits for left singleton, select best alignments
		read_best_alignments(curr_left_hit_group, grade, left_best_hits, *gtf_junctions);
		update_junctions(left_best_hits, *junctions);
		update_fusions(left_best_hits, *rt, *fusions);
	      }
	    else
	      {
		HitsForRead left_best_hits;
		HitsForRead right_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
		right_best_hits.insert_id = curr_right_obs_order;
		
		// daehwan - apply gtf_junctions here, too!

		InsertAlignmentGrade grade;
		pair_best_alignments(curr_left_hit_group,
				     curr_right_hit_group,
				     grade,
				     left_best_hits,
				     right_best_hits);
		update_junctions(left_best_hits, *junctions);
		update_junctions(right_best_hits, *junctions);
		update_fusions(left_best_hits, *rt, *fusions);
		update_fusions(right_best_hits, *rt, *fusions);
	      }
            
	    l_hs.next_read_hits(curr_left_hit_group);
	    curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
            
	    r_hs.next_read_hits(curr_right_hit_group);
	    curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	  }
      }

    for (size_t i = 0; i < hit_factories.size(); ++i)
      delete hit_factories[i];
    
    hit_factories.clear();
  }

  string left_map_fname;
  string right_map_fname;

  JunctionSet* junctions;
  JunctionSet* gtf_junctions;
  FusionSet* fusions;
  RefSequenceTable* rt;

  uint64_t begin_id;
  uint64_t end_id;
  int64_t left_map_offset;
  int64_t right_map_offset;
};

struct ReportWorker
{
  void operator()()
  {
    ReadTable it;
    GBamWriter bam_writer(bam_output_fname.c_str(), sam_header.c_str());

    ReadStream left_reads_file(left_reads_fname);
    if (left_reads_file.file() == NULL)
      err_die("Error: cannot open %s for reading\n", left_reads_fname.c_str());

    if (left_reads_offset > 0)
      left_reads_file.seek(left_reads_offset);
    
    if (!zpacker.empty()) left_um_fname+=".z";
    FZPipe left_um_out;
    if (left_um_out.openWrite(left_um_fname.c_str(), zpacker)==NULL)
      err_die("Error: cannot open file %s for writing!\n",left_um_fname.c_str());
    
    ReadStream right_reads_file(right_reads_fname);
    if (right_reads_offset > 0)
      right_reads_file.seek(right_reads_offset);
    
    FZPipe right_um_out;
    if (right_reads_fname != "")
      {
	if (!zpacker.empty()) right_um_fname+=".z";
	if (right_um_out.openWrite(right_um_fname.c_str(), zpacker)==NULL)
	  err_die("Error: cannot open file %s for writing!\n",right_um_fname.c_str());
      }

    vector<BAMHitFactory*> hit_factories;
    hit_factories.push_back(new BAMHitFactory(it, *rt));
    HitStream left_hs(left_map_fname, hit_factories.back(), false, true, true, true);
    if (left_map_offset > 0)
      left_hs.seek(left_map_offset);

    hit_factories.push_back(new BAMHitFactory(it, *rt));
    HitStream right_hs(right_map_fname, hit_factories.back(), false, true, true, true);
    if (right_map_offset > 0)
      right_hs.seek(right_map_offset);
    
    HitsForRead curr_left_hit_group;
    HitsForRead curr_right_hit_group;
    
    left_hs.next_read_hits(curr_left_hit_group);
    right_hs.next_read_hits(curr_right_hit_group);
    
    uint64_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
    uint64_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);

    // While we still have unreported hits...
    Read l_read;
    Read r_read;
    while((curr_left_obs_order != VMAXINT32 || curr_right_obs_order != VMAXINT32) &&
	  (curr_left_obs_order < end_id || curr_right_obs_order < end_id))
      {
	// Chew up left singletons (pairs with right reads unmapped)
	while (curr_left_obs_order < curr_right_obs_order &&
	       curr_left_obs_order < end_id &&
	       curr_left_obs_order != VMAXINT32)
	  {
	    HitsForRead best_hits;
	    best_hits.insert_id = curr_left_obs_order;
	    FragmentAlignmentGrade grade;
	    bool got_read=left_reads_file.getRead(curr_left_obs_order, l_read,
						  reads_format, false, left_um_out.file,
						  false, begin_id, end_id);
	    assert(got_read);

	    if (right_reads_file.file()) {
	      fprintf(left_um_out.file, "@%s #MAPPED#\n%s\n+\n%s\n", l_read.alt_name.c_str(),
		      l_read.seq.c_str(), l_read.qual.c_str());
	      got_read=right_reads_file.getRead(curr_left_obs_order, r_read,
						reads_format, false, right_um_out.file,
						true, begin_id, end_id);
	      assert(got_read);
	    }
    
	    exclude_hits_on_filtered_junctions(*junctions, curr_left_hit_group);

	    // Process hits for left singleton, select best alignments
	    read_best_alignments(curr_left_hit_group, grade, best_hits, *gtf_junctions);

	    if (best_hits.hits.size()>0 && best_hits.hits.size() <= max_multihits)
	      {
		update_junctions(best_hits, *final_junctions);
		update_insertions_and_deletions(best_hits, *final_insertions, *final_deletions);
		update_fusions(best_hits, *rt, *final_fusions, *fusions);

		print_sam_for_single(*rt,
				     best_hits,
				     grade,
				     (right_map_fname.empty() ? FRAG_UNPAIRED : FRAG_LEFT),
				     l_read,
				     bam_writer,
				     left_um_out.file);
	      }
	    // Get next hit group
	    left_hs.next_read_hits(curr_left_hit_group);
	    curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	  } //left singletons 
	
	// Chew up right singletons
	while (curr_left_obs_order > curr_right_obs_order &&
	       curr_right_obs_order < end_id &&
	       curr_right_obs_order != VMAXINT32)
	  {
	    HitsForRead best_hits;
	    best_hits.insert_id = curr_right_obs_order;
	    FragmentAlignmentGrade grade;

	    if (curr_right_obs_order >=  begin_id)
	      {
		bool got_read=right_reads_file.getRead(curr_right_obs_order, r_read,
						       reads_format, false, right_um_out.file,
						       false, begin_id, end_id);

		assert(got_read);

		fprintf(right_um_out.file, "@%s #MAPPED#\n%s\n+\n%s\n", r_read.alt_name.c_str(),
			r_read.seq.c_str(), r_read.qual.c_str());
		got_read=left_reads_file.getRead(curr_right_obs_order, l_read,
						 reads_format, false, left_um_out.file,
						 true, begin_id, end_id);
		assert(got_read);

		exclude_hits_on_filtered_junctions(*junctions, curr_right_hit_group);

		// Process hit for right singleton, select best alignments
		read_best_alignments(curr_right_hit_group, grade, best_hits, *gtf_junctions);
		if (best_hits.hits.size()>0 && best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(best_hits, *final_junctions);
		    update_insertions_and_deletions(best_hits, *final_insertions, *final_deletions);
		    update_fusions(best_hits, *rt, *final_fusions, *fusions);
		    
		    print_sam_for_single(*rt,
					 best_hits,
					 grade,
					 FRAG_RIGHT,
					 r_read,
					 bam_writer,
					 right_um_out.file);
		  }
	      }

	    // Get next hit group
	    right_hs.next_read_hits(curr_right_hit_group);
	    curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	  }
	
	// Since we have both left hits and right hits for this insert,
	// Find the best pairing and print both
	while (curr_left_obs_order == curr_right_obs_order &&
	       curr_left_obs_order < end_id &&
	       curr_left_obs_order != VMAXINT32)
	  {
	    exclude_hits_on_filtered_junctions(*junctions, curr_left_hit_group);
	    exclude_hits_on_filtered_junctions(*junctions, curr_right_hit_group);
	    if (curr_left_hit_group.hits.empty())
	      {   //only right read mapped
                //write it in the mapped file with the #MAPPED# flag
		
		bool got_read=right_reads_file.getRead(curr_left_obs_order, r_read,
						       reads_format, false, right_um_out.file,
						       false, begin_id, end_id);
		assert(got_read);
		fprintf(right_um_out.file, "@%s #MAPPED#\n%s\n+\n%s\n", r_read.alt_name.c_str(),
			r_read.seq.c_str(), r_read.qual.c_str());
		HitsForRead right_best_hits;
		right_best_hits.insert_id = curr_right_obs_order;
		
		FragmentAlignmentGrade grade;
		read_best_alignments(curr_right_hit_group, grade, right_best_hits, *gtf_junctions);
		
		if (right_best_hits.hits.size()>0 && right_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(right_best_hits, *final_junctions);
		    update_insertions_and_deletions(right_best_hits, *final_insertions, *final_deletions);
		    update_fusions(right_best_hits, *rt, *final_fusions, *fusions);
		    
		    print_sam_for_single(*rt,
					 right_best_hits,
					 grade,
					 FRAG_RIGHT,
					 r_read,
					 bam_writer,
					 right_um_out.file);
		  }
	      }
	    else if (curr_right_hit_group.hits.empty())
	      {
		HitsForRead left_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
		//only left read mapped
		bool got_read=left_reads_file.getRead(curr_left_obs_order, l_read,
						      reads_format, false, left_um_out.file,
						      false, begin_id, end_id);
		assert(got_read);
		fprintf(left_um_out.file, "@%s #MAPPED#\n%s\n+\n%s\n", l_read.alt_name.c_str(),
			l_read.seq.c_str(), l_read.qual.c_str());
		FragmentAlignmentGrade grade;
		// Process hits for left singleton, select best alignments
		read_best_alignments(curr_left_hit_group, grade, left_best_hits, *gtf_junctions);
		if (left_best_hits.hits.size()>0 && left_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(left_best_hits, *final_junctions);
		    update_insertions_and_deletions(left_best_hits, *final_insertions, *final_deletions);
		    update_fusions(left_best_hits, *rt, *final_fusions, *fusions);
		    
		    print_sam_for_single(*rt,
					 left_best_hits,
					 grade,
					 FRAG_LEFT,
					 l_read,
					 bam_writer,
					 left_um_out.file);
		  }
	      }
	    else
	      {   //hits for both left and right reads
		HitsForRead left_best_hits;
		HitsForRead right_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
		right_best_hits.insert_id = curr_right_obs_order;

		InsertAlignmentGrade grade;
		pair_best_alignments(curr_left_hit_group,
				     curr_right_hit_group,
				     grade,
				     left_best_hits,
				     right_best_hits);
		
		if (left_best_hits.hits.size()>0 && left_best_hits.hits.size() <= max_multihits &&
		    right_best_hits.hits.size()>0 && right_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(left_best_hits, *final_junctions);
		    update_junctions(right_best_hits, *final_junctions);
		    update_insertions_and_deletions(left_best_hits, *final_insertions, *final_deletions);
		    update_insertions_and_deletions(right_best_hits, *final_insertions, *final_deletions);

		    pair_support(left_best_hits, right_best_hits, *final_fusions, *fusions);
		    update_fusions(left_best_hits, *rt, *final_fusions, *fusions);
		    update_fusions(right_best_hits, *rt, *final_fusions, *fusions);

		    print_sam_for_pair(*rt,
				       left_best_hits,
				       right_best_hits,
				       grade,
				       left_reads_file,
				       right_reads_file,
				       bam_writer,
				       left_um_out.file,
				       right_um_out.file,
				       begin_id,
				       end_id);
		  }
	      }
	    
	    left_hs.next_read_hits(curr_left_hit_group);
	    curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	    
	    right_hs.next_read_hits(curr_right_hit_group);
	    curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	  }
      } //while we still have unreported hits..
    //print the remaining unmapped reads at the end of each reads' stream
    left_reads_file.getRead(VMAXINT32, l_read,
			    reads_format, false, left_um_out.file,
			    false, begin_id, end_id);
    if (right_reads_file.file())
      right_reads_file.getRead(VMAXINT32, r_read,
			       reads_format, false, right_um_out.file,
			       false, begin_id, end_id);

    left_um_out.close();
    right_um_out.close();

    for (size_t i = 0; i < hit_factories.size(); ++i)
      delete hit_factories[i];

    hit_factories.clear();
  }

  string bam_output_fname;
  string sam_header_fname;

  string left_reads_fname;
  string left_map_fname;
  string right_reads_fname;
  string right_map_fname;

  string left_um_fname;
  string right_um_fname;

  JunctionSet* gtf_junctions;
  JunctionSet* junctions;
  FusionSet* fusions;
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
};

void driver(const string& bam_output_fname,
	    istream& ref_stream,
	    const string& left_map_fname,
	    const string& left_reads_fname,
	    const string& right_map_fname,
	    const string& right_reads_fname,
	    FILE* junctions_out,
	    FILE* insertions_out,
	    FILE* deletions_out,
	    FILE* fusions_out)
{
  if (!parallel)
    num_threads = 1;

  RefSequenceTable rt(sam_header, true);

  if (fusion_search)
    get_seqs(ref_stream, rt, true);

  srandom(1);
  
  JunctionSet gtf_junctions;
  if (!gtf_juncs.empty())
    {
      char splice_buf[2048];
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

              // add 1 to left_coord to meet BED format
              gtf_junctions.insert(make_pair<Junction, JunctionStats>(Junction(ref_id, left_coord + 1, right_coord, antisense), JunctionStats()));
            }
        }
      fprintf(stderr, "Loaded %d GFF junctions from %s.\n", (int)(gtf_junctions.size()), gtf_juncs.c_str());
    }

  vector<uint64_t> read_ids;
  vector<vector<int64_t> > offsets;
  if (num_threads > 1)
    {
      vector<string> fnames;
      if (right_map_fname != "")
	{
	  fnames.push_back(right_reads_fname);
	  fnames.push_back(right_map_fname);
	}
      fnames.push_back(left_reads_fname);
      fnames.push_back(left_map_fname);
      bool enough_data = calculate_offsets(fnames, read_ids, offsets);
      if (!enough_data)
	num_threads = 1;
    }

  JunctionSet vjunctions[num_threads];
  FusionSet vfusions[num_threads];
  vector<boost::thread*> threads;
  for (int i = 0; i < num_threads; ++i)
    {
      ConsensusEventsWorker worker;

      worker.left_map_fname = left_map_fname;
      worker.right_map_fname = right_map_fname;
      worker.junctions = &vjunctions[i];
      worker.fusions = &vfusions[i];
      worker.gtf_junctions = &gtf_junctions;
      worker.rt = &rt;

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

      if (num_threads > 1)
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

  JunctionSet junctions = vjunctions[0];
  FusionSet fusions = vfusions[0];
  for (size_t i = 1; i < num_threads; ++i)
    {
      merge_with(junctions, vjunctions[i]);
      vjunctions[i].clear();

      merge_with(fusions, vfusions[i]);
      vfusions[i].clear();
    }

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

  for (int i = 0; i < num_threads; ++i)
    {
      ReportWorker worker;

      char filename[1024] = {0};
      sprintf(filename, "%s%d.bam", bam_output_fname.c_str(), i);
      worker.bam_output_fname = filename;
      worker.sam_header_fname = sam_header;
      sprintf(filename, "%s/unmapped_left%d.fq", output_dir.c_str(), i);
      worker.left_um_fname = filename;
      
      if (right_reads_fname != "")
	{
	  sprintf(filename, "%s/unmapped_right%d.fq", output_dir.c_str(), i);
	  worker.right_um_fname = filename;
	}
      
      worker.left_reads_fname = left_reads_fname;
      worker.left_map_fname = left_map_fname;
      worker.right_reads_fname = right_reads_fname;
      worker.right_map_fname = right_map_fname;

      worker.gtf_junctions = &gtf_junctions;
      worker.junctions = &junctions;
      worker.fusions = &fusions;
      worker.final_junctions = &vfinal_junctions[i];
      worker.final_insertions = &vfinal_insertions[i];
      worker.final_deletions = &vfinal_deletions[i];
      worker.final_fusions = &vfinal_fusions[i];
      worker.rt = &rt;

      worker.right_reads_offset = 0;
      worker.right_map_offset = 0;
      
      if (i == 0)
	{
	  worker.begin_id = 0;
	  worker.left_reads_offset = 0;
	  worker.left_map_offset = 0;
	}
      else
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

      if (num_threads > 1)
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

  JunctionSet final_junctions = vfinal_junctions[0];
  InsertionSet final_insertions = vfinal_insertions[0];
  DeletionSet final_deletions = vfinal_deletions[0];
  FusionSet final_fusions = vfinal_fusions[0];
  for (size_t i = 1; i < num_threads; ++i)
    {
      merge_with(final_junctions, vfinal_junctions[i]);
      vfinal_junctions[i].clear();

      merge_with(final_insertions, vfinal_insertions[i]);
      vfinal_insertions[i].clear();
      
      merge_with(final_deletions, vfinal_deletions[i]);
      vfinal_deletions[i].clear();

      merge_with(final_fusions, vfinal_fusions[i]);
      vfinal_fusions[i].clear();
    }

  fprintf (stderr, "Reporting final accepted alignments...");
  fprintf (stderr, "done.\n");
  
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
  
  reads_format = FASTQ;
  
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
  
  string left_map_filename = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_reads_filename = argv[optind++];
  string unzcmd=getUnpackCmd(left_reads_filename, false);
  
  string right_map_filename;
  string right_reads_filename;
  
  if (optind < argc)
    {
      right_map_filename = argv[optind++];
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
	 left_map_filename,
	 left_reads_filename,
	 right_map_filename,
	 right_reads_filename,
	 junctions_file,
	 insertions_file,
	 deletions_file,
	 fusions_file);
  
  return 0;
}
