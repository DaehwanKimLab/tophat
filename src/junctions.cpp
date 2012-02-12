/*
 *  junctions.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include "common.h"
#include "junctions.h"
#include "bwt_map.h"

void junctions_from_spliced_hit(const BowtieHit& h, 
				vector<pair<Junction, JunctionStats> >& new_juncs)
{
  const vector<CigarOp>& cigar = h.cigar();
  int j = h.left();

  bool bSawFusion = false;
  for (size_t c = 0 ; c < cigar.size(); ++c)
    {
      Junction junc;
      JunctionStats stats;

      int opcode = cigar[c].opcode;
      int length = cigar[c].length;
      
      switch(opcode)
	{
	case REF_SKIP:
	case rEF_SKIP:
	  if (bSawFusion)
	    junc.refid = h.ref_id2();
	  else
	    junc.refid = h.ref_id();

	  // daehwan - we need to consider indels very next to REF_SKIP,
	  // which is possible due to Bowtie2
	  assert (c > 0 && c < cigar.size() - 1); 
	  assert (cigar[c - 1].length);
	  assert (cigar[c + 1].length);

	  if (opcode == REF_SKIP)
	    {
	      junc.left = j - 1;
	      junc.right = j + length;
	      stats.left_extent = cigar[c - 1].length;
	      stats.right_extent = cigar[c + 1].length;
	      j += length;
	    }
	  else
	    {
	      junc.right = j + 1;
	      junc.left = j - length;
	      stats.right_extent = cigar[c - 1].length;
	      stats.left_extent = cigar[c + 1].length;
	      j -= length;
	    }
	  
	  junc.antisense = h.antisense_splice();

	  /*
	   * Note that in valid_hit() in tophat_report.cpp
	   * we have tried to ensure that the REF_SKIP operator
	   * will only be surrounded by match operators
	   */	

	  stats.min_splice_mms = h.splice_mms();
	  stats.supporting_hits++;
	  new_juncs.push_back(make_pair(junc, stats));
	  break;
	case MATCH:
	case DEL:
	  j += cigar[c].length;
	  break;
	case mATCH:
	case dEL:
	  j -= cigar[c].length;
	  break;
	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	  j = cigar[c].length;
	  bSawFusion = true;
	  break;
	default:
	  break;
	}
    }
}

void print_junction(FILE* junctions_out, 
		    const char* name, 
		    const Junction& j, 
		    const JunctionStats& s, 
		    uint64_t junc_id)
{
  int left_plus_one = j.left + 1;
  fprintf(junctions_out,
	  "%s\t%d\t%d\tJUNC%08d\t%d\t%c\t%d\t%d\t255,0,0\t2\t%d,%d\t0,%d\n",
	  name,
	  left_plus_one - s.left_extent,
	  j.right + s.right_extent,
	  (int)junc_id,
	  s.supporting_hits,
	  j.antisense ? '-' : '+',
	  left_plus_one - s.left_extent,
	  j.right + s.right_extent,
	  s.left_extent,
	  s.right_extent,
	  j.right - (left_plus_one - s.left_extent));
}

void junctions_from_alignment(const BowtieHit& spliced_alignment,
			      JunctionSet& junctions)
{
  vector<pair<Junction, JunctionStats> > juncs;
  junctions_from_spliced_hit(spliced_alignment, juncs);

  for (size_t i = 0; i < juncs.size(); ++i)
    {
      pair<Junction, JunctionStats>& junc = juncs[i];
      JunctionSet::iterator itr = junctions.find(junc.first);
      
      if (itr != junctions.end())
	{
	  JunctionStats& j = itr->second;
	  j.merge_with(junc.second);
	}
      else
	{
	  assert(junc.first.refid != VMAXINT32);
	  junctions[junc.first] = junc.second;
	}
    }
}

#if !NDEBUG
void validate_junctions(const JunctionSet& junctions)
{
  uint32_t invalid_juncs = 0;
  for (JunctionSet::const_iterator i = junctions.begin();
       i != junctions.end();
       ++i)
    {
      if (!i->first.valid())
	invalid_juncs++;
    }
  fprintf(stderr, "Found %d invalid junctions\n", invalid_juncs);
}
#endif

int rejected = 0;
int rejected_spliced = 0;
int total_spliced = 0;
int total = 0;

/*
void junctions_from_alignments(HitTable& hits,
JunctionSet& junctions)
{
	for (HitTable::iterator ci = hits.begin();
		 ci != hits.end();
		 ++ci)
	{
		HitList& rh = ci->second;
		if (rh.size() == 0)
			continue;
		for (size_t i = 0; i < rh.size(); ++i)
		{
			BowtieHit& bh = rh[i];
			AlignStatus s = status(&bh);
			total++;
			if (s == SPLICED)
				total_spliced++;
			
			if (s == SPLICED)
			{
				junctions_from_alignment(bh, junctions); 
			}
		}
	}
}
*/

bool accept_if_valid(const Junction& j, JunctionStats& s)
{
  if (min(s.left_extent, s.right_extent) < min_anchor_len)
    {
      s.accepted = false;
      return false;
    }
  
  if (s.min_splice_mms > max_splice_mismatches)
    {
      s.accepted = false;
      return false;
    }
  
  //	uint32_t junc_doc = 0;
  //	uint8_t extent = 0;
  //	if (s.left_exon_doc > s.right_exon_doc)
  //	{
  //		junc_doc = s.left_exon_doc;
  //		extent = s.left_extent;
  //	}
  //	else
  //	{
  //		junc_doc = s.right_exon_doc;
  //		extent = s.right_extent;
  //	}
  //	
  //	double avg_junc_doc = junc_doc / (double)(extent);
  
  //if (avg_junc_doc / (float) s.num_reads > 100.0)
  //	if (s.supporting_hits / avg_junc_doc < min_isoform_fraction)
  //	{
  //		s.accepted = false;
  //	}
  //	else
  {
    //fprintf (stderr, "Junction size = %d\n, support = %d", (int)j.right - (int)j.left, (int)s.supporting_hits.size()  );
    if ((int)j.right - (int)j.left > 50000)
      {
	s.accepted = (s.supporting_hits >= 2 && 
		      min(s.left_extent, s.right_extent) > 12);
      }
    else
      {
	s.accepted = true;
      }
  }
  return s.accepted;
}

void knockout_shadow_junctions(JunctionSet& junctions)
{
  vector<uint32_t> ref_ids;
  
  for (JunctionSet::iterator i = junctions.begin(); i != junctions.end(); ++i)
    {
      ref_ids.push_back(i->first.refid);
    }
  
  sort(ref_ids.begin(), ref_ids.end());
  vector<uint32_t>::iterator new_end = unique(ref_ids.begin(), ref_ids.end());
  ref_ids.erase(new_end, ref_ids.end());
  
  for(size_t i = 0; i < ref_ids.size(); ++i)
    {
      uint32_t refid = ref_ids[i];
		
      Junction dummy_left(refid, 0, 0, true);
      Junction dummy_right(refid, VMAXINT32, VMAXINT32, true);
      
      pair<JunctionSet::iterator, JunctionSet::iterator> r;
      r.first = junctions.lower_bound(dummy_left);
      r.second = junctions.upper_bound(dummy_right);
      
      JunctionSet::iterator itr = r.first;
      
      while(itr != r.second && itr != junctions.end())
      {
        if (itr->second.accepted && !itr->second.gtf_match)
          {
            Junction fuzzy_left = itr->first;
            Junction fuzzy_right = itr->first;
            fuzzy_left.left -= min_anchor_len;
            fuzzy_right.right += min_anchor_len;
            fuzzy_left.antisense = !itr->first.antisense;
            fuzzy_right.antisense = !itr->first.antisense;

            pair<JunctionSet::iterator, JunctionSet::iterator> s;
            s.first = junctions.lower_bound(fuzzy_left);
            s.second = junctions.upper_bound(fuzzy_right);
            JunctionSet::iterator itr2 = s.first;

            int junc_support = itr->second.supporting_hits;

            while(itr2 != s.second && itr2 != junctions.end())
              {
                int left_diff = itr->first.left - itr2->first.left;
                int right_diff = itr->first.right - itr2->first.right;
                if (itr != itr2 &&
                    itr->first.antisense != itr2->first.antisense &&
                    (left_diff < min_anchor_len || right_diff < min_anchor_len))
                  {
                    if (junc_support < itr2->second.supporting_hits)
                      itr->second.accepted = false;
                  }
                ++itr2;
              }
          }
        ++itr;
      }
    }
}

void filter_junctions(JunctionSet& junctions, const JunctionSet& gtf_junctions)
{
  for (JunctionSet::iterator i = junctions.begin(); i != junctions.end(); ++i)
    {
      if (gtf_junctions.find(i->first) == gtf_junctions.end())
        accept_if_valid(i->first, i->second);
      else {//automatically accept junctions matching GTF
        i->second.accepted = true;
        i->second.gtf_match = true;
        }
    }

  knockout_shadow_junctions(junctions);
}

void accept_all_junctions(JunctionSet& junctions,
			  const uint32_t refid)
{
  fprintf(stderr, "Accepting all junctions\n");
  for (JunctionSet::iterator itr = junctions.begin(); itr != junctions.end(); ++itr)
    {
      itr->second.accepted = true;
    }
}


void print_junctions(FILE* junctions_out, 
		     const JunctionSet& junctions,
		     RefSequenceTable& ref_sequences)
{
  uint64_t junc_id = 1;
  fprintf(junctions_out, "track name=junctions description=\"TopHat junctions\"\n");
  for (JunctionSet::const_iterator i = junctions.begin();
       i != junctions.end();
       ++i)
    {
      const pair<Junction, JunctionStats>& j_itr = *i; 
      const Junction& j = j_itr.first;
      const JunctionStats& s = j_itr.second;			
      
      assert(ref_sequences.get_name(j.refid));
      //fprintf(stdout,"%d\t%d\t%d\t%c\n", j.refid, j.left, j.right, j.antisense ? '-' : '+');
      print_junction(junctions_out, 
		     ref_sequences.get_name(j.refid),
		     j,
		     s, 
		     junc_id++);
    }
  //fprintf(stderr, "Rejected %d / %d alignments, %d / %d spliced\n", rejected, total, rejected_spliced, total_spliced);
}

// Extracts junctions from all the SAM hits (based on REF_SKIPs) in the hit file
// resets the stream when finished.
void get_junctions_from_hits(HitStream& hit_stream, 
			     ReadTable& it, 
			     JunctionSet& junctions)
{
  HitsForRead curr_hit_group;
  hit_stream.next_read_hits(curr_hit_group);
  
  uint32_t curr_obs_order = it.observation_order(curr_hit_group.insert_id);
  
  while(curr_obs_order != VMAXINT32)
    {
      for (size_t i = 0; i < curr_hit_group.hits.size(); ++i)
	{
	  BowtieHit& bh = curr_hit_group.hits[i];
	  if (!bh.contiguous())
	    {
	      junctions_from_alignment(bh, junctions);
	    }
	  hit_stream.next_read_hits(curr_hit_group);
	  curr_obs_order = it.observation_order(curr_hit_group.insert_id);
	}
    }
  
  hit_stream.reset();
}

void merge_with(JunctionSet& juncs, const JunctionSet& other_juncs)
{
  for (JunctionSet::const_iterator junc = other_juncs.begin(); junc != other_juncs.end(); ++junc)
    {
      JunctionSet::iterator itr = juncs.find(junc->first);
      if (itr != juncs.end())
	{
	  JunctionStats& curr  = itr->second;
	  curr.merge_with(junc->second);
	}
      else
	{
	  juncs[junc->first] = junc->second;
	}
    }
}
