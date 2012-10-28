/*
 *  align_status.cpp
 *  TopHat
 *
 *  Created by Ryan Kelley on 11/09/2010.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "reads.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "fusions.h"
#include "coverage.h"
#include "align_status.h"

using namespace std;

AlignStatus::AlignStatus()
{
  _alignment_score = std::numeric_limits<int>::min();
} 

/**
 * Parse the cigar string of a BowtieHit in order to determine the alignment status.
 */
AlignStatus::AlignStatus(const BowtieHit& bh,
			 const JunctionSet& gtf_junctions,
			 const JunctionSet& junctions,
			 const InsertionSet& insertions,
			 const DeletionSet& deletions,
			 const FusionSet& fusions,
			 const Coverage& coverage)
{
  // it seems like we need to work on this more
  // daehwan - it doesn't seem to work.
  const bool recalculate_indel_score = false;
  
  const vector<CigarOp>& cigar = bh.cigar();
  _alignment_score = bh.alignment_score();

  const int read_len = bh.read_len();
  const int min_extent = min(read_len / 4, 10);
  bool recalculate_score = !junctions.empty();

  int j = bh.left();
  int r = 0;
  RefID ref_id = bh.ref_id();
  for (size_t c = 0 ; c < cigar.size(); ++c)
    {
      int opcode = cigar[c].opcode;
      int length = cigar[c].length;
      switch(opcode)
	{
	case REF_SKIP:
	case rEF_SKIP:
	  {
	    Junction junc;
	    junc.refid = bh.ref_id();
	    
	    if (opcode == REF_SKIP)
	      {
		junc.left = j - 1;
		junc.right = j + length;
		j += length;
	      }
	    else
	      {
		junc.right = j + 1;
		junc.left = j - length;
		j -= length;
	      }

	    if (recalculate_score)
	      {
		junc.antisense = bh.antisense_splice();
		if (gtf_junctions.find(junc) == gtf_junctions.end())
		  {
		    JunctionSet::const_iterator itr = junctions.find(junc);
		    if (itr == junctions.end())
		      {
			_alignment_score -= bowtie2_max_penalty;
		      }
		    else
		      {
			const int left_cov = coverage.get_coverage(ref_id, junc.left + 1);
			const int right_cov = coverage.get_coverage(ref_id, junc.right - 1);
			const int avg_cov = (left_cov + right_cov) / 2;

			int penalty = bowtie2_max_penalty + 2;
			const int supporting_hits = itr->second.supporting_hits;
			const int left_extent = itr->second.left_extent;
			const int right_extent = itr->second.right_extent;
			float extent_penalty = 0.0f;
			if (left_extent < min_extent || right_extent < min_extent)
			  extent_penalty = 0.5f;
			
			if (supporting_hits >= 5)
			  penalty *= min((float)avg_cov/supporting_hits + extent_penalty, 1.f);

			// daehwan - check this out
			// add two points to prefer junction alignments to other that may span one side of split site.
			// penalty -= 2;

			if (itr->second.gtf_match)
			  penalty -= bowtie2_max_penalty;

			int prev_alignment_score = _alignment_score;
			_alignment_score -= penalty;

			// daehwan - for debugging purposes
			if (bh.insert_id() == 325708 && false)
			  {
			    fprintf(stderr, "junc(%d:%d-%d) %d / (%d + %d) = %d => %d\n",
				    junc.refid, junc.left, junc.right,
				    itr->second.supporting_hits, left_cov, right_cov,
				    prev_alignment_score, _alignment_score);
			    fprintf(stderr, "\textent: %d-%d\n",
				left_extent, right_extent);
			  }
		      }
		  }
		else
		  {
		    _alignment_score += 2;
		  }
	      }
	  }
	  break;
	  
	case MATCH:
	case mATCH:
	  {
	    if (opcode == MATCH)
	      j += length;
	    else
	      j -= length;

	    r += length;
	  }
	  break;
	  
	case DEL:
	case dEL:
	  {
	    Deletion deletion;
	    deletion.refid = bh.ref_id();
	    if (opcode == DEL)
	      {
		deletion.left = j - 1;
		deletion.right = j + length;
		j += length;
	      }
	    else
	      {
		deletion.right = j + 1;
		deletion.left = j - length;
		j -= length;
	      }

	    if (recalculate_score && recalculate_indel_score)
	      {
		DeletionSet::const_iterator itr = deletions.find(deletion);
		if (itr != deletions.end())
		  {
		    const int left_cov = coverage.get_coverage(ref_id, deletion.left + 1);
		    const int right_cov = (length == 1 ? left_cov : coverage.get_coverage(ref_id, deletion.right - 1));
		    const int avg_cov = (left_cov + right_cov) / 2;
		    const int del_penalty = bowtie2_ref_gap_open + bowtie2_ref_gap_cont * length;
		    int addition = del_penalty;
		    
		    const int supporting_hits = itr->second.supporting_hits;
		    const int left_extent = itr->second.left_extent;
		    const int right_extent = itr->second.right_extent;
		    int penalty = 0;
		    if (left_extent < min_extent || right_extent < min_extent)
			  penalty = del_penalty * 0.5f;
			
		    if (avg_cov > 0 && supporting_hits >= 5)
		      addition *= min((float)supporting_hits/avg_cov, 1.f);
		    else
		      addition = 0;

		    addition -= penalty;
		    if (addition < 0)
		      addition = 0;

		    int prev_alignment_score = _alignment_score;
		    _alignment_score += addition;
		    _alignment_score = min(0, _alignment_score);

		    // daehwan - for debug purposes
		    if (bh.insert_id() == 325708 && false)
		      {
			fprintf(stderr, "del(%d:%d-%d) %d / (%d + %d) = %d => %d\n",
				deletion.refid, deletion.left, deletion.right,
				supporting_hits, left_cov, right_cov,
				prev_alignment_score, _alignment_score);
			fprintf(stderr, "\textent: %d-%d\n",
				left_extent, right_extent);
		      }
		  }
	      }
	  }
	  break;

	case INS:
	case iNS:
	  {
	    if (recalculate_score && recalculate_indel_score)
	      {
		string seq = bh.seq().substr(r, length);
		Insertion ins(ref_id, j, seq);
		InsertionSet::const_iterator itr = insertions.find(ins);
		if (itr != insertions.end())
		  {
		    const int supporting_hits = itr->second.supporting_hits;
		    const int left_extent = itr->second.left_extent;
		    const int right_extent = itr->second.right_extent;
		    
		    const int left_cov = coverage.get_coverage(ref_id, j);
		    const int right_cov = coverage.get_coverage(ref_id, j + (opcode == INS ? 1 : -1));
		    const int avg_cov = (left_cov + right_cov) / 2 - supporting_hits;
		    const int ins_penalty = bowtie2_read_gap_open + bowtie2_read_gap_cont * length;
		    int addition = ins_penalty;
		    int extent_penalty = 0.0f;
		    if (left_extent < min_extent || right_extent < min_extent)
		      extent_penalty = ins_penalty * 0.5f;
		    
		    if (avg_cov > 0 && supporting_hits >= 5)
		      addition *= min((float)supporting_hits/avg_cov, 1.f);
		    else
		      addition = 0;

		    addition -= extent_penalty;
		    if (addition < 0)
		      addition = 0;

		    // int prev_alignment_score = _alignment_score;
		    _alignment_score += addition;
		    _alignment_score = min(0, _alignment_score);
		    
		    /*
		      fprintf(stderr, "ins(%d:%d:%s) %d / (%d - %d) = %d => %d (%d)\n",
		      ref_id, ins.left, seq.c_str(),
		      supporting_hits, avg_cov, supporting_hits,
		      prev_alignment_score, _alignment_score, ins_penalty);
		      fprintf(stderr, "\textent: %d-%d\n",
		      left_extent, right_extent);
		    */
		  }
	      }

	    r += length;
	  }
	  break;
	  
	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  // daehwan - implement this later
	  j = length;
	  ref_id = bh.ref_id2();
	  break;
	  
	default:
	  break;
	}
    }
}

/**
 * Establish an ordering on alignments.
 * Prefer aligned reads over unaligned reads
 * Within the space of aligned reads
 * prefer splice-free reads over splice reads, and
 * indel-free reads over indel reads.
 * If a read can either be indel-free or splice-free,
 * prefer the indel-free alignment
 */
bool AlignStatus::operator<(const AlignStatus& rhs) const
{
  if (_alignment_score != rhs._alignment_score)
    return _alignment_score > rhs._alignment_score;

  return false;
}

/**
 * Alignments are only equal if all fields are identical.
 */
bool AlignStatus::operator==(const AlignStatus& rhs) const
{
  return _alignment_score == rhs._alignment_score;
}

bool AlignStatus::operator!=(const AlignStatus& rhs) const
{
  return _alignment_score != rhs._alignment_score;
}
