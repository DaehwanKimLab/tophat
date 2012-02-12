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
			 const Coverage& coverage) {
  const vector<CigarOp>& cigar = bh.cigar();
  _alignment_score = bh.alignment_score();

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
		  
		  int penalty = bowtie2_max_penalty;
		  const int supporting_hits = itr->second.supporting_hits;
		  if (supporting_hits > 0)
		    penalty *= min((float)avg_cov/supporting_hits, 1.f);

		  _alignment_score -= penalty;
		  _alignment_score = min(0, _alignment_score);
		  
		  /*
		  fprintf(stderr, "junc(%d:%d-%d) %d / (%d + %d) = %d => %d\n",
			  junc.refid, junc.left, junc.right,
			  itr->second.supporting_hits, left_cov, right_cov,
			  _alignment_score + penalty, _alignment_score);
		  */
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
	    Junction junc;
	    junc.refid = bh.ref_id();
	    if (opcode == DEL)
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
	    
	    DeletionSet::const_iterator itr = deletions.find(junc);
	    if (itr != deletions.end())
	      {
		const int left_cov = coverage.get_coverage(ref_id, junc.left + 1);
		const int right_cov = (length == 1 ? left_cov : coverage.get_coverage(ref_id, junc.right - 1));
		const int avg_cov = (left_cov + right_cov) / 2;
		const int del_penalty = bowtie2_ref_gap_open + bowtie2_ref_gap_cont * length;
		int addition = del_penalty;

		const int supporting_hits = itr->second;
		if (avg_cov > 0)
		  addition *= min((float)supporting_hits/avg_cov * 2.f, 1.f);

		_alignment_score += addition;
		_alignment_score = min(0, _alignment_score);
		
		/*
		fprintf(stderr, "del(%d:%d-%d) %d / (%d + %d) = %d => %d (%d)\n",
			junc.refid, junc.left, junc.right,
			itr->second, left_cov, right_cov,
			_alignment_score - addition, _alignment_score, del_penalty);
		*/
	      }
	  }
	  break;

	case INS:
	case iNS:
	  {
	    string seq = bh.seq().substr(r, length);
	    Insertion ins(ref_id, j, seq);
	    InsertionSet::const_iterator itr = insertions.find(ins);
	    if (itr != insertions.end())
	      {
		const int cov = coverage.get_coverage(ref_id, j);
		const int ins_penalty = bowtie2_read_gap_open + bowtie2_read_gap_cont * length;
		int addition = ins_penalty;

		const int supporting_hits = itr->second;
		if (cov > 0)
		  addition *= min((float)supporting_hits/cov * 2.f, 1.f);

		_alignment_score += addition;
		_alignment_score = min(0, _alignment_score);
		
		/*
		fprintf(stderr, "ins(%d:%d:%s) %d / %d = %d => %d (%d)\n",
			ref_id, ins.left, seq.c_str(),
			itr->second, cov,
			_alignment_score - addition, _alignment_score, ins_penalty);
		*/
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
