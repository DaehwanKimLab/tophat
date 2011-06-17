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
#include "align_status.h"

using namespace std;

AlignStatus::AlignStatus()
{
  _aligned = false;
  _indelFreeAlignment = false;
  _unannotatedSpliceFreeAlignment = false;
} 

/**
 * Parse the cigar string of a BowtieHit in order to determine the alignment status.
 */
AlignStatus::AlignStatus(const BowtieHit& bh, const JunctionSet& gtf_junctions)
{
  const vector<CigarOp>& cigar = bh.cigar();
  _aligned = cigar.size() > 0;
  _indelFreeAlignment = true;
  _unannotatedSpliceFreeAlignment = true;

  int j = bh.left();
  for (size_t c = 0 ; c < cigar.size(); ++c)
    {
      Junction junc;
      switch(cigar[c].opcode)
	{
	case REF_SKIP:
	  junc.refid = bh.ref_id();
	  junc.left = j;
	  junc.right = junc.left + cigar[c].length;
	  junc.antisense = bh.antisense_splice();
	  j += cigar[c].length;

	  if (gtf_junctions.find(junc) == gtf_junctions.end())
	    _unannotatedSpliceFreeAlignment = false;
	  break;	  
	case MATCH:
	  j += cigar[c].length;
	  break;
	case DEL:
	  j += cigar[c].length;
	  _indelFreeAlignment = false;
	  break;
	case INS:
	  _indelFreeAlignment = false;
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
	int lhs_value = _aligned ? 1 : 0;
	lhs_value += _indelFreeAlignment ? 4 : 0;
	lhs_value += _unannotatedSpliceFreeAlignment ? 2 : 0;

	int rhs_value = rhs._aligned ? 1 : 0;
	rhs_value += rhs._indelFreeAlignment ? 4 : 0;
	rhs_value += rhs._unannotatedSpliceFreeAlignment ? 2 : 0;
	return lhs_value < rhs_value;
}

/**
 * Alignments are only equal if all fields are identical.
 */
bool AlignStatus::operator==(const AlignStatus& rhs) const
{
	return ((_aligned == rhs._aligned) &&
		(_indelFreeAlignment == rhs._indelFreeAlignment) &&
		(_unannotatedSpliceFreeAlignment == rhs._unannotatedSpliceFreeAlignment));
}

bool AlignStatus::operator!=(const AlignStatus& rhs) const
{
	return !((*this) == rhs);
}



