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

AlignStatus::AlignStatus(){
	_aligned = 0;
	_indelFreeAlignment = 0;
	_spliceFreeAlignment = 0;
} 

/**
 * Parse the cigar string of a BowtieHit in order to determine the alignment status.
 */
AlignStatus::AlignStatus(const BowtieHit& bh){
	const vector<CigarOp>& cigar = bh.cigar();
	_aligned = 0;
	_indelFreeAlignment = 0;
	_spliceFreeAlignment = 0;
	if(cigar.size() > 0){
		_aligned = 1;
		_indelFreeAlignment = 1;
		_spliceFreeAlignment = 1;
		vector<CigarOp>::const_iterator it;
		for(it = cigar.begin(); it < cigar.end(); ++it){
			if(it->opcode == INS || it->opcode == DEL){
				_indelFreeAlignment = 0;
			}
			if(it->opcode == REF_SKIP){
				_spliceFreeAlignment = 0;
			}
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
	lhs_value += _spliceFreeAlignment ? 2 : 0;

	int rhs_value = rhs._aligned ? 1 : 0;
	rhs_value += rhs._indelFreeAlignment ? 4 : 0;
	rhs_value += rhs._spliceFreeAlignment ? 2 : 0;
	return lhs_value < rhs_value;
}

/**
 * Alignments are only equal if all fields are identical.
 */
bool AlignStatus::operator==(const AlignStatus& rhs) const
{
	return ((_aligned == rhs._aligned) && (_indelFreeAlignment == rhs._indelFreeAlignment) && (_spliceFreeAlignment == rhs._spliceFreeAlignment));
}

bool AlignStatus::operator!=(const AlignStatus& rhs) const
{
	return !((*this) == rhs);
}



