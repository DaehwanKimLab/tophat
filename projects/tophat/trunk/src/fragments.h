#ifndef FRAGMENTS_H
#define FRAGMENTS_H
/*
 *  fragments.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include "bwt_map.h"

typedef BowtieHit FragmentAlignment;

struct FragmentAlignmentGrade
{
	FragmentAlignmentGrade() : 
	status(UNALIGNED) {}
	
	FragmentAlignmentGrade(const BowtieHit& h1) 
	{
		if (h1.splice_pos_left != -1)
		{
			status = SPLICED;
		}
		else
		{
			status = CONTIGUOUS;
		}
	}
	
	FragmentAlignmentGrade& operator=(const FragmentAlignmentGrade& rhs)
	{
		status = rhs.status;
		return *this;
	}
	
	// Returns true if rhs is a "happier" alignment for the ends of this insert
	// than this InsertStatus.
	bool operator<(const FragmentAlignmentGrade& rhs)
	{
		return status < rhs.status;
	}
	
	uint8_t status;
};

typedef vector<pair<FragmentAlignmentGrade, vector<FragmentAlignment*> > > BestFragmentAlignmentTable;

void best_fragment_mappings(uint32_t refid,
							const string& name,
							HitList& hits1_in_ref,
							BestFragmentAlignmentTable& best_status_for_fragments);

void accept_valid_hits(BestFragmentAlignmentTable& best_status_for_fragments);
void accept_unique_hits(BestFragmentAlignmentTable& best_status_for_fragments);

#endif
