#ifndef INSERTS_H
#define INSERTS_H
/*
 *  inserts.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include "bwt_map.h"

struct InsertAlignment
{
	InsertAlignment(uint32_t _refid, 
					BowtieHit* _left_alignment, 
					BowtieHit* _right_alignment) : 
	refid(_refid), 
	left_alignment(_left_alignment),
	right_alignment(_right_alignment) {}
	
	uint32_t refid;
	BowtieHit* left_alignment;
	BowtieHit* right_alignment;
};

pair<int, int> pair_distances(const BowtieHit& h1, const BowtieHit& h2);


struct InsertAlignmentGrade
{
	InsertAlignmentGrade() : 
	too_close(false),
	too_far(false),
	one_spliced(false),
	both_spliced(false),
	one_mapped(false),
	both_mapped(false),
	opposite_strands(false),
	consistent_splices(false) {}
	
	InsertAlignmentGrade(const BowtieHit& h1) : 
	too_close(false),
	too_far(false),
	one_spliced(false),
	both_spliced(false),
	one_mapped(false),
	both_mapped(false),
	opposite_strands(false),
	consistent_splices(false)
	{
		if (h1.splice_pos_left != -1)
			one_spliced = true;
		one_mapped = true;
	}
	
	InsertAlignmentGrade(const BowtieHit& h1, 
						 const BowtieHit& h2, 
						 int max_inner_distance, 
						 int min_inner_distance)
	{
		pair<int, int> distances = pair_distances(h1,h2);
		int inner_dist = distances.second;
		
		both_mapped = true;
		
		one_spliced = ((h1.splice_pos_left == -1 && h2.splice_pos_right != -1) ||
					   (h1.splice_pos_left != -1 && h2.splice_pos_right == -1));
		
		both_spliced = (h1.splice_pos_left != -1 && h2.splice_pos_right != -1);
		
		
		too_far = (inner_dist > max_inner_distance);
		
		too_close = (inner_dist < min_inner_distance);
		
		opposite_strands = (h1.antisense_aln != h2.antisense_aln);
		
		consistent_splices = (h1.splice_pos_left != -1 && h2.splice_pos_left != -1 &&
							  h1.antisense_splice != h2.antisense_splice);
		assert(!(too_far && too_close));
	}
	
	InsertAlignmentGrade& operator=(const InsertAlignmentGrade& rhs)
	{
		too_close = rhs.too_close;
		too_far = rhs.too_far;
		one_spliced = rhs.one_spliced;
		both_spliced = rhs.both_spliced;
		one_mapped = rhs.one_mapped;
		both_mapped = rhs.both_mapped;
		opposite_strands = rhs.opposite_strands;
		consistent_splices = rhs.consistent_splices;
		return *this;
	}
	
	// Returns true if rhs is a "happier" alignment for the ends of this insert
	// than this InsertStatus.
	bool operator<(const InsertAlignmentGrade& rhs)
	{
		assert(!(both_mapped && one_mapped));
		// We always prefer a insert alignment with both ends mapped than a
		// singleton
		if (both_mapped != rhs.both_mapped)
		{
			return both_mapped < rhs.both_mapped;
		}
		else
		{
			
			// Prefer an alignment with one end contiguously mapped over a pair 
			// of spliced alignments
			if (both_spliced != rhs.both_spliced)
				return rhs.both_spliced < both_spliced;
			
			// Prefer a pair of contiguously aligned ends to a single 
			// contiguously aligned end and a spliced end
			if (one_spliced != rhs.one_spliced)
				return rhs.one_spliced < one_spliced;
			
			// Prefer a pair that is too close or perfect to one that is too far
			if (too_far && !rhs.too_far)
				return true;
			
			// Prefer a pair that is perfect to one that is too close
			if (too_close && !(rhs.too_close || rhs.too_far))
				return true;
			
			return false;
		}
		
		// We prefer a singleton mapping to an insert with neither end mapped
		if (one_mapped != rhs.one_mapped)
		{
			return one_mapped < rhs.one_mapped;
			
		}
		else
		{
			return rhs.one_spliced < one_spliced;
		}
	}
	
	bool too_close;
	bool too_far;
	
	bool one_spliced;
	bool both_spliced;
	
	bool one_mapped;
	bool both_mapped;
	
	bool opposite_strands;
	bool consistent_splices;
};		

typedef vector<pair<InsertAlignmentGrade, vector<InsertAlignment> > > BestInsertAlignmentTable;
void accept_valid_hits(BestInsertAlignmentTable& best_status_for_inserts);

void best_insert_mappings(uint32_t refid,
						  const string& name,
						  HitList& hits1_in_ref,
						  HitList& hits2_in_ref,
						  BestInsertAlignmentTable& best_insert_alignments);



#endif

