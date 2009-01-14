/*
 *  inserts.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <map>
#include <set>
#include "bwt_map.h"
#include "common.h"
#include "inserts.h"

using namespace std;

pair<int, int> pair_distances(const BowtieHit& h1, 
							  const BowtieHit& h2)
{
	int minor_hit_start, major_hit_start;
	int minor_hit_end, major_hit_end;
	if (h1.left < h2.left)
	{
		minor_hit_start = (int)h1.left;
		minor_hit_end = (int)h1.right;
		major_hit_start = (int)h2.left;
		major_hit_end = (int)h2.right;
	}
	else
	{
		minor_hit_start = (int)h2.left;
		minor_hit_end = (int)h2.right;
		major_hit_start = (int)h1.left;
		major_hit_end = (int)h1.right;
	}
	
	int inner_dist = major_hit_start - minor_hit_end;
	int outer_dist = major_hit_end - minor_hit_start;
	return make_pair(outer_dist, inner_dist);
}

void best_insert_mappings(uint32_t refid,
						  const string& name,
						  HitList& hits1_in_ref,
						  HitList& hits2_in_ref,
						  BestInsertAlignmentTable& best_status_for_inserts)
{	
	std::set<size_t> marked;
	HitList::iterator last_good = hits2_in_ref.begin();
	
	for (size_t i = 0; i < hits1_in_ref.size(); ++i)
	{
		BowtieHit& h1 = hits1_in_ref[i];
		pair<HitList::iterator, HitList::iterator> range_pair;
		range_pair = equal_range(last_good, hits2_in_ref.end(),
								 h1, hit_insert_id_lt);
		bool found_hit = false;
		if (range_pair.first != range_pair.second)
			last_good = range_pair.first;
		
		for (HitList::iterator f = range_pair.first;
			 f != range_pair.second;
			 ++f)
		{
			BowtieHit& h2 = *f;
			
			if (h1.insert_id == h2.insert_id)
			{
				uint32_t insert_id = h2.insert_id;
				
				int min_mate_inner_dist = insert_len - h1.read_len() - 
				h2.read_len() - insert_len_std_dev;
				if (max_mate_inner_dist == -1)
				{
					max_mate_inner_dist = insert_len - h1.read_len() - 
					h2.read_len() + insert_len_std_dev;
				}
				InsertAlignmentGrade s(h1, h2, max_mate_inner_dist, min_mate_inner_dist);
				
				pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
				= best_status_for_inserts[insert_id];
				InsertAlignmentGrade& current = insert_best.first;
				// Is the new status better than the current best one?
				if (current < s)
				{
					insert_best.second.clear();
					current = s;
					insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
				}
				else if (! (s < current))
				{
					insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
				}
				
				marked.insert(f - hits2_in_ref.begin());
				found_hit = true;
			}
			
		}
		if (!found_hit)
		{
			pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
			= best_status_for_inserts[h1.insert_id];
			InsertAlignmentGrade& current = insert_best.first;	
			
			InsertAlignmentGrade s(h1);
			
			if (current < s)
			{
				insert_best.second.clear();
				current = s;
				insert_best.second.push_back(InsertAlignment(refid, &h1, NULL));
			}
			else if (! (s < current))
			{
				insert_best.second.push_back(InsertAlignment(refid, &h1, NULL));
			}
			
		}
	}
	
	for (size_t i = 0; i < hits2_in_ref.size(); ++i)
	{
		BowtieHit& h2 = hits2_in_ref[i];
		
		pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
		= best_status_for_inserts[h2.insert_id];
		InsertAlignmentGrade& current = insert_best.first;	
		
		InsertAlignmentGrade s(h2);
		// Did we include h2 as part of a pairing already, or is this first time
		// we've seen it?  If so, it's a singleton.
		if (marked.find(i) == marked.end())
		{
			if (current < s)
			{
				insert_best.second.clear();
				current = s;
				insert_best.second.push_back(InsertAlignment(refid, NULL, &h2));
			}
			else if (! (s < current))
			{
				insert_best.second.push_back(InsertAlignment(refid, NULL, &h2));
			}
		}
	}	
	
}

int long_spliced = 0;
int short_spliced = 0;
int singleton_splices = 0;

bool valid_insert_alignment(const InsertAlignmentGrade& g, const InsertAlignment& a)
{
	if (!a.left_alignment || !a.right_alignment)
		return false;
	if (g.both_mapped)
	{ 
		// Take all the contiguously mapped pairs
		if (!g.one_spliced && !g.both_spliced)
			return true;
		
		// Take the pairs that include one or more spliced reads as long as 
		// the inner dist isn't too big
		if (g.one_spliced || g.both_spliced)
		{
			if (g.too_far || g.too_close)
				return false;
		}
		return true;
	}
	return false;
}

void accept_valid_hits(BestInsertAlignmentTable& best_status_for_inserts)
{
	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
	{
		const pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
		= best_status_for_inserts[i];
		for (size_t j = 0; j < insert_best.second.size(); ++j)
		{
			const InsertAlignment& a = insert_best.second[j];
			bool valid = valid_insert_alignment(insert_best.first, a);
			if (a.left_alignment)
			{
				a.left_alignment->accepted = valid;
			}
			if (a.right_alignment)
			{
				a.right_alignment->accepted = valid;
			}
		}
	}
}