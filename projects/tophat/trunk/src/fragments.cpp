/*
 *  fragments.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */
#include <cassert>
#include "bwt_map.h"
#include "fragments.h"

void best_fragment_mappings(uint32_t refid,
							const string& name,
							HitList& hits_in_ref,
							BestFragmentAlignmentTable& best_status_for_fragments)
{	
	for (size_t i = 0; i < hits_in_ref.size(); ++i)
	{
		BowtieHit& h1 = hits_in_ref[i];
		
		uint32_t fragment_id = h1.insert_id;
		
		FragmentAlignmentGrade s(h1);
		
		pair<FragmentAlignmentGrade, vector<FragmentAlignment*> >& fragment_best
		= best_status_for_fragments[fragment_id];
		FragmentAlignmentGrade& current = fragment_best.first;
		// Is the new status better than the current best one?
		if (current < s)
		{
			fragment_best.second.clear();
			current = s;
			fragment_best.second.push_back(&h1);
		}
		else if (! (s < current)) // is it just as good?
		{
			fragment_best.second.push_back(&h1);
		}
	}
}

bool valid_fragment_alignment(const FragmentAlignmentGrade& g, const FragmentAlignment& a)
{
	// stub
	return true;
}

void accept_unique_hits(BestFragmentAlignmentTable& best_status_for_fragments)
{
	for (size_t i = 0; i < best_status_for_fragments.size(); ++i)
	{
		const pair<FragmentAlignmentGrade, vector<FragmentAlignment*> >& fragment_best
		= best_status_for_fragments[i];
		if (fragment_best.first.status == SPLICED)
		{
			//			if (fragment_best.second.size() == 1)
			//				fragment_best.second[0].alignment->accepted = true;
			//			else
			//				fragment_best.second[0].alignment->accepted = false;
			for (size_t j = 0; j < fragment_best.second.size(); ++j)
			{
				FragmentAlignment& a = *(fragment_best.second[j]);
				a.accepted = true;
			}
			for (size_t j = 0; j < fragment_best.second.size(); ++j)
			{
				FragmentAlignment& a = *(fragment_best.second[j]);
				for (size_t k = 0; k < fragment_best.second.size(); ++k)
				{
					FragmentAlignment& b = *(fragment_best.second[k]);
					if (k != j && 
						(a.left == b.left || 
						 a.right == b.right))
					{
						a.accepted = false;
						b.accepted = false;
					}
				}
			}
		}
		else
		{
			for (size_t j = 0; j < fragment_best.second.size(); ++j)
			{
				FragmentAlignment& a = *(fragment_best.second[j]);
				assert (a.splice_pos_left == -1);
				a.accepted = true;
			}
		}
	}
}

void accept_valid_hits(BestFragmentAlignmentTable& best_status_for_fragments)
{
	for (size_t i = 0; i < best_status_for_fragments.size(); ++i)
	{
		const pair<FragmentAlignmentGrade, vector<FragmentAlignment*> >& fragment_best
		= best_status_for_fragments[i];
		
		for (size_t j = 0; j < fragment_best.second.size(); ++j)
		{
			FragmentAlignment& a = *(fragment_best.second[j]);
			bool valid = valid_fragment_alignment(fragment_best.first, a);
			a.accepted = valid;
		}
	}
}


