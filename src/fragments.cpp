/*
 *  fragments.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include "bwt_map.h"
#include "fragments.h"

void best_fragment_mappings(uint64_t refid,
			    const string& name,
			    HitList& hits_in_ref,
			    ReadTable& it, 
			    BestFragmentAlignmentTable& best_status_for_fragments)
{
  return;
#if 0
  for (size_t i = 0; i < hits_in_ref.size(); ++i)
    {
      BowtieHit& h1 = hits_in_ref[i];
      uint64_t fragment_id = h1.insert_id();
      uint32_t obs_order = it.observation_order(fragment_id);

      JunctionSet dummy;
      FragmentAlignmentGrade s(h1, dummy);
      
      pair<FragmentAlignmentGrade, vector<FragmentAlignment*> >& fragment_best
        = best_status_for_fragments[obs_order];
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
#endif
}

bool valid_fragment_alignment(const FragmentAlignmentGrade& g, const FragmentAlignment& a)
{
  // stub
  return true;
}
