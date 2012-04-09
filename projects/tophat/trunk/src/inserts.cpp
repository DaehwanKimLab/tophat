/*
 *  inserts.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cassert>
#include <map>
#include <set>
#include "bwt_map.h"
#include "common.h"
#include "inserts.h"
using namespace std;

bool InsertAlignmentGrade::operator<(const InsertAlignmentGrade& rhs)
{
  if (fusion && !rhs.fusion)
    return true;

  if (!fusion && rhs.fusion)
    return false;

  // We always prefer a insert alignment with both ends mapped than a
  // singleton
  if (num_mapped != rhs.num_mapped)
    {
      return num_mapped < rhs.num_mapped;
    }
  else if (num_mapped == 2)
    {
      // daehwan - I'm testing this!
#if 1
      if (alignment_score != rhs.alignment_score)
	return alignment_score < rhs.alignment_score;

      return false;
#else
      // if significant difference in their inner mate distances 
      if (abs(rhs.inner_dist - inner_dist) >= 30)
	{
	  // Prefer a pair that is too close or perfect to one that is too far
	  if (too_far && !rhs.too_far)
	    return true;
	  
	  // Prefer a pair that is perfect to one that is too close
	  if (too_close && !(rhs.too_close || rhs.too_far))
	    return true;
	  
	  // Prefer closer mates
	  if (rhs.inner_dist < inner_dist)
	    return true;
	}

      if (edit_dist != rhs.edit_dist)
	return rhs.edit_dist < edit_dist;

      // daehwan - do somethings here
      if (!bowtie2)
	{
	  // Prefer shorter introns
	  if (longest_ref_skip != rhs.longest_ref_skip)
	    return rhs.longest_ref_skip < longest_ref_skip;
	}
#endif
	  
      return false;
    }
  else
    {
      // We prefer a singleton mapping to an insert with neither end mapped
      if (num_mapped != rhs.num_mapped)
	{
	  return num_mapped < rhs.num_mapped; // if RHS has MORE READS, RHS is BETTER (lhs < rhs)
	}
      else
	{
	  if (rhs.num_spliced != num_spliced)
	    return rhs.num_spliced < num_spliced;// if RHS is LESS SPLICED, RHS is BETTER (lhs < rhs)
	  
	  // Prefer shorter introns
	  if (longest_ref_skip != rhs.longest_ref_skip)
	    return rhs.longest_ref_skip < longest_ref_skip; // if RHS intron is SHORTER, RHS is BETTER (lhs < rhs)
	  
	  return rhs.edit_dist < edit_dist; // if RHS edit is LOWER, RHS is BETTER (lhs < rhs)
	}
    }
  return false;
}

bool gap_lt(const pair<int, int>& lhs, const pair<int, int>& rhs)
{
  return abs(lhs.second - lhs.first) < abs(rhs.second - rhs.first);
}

pair<int, int> pair_distances(const BowtieHit& h1, const BowtieHit& h2)
{
  int minor_hit_start, major_hit_start;
  int minor_hit_end, major_hit_end;
  if (h1.left() < h2.left())
    {
      minor_hit_start = (int)h1.left();
      minor_hit_end = (int)h1.right();
      major_hit_start = (int)h2.left();
      major_hit_end = (int)h2.right();
    }
  else
    {
      minor_hit_start = (int)h2.left();
      minor_hit_end = (int)h2.right();
      major_hit_start = (int)h1.left();
      major_hit_end = (int)h1.right();
    }
  
  int inner_dist = major_hit_start - minor_hit_end;
  int outer_dist = major_hit_end - minor_hit_start;
  return make_pair(outer_dist, inner_dist);
}

void best_insert_mappings(uint64_t refid,
			  ReadTable& it,
			  HitList& hits1_in_ref,
			  HitList& hits2_in_ref,
			  BestInsertAlignmentTable& best_status_for_inserts,
			  bool prefer_shorter_pairs)
{	
  long chucked_for_shorter_pair = 0;
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
      
      uint32_t obs_order = it.observation_order(h1.insert_id());
      
      for (HitList::iterator f = range_pair.first;
	   f != range_pair.second;
	   ++f)
	{
	  BowtieHit& h2 = *f;
	  
	  if (h1.insert_id() == h2.insert_id())
	    {
	      InsertAlignmentGrade s(h1, h2);
	      
	      pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
		= best_status_for_inserts[obs_order];
	      InsertAlignmentGrade& current = insert_best.first;
	      // Is the new status better than the current best one?
	      if (current < s)
		{
		  insert_best.second.clear();
		  current = s;
		  insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
		}
	      else if (!(s < current))
		{
		  if (prefer_shorter_pairs && current.num_mapped == 2)
		    {
		      pair<int, int> dc = pair_distances(*(insert_best.second[0].left_alignment), *(insert_best.second[0].right_alignment));
		      pair<int, int> ds = pair_distances(h1,h2);
		      if (ds.second < dc.second)
			{
			  chucked_for_shorter_pair += insert_best.second.size();
			  insert_best.second.clear();
			  current = s;
			  insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
			}
		    }
		  else
		    {
		      insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
		    }
		}
	      
	      marked.insert(f - hits2_in_ref.begin());
	      found_hit = true;
	    }
	  
	}
      if (!found_hit)
	{
	  pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
	    = best_status_for_inserts[obs_order];
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
      uint32_t obs_order = it.observation_order(h2.insert_id());
      pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
	= best_status_for_inserts[obs_order];
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
  fprintf(stderr, "Chucked %ld pairs for shorter pairing of same mates\n", chucked_for_shorter_pair);
}

int long_spliced = 0;
int short_spliced = 0;
int singleton_splices = 0;

bool valid_insert_alignment(const InsertAlignmentGrade& g, const InsertAlignment& a)
{
  if (!a.left_alignment || !a.right_alignment)
    return false;
  if (g.num_mapped == 2)
    { 
      // Take all the contiguously mapped pairs
      if (g.num_spliced == 0)
	return true;
      
      // Take the pairs that include one or more spliced reads as long as 
      // the inner dist isn't too big
      //		if (g.one_spliced || g.both_spliced)
      //		{
      //			if (g.too_far || g.too_close)
      //				return false;
      //		}
      return true;
    }
  return false;
}


void insert_best_pairings(RefSequenceTable& rt,
			  ReadTable& it,
                          HitTable& hits1,
                          HitTable& hits2,
                          BestInsertAlignmentTable& best_pairings,
			  bool prefer_shorter_pairs)
{
  for(RefSequenceTable::const_iterator ci = rt.begin();
      ci != rt.end();
      ++ci)
    {
      
      // Tracks the number of singleton ALIGNMENTS, not the number of singleton
      // READS in each Bowtie map.
      vector<size_t> map1_singletons;
      vector<size_t> map2_singletons;
      vector<pair<size_t, size_t> > happy_mates;
      
      uint64_t ref_id = ci->second.observation_order;
      HitList* hits1_in_ref = hits1.get_hits(ref_id);
      HitList* hits2_in_ref = hits2.get_hits(ref_id);
      
      if (!hits1_in_ref || !hits2_in_ref)
	continue;
      
      //if (verbose)
      //    fprintf(stderr, "Looking for best insert mappings in %s\n", name.c_str());
      
      best_insert_mappings(ref_id,
			   it,
			   *hits1_in_ref,
			   *hits2_in_ref,
			   best_pairings, 
			   prefer_shorter_pairs);
    }
}
