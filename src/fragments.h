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
#include "align_status.h"

typedef BowtieHit FragmentAlignment;

struct FragmentAlignmentGrade
{
  FragmentAlignmentGrade() 
  {
    num_alignments = 0;
    status = AlignStatus();
  }
  
  FragmentAlignmentGrade(const BowtieHit& h1,
			 const JunctionSet& gtf_junctions,
			 const JunctionSet& junctions,
			 const InsertionSet& insertions,
			 const DeletionSet& deletions,
			 const FusionSet& fusions,
			 const Coverage& coverage) 
  {
    status = AlignStatus(h1,
			 gtf_junctions,
			 junctions,
			 insertions,
			 deletions,
			 fusions,
			 coverage);
    num_alignments = 1;
  }
  
  FragmentAlignmentGrade& operator=(const FragmentAlignmentGrade& rhs)
  {
    status = rhs.status;
    num_alignments = rhs.num_alignments;
    
    return *this;
  }
    
  // Returns true if rhs is a "happier" alignment for the ends of this insert
  // than this InsertStatus.
    
  bool operator<(const FragmentAlignmentGrade& rhs)
  {
    return status < rhs.status;
  }
  
  AlignStatus status;
  int num_alignments; // number of equally good alignments for this fragment 
};

typedef vector<pair<FragmentAlignmentGrade, vector<FragmentAlignment*> > > BestFragmentAlignmentTable;

void best_fragment_mappings(uint64_t refid,
			    const string& name,
			    HitList& hits1_in_ref,
			    ReadTable& it,
			    BestFragmentAlignmentTable& best_status_for_fragments);

void accept_best_hits(BestFragmentAlignmentTable& best_status_for_fragments);
void accept_unique_hits(BestFragmentAlignmentTable& best_status_for_fragments);

#endif
