#ifndef DELETIONS_H
#define DELETIONS_H
/*
 *  deletions.h
 *  TopHat
 *
 *  Created by Ryan Kelley on 11/09/2010 
 *
 */

#include <cstdio>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cstring>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>

#include "bwt_map.h"
#include "junctions.h"
using namespace std;

typedef Junction Deletion;
struct DeletionStats
{
DeletionStats() : left_extent(0), right_extent(0), supporting_hits(0) {}

  DeletionStats& merge_with(const DeletionStats& other)
  {
    if (this == &other)
      return *this;

    left_extent = max(left_extent, other.left_extent);
    right_extent = max(right_extent, other.right_extent);
    supporting_hits += other.supporting_hits;

    return *this;
  }
  
  int left_extent;
  int right_extent;
  int supporting_hits;
};

typedef std::map<Deletion, DeletionStats> DeletionSet;

void deletions_from_alignment(const BowtieHit& spliced_alignment, DeletionSet& junctions);
void deletions_from_spliced_hit(const BowtieHit& bh, vector<pair<Deletion, DeletionStats> >& deletions);
void print_deletions(FILE* deletions_out, const DeletionSet& deletions, RefSequenceTable& ref_sequences);
void merge_with(DeletionSet& deletions, const DeletionSet& other);

#endif
