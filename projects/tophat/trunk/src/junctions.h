#ifndef JUNCTIONS_H
#define JUNCTIONS_H
/*
 *  junctions.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/22/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
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

using namespace std;

struct Junction
{
  Junction (uint32_t ref, uint32_t l, uint32_t r, bool a, int sc = 0)
  : refid(ref), left(l), right(r), antisense(a), skip_count(sc){}
Junction() : refid(0), left(0), right(0), antisense(false), skip_count(0) {}
  uint32_t refid;
  uint32_t left;
  uint32_t right;
  bool antisense;
  
  int skip_count;
  
  bool operator<(const Junction& rhs) const
  {
    if (refid < rhs.refid)
      return true;
    else if (refid > rhs.refid)
      return false;
    
    if (left < rhs.left)
      return true;
    else if (left > rhs.left)
      return false;
    
    if (right < rhs.right)
      return true;
    else if (right > rhs.right)
      return false;
    
    return antisense < rhs.antisense;
  }
  
  bool operator==(const Junction& rhs) const
  {
    return  (refid == rhs.refid && left == rhs.left && right == rhs.right && antisense == rhs.antisense);
  }
  
#if !NDEBUG
  bool valid() const
  {
    return refid != VMAXINT32 && left < right && (left != right);
  }
#endif
};

struct skip_count_lt
{
  bool operator()(const Junction& lhs, const Junction& rhs)
  {
    if (lhs.skip_count != rhs.skip_count)
      return lhs.skip_count < rhs.skip_count;
    return lhs < rhs;
  }
};

struct JunctionStats
{
JunctionStats() : left_extent(0), right_extent(0), left_exon_doc(0), right_exon_doc(0),
    min_splice_mms(0), supporting_hits(0), gtf_match(false), accepted(false) {}

  JunctionStats& merge_with(const JunctionStats& other)
  {
    if (this == &other)
      return *this;

    left_extent = max(left_extent, other.left_extent);
    right_extent = max(right_extent, other.right_extent);
    min_splice_mms = min(min_splice_mms, other.min_splice_mms);
    supporting_hits += other.supporting_hits;

    return *this;
  }
  
  int left_extent;
  int right_extent;
  int left_exon_doc;
  int right_exon_doc;
  int min_splice_mms;
  int supporting_hits;
  bool gtf_match;
  bool accepted;
};

typedef std::map<Junction, JunctionStats> JunctionSet;

// This routine DOES NOT set the real refid!  
pair<Junction, JunctionStats> junction_from_spliced_hit(const BowtieHit& h);

void print_junction(FILE* junctions_out, 
		    const string& name, 
		    const Junction& j, 
		    const JunctionStats& s, 
		    uint32_t junc_id);


void junctions_from_alignment(const BowtieHit& spliced_alignment,
			      JunctionSet& junctions);

void accept_valid_junctions(JunctionSet& junctions,
			    const uint32_t refid,
			    const vector<unsigned short>& DoC,
			    double min_isoform_fraction);

void accept_all_junctions(JunctionSet& junctions,
			  const uint32_t refid);

void print_junctions(FILE* junctions_out, 
		     const JunctionSet& junctions,
		     RefSequenceTable& ref_sequences);

bool accept_if_valid(const Junction& j, JunctionStats& s);

void filter_junctions(JunctionSet& junctions, const JunctionSet& gtf_junctions);

void get_junctions_from_hits(HitStream& hit_stream, 
			     ReadTable& it, 
			     JunctionSet& junctions);

void merge_with(JunctionSet& juncs, const JunctionSet& other_juncs);

#endif
