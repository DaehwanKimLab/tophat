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
InsertAlignment(uint64_t _refid, 
		BowtieHit* _left_alignment, 
		BowtieHit* _right_alignment) : 
  refid(_refid), 
    left_alignment(_left_alignment),
    right_alignment(_right_alignment) {}
  
  uint64_t refid;
  BowtieHit* left_alignment;
  BowtieHit* right_alignment;
};

pair<int, int> pair_distances(const BowtieHit& h1, const BowtieHit& h2);

bool gap_lt(const pair<int, int>& lhs, const pair<int, int>& rhs);

struct InsertAlignmentGrade
{
InsertAlignmentGrade() : 
  too_close(false),
    too_far(false),
    num_spliced(0),
    num_mapped(0),
    opposite_strands(false),
    consistent_splices(false),
    longest_ref_skip(0x7FFFFu),
    edit_dist(0x1F),
    fusion(true),
    inner_dist(99999999),
    alignment_score(std::numeric_limits<int>::min())
  {}

InsertAlignmentGrade(const BowtieHit& h1, bool fusion = false) : 
  too_close(false),
    too_far(false),
    num_spliced(0),
    num_mapped(0),
    opposite_strands(false),
    consistent_splices(false),
    edit_dist(0x1F),
    fusion(fusion),
    num_alignments(0),
    inner_dist(99999999),
    alignment_score(std::numeric_limits<int>::min())
  {
    if (!h1.contiguous())
      num_spliced++;
    
    num_mapped = 1;
    
    longest_ref_skip = min(0x3FFFFu, (unsigned int)get_longest_ref_skip(h1) / 100);
    edit_dist = h1.edit_dist();
    num_alignments = 1;
  }
  
InsertAlignmentGrade(const BowtieHit& h1, 
		     const BowtieHit& h2, 
		     bool fusion = false) :
  too_close(false),
    too_far(false),
    num_spliced(0),
    num_mapped(0),
    opposite_strands(false),
    consistent_splices(false),
    edit_dist(0x1F),
    fusion(fusion),
    num_alignments(0),
    alignment_score(std::numeric_limits<int>::min())
  {
    int min_inner_distance = inner_dist_mean - inner_dist_std_dev;
    int max_inner_distance = inner_dist_mean + inner_dist_std_dev;
  
    pair<int, int> distances = pair_distances(h1, h2);
    inner_dist = distances.second;
    
    num_mapped = 2;
    
    if (!h1.contiguous())
      num_spliced++;
    if (!h2.contiguous())
      num_spliced++;
    
    too_far = (inner_dist > max_inner_distance);
    too_close = (inner_dist < min_inner_distance);
    opposite_strands = (h1.antisense_align() != h2.antisense_align());
    consistent_splices = (num_spliced == 2 && h1.antisense_splice() == h2.antisense_splice());
    
    uint32_t ls = max(get_longest_ref_skip(h1), get_longest_ref_skip(h2));
    longest_ref_skip = min (ls / 100, 0x3FFFFu);
    edit_dist = h1.edit_dist() + h2.edit_dist();
    num_alignments = 1;
    assert(!(too_far && too_close));

    //
    static const int penalty_for_long_inner_dist = bowtie2_max_penalty;
    alignment_score = h1.alignment_score() + h2.alignment_score();
    if (!fusion)
      {
	if (too_far)
	  {
	    int penalty = penalty_for_long_inner_dist;
	    if (inner_dist - max_inner_distance < inner_dist_std_dev)
	      {
		penalty = penalty_for_long_inner_dist / 2;
	      }
	      
	    alignment_score -= penalty;
	  }
	else if (too_close)
	  {
	    int penalty = min(penalty_for_long_inner_dist/2,
			      (min_inner_distance - inner_dist) / inner_dist_std_dev + 1);
	    alignment_score -= penalty;
	  }
	
	static const int penalty_for_same_strand = bowtie2_max_penalty;
	if (!opposite_strands)
	  alignment_score -= penalty_for_same_strand;
      }
  }
  
  InsertAlignmentGrade& operator=(const InsertAlignmentGrade& rhs)
  {
    too_close = rhs.too_close;
    too_far = rhs.too_far;
    num_spliced = rhs.num_spliced;
    num_mapped = rhs.num_mapped;
    opposite_strands = rhs.opposite_strands;
    consistent_splices = rhs.consistent_splices;
    longest_ref_skip = rhs.longest_ref_skip;
    edit_dist = rhs.edit_dist;
    fusion = rhs.fusion;
    num_alignments = rhs.num_alignments;
    inner_dist = rhs.inner_dist;
    alignment_score = rhs.alignment_score;
    return *this;
  }
    
    static int get_longest_ref_skip(const BowtieHit& h1)
  {
    vector<pair<int, int> > gaps;
    h1.gaps(gaps);
    if (gaps.empty())
      return 0;
    
    vector<pair<int, int> >::iterator max_itr = max_element(gaps.begin(), gaps.end(), gap_lt);
    return abs(max_itr->second - max_itr->first);
  }
  
  // Returns true if rhs is a "happier" alignment for the ends of this insert
  // than this InsertStatus.
  bool operator<(const InsertAlignmentGrade& rhs);
  
  bool happy() const
  {
    return num_mapped == 2 && opposite_strands && (num_spliced != 2 || consistent_splices) && !too_far;
  }

  int align_score() { return alignment_score; }

  bool is_fusion() const { return fusion; }
  
  bool too_close;
  bool too_far;
  
  uint8_t num_spliced;
  uint8_t num_mapped;
  
  bool opposite_strands;
  bool consistent_splices;
  uint32_t longest_ref_skip; // in 100s of bp
  unsigned char edit_dist;
  bool fusion;
  int num_alignments; // number of equally good alignments for the insert 
  int inner_dist; // distance between inner edges of mates

  int alignment_score;
};		

typedef vector<pair<InsertAlignmentGrade, vector<InsertAlignment> > > BestInsertAlignmentTable;
void accept_valid_hits(BestInsertAlignmentTable& best_status_for_inserts);
void accept_all_best_hits(BestInsertAlignmentTable& best_status_for_inserts);

void best_insert_mappings(uint64_t refid,
			  ReadTable& it,
			  HitList& hits1_in_ref,
			  HitList& hits2_in_ref,
			  BestInsertAlignmentTable& best_insert_alignments,
			  bool prefer_shorter_pairs = false);


void insert_best_pairings(RefSequenceTable& rt,
			  ReadTable& it,
                          HitTable& hits1,
                          HitTable& hits2,
                          BestInsertAlignmentTable& best_pairings, 
			  bool prefer_shorter_pairs = false);
#endif
