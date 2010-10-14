#ifndef CLOSURES_H
#define CLOSURES_H

/*
 *  closures.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/14/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <set>

#include "common.h"
#include "junctions.h"

struct SpliceMotif
{
	SpliceMotif(uint32_t p, bool l) : pos(p), left_motif(l) {}
	int pos : 31;
	bool left_motif : 1;
	
	bool operator<(const SpliceMotif& rhs) const
	{
		return this->pos < rhs.pos;
	}
	
	bool operator<=(const SpliceMotif& rhs) const
	{
		return this->pos <= rhs.pos;
	}
};

// This simply finds all donor and acceptors in a given string,
// regardless of whether they are within or near an initially mapped 
// region.
template<typename TStr>
class SimpleIntronFinder
{
	vector<SpliceMotif> _hits;

public:
	// This should take three TStr, and we should be using some search function
	// of seqan to look up the hits.  Don't ask.
	SimpleIntronFinder(TStr& ref_str, const string& left_motif, const string& right_motif)
	{
		_hits.clear();
		
		for (size_t i = 0; i < length(ref_str) - 1; ++i)
		{
			// Look at a slice of the reference without creating a copy.
			seqan::Infix<RefSequenceTable::Sequence>::Type curr = seqan::infix(ref_str,i, i + 2);
			if (curr == left_motif)
				_hits.push_back(SpliceMotif(i, true));
			else if (curr == right_motif)
				_hits.push_back(SpliceMotif(i, false));
		}
	}
	
	const vector<size_t>& hits() const { return _hits; }
};

uint32_t searched = 0;
uint32_t closed = 0;

typedef std::set<Junction, skip_count_lt > ClosureJunctionSet;

template<typename TStr, typename IntronFinder>
class JunctionFinder
{
public:
	JunctionFinder(const IntronFinder& intron_finder,
				   uint32_t inner_dist_mean,
				   uint32_t inner_dist_std_dev,
				   uint32_t min_intron_length,
				   uint32_t min_exon_length,
				   uint32_t bowtie_overlap_padding = 0, 
				   uint32_t max_gap = 0,
				   uint32_t max_juncs = 7500000 /*this is the max juncs from a single strand*/) :
	_inner_dist_mean(inner_dist_mean),
	_inner_dist_std_dev(inner_dist_std_dev),
	_min_intron_len(min_intron_length),
	_min_exon_len(min_exon_length),
	_motif_hits(intron_finder.hits()),
	_padding(bowtie_overlap_padding),
	_max_gap(max_gap),
	_max_juncs(max_juncs)
	{
		_left_memo.resize(_motif_hits.size());
		_right_memo.resize(_motif_hits.size());
	}
	
	void possible_junctions(ClosureJunctionSet& potential_splices,
							const BowtieHit& h1, 
							const BowtieHit& h2)
	{
		assert (h1.ref_id() == h2.ref_id());
	
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
		
		if (inner_dist > (int)max_closure_intron_length || 
			inner_dist < (int)min_closure_intron_length)
			return;
		// Bowtie allows some mismatches, so we need to allow for the fact that
		// the donor/acceptor is actually "inside" the alignment for one or both
		// ends of the insert.  That is, intron may be longer than the pair's 
		// inner distance in the genomic coordinate space.  This variable
		// padding number accounts for this fact in a silly way that could be 
		// improved.
		
		int expected_inner_dist = _inner_dist_mean;
		
		if (abs(inner_dist - expected_inner_dist) <= (int)_inner_dist_std_dev)
			return;
		
		++searched;
		
		if (search(potential_splices, h1.ref_id(), minor_hit_end, major_hit_start, expected_inner_dist))
			++closed;
	}
		
private:
	
	vector<SpliceMotif>::const_iterator _motif_upper_bound;
	vector<SpliceMotif>::const_iterator _motif_lower_bound;
	
	/* The following computation identifies possible splice sites.
	 
	 |||||||||||||||-----------------------------------------||||||||||||||||
	 read_len                  inner_dist                      read_len    
	 Where r is the length of a read and x is the *internal* distance between 
	 mates in the *genomic* coordinate space.  We first check to see whether
	 x is too long to be an insert from a single contiguous substring of the 
	 genome.  That is, if x is longer than you'd expect from a unspliced 
	 molecule, than the insert likely spans a splice junction.
	 
	 The expected value of x when the insert doesn't span a splice is the
	 library insert size mean (I) minus twice the read length, plus or minus 
	 the library insert size standard deviation.  If x is longer than that,
	 the insert probably spans a splice.
	 
	 For an insert that spans a splice junction, that junction must fall between
	 the ends of the mates.  Let the right side of the alignment of the left
	 read be called minor_hit_end, and the left side of the alignment of the right 
	 read be major_hit_start.  Then an insert's inner_dist = major_hit_start - 
	 minor_hit_end.  Let the expected inner distance for a insert that doesn't 
	 cross a junction be Insert_mean - 2* read_len. Then a splice-crossing insert 
	 must have inner_dist >= expected_inner_dist.  
	 
	 Let the actual splice position = (splice_left, splice_right). Then 
	 (splice_left - minor_hit_end) + (major_hit_start - splice_right) = 
	 expected_inner_dist +/- std_dev.
	 */
	
	/*
	 
	 The search for possible junctions works as follows.  Given two genomic 
	 coordinates, start, and ref_target, and an expected inner distance between
	 their projections in the transcriptome coordinate space, we can enumerate 
	 paths through the genome that represent possible spliced transcript 
	 segments and connect start and ref_target.  That is, assuming start and 
	 ref_target are both in exonic regions, we can find all possible 
	 genomic subsequences formed by concatenating segments flanked by an "AG" on 
	 the left and a "GT" on the right.  The resulting concatentation is a 
	 possible 'closure' if its length is approximately equal to the expected
	 inner distance between mate pairs for this library.
	 
	 To actually perform the search, this class uses two recursive functions
	 that call each other in an alternating fashion.  search_hop_left_motif takes
	 a genomic coordinate lstart and finds paths to ref_target that start with 
	 a "GT" that occurs in (lstart, ref_target].  search_hop_right_motif takes a 
	 genomic coordinate rstart and finds paths to ref_target that start with an 
	 "AG" that occurs in (rstart, ref_target].  A valid path from start to 
	 ref_target must have a "GT" on the left and an "AG" on the right, so valid 
	 paths may only be found in a call to search_hop_right_motif. 
	 
	 */
	
	bool search(ClosureJunctionSet& splices,
				uint32_t ref_id,
				uint32_t start,
				uint32_t ref_target,
				int expected_inner_dist)
	{
		uint32_t left_start_pos = max((int)start - (int)_padding, 0);
		uint32_t right_end_pos = ref_target + _padding;
		
		SpliceMotif left_start(left_start_pos, true);
		SpliceMotif right_end(right_end_pos, false);
		
		
		_motif_lower_bound = lower_bound(_motif_hits.begin(), 
										 _motif_hits.end(), 
										 left_start);
		
		_motif_upper_bound = upper_bound(_motif_hits.begin(), 
										 _motif_hits.end(), 
										 right_end);
		
		//fprintf(stderr, "Curr junctions = %d\n", splices.size());
		
		/* TODO: consider some kind of long-term memo.  Closures found at the
		 top level of recursion here can be saved and potentially reused in 
		 future searches on overlapping intervals.
		 */
		bool c = false;
		if (search_hop_left_motif(splices,
								  ref_id,
								  _motif_lower_bound, 
								  right_end_pos, 
								  expected_inner_dist, 
								  0,
								  true) != -1)
			c = true;
		
		/* Clear left short-term memo */
		vector<vector<int> >::iterator _lm_begin = _left_memo.begin() + 
		(_motif_lower_bound - _motif_hits.begin());
		
		vector<vector<int> >::iterator _lm_end = _left_memo.begin() + 
		(_motif_upper_bound - _motif_hits.begin()) + 1;
		
		_lm_end = min(_lm_end, _left_memo.end());
		
		/* Clear right short-term memo */
		vector<vector<int> >::iterator _rm_begin = _right_memo.begin() + 
		(_motif_lower_bound - _motif_hits.begin());
		
		vector<vector<int> >::iterator _rm_end = _right_memo.begin() + 
		(_motif_upper_bound - _motif_hits.begin()) + 1;
		
		_rm_end = min(_rm_end, _right_memo.end());
		
		if (_lm_begin < _lm_end)
			fill(_lm_begin, _lm_end, vector<int>());
		
		if (_rm_begin < _rm_end)
			fill(_rm_begin, _rm_end, vector<int>());
		
		while(splices.size() > _max_juncs)
		{
			splices.erase(*splices.rbegin());
		}
		
		return c;
	}
	
	inline int _abs(int x)
	{
		return (x > 0) ? x : -x;
	}
	
	inline int check_memo(const vector<vector<int> >& memo, 
						  size_t motif_pos, 
						  int curr_path_distance, 
						  int target_path_distance,
						  int tolerance)
	{
		const vector<int>& distances = memo[motif_pos];
		for (size_t i = 0; i < distances.size(); ++i)
		{
			int dist = distances[i];
			if (_abs(curr_path_distance + dist - target_path_distance) < tolerance)
				return dist;
		}
		return -1;
	}
	
	int search_hop_left_motif(ClosureJunctionSet& splices,
							  uint32_t ref_id,
							  vector<SpliceMotif>::const_iterator start, 
							  int ref_target,
							  int expected_inner_dist,
							  uint32_t curr_path_length,
							  bool initial_hop = false)
	{
		int in_valid_path = -1;
//		SpliceMotif left_bound(start + (initial_hop ? 0 : _min_exon_len), true);
//		SpliceMotif right_bound(start + (expected_inner_dist - curr_path_length) + (2 * _padding) + _max_gap, true);
		
		vector<SpliceMotif>::const_iterator lb_left = start;
		
		while (lb_left <= _motif_upper_bound && 
			   lb_left < _motif_hits.end() &&
			   (lb_left->pos - start->pos < (initial_hop ? 0 :(int) _min_exon_len)  || !lb_left->left_motif))
		{
			lb_left++;
		}
		
		
		for (vector<SpliceMotif>::const_iterator li = lb_left; 
			 li != _motif_hits.end(); 
			 ++li)
		{
			if (!li->left_motif) //skip over the right motifs
				continue;
			//assert (*li >= SpliceMotif(start, true) || initial_hop);
			int dist = (int)li->pos - (int)start->pos;
			int next_path_len = curr_path_length + dist;
			
			int tolerance = _inner_dist_std_dev +  2 * _padding + _max_gap;
			
			if (next_path_len - expected_inner_dist > tolerance)
			{
				return in_valid_path;
			}
			
			int next = li->pos;
			//assert(next + _min_intron_len <= ref_target);
			if(next + (int)_min_intron_len <= ref_target)
			{
				int ret = -1;
				ret = check_memo(_left_memo, 
								 li - _motif_hits.begin(), 
								 curr_path_length + dist, 
								 expected_inner_dist, 
								 tolerance);
				if (ret != -1)
				{
					in_valid_path = ret + dist;
					_left_memo[li - _motif_hits.begin()].push_back(in_valid_path);
				}
				else
				{
					ret = search_hop_right_motif(splices,
												 ref_id,
												 li,
												 ref_target,
												 expected_inner_dist,
												 next_path_len);
					if (ret != -1)
					{
						in_valid_path = ret + dist;
						_left_memo[li - _motif_hits.begin()].push_back(in_valid_path);
					}
				}
			}
			
		}
		
		return in_valid_path;
		
	}
	
	int search_hop_right_motif(ClosureJunctionSet& splices,
							   uint32_t ref_id,
							   vector<SpliceMotif>::const_iterator start, 
							   int ref_target,
							   int expected_inner_dist,
							   uint32_t curr_path_length)
	{
		int in_valid_path = -1;
		
//		SpliceMotif left_start(start + _min_intron_len - 2, false);
//		SpliceMotif right_start(ref_target, false);
		if (ref_target < start->pos)
		{
			return in_valid_path;
		}
		vector<SpliceMotif>::const_iterator lb_right = start;
		while (lb_right <= _motif_upper_bound && 
			   lb_right < _motif_hits.end() &&
			   (lb_right->pos - start->pos < (int)_min_intron_len - 2 || lb_right->left_motif))
		{
			lb_right++;
		}
			
		for (vector<SpliceMotif>::const_iterator ri = lb_right; 
			 ri != _motif_hits.end(); 
			 ++ri)
		{
			if (ri->left_motif) //skip over the left motifs
				continue;
			if (ri->pos > ref_target)
			{
				return in_valid_path;
			}
			int next = ri->pos + 2;
			int final_hop_dist = (int)ref_target - (int)next;
			int final_path_len = curr_path_length + final_hop_dist;
			
			int tolerance = _inner_dist_std_dev +  2 * _padding + _max_gap;
			
			if (abs(expected_inner_dist - final_path_len) <= tolerance)
			{
				splices.insert(Junction(ref_id,
										start->pos - 1,
										next,
										false,
										start->pos - 1 - next));
				in_valid_path = final_hop_dist;
				_right_memo[ri - _motif_hits.begin()].push_back(final_hop_dist);
			}
			else 
			{
				int ret = -1;
				ret = check_memo(_right_memo, 
								 ri - _motif_hits.begin(), 
								 curr_path_length, 
								 expected_inner_dist, 
								 tolerance);
				if (ret != -1)
				{
					splices.insert(Junction(ref_id,
											start->pos - 1,
											next,
											false,
											start->pos - 1 - next));
					
					in_valid_path = ret;
					_right_memo[ri - _motif_hits.begin()].push_back(ret);
				}
				else
				{
					
					ret = search_hop_left_motif(splices,
												ref_id,
												ri,
												ref_target,
												expected_inner_dist,
												curr_path_length);
					if (ret != -1)
					{
						splices.insert(Junction(ref_id,
												start->pos - 1,
												next,
												false,
												start->pos - 1 - next));
						
						in_valid_path = ret;
						_right_memo[ri - _motif_hits.begin()].push_back(in_valid_path);
					}
				}
			}
		}
		
		return in_valid_path;
	}
	
	uint32_t _inner_dist_mean;
	uint32_t _inner_dist_std_dev;
	uint32_t _min_intron_len;
	uint32_t _min_exon_len;
	
	vector<SpliceMotif> _motif_hits;
	uint32_t _padding;
	uint32_t _max_gap;
	uint32_t _max_juncs;
	vector<vector<int> > _left_memo;
	vector<vector<int> > _right_memo;
};

typedef std::set<pair<size_t, size_t> > CoordSet;
void check_mates(const HitList& hits1_in_ref,
				 const HitList& hits2_in_ref,
				 vector<pair<size_t, size_t> >& happy_mates,
				 vector<size_t>& map1_singletons,
				 vector<size_t>& map2_singletons);

#endif
