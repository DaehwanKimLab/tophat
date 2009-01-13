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
//#include <map>
//#include <algorithm>
#include <set>
//#include <stdexcept>
#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>

#include "bwt_map.h"

using namespace std;

typedef std::set<pair<size_t, size_t> > SpliceSet;

// This simply finds all donor and acceptors in a given string,
// regardless of whether they are within or near an initially mapped 
// region.
template<typename TStr>
class SimpleIntronFinder
{
	vector<size_t> _l_hits;
	vector<size_t> _r_hits;
public:
	// This should take three TStr, and we should be using some search function
	// of seqan to look up the hits.  Don't ask.
	SimpleIntronFinder(TStr& ref_str, const string& left_motif, const string& right_motif)
	{
		_l_hits.clear();
		_r_hits.clear();
		
		for (size_t i = 0; i < length(ref_str) - 1; ++i)
		{
			// Look at a slice of the reference without creating a copy.
			seqan::Infix<seqan::String<seqan::Dna5, seqan::Alloc<> > >::Type curr = seqan::infix(ref_str,i, i + 2);
			if (curr == left_motif)
				_l_hits.push_back(i);
			else if (curr == right_motif)
				_r_hits.push_back(i);
		}
	}
	
	const vector<size_t>& l_hits() const { return _l_hits; }
	const vector<size_t>& r_hits() const { return _r_hits; }
};

template<typename TStr>
class MappedIntronFinder
{	
	vector<size_t> _l_hits;
	vector<size_t> _r_hits;
	
	public:
	// This should take three TStr, and we should be using some search function
	// of seqan to look up the hits.  Don't ask.
	MappedIntronFinder(TStr& ref_str, 
					   const vector<const HitList*>& all_hits,
					   const string& left_motif, 
					   const string& right_motif,
					   uint32_t extension)
	{
		_l_hits.clear();
		_r_hits.clear();
		
		set<size_t> l_tmp;
		set<size_t> r_tmp;
		fprintf (stderr, "\tLooking for motifs...");
		for (size_t i = 0; i < length(ref_str) - 1; ++i)
		{
			seqan::Infix<seqan::String<seqan::Dna5, 
									   seqan::Alloc<> > >::Type curr
				= seqan::infix(ref_str,i, i + 2);
			if (curr == left_motif)
				l_tmp.insert(i);
			else if (curr == right_motif)
				r_tmp.insert(i);
		}
		fprintf (stderr, "done\n");
		fprintf (stderr, "\tFiltering motifs through the map...");
		set<size_t> l_keep;
		set<size_t> r_keep;
		
		for (vector<const HitList*>::const_iterator hi = all_hits.begin();
			 hi != all_hits.end(); 
			 ++hi)
		{
			const HitList* hl = *hi;
			assert (hl);
			for (HitList::const_iterator ci = hl->begin();
				ci != hl->end();
				++ci)
			{
				size_t left = max(0,  (int)ci->left - (int)extension);
				size_t right = (int)ci->right + (int)extension;

				set<size_t>::iterator lb = l_tmp.lower_bound(left);
				set<size_t>::iterator ub = l_tmp.upper_bound(right);
				l_keep.insert(lb, ub);
				
				lb = r_tmp.lower_bound(left);
				ub = r_tmp.upper_bound(right);
				r_keep.insert(lb, ub);
			}
		}
		fprintf (stderr, "done\n");
		_l_hits.resize(l_keep.size());
		copy(l_keep.begin(), l_keep.end(), _l_hits.begin());
		_r_hits.resize(r_keep.size());
		copy(r_keep.begin(), r_keep.end(), _r_hits.begin());
		
	}
	
	const vector<size_t>& l_hits() const { return _l_hits; }
	const vector<size_t>& r_hits() const { return _r_hits; }
};

template<typename TStr, typename IntronFinder>
class JunctionFinder
{
public:
	JunctionFinder(const TStr& ref_str,
				   const IntronFinder& intron_finder,
				   uint32_t insert_length,
				   uint32_t insert_length_std_dev,
				   uint32_t min_intron_length,
				   uint32_t min_exon_length,
				   uint32_t bowtie_overlap_padding = 0, 
				   uint32_t max_gap = 0) :
		_ref_str(ref_str),
		_insert_len(insert_length),
		_std_dev(insert_length_std_dev),
		_min_intron_len(min_intron_length),
		_min_exon_len(min_exon_length),
		_l_motif_hits(intron_finder.l_hits()),
		_r_motif_hits(intron_finder.r_hits()),
		_padding(bowtie_overlap_padding),
		_max_gap(max_gap)
	{
		_left_memo.resize(_l_motif_hits.size());
		_right_memo.resize(_r_motif_hits.size());
	}
	
	void possible_junctions(std::set<pair<size_t, size_t> >& potential_splices,
							const vector<pair<size_t, size_t> >& happy_mates,
							const HitList& hits1,
							const HitList& hits2)
	{

		fprintf(stderr,"\n");
		for (size_t i = 0; i < happy_mates.size(); ++i)
		{
			const pair<size_t, size_t>& p = happy_mates[i];
			const BowtieHit& h1 = hits1[p.first];
			const BowtieHit& h2 = hits2[p.second];
			
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
			
			// Bowtie allows some mismatches, so we need to allow for the fact that
			// the donor/acceptor is actually "inside" the alignment for one or both
			// ends of the insert.  That is, intron may be longer than the pair's 
			// inner distance in the genomic coordinate space.  This variable
			// padding number accounts for this fact in a silly way that could be 
			// improved.
			
			int expected_inner_dist = _insert_len - h1.read_len() - h2.read_len();
			
			if (i % 1000 == 0)
			{
				fprintf(stderr, "\tProcessing mate %d of %d\n", (int)i, (int)happy_mates.size());
			}
			
			if (abs(inner_dist - expected_inner_dist) <= (int)_std_dev)
				continue;
			//fprintf(stderr, "\tProcessing mate %d\n", i);

			search(potential_splices, minor_hit_end, major_hit_start, expected_inner_dist);
		}
	}
	
private:

	vector<size_t>::const_iterator _right_motif_upper_bound;
	vector<size_t>::const_iterator _right_motif_lower_bound;
	
	vector<size_t>::const_iterator _left_motif_upper_bound;
	vector<size_t>::const_iterator _left_motif_lower_bound;
	//vector<size_t> _right_motif_upper_bound;
	
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
	
	void search(std::set<pair<size_t, size_t> >& splices,
				uint32_t start,
				uint32_t ref_target,
				int expected_inner_dist)
	{
		uint32_t left_start = max((int)start - (int)_padding, 0);
		uint32_t right_end = ref_target + _padding;
		
		
		_right_motif_lower_bound = lower_bound(_r_motif_hits.begin(), 
											   _r_motif_hits.end(), 
											   left_start);
		
		_right_motif_upper_bound = upper_bound(_r_motif_hits.begin(), 
											   _r_motif_hits.end(), 
											   right_end);
		
		_left_motif_lower_bound = lower_bound(_l_motif_hits.begin(), 
											  _l_motif_hits.end(), 
											  left_start);
		
		_left_motif_upper_bound = upper_bound(_l_motif_hits.begin(), 
											  _l_motif_hits.end(), 
											  right_end);
		
		
		//fprintf(stderr, "\t[%d,%d] : %d - ", _left_motif_upper_bound - _left_motif_lower_bound, _right_motif_upper_bound - _right_motif_lower_bound, right_end - left_start);
		_l_motif_calls = 0;
		_r_motif_calls = 0;
		
		/* TODO: consider some kind of long-term memo.  Closures found at the
		   top level of recursion here can be saved and potentially reused in 
		   future searches on overlapping intervals.
		 */
		search_hop_left_motif(splices, 
							  left_start, 
							  right_end, 
							  expected_inner_dist, 
							  0, 
							  true);
		
		/* Clear short-term memo */
		vector<vector<int> >::iterator _lm_begin = _left_memo.begin() + 
			(_left_motif_lower_bound - _l_motif_hits.begin());
		vector<vector<int> >::iterator _lm_end = _left_memo.begin() + 
			(_left_motif_upper_bound - _l_motif_hits.begin()) + 1;
		_lm_end = min(_lm_end, _left_memo.end());
		
		vector<vector<int> >::iterator _rm_begin = _right_memo.begin() + 
			(_right_motif_lower_bound - _r_motif_hits.begin());
		vector<vector<int> >::iterator _rm_end = _right_memo.begin() + 
			(_right_motif_upper_bound - _r_motif_hits.begin()) + 1;
		_rm_end = min(_rm_end, _right_memo.end());
		
		fill(_lm_begin, _lm_end, vector<int>());
		fill(_rm_begin, _rm_end, vector<int>());
		
#if !NDEBUG
		for (size_t i = 0; i < _left_memo.size(); ++i)
		{
			assert(_left_memo[i].size() == 0);
		}
		for (size_t i = 0; i < _right_memo.size(); ++i)
		{
			assert(_right_memo[i].size() == 0);
		}
#endif
		//fprintf(stderr, "(%d,%d)\n", _l_motif_calls, _r_motif_calls);
	}
	
//	void memo_left(size_t left_motif_pos, int path_distance)
//	{
//		_left_memo[li - _l_motif_hits.begin()].push_back(path_distance);
//	}
//	
//	void memo_right(size_t right_motif_pos, int path_distance)
//	{
//		_right_memo[right_motif_pos].push_back(path_distance);
//	}
	
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
	
	int search_hop_left_motif(std::set<pair<size_t, size_t> >& splices,
							   size_t start, 
							   uint32_t ref_target,
							   int expected_inner_dist,
							   uint32_t curr_path_length,
							   bool initial_hop = false)
	{
		int in_valid_path = -1;
		size_t left_bound = start + (initial_hop ? 0 : _min_exon_len);
		_l_motif_calls++;
		
		size_t right_bound = start + (expected_inner_dist - curr_path_length) + (2 * _padding) + _max_gap;
		
		vector<size_t>::const_iterator lb_left = lower_bound(_left_motif_lower_bound, 
															 _left_motif_upper_bound, 
															 left_bound);
		
		while (lb_left <= _left_motif_upper_bound && 
			   lb_left < _l_motif_hits.end() &&
			   *lb_left < start)
		{
			lb_left++;
		}
		
		vector<size_t>::const_iterator ub_left = lb_left;
		while (ub_left <= _left_motif_upper_bound && 
			   ub_left < _l_motif_hits.end() &&
			   *ub_left < right_bound)
		{
			ub_left++;
		}
		
		for (vector<size_t>::const_iterator li = lb_left; 
			 li <= ub_left && li != _l_motif_hits.end(); 
			 ++li)
		{
			assert (*li >= start || initial_hop);
			int dist = (*li - start);
			int next_path_len = curr_path_length + dist;
			if (next_path_len - expected_inner_dist >  (int)(_std_dev +  2 * _padding + _max_gap))
				return in_valid_path;
			int tolerance = _std_dev +  2 * _padding + _max_gap;
			
			size_t next = *li;
			//assert(next + _min_intron_len <= ref_target);
			if(next + _min_intron_len <= ref_target)
			{
				int ret = -1;
				ret = check_memo(_left_memo, 
								 li - _l_motif_hits.begin(), 
								 curr_path_length + dist, 
								 expected_inner_dist, 
								 tolerance);
				if (ret != -1)
				{
					in_valid_path = ret + dist;
					_left_memo[li - _l_motif_hits.begin()].push_back(in_valid_path);
				}
				else
				{
				
					ret = search_hop_right_motif(splices,
													 next,
													 ref_target,
													 expected_inner_dist,
													 next_path_len);
					if (ret != -1)
					{
						in_valid_path = ret + dist;
						_left_memo[li - _l_motif_hits.begin()].push_back(in_valid_path);
					}
				}
			}
			
		}
		
		return in_valid_path;
		
	}
	
	int search_hop_right_motif(std::set<pair<size_t, size_t> >& splices,
							    size_t start, 
							    uint32_t ref_target,
							    int expected_inner_dist,
							    uint32_t curr_path_length)
	{
		
		_r_motif_calls++;	
		int in_valid_path = -1;
		
		size_t left_start = start + _min_intron_len - 2;
		size_t right_start = ref_target;
		if (left_start > right_start)
			return in_valid_path;
		
		vector<size_t>::const_iterator lb_right = lower_bound(_right_motif_lower_bound, 
															  _right_motif_upper_bound, 
															  left_start);
		
		vector<size_t>::const_iterator ub_right = _right_motif_upper_bound;
		
		for (vector<size_t>::const_iterator ri = lb_right; 
			 ri <= ub_right && ri != _r_motif_hits.end(); 
			 ++ri)
		{
			size_t next = *ri + 2;
			int final_hop_dist = (int)ref_target - next;
			int final_path_len = curr_path_length + final_hop_dist;
			
			//fprintf (stderr, "%d, %d\n", start, next);
			int tolerance = _std_dev +  2 * _padding + _max_gap;
			
			if (abs(expected_inner_dist - final_path_len) <= tolerance)
			{
				pair<size_t, size_t> p = make_pair(start - 1, next);
				splices.insert(p);
				in_valid_path = final_hop_dist;
				_right_memo[ri - _r_motif_hits.begin()].push_back(final_hop_dist);
			}
			else 
			{
				if (next + _min_intron_len + _min_exon_len > ref_target)
					return in_valid_path;
				int ret = -1;
				ret = check_memo(_right_memo, 
								 ri - _r_motif_hits.begin(), 
								 curr_path_length, 
								 expected_inner_dist, 
								 tolerance);
				if (ret != -1)
				{
					pair<size_t, size_t> p = make_pair(start - 1, *ri + 2);
					splices.insert(p);
					in_valid_path = ret;
					_right_memo[ri - _r_motif_hits.begin()].push_back(ret);
				}
				else
				{
					
					 ret = search_hop_left_motif(splices,
													next,
													ref_target,
													expected_inner_dist,
													curr_path_length);
					if (ret != -1)
					{
						pair<size_t, size_t> p = make_pair(start - 1, *ri + 2);
						splices.insert(p);
						in_valid_path = ret;
						_right_memo[ri - _r_motif_hits.begin()].push_back(in_valid_path);
					}
				}
			}
			
		}
		
		return in_valid_path;
	}
	
	TStr _ref_str;
	uint32_t _insert_len;
	uint32_t _std_dev;
	uint32_t _min_intron_len;
	uint32_t _min_exon_len;
	
	const vector<size_t>& _l_motif_hits;
	const vector<size_t>& _r_motif_hits;
	
	uint32_t _padding;
	uint32_t _max_gap;
	uint32_t _r_motif_calls;
	uint32_t _l_motif_calls;
	
	vector<vector<int> > _left_memo;
	vector<vector<int> > _right_memo;
};

struct Junction
{
	
	Junction (uint32_t ref, uint32_t l, uint32_t r, bool a)
	: refid(ref), left(l), right(r), antisense(a) {}
	Junction() {}
	uint32_t refid;
	uint32_t left;
	uint32_t right;
	bool antisense;
	
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
	
#if !NDEBUG
	bool valid() const
	{
		return refid != 0xFFFFFFFF && left < right && (left != right);
	}
#endif
};

struct JunctionStats
{
	uint8_t left_extent;
	uint8_t right_extent;
	uint16_t num_reads : 15;
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

void junction_from_alignment(const BowtieHit& spliced_alignment,
							 uint32_t refid,
							 JunctionSet& junctions);

void junctions_from_alignments(const HitTable& hits,
							   JunctionSet& junctions);

void accept_valid_junctions(JunctionSet& junctions,
							const uint32_t refid,
							const vector<short>& DoC,
							float min_isoform_fraction);

void accept_all_junctions(JunctionSet& junctions,
							const uint32_t refid);

void print_junctions(FILE* junctions_out, 
					 const JunctionSet& junctions,
					 SequenceTable& ref_sequences);

#endif

