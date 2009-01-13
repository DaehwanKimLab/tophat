#ifndef BWT_MAP_H
#define BWT_MAP_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>

using namespace std;

/*
 *  bwt_map.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

/*  Stores the information from a single record of the bowtie map. A given read
    may have many of these.  Reads up to 255bp are supported.
 
	TODO: this struct takes up 20 bytes, but it could be as small as 16 and 
	still be word aligned.  
*/
struct BowtieHit
{
	BowtieHit(uint32_t _insert_id, uint32_t _left, uint32_t read_len, bool _antisense) : 
		insert_id(_insert_id), 
		left(_left), 
		right(_left + read_len),
		splice_pos_left(-1),
		splice_pos_right(-1),
		antisense_aln(_antisense) {}
	
	BowtieHit(uint32_t _insert_id, 
			  uint32_t _left, 
			  uint32_t _right, 
			  uint32_t _sp_left,
			  uint32_t _sp_right,
			  uint32_t read_len, 
			  bool _antisense_aln,
			  bool _antisense_splice) : 
		insert_id(_insert_id), 
		left(_left), 
		right(_right),
		splice_pos_left(_sp_left),
		splice_pos_right(_sp_right),
		antisense_aln(_antisense_aln),
		antisense_splice(_antisense_splice), 
		accepted(false) {}
	
	uint8_t read_len() const
	{
		if (splice_pos_left == -1 || splice_pos_right == -1)
		{
			return right - left;
		}
		else
		{
			return splice_pos_left + splice_pos_right;
		}
	}
	
	uint32_t insert_id;   // Id of the sequencing insert
	uint32_t left;        // Position in the reference of the left side of the alignment
	uint32_t right;       // Position in the reference of the right side of the alignment
	int8_t splice_pos_left; // Offset from left where the splice begins, or -1 for unspliced alignments (ADD to left)
	int8_t splice_pos_right;// Offset from right where the splice begins, or -1 for unspliced alignments (SUBTRACT from right)
	bool antisense_aln : 1;       // Whether the alignment is to the reverse strand
	bool antisense_splice : 1;    // Whether the junction spanned is on the reverse strand
	bool accepted : 1;
	uint16_t packing : 13;
};

class SequenceTable
{
public:
	
	typedef map<string, uint32_t> TableType;
	typedef map<uint32_t, string> InvertedTableType;
	typedef TableType::iterator iterator;
	typedef TableType::const_iterator const_iterator;
	
	SequenceTable() : _next_id(0) {}
	
	uint32_t get_id(const string& seq_name)
	{
		pair<TableType::iterator, bool> ret = 
		_sequences_by_name.insert(make_pair(seq_name, _next_id));
		if (ret.second == true)
		{
			++_next_id;
		}
		return ret.first->second;
	}
	
	// You must call invert() before using this function
	const string* get_name(uint32_t ID) const
	{
		InvertedTableType::const_iterator i = _sequences_by_id.find(ID);
		if (i != _sequences_by_id.end())
			return &(i->second);
		return NULL;
	}
	
	iterator begin() { return _sequences_by_name.begin(); }
	iterator end() { return _sequences_by_name.end(); }
	
	const_iterator begin() const { return _sequences_by_name.begin(); }
	const_iterator end() const { return _sequences_by_name.end(); }
	
	size_t size() const { return _sequences_by_name.size(); }
	
	void invert()
	{
		for (TableType::const_iterator i = _sequences_by_name.begin();
			 i != _sequences_by_name.end();
			 ++i)
		{
			_sequences_by_id[i->second] = i->first;
		}
	}
	
private:
	
	TableType _sequences_by_name;
	uint32_t _next_id;
	
	InvertedTableType _sequences_by_id;
};

bool hit_insert_id_lt(const BowtieHit& h1, const BowtieHit& h2);

typedef vector<BowtieHit> HitList;

/* This class stores all the hits from a Bowtie map */
class HitTable
{

	
public:

	typedef map<uint32_t, HitList> RefHits;
	typedef RefHits::const_iterator const_iterator;
	typedef RefHits::iterator iterator;
	
	HitTable(SequenceTable& insert_table, 
			 SequenceTable& reference_table) : 
	_insert_table(insert_table), _ref_table(reference_table) {}
	
	const_iterator begin() const { return _hits_for_ref.begin(); }
	const_iterator end() const { return _hits_for_ref.end(); }
	
	iterator begin() { return _hits_for_ref.begin(); }
	iterator end() { return _hits_for_ref.end(); }
	
	
	void add_spliced_hit(const string& insert_name, 
						 const string& ref_name,
						 uint32_t left,
						 uint32_t right,
						 char splice_pos_left,
						 char splice_pos_right,
						 uint32_t read_len,
						 bool antisense_aln,
						 bool antisense_splice)
	{
		uint32_t insert_id = _insert_table.get_id(insert_name);
		uint32_t reference_id = _ref_table.get_id(ref_name);
		
		pair<RefHits::iterator, bool> ret = 
		_hits_for_ref.insert(make_pair(reference_id, HitList()));
		
		BowtieHit bh = BowtieHit(insert_id, 
								 left, 
								 right, 
								 splice_pos_left, 
								 splice_pos_right, 
								 read_len, 
								 antisense_aln,
								 antisense_splice);
		
		(*(ret.first)).second.push_back(bh);
	}
	
	void add_hit(const string& insert_name, 
				 const string& ref_name,
				 uint32_t left,
				 uint32_t read_len,
				 bool antisense)
	{
		uint32_t insert_id = _insert_table.get_id(insert_name);
		uint32_t reference_id = _ref_table.get_id(ref_name);
		
		pair<RefHits::iterator, bool> ret = 
		_hits_for_ref.insert(make_pair(reference_id, HitList()));
		
		(*(ret.first)).second.push_back(BowtieHit(insert_id, left, read_len, antisense));
	}
	
	void finalize()
	{
		for (RefHits::iterator i = _hits_for_ref.begin();
			 i != _hits_for_ref.end();
			 ++i)
		{
			sort(i->second.begin(), i->second.end(), hit_insert_id_lt);
		}
	}
	
	HitList* get_hits(uint32_t ref_id)
	{
		RefHits::iterator i = _hits_for_ref.find(ref_id);
		if (i == _hits_for_ref.end())
			return NULL;
		else
			return &(i->second);
	}
	
	
	
private:
	SequenceTable& _insert_table;
	SequenceTable& _ref_table;
	RefHits _hits_for_ref;
};

typedef uint32_t MateStatusMask;

void get_mapped_reads(FILE* bwtf, HitTable& hits, bool spliced, bool verbose = false);
pair<int, int> pair_distances(const BowtieHit& h1, const BowtieHit& h2);

bool left_status_better(MateStatusMask left, MateStatusMask right);
bool status_equivalent(MateStatusMask left, MateStatusMask right);
typedef uint32_t MateStatusMask;

enum AlignStatus {UNALIGNED, SPLICED, CONTIGUOUS};

struct FragmentAlignment
{
	FragmentAlignment(uint32_t _refid, 
					BowtieHit* _alignment) : 
	refid(_refid), 
	alignment(_alignment) {}
	
	uint32_t refid;
	BowtieHit* alignment;
};

struct FragmentAlignmentGrade
{
	FragmentAlignmentGrade() : 
		status(UNALIGNED) {}
	
	FragmentAlignmentGrade(const BowtieHit& h1) 
	{
		if (h1.splice_pos_left != -1)
		{
			status = SPLICED;
		}
		else
		{
			status = CONTIGUOUS;
		}
	}
	
	FragmentAlignmentGrade& operator=(const FragmentAlignmentGrade& rhs)
	{
		status = rhs.status;
		return *this;
	}
	
	// Returns true if rhs is a "happier" alignment for the ends of this insert
	// than this InsertStatus.
	bool operator<(const FragmentAlignmentGrade& rhs)
	{
		return status < rhs.status;
	}
	
	uint8_t status;
};

struct InsertAlignment
{
	InsertAlignment(uint32_t _refid, 
					BowtieHit* _left_alignment, 
					BowtieHit* _right_alignment) : 
		refid(_refid), 
		left_alignment(_left_alignment),
		right_alignment(_right_alignment) {}
	
	uint32_t refid;
	BowtieHit* left_alignment;
	BowtieHit* right_alignment;
};

AlignStatus status(const BowtieHit* align);

struct InsertAlignmentGrade
{
	InsertAlignmentGrade() : 
		too_close(false),
		too_far(false),
		one_spliced(false),
		both_spliced(false),
		one_mapped(false),
		both_mapped(false),
		opposite_strands(false),
		consistent_splices(false) {}
	
	InsertAlignmentGrade(const BowtieHit& h1) : 
		too_close(false),
		too_far(false),
		one_spliced(false),
		both_spliced(false),
		one_mapped(false),
		both_mapped(false),
		opposite_strands(false),
		consistent_splices(false)
	{
		if (h1.splice_pos_left != -1)
			one_spliced = true;
		one_mapped = true;
	}
	
	InsertAlignmentGrade(const BowtieHit& h1, 
				 const BowtieHit& h2, 
				 int max_inner_distance, 
				 int min_inner_distance)
	{
		pair<int, int> distances = pair_distances(h1,h2);
		int inner_dist = distances.second;
		
		both_mapped = true;
		
		one_spliced = ((h1.splice_pos_left == -1 && h2.splice_pos_right != -1) ||
					   (h1.splice_pos_left != -1 && h2.splice_pos_right == -1));
		
		both_spliced = (h1.splice_pos_left != -1 && h2.splice_pos_right != -1);
		

		too_far = (inner_dist > max_inner_distance);

		too_close = (inner_dist < min_inner_distance);
		
		opposite_strands = (h1.antisense_aln != h2.antisense_aln);
		
		consistent_splices = (h1.splice_pos_left != -1 && h2.splice_pos_left != -1 &&
							  h1.antisense_splice != h2.antisense_splice);
		assert(!(too_far && too_close));
	}
	
	InsertAlignmentGrade& operator=(const InsertAlignmentGrade& rhs)
	{
		too_close = rhs.too_close;
		too_far = rhs.too_far;
		one_spliced = rhs.one_spliced;
		both_spliced = rhs.both_spliced;
		one_mapped = rhs.one_mapped;
		both_mapped = rhs.both_mapped;
		opposite_strands = rhs.opposite_strands;
		consistent_splices = rhs.consistent_splices;
		return *this;
	}
	
	// Returns true if rhs is a "happier" alignment for the ends of this insert
	// than this InsertStatus.
	bool operator<(const InsertAlignmentGrade& rhs)
	{
		assert(!(both_mapped && one_mapped));
		// We always prefer a insert alignment with both ends mapped than a
		// singleton
		if (both_mapped != rhs.both_mapped)
		{
			return both_mapped < rhs.both_mapped;
		}
		else
		{
			
			// Prefer an alignment with one end contiguously mapped over a pair 
			// of spliced alignments
			if (both_spliced != rhs.both_spliced)
				return rhs.both_spliced < both_spliced;
			
			// Prefer a pair of contiguously aligned ends to a single 
			// contiguously aligned end and a spliced end
			if (one_spliced != rhs.one_spliced)
				return rhs.one_spliced < one_spliced;
			
			// Prefer a pair that is too close or perfect to one that is too far
			if (too_far && !rhs.too_far)
				return true;
			
			// Prefer a pair that is perfect to one that is too close
			if (too_close && !(rhs.too_close || rhs.too_far))
				return true;
			
			return false;
		}
		
		// We prefer a singleton mapping to an insert with neither end mapped
		if (one_mapped != rhs.one_mapped)
		{
			return one_mapped < rhs.one_mapped;
		
		}
		else
		{
			return rhs.one_spliced < one_spliced;
		}
	}
	
	bool too_close;
	bool too_far;
	
	bool one_spliced;
	bool both_spliced;
	
	bool one_mapped;
	bool both_mapped;
	
	bool opposite_strands;
	bool consistent_splices;
};		

typedef vector<pair<InsertAlignmentGrade, vector<InsertAlignment> > > BestInsertAlignmentTable;
//vector<InsertAlignment> best_alignments;

typedef vector<pair<FragmentAlignmentGrade, vector<FragmentAlignment> > > BestFragmentAlignmentTable;


void best_insert_mappings(uint32_t refid,
						  const string& name,
						  HitList& hits1_in_ref,
						  HitList& hits2_in_ref,
						  BestInsertAlignmentTable& best_insert_alignments);

void best_fragment_mappings(uint32_t refid,
							const string& name,
							HitList& hits1_in_ref,
							BestFragmentAlignmentTable& best_status_for_fragments);

//typedef vector<InsertStatus> InsertStatusTable;

void add_hits_to_coverage(const HitList& hits, vector<short>& DoC);

void accept_all_hits(HitTable& hits);

void accept_valid_hits(BestInsertAlignmentTable& best_status_for_inserts);
void accept_valid_hits(BestFragmentAlignmentTable& best_status_for_fragments);
void accept_unique_hits(BestFragmentAlignmentTable& best_status_for_fragments);
#endif
