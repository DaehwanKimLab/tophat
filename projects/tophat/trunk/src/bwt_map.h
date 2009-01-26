#ifndef BWT_MAP_H
#define BWT_MAP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
	BowtieHit(uint16_t _ref_id,
			  uint32_t _insert_id, 
			  uint32_t _left, 
			  uint32_t read_len, 
			  bool _antisense) :
		ref_id(_ref_id),
		antisense_aln(_antisense),
		insert_id(_insert_id), 
		left(_left), 
		antisense_splice(false),
		right(_left + read_len),
		accepted(false),
		splice_pos_left(-1),
		splice_pos_right(-1) {}
	
	BowtieHit(uint16_t _ref_id,
			  uint32_t _insert_id, 
			  uint32_t _left, 
			  uint32_t _right, 
			  uint32_t _sp_left,
			  uint32_t _sp_right,
			  uint32_t read_len, 
			  bool _antisense_aln,
			  bool _antisense_splice) : 
		ref_id(_ref_id),
		antisense_aln(_antisense_aln),	
		insert_id(_insert_id), 	
		left(_left), 
		antisense_splice(_antisense_splice),
		right(_right),
		accepted(false),
		splice_pos_left(_sp_left),
		splice_pos_right(_sp_right) {}
	
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
	
	bool operator==(const BowtieHit& rhs)
	{
	    return (insert_id == rhs.insert_id &&
	            ref_id == rhs.ref_id &&
	            antisense_aln == rhs.antisense_aln &&
	            left == rhs.left && 
	            right == rhs.right &&
	            antisense_splice == rhs.antisense_splice &&
	            /* DO NOT USE ACCEPTED IN COMPARISON */
	            splice_pos_left == rhs.splice_pos_left &&
                splice_pos_right == rhs.splice_pos_right);
    }
	
	uint16_t ref_id : 15;
	bool antisense_aln : 1;       // Whether the alignment is to the reverse strand
	uint32_t insert_id;   // Id of the sequencing insert
	
	uint32_t left : 31;        // Position in the reference of the left side of the alignment
	bool antisense_splice : 1;    // Whether the junction spanned is on the reverse strand
	uint32_t right : 31;       // Position in the reference of the right side of the alignment
	bool accepted : 1;
	int8_t splice_pos_left; // Offset from left where the splice begins, or -1 for unspliced alignments (ADD to left)
	int8_t splice_pos_right;// Offset from right where the splice begins, or -1 for unspliced alignments (SUBTRACT from right)
	
	//uint16_t packing : 13;
} __attribute__((packed));

class SequenceTable
{
public:
	
	typedef map<string, uint32_t> TableType;
	typedef map<uint32_t, string> InvertedTableType;
	typedef TableType::iterator iterator;
	typedef TableType::const_iterator const_iterator;
	
	SequenceTable(bool keep_names) : _next_id(0), _keep_names(keep_names) {}
	
	uint32_t get_id(const string& seq_name)
	{
		if (_keep_names)
		{
			pair<TableType::iterator, bool> ret = 
			_sequences_by_name.insert(make_pair(seq_name, _next_id));
			if (ret.second == true)
			{
				++_next_id;
			}
			return ret.first->second;
		}
		else
		{
			return hash_string(seq_name.c_str());
		}
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
	
	void clear()
	{
		_sequences_by_name.clear();
		_sequences_by_id.clear();
	}
	
private:
	
	inline uint32_t hash_string(const char* __s)
	{
		unsigned long __h = 0;
		for ( ; *__s; ++__s)
			__h = __h*5 + *__s;
		return uint32_t(__h);
	}
	
	TableType _sequences_by_name;
	uint32_t _next_id;
	bool _keep_names;
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
	_insert_table(insert_table), _ref_table(reference_table), _total_hits(0) {}
	
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
                         bool antisense_splice);
	
	void add_hit(const string& insert_name, 
				 const string& ref_name,
				 uint32_t left,
				 uint32_t read_len,
                 bool antisense);
	
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
	
    uint32_t total_hits() const { return _total_hits; }
	
private:
	SequenceTable& _insert_table;
	SequenceTable& _ref_table;
	RefHits _hits_for_ref;
    uint32_t _total_hits;
};

typedef uint32_t MateStatusMask;

void get_mapped_reads(FILE* bwtf, 
					  HitTable& hits, 
					  bool spliced, 
					  bool verbose = false);


//bool left_status_better(MateStatusMask left, MateStatusMask right);
//bool status_equivalent(MateStatusMask left, MateStatusMask right);
typedef uint32_t MateStatusMask;

enum AlignStatus {UNALIGNED, SPLICED, CONTIGUOUS};
AlignStatus status(const BowtieHit* align);

void add_hits_to_coverage(const HitList& hits, vector<short>& DoC);

void accept_all_hits(HitTable& hits);

#endif
