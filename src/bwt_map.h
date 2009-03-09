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
#include <seqan/sequence.h>

using namespace std;

/*
 *  bwt_map.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */


enum CigarOpCode { MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD };

struct CigarOp
{
	CigarOp(CigarOpCode o, uint32_t l) : opcode(o), length(l) {}
	CigarOpCode opcode : 3;
	uint32_t length : 29;
	
	bool operator==(const CigarOp& rhs) const { return opcode == rhs.opcode && length == rhs.length; }
	
};

typedef uint32_t ReadID;

/*  Stores the information from a single record of the bowtie map. A given read
    may have many of these.  Reads up to 255bp are supported. 
*/
struct BowtieHit
{
	BowtieHit() : 
		_ref_id(0),
		_insert_id(0),
		_accepted(false) {}
	
	BowtieHit(uint32_t ref_id,
			  ReadID insert_id, 
			  int left, 
			  int read_len, 
			  bool antisense) :
		_ref_id(ref_id),
		_insert_id(insert_id), 
		_left(left), 
		_cigar(vector<CigarOp>(1,CigarOp(MATCH,read_len))),
		_antisense_splice(false),
		_antisense_aln(antisense),
		_accepted(false)
	{
		assert(_cigar.capacity() == _cigar.size());
	}
	
	BowtieHit(uint32_t ref_id,
			  ReadID insert_id, 
			  int left,  
			  const vector<CigarOp>& cigar,
			  bool antisense_aln,
			  bool antisense_splice) : 
		_ref_id(ref_id),
		_insert_id(insert_id), 	
		_left(left),
		_cigar(cigar),
		_antisense_splice(antisense_splice),
		_antisense_aln(antisense_aln),
		_accepted(false)
	{
		assert(_cigar.capacity() == _cigar.size());
	}
	
	uint8_t read_len() const
	{
		uint32_t len = 0;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			switch(op.opcode)
			{
				case MATCH:
				case INS:
				case SOFT_CLIP:
					len += op.length;
					break;
				default:
					break;
			}
		}
		
		return len;
	}
	
	bool operator==(const BowtieHit& rhs) const
	{
	    return (_insert_id == rhs._insert_id &&
	            _ref_id == rhs._ref_id &&
	            _antisense_aln == rhs._antisense_aln &&
	            _left == rhs._left && 
	            _antisense_splice == rhs._antisense_splice &&
	            /* DO NOT USE ACCEPTED IN COMPARISON */
	            _cigar == rhs._cigar);
    }
	
	uint32_t ref_id() const				{ return _ref_id;			}
	ReadID insert_id() const			{ return _insert_id;		}
	int left() const				{ return _left;				}
	int right() const	
	{
		int r = _left;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			
			switch(op.opcode)
			{
				case MATCH:
				case REF_SKIP:
				case DEL:
					r += op.length;
					break;
				default:
					break;
			}
		}
		return r;			
	}
	bool antisense_splice()	const	{ return _antisense_splice; }
	bool antisense_align() const	{ return _antisense_aln;	}
	bool accepted()	const			{ return _accepted;			}
	
	void accepted(bool accept)		{ _accepted = accept; }
	
	// For convenience, if you just want a copy of the gap intervals
	// for this hit.
	void gaps(vector<pair<int,int> >& gaps_out) const
	{
		gaps_out.clear();
		int pos = _left;
		for (size_t i = 0; i < _cigar.size(); ++i)
		{
			const CigarOp& op = _cigar[i];
			
			switch(op.opcode)
			{
				case REF_SKIP:
					gaps_out.push_back(make_pair(pos, pos + op.length - 1));
					pos += op.length;
					break;
				case MATCH:
					pos += op.length;
					break;
				default:
					break;
			}
		}
	}
	
	const vector<CigarOp>& cigar() const { return _cigar; }
	
	bool contiguous() const 
	{ 
		return _cigar.size() == 1 && _cigar[0].opcode == MATCH;
	}
	
private:
	
	uint32_t _ref_id;
	ReadID _insert_id;   // Id of the sequencing insert
	int _left;        // Position in the reference of the left side of the alignment
	
	vector<CigarOp> _cigar;
	
	bool _antisense_splice;    // Whether the junction spanned is on the reverse strand
	bool _antisense_aln;       // Whether the alignment is to the reverse strand
	bool _accepted;
};

class ReadTable
{
public:
	
	typedef seqan::Dna5String Sequence;
	
	ReadTable() : 
		_next_id(1) {}
	
	// This function should NEVER return zero
	ReadID get_id(const string& name)
	{
		uint32_t _id = atoi(name.c_str());
		assert(_id);
		_next_id = max(_next_id, (size_t)_id);
		return _id;
	}
	
	
	uint32_t observation_order(ReadID ID)
	{
		if (ID == 0)
			return 0xFFFFFFFF;
		return ID;
	}
	
	
	size_t size() const { return _next_id - 1; }
	
private:
	size_t _next_id;
};

class RefSequenceTable
{
public:
	
	typedef seqan::Dna5String Sequence;
	
	struct SequenceInfo
	{
		SequenceInfo(uint32_t _order, 
					 char* _name, 
					 Sequence* _seq) :
		observation_order(_order),
		name(_name),
		seq(_seq) {}
		uint32_t observation_order;
		char* name;
		Sequence* seq;
	};
	
	typedef map<string, uint64_t> IDTable;
	typedef map<uint32_t, SequenceInfo> InvertedIDTable;
	typedef InvertedIDTable::iterator iterator;
	typedef InvertedIDTable::const_iterator const_iterator;
	
	RefSequenceTable(bool keep_names, bool keep_seqs = false) : 
	_next_id(1), 
	_keep_names(keep_names) {}
	
	// This function should NEVER return zero
	uint32_t get_id(const string& name,
					Sequence* seq)
	{
		uint32_t _id = hash_string(name.c_str());
		pair<InvertedIDTable::iterator, bool> ret = 
		_by_id.insert(make_pair(_id, SequenceInfo(_next_id, NULL, NULL)));
		if (ret.second == true)
		{			
			char* _name = NULL;
			if (_keep_names)
				_name = strdup(name.c_str());
			ret.first->second.name = _name;
			ret.first->second.seq	= seq;
			++_next_id;
		}
		assert (_id);
		return _id;
	}
	
	// You must call invert() before using this function
	const char* get_name(uint32_t ID)
	{
		InvertedIDTable::iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.name;
		else
			return NULL;
	}

	Sequence* get_seq(uint32_t ID) 
	{
		InvertedIDTable::iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.seq;
		else
			return NULL;
	}

	const SequenceInfo* get_info(uint32_t ID)
	{
		
		InvertedIDTable::iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
		{
			return &(itr->second);
		}
		else
			return NULL;
	}
	
	int observation_order(uint32_t ID)
	{
		InvertedIDTable::iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
		{
			return itr->second.observation_order;
		}
		else
			return -1;
	}
	
	iterator begin() { return _by_id.begin(); }
	iterator end() { return _by_id.end(); }
	
	const_iterator begin() const { return _by_id.begin(); }
	const_iterator end() const { return _by_id.end(); }
	
	size_t size() const { return _by_id.size(); }
	
	void clear()
	{
		//_by_name.clear();
		_by_id.clear();
	}
	
private:
	
	// This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
	inline uint32_t hash_string(const char* __s)
	{
		uint32_t hash = 0x811c9dc5;
		for ( ; *__s; ++__s)
		{
			hash *= 16777619;
			hash ^= *__s;
		}
		return hash;
	}
	
	//IDTable _by_name;
	uint32_t _next_id;
	bool _keep_names;
	InvertedIDTable _by_id;
};

bool hit_insert_id_lt(const BowtieHit& h1, const BowtieHit& h2);

typedef vector<BowtieHit> HitList;

/* This class stores all the hits from a Bowtie map */

// TODO: HitTable also acts like a factory for BowtieHits.  This is a poor
// design, and should be refactored, decoupling its container and factory
// roles.
class HitTable
{
public:

	typedef map<uint64_t, HitList> RefHits;
	typedef RefHits::const_iterator const_iterator;
	typedef RefHits::iterator iterator;
	
	HitTable() :  _total_hits(0) {}
	
	const_iterator begin() const { return _hits_for_ref.begin(); }
	const_iterator end() const { return _hits_for_ref.end(); }
	
	iterator begin() { return _hits_for_ref.begin(); }
	iterator end() { return _hits_for_ref.end(); }
	
	void add_hit(const BowtieHit& bh, bool check_uniqueness);
	
	void finalize()
	{
		for (RefHits::iterator i = _hits_for_ref.begin();
			 i != _hits_for_ref.end();
			 ++i)
		{
			sort(i->second.begin(), i->second.end(), hit_insert_id_lt);
		}
	}
	
	HitList* get_hits(uint64_t ref_id)
	{
		RefHits::iterator i = _hits_for_ref.find(ref_id);
		if (i == _hits_for_ref.end())
			return NULL;
		else
			return &(i->second);
	}
	
    uint32_t total_hits() const { return _total_hits; }
	
private:
	RefHits _hits_for_ref;
    uint32_t _total_hits;
};


class HitFactory
{
public:
	HitFactory(ReadTable& insert_table, 
			   RefSequenceTable& reference_table) : 
	_insert_table(insert_table), _ref_table(reference_table) {}
	
	virtual ~HitFactory() {}
	
	BowtieHit create_hit(const string& insert_name, 
						 const string& ref_name,
						 int left,
						 const vector<CigarOp>& cigar,
						 bool antisense_aln,
						 bool antisense_splice);
	
	BowtieHit create_hit(const string& insert_name, 
						 const string& ref_name,
						 uint32_t left,
						 uint32_t read_len,
						 bool antisense_aln);
	
	virtual bool get_hit_from_buf(const char* bwt_buf, 
								  BowtieHit& bh,
								  bool strip_slash,
								  char* name_out = NULL,
								  char* name_tags = NULL) = 0;
	
private:
	ReadTable& _insert_table;
	RefSequenceTable& _ref_table;
};

class BowtieHitFactory : public HitFactory
{
public:
	BowtieHitFactory(ReadTable& insert_table, 
					 RefSequenceTable& reference_table) : 
		HitFactory(insert_table, reference_table) {}
	//virtual ~BowtieHitFactory() {}
	
	bool get_hit_from_buf(const char* bwt_buf, 
						  BowtieHit& bh,
						  bool strip_slash,
						  char* name_out = NULL,
						  char* name_tags = NULL);
};

class SplicedBowtieHitFactory : public HitFactory
{
public:
	SplicedBowtieHitFactory(ReadTable& insert_table, 
							RefSequenceTable& reference_table,
							bool filter_anchor_mismatches) : 
	HitFactory(insert_table, reference_table),
	_filter_anchor_mismatches(filter_anchor_mismatches) {}
	
	//virtual ~SplicedBowtieHitFactory() {}
	
	bool get_hit_from_buf(const char* bwt_buf, 
						  BowtieHit& bh,
						  bool strip_slash,
						  char* name_out = NULL,
						  char* name_tags = NULL);
private:
	bool _filter_anchor_mismatches;
};

class SAMHitFactory : public HitFactory
{
public:
	SAMHitFactory(ReadTable& insert_table, 
				  RefSequenceTable& reference_table) : 
	HitFactory(insert_table, reference_table) {}
	
	bool get_hit_from_buf(const char* bwt_buf, 
						  BowtieHit& bh,
						  bool strip_slash,
						  char* name_out = NULL,
						  char* name_tags = NULL);
};

struct HitsForRead
{
	HitsForRead() : insert_id(0) {}
	uint64_t insert_id;
	vector<BowtieHit> hits;
	
//	HitsForRead& operator=(const HitsForRead& rhs)
//	{
//		insert_id = rhs.insert_id;
//		hits = rhs.hits;
//		return *this;
//	}
};

class HitStream
{
	
public:
	HitStream(FILE* hit_file, HitFactory* hit_factory, bool spliced, bool strip_slash) : 
	_hit_file(hit_file),
	_factory(hit_factory),
	_spliced(spliced),
	_strip_slash(strip_slash),
	buffered_hit(BowtieHit())
	{
		// Prime the stream by reading a single hit into the buffered_hit
		HitsForRead dummy = HitsForRead();
		next_read_hits(dummy);
	}
	
	bool next_read_hits(HitsForRead& hits_for_read)
	{
		hits_for_read.hits.clear();
		hits_for_read.insert_id = 0; 
		
		if (feof(_hit_file))
		{
			buffered_hit = BowtieHit();
			return false;
		}
		char bwt_buf[2048];
		
		hits_for_read.insert_id = buffered_hit.insert_id();
		hits_for_read.hits.push_back(buffered_hit);
		
		while (fgets(bwt_buf, 2048, _hit_file))
		{
			// Chomp the newline
			char* nl = strrchr(bwt_buf, '\n');
			if (nl) *nl = 0;
			string clean_buf = bwt_buf;
			// Get a new record from the tab-delimited Bowtie map
			BowtieHit bh;
			
			if (!_factory->get_hit_from_buf(bwt_buf, bh, _strip_slash))
				continue;
			if (bh.insert_id() == hits_for_read.insert_id)
			{
				hits_for_read.hits.push_back(bh);
			}
			else
			{
				buffered_hit = bh;
				break;
			}
		}	
		if (feof(_hit_file))
			buffered_hit = BowtieHit();
		
		return (!hits_for_read.hits.empty() && hits_for_read.insert_id != 0);
	}
	
	uint64_t next_group_id() const 
	{
		return buffered_hit.insert_id();
	}
	
private:
	FILE* _hit_file;
	HitFactory* _factory;
	bool _spliced;
	bool _strip_slash;
	BowtieHit buffered_hit;
};

typedef uint32_t MateStatusMask;

void get_mapped_reads(FILE* bwtf, 
					  HitTable& hits,
					  HitFactory& hit_factory,
					  bool strip_slash,
					  bool verbose = false);

//bool left_status_better(MateStatusMask left, MateStatusMask right);
//bool status_equivalent(MateStatusMask left, MateStatusMask right);
typedef uint32_t MateStatusMask;

enum AlignStatus {UNALIGNED, SPLICED, CONTIGUOUS};
AlignStatus status(const BowtieHit* align);

void add_hits_to_coverage(const HitList& hits, vector<short>& DoC);

void accept_all_hits(HitTable& hits);

void print_hit(FILE* fout, 
			   const char* read_name,
			   const BowtieHit& bh,
			   const char* ref_name,
			   const char* sequence,
			   const char* qualities);

#endif
