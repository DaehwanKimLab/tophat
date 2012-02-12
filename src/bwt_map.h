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
#include <cassert>
#include <cstring>
#include <algorithm>
#include <seqan/sequence.h>

#include <bam/sam.h>
using namespace std;

#include "common.h"
#include "reads.h"

#define _FBUF_SIZE 10*1024

/*
 *  bwt_map.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

enum CigarOpCode
  {
    CIGAR_NOTHING = 0,
    FUSION_NOTHING = 0,
    MATCH,
    mATCH,
    INS,
    iNS,
    DEL,
    dEL,
    FUSION_FF,
    FUSION_FR,
    FUSION_RF,
    FUSION_RR,
    REF_SKIP,
    rEF_SKIP,
    SOFT_CLIP,
    HARD_CLIP,
    PAD
  };

struct CigarOp
{
	CigarOp(CigarOpCode o, uint32_t l) : opcode(o), length(l) {}
	CigarOpCode opcode;
	uint32_t length;
	
	bool operator==(const CigarOp& rhs) const 
	{ 
		return opcode == rhs.opcode && length == rhs.length; 
	}
	
};

typedef uint32_t ReadID;
typedef uint32_t RefID;

class RefSequenceTable;

/*  Stores the information from a single record of the bowtie map. A given read
    may have many of these.  Reads up to 255bp are supported. 
*/
struct BowtieHit
{
BowtieHit() : 
  _ref_id(0),
    _ref_id2(0),
    _insert_id(0),
    _left(0),
    _antisense_splice(false),
    _antisense_aln(false),
    _edit_dist(0),
    _splice_mms(0),
    _alignment_score(0),
    _end(false){}
  
BowtieHit(uint32_t ref_id,
	  uint32_t ref_id2,
	  ReadID insert_id, 
	  int left, 
	  int read_len, 
	  bool antisense,
	  unsigned char edit_dist,
	  bool end) :
  _ref_id(ref_id),
    _ref_id2(ref_id2),
    _insert_id(insert_id), 
    _left(left), 
    _cigar(vector<CigarOp>(1,CigarOp(MATCH,read_len))),
    _antisense_splice(false),
    _antisense_aln(antisense),
    _edit_dist(edit_dist),
    _splice_mms(0),
    _alignment_score(0),
    _end(end)
  {
    assert(_cigar.capacity() == _cigar.size());
  }
  
BowtieHit(uint32_t ref_id,
	  uint32_t ref_id2,
	  ReadID insert_id, 
	  int left,  
	  const vector<CigarOp>& cigar,
	  bool antisense_aln,
	  bool antisense_splice,
	  unsigned char edit_dist,
	  unsigned char splice_mms,
	  bool end) : 
  _ref_id(ref_id),
    _ref_id2(ref_id2),
    _insert_id(insert_id), 	
    _left(left),
    _cigar(cigar),
    _antisense_splice(antisense_splice),
    _antisense_aln(antisense_aln),
    _edit_dist(edit_dist),
    _splice_mms(splice_mms),
    _alignment_score(0),
    _end(end)
  {
    assert(_cigar.capacity() == _cigar.size());
  }

  int read_len() const
  {
    uint32_t len = 0;
    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	switch(op.opcode)
	  {
	  case MATCH:
	  case mATCH:
	  case INS:
	  case iNS:
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
	    _ref_id2 == rhs._ref_id2 &&
	    _antisense_aln == rhs._antisense_aln &&
	    _left == rhs._left && 
	    _antisense_splice == rhs._antisense_splice &&
	    _edit_dist == rhs._edit_dist && 
	    /* DO NOT USE ACCEPTED IN COMPARISON */
	    _cigar == rhs._cigar);
  }
  
  bool operator<(const BowtieHit& rhs) const 
  {
    if (_insert_id != rhs._insert_id)
      return _insert_id < rhs._insert_id;
    if (_ref_id != rhs._ref_id)
      return _ref_id < rhs._ref_id;
    if (_ref_id2 != rhs._ref_id2)
      return _ref_id2 < rhs._ref_id2;
    if (_left != rhs._left)
      return _left < rhs._left;
    if (_antisense_aln != rhs._antisense_aln)
      return _antisense_aln < rhs._antisense_aln;
    if (_edit_dist != rhs._edit_dist)
      return _edit_dist < rhs._edit_dist;
    if (_cigar != rhs._cigar)
      {
	if (_cigar.size() != rhs._cigar.size())
	  return _cigar.size() < rhs._cigar.size();
	for (size_t i = 0; i < _cigar.size(); ++i)
	  {
	    if (!(_cigar[i] == rhs._cigar[i]))
	      return (_cigar[i].opcode < rhs._cigar[i].opcode || (_cigar[i].opcode == rhs._cigar[i].opcode && _cigar[i].length < rhs._cigar[i].length)); 
	  }
      }
    return false;
  }
  
  uint32_t ref_id() const			{ return _ref_id;			}
  uint32_t ref_id2() const			{ return _ref_id2;			}
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
	  case mATCH:
	  case rEF_SKIP:
	  case dEL:
	    r -= op.length;
	    break;
	  case FUSION_FF:
	  case FUSION_FR:
	  case FUSION_RF:
	  case FUSION_RR:
	    r = op.length;
	    break;
	  default:
	    break;
	  }
      }
    return r;			
  }

  bool is_spliced() const
  {
    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	
	if (op.opcode == REF_SKIP || op.opcode == rEF_SKIP)
	  return true;
      }
    return false;
  }

  CigarOpCode fusion_opcode() const
  {
    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	
	if (op.opcode == FUSION_FF || op.opcode == FUSION_FR || op.opcode == FUSION_RF || op.opcode == FUSION_RR)
	  return op.opcode;
      }
    return FUSION_NOTHING;
  }

  /*
   * checks whether its coordinate is increasing or decreasing
   *  before its fusion or until the end.
   */
  bool is_forwarding_left() const
  {
    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	
	if (op.opcode == MATCH || op.opcode == REF_SKIP || op.opcode == INS || op.opcode == DEL)
	  return true;

	if (op.opcode == mATCH || op.opcode == rEF_SKIP || op.opcode == iNS || op.opcode == dEL)
	  return false;
	
	if (op.opcode == FUSION_FF || op.opcode == FUSION_FR || op.opcode == FUSION_RF || op.opcode == FUSION_RR)
	  break;
      }
    
    return true;
  }

  /*
   * checks whether its coordinate is increasing or decreasing
   *  before its fusion or until the end.
   */
  bool is_forwarding_right() const
  {
    for (int i = _cigar.size() - 1; i >= 0; --i)
      {
	const CigarOp& op = _cigar[i];
	
	if (op.opcode == MATCH || op.opcode == REF_SKIP || op.opcode == INS || op.opcode == DEL)
	  return true;

	if (op.opcode == mATCH || op.opcode == rEF_SKIP || op.opcode == iNS || op.opcode == dEL)
	  return false;
	
	if (op.opcode == FUSION_FF || op.opcode == FUSION_FR || op.opcode == FUSION_RF || op.opcode == FUSION_RR)
	  break;
      }
    
    return true;
  }

  bool antisense_splice() const { return _antisense_splice; }
  bool antisense_align() const { return _antisense_aln;	}
  void antisense_align(bool antisense_align) { _antisense_aln = antisense_align;	}
  
  bool antisense_align2() const
  {
    /*
     * antisense_splice is also used to indicate whether fusion is ff, fr, or rf.
     */
    CigarOpCode fusionOpCode = fusion_opcode();
    if (fusionOpCode == FUSION_NOTHING || fusionOpCode == FUSION_FF || fusionOpCode == FUSION_RR)
	return antisense_align();
	
    return !antisense_align();
  }

  BowtieHit reverse() const
  {
    BowtieHit result;
    result._ref_id = _ref_id2;
    result._ref_id2 = _ref_id;
    result._insert_id = _insert_id;

    uint32_t right, fusion_pos;
    right = fusion_pos = _left;

    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	
	switch(op.opcode)
	  {
	  case MATCH:
	  case REF_SKIP:
	  case DEL:
	    right += op.length;
	    break;
	  case mATCH:
	  case rEF_SKIP:
	  case dEL:
	    right -= op.length;
	    break;
	  case FUSION_FF:
	  case FUSION_FR:
	  case FUSION_RF:
	  case FUSION_RR:
	    fusion_pos = right;
	    right = op.length;
	    break;
	  default:
	    break;
	  }
      }

    if (is_forwarding_left())
      fusion_pos -= 1;
    else
      fusion_pos += 1;

    CigarOpCode fusionOpCode = fusion_opcode();
    if (fusionOpCode == FUSION_NOTHING || fusionOpCode == FUSION_FF || fusionOpCode == FUSION_RR)
      {
	if (is_forwarding_left())
	  result._left = right - 1;
	else
	  result._left = right + 1;
      }
    else
      {
	if (fusionOpCode == FUSION_FR)
	  result._left = right + 1;
	else
	  result._left = right - 1;
      }

    result._cigar.clear();
    for (int i = _cigar.size() - 1; i >= 0; --i)
      {
	CigarOp cigar = _cigar[i];

	switch(cigar.opcode)
	  {
	  case MATCH:
	    cigar.opcode = mATCH; break;
	  case mATCH:
	    cigar.opcode = MATCH; break;
	  case INS:
	    cigar.opcode = iNS; break;
	  case iNS:
	    cigar.opcode = INS; break;
	  case DEL:
	    cigar.opcode = dEL; break;
	  case dEL:
	    cigar.opcode = DEL; break;
	  case REF_SKIP:
	    cigar.opcode = rEF_SKIP; break;
	  case rEF_SKIP:
	    cigar.opcode = REF_SKIP; break;
	  case FUSION_FF:
	  case FUSION_FR:
	  case FUSION_RF:
	  case FUSION_RR:
	    cigar.length = fusion_pos; break;
	  default:
	    break;
	  }

	result._cigar.push_back(cigar);
      }

    if (fusionOpCode == FUSION_FR || fusionOpCode == FUSION_RF)
      result._antisense_aln = !_antisense_aln;
    else
      result._antisense_aln = _antisense_aln;
    
    result._antisense_splice = _antisense_splice;
    result._edit_dist = _edit_dist;
    result._splice_mms = _splice_mms;
    result._end = _end;

    result._seq = _seq;
    reverse_complement(result._seq);
    result._qual = _qual;
    ::reverse(result._qual.begin(), result._qual.end());
    
    return result;
  }
  
  unsigned char edit_dist() const	{ return _edit_dist;		}
  unsigned char splice_mms() const	{ return _splice_mms;		}

  int alignment_score() const           { return _alignment_score;      }
  void alignment_score(int as)          { _alignment_score = as;        }
  
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
	  case rEF_SKIP:
	    gaps_out.push_back(make_pair(pos, pos - op.length + 1));
	    pos -= op.length;
	    break;
	  case MATCH:
	  case DEL:
	    pos += op.length;
	    break;
	  case mATCH:
	  case dEL:
	    pos -= op.length;
	    break;
	  case FUSION_FF:
	  case FUSION_FR:
	  case FUSION_RF:
	  case FUSION_RR:
	    pos = op.length;
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
  
  const string& hitfile_rec() const { return _hitfile_rec; }
  void hitfile_rec(const string& rec) { _hitfile_rec = rec; }
  
  const string& seq() const { return _seq; }
  void seq(const string& seq) { _seq = seq; }
  
  const string& qual() const { return _qual; }
  void qual(const string& qual) { _qual = qual; }
  
  bool end() const { return _end; }
  void end(bool end) { _end = end; }

  // this is for debugging purpose
  bool check_editdist_consistency(const RefSequenceTable& rt, bool bDebug = false);
  
private:
  uint32_t _ref_id;
  uint32_t _ref_id2;
  ReadID _insert_id;   // Id of the sequencing insert
  int _left;        // Position in the reference of the left side of the alignment
  
  vector<CigarOp> _cigar;
  
  bool _antisense_splice;    // Whether the junction spanned is on the reverse strand
  bool _antisense_aln;       // Whether the alignment is to the reverse strand
  
  unsigned char _edit_dist;  // Total mismatches (note this is not including insertions or deletions as mismatches, ie, not equivalent to NM field of a SAM record)
  unsigned char _splice_mms; // Mismatches within min_anchor_len of a splice junction
  string _hitfile_rec; // Points to the buffer for the record from which this hit came
  string _seq;
  string _qual;

  int _alignment_score; // Bowtie2 outputs AS (alignment score) in SAM, TopHat2 uses the value when selecting the best alignments.
  bool _end; // Whether this segment is the last one of the read it belongs to
};

class ReadTable
{
public:
	
	ReadTable() : 
		_next_id(1) {}
	
	// This function should NEVER return zero
	ReadID get_id(const string& name)
	{
		uint32_t _id = atoi(name.c_str());
		//assert(_id);
		_next_id = max(_next_id, (size_t)_id);
		return _id;
	}
	
	
	uint32_t observation_order(ReadID ID)
	{
		if (ID == 0)
			return VMAXINT32;
		return ID;
	}
	
	
	size_t size() const { return _next_id; }
	
private:
	size_t _next_id;
};

inline bool REFID_Less(uint32_t ref_id1, uint32_t ref_id2)
{
  return false;
}

inline bool REFID_Equal(uint32_t ref_id1, uint32_t ref_id2)
{
  return ref_id1 == ref_id2;
}

#if 0
class RefSequenceTable
{
public:
	
	typedef seqan::String<seqan::Dna5, seqan::Packed<seqan::Alloc<> > > Sequence;
	
	struct SequenceInfo
	{
		SequenceInfo(uint32_t _order, 
					 char* _name, 
					 Sequence* _seq, uint32_t _len) :
            observation_order(_order),
            name(_name),
            seq(_seq),
            len(_len) {}
        
		uint32_t observation_order;
		char* name;
		Sequence* seq;
        uint32_t len;
	};
	
	typedef map<string, uint64_t> IDTable;
	typedef map<uint32_t, SequenceInfo> InvertedIDTable;
	typedef InvertedIDTable::iterator iterator;
	typedef InvertedIDTable::const_iterator const_iterator;
	
	RefSequenceTable(bool keep_names, bool keep_seqs = false) : 
	_next_id(1), 
	_keep_names(keep_names) {}
    
    RefSequenceTable(const string& sam_header_filename, 
                     bool keep_names, 
                     bool keep_seqs = false) : 
        _next_id(1), 
        _keep_names(keep_names) 
    {
        if (sam_header_filename != "")
        {
            samfile_t* fh = samopen(sam_header_filename.c_str(), "r", 0);
            if (fh == 0) {
                fprintf(stderr, "Failed to open SAM header file %s\n", sam_header_filename.c_str());
                exit(1);
            }
            
            for (size_t i = 0; i < (size_t)fh->header->n_targets; ++i)
            {
                const char* name = fh->header->target_name[i];
                uint32_t len  = fh->header->target_len[i];
                get_id(name, NULL, len);
                //fprintf(stderr, "SQ: %s - %d\n", name, len);
            }
        }
    }
	
	// This function should NEVER return zero
	uint32_t get_id(const string& name,
			Sequence* seq = NULL,
			uint32_t len = 0)
	{
		uint32_t _id = hash_string(name.c_str());
		pair<InvertedIDTable::iterator, bool> ret = 
		_by_id.insert(make_pair(_id, SequenceInfo(_next_id, NULL, NULL, 0)));
		if (ret.second == true)
		{			
			char* _name = NULL;
			if (_keep_names)
				_name = strdup(name.c_str());
			ret.first->second.name  = _name;
			ret.first->second.seq	= seq;
            ret.first->second.len   = len;
			++_next_id;
		}
		assert (_id);
		return _id;
	}
	
	const char* get_name(uint32_t ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.name;
		else
			return NULL;
	}
    
    uint32_t get_len(uint32_t ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.len;
		else
			return 0;
	}
	
	Sequence* get_seq(uint32_t ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
			return itr->second.seq;
		else
			return NULL;
	}
	
	const SequenceInfo* get_info(uint32_t ID) const
	{
		
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
		if (itr != _by_id.end())
		{
			return &(itr->second);
		}
		else
			return NULL;
	}
	
	int observation_order(uint32_t ID) const
	{
		InvertedIDTable::const_iterator itr = _by_id.find(ID);
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

	// daehwan
	// This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
	static inline uint32_t hash_string(const char* __s)
	{
		uint32_t hash = 0x811c9dc5;
		for ( ; *__s; ++__s)
		{
			hash *= 16777619;
			hash ^= *__s;
		}
		return hash;
	}
	
private:
	
	//IDTable _by_name;
	uint32_t _next_id;
	bool _keep_names;
	InvertedIDTable _by_id;
};

#else

class RefSequenceTable
{
 public:
  typedef seqan::String<seqan::Dna5, seqan::Packed<seqan::Alloc<> > > Sequence;
  
  struct SequenceInfo
  {
  SequenceInfo(uint32_t _order, 
	       char* _name, 
	       Sequence* _seq, uint32_t _len) :
    observation_order(_order),
      name(_name),
      seq(_seq),
      len(_len) {}
    
    uint32_t observation_order;
    char* name;
    Sequence* seq;
    uint32_t len;
  };
  
  typedef map<string, SequenceInfo> IDTable;
  typedef IDTable::iterator iterator;
  typedef IDTable::const_iterator const_iterator;
  
 RefSequenceTable(bool keep_names) : 
  _next_id(1), 
    _keep_names(keep_names) {}
  
 RefSequenceTable(const string& sam_header_filename, 
		  bool keep_names) : 
  _next_id(1), 
    _keep_names(keep_names) 
    {
      if (sam_header_filename != "")
        {
	  samfile_t* fh = samopen(sam_header_filename.c_str(), "r", 0);
	  if (fh == 0)
	    {
	      fprintf(stderr, "Failed to open SAM header file %s\n", sam_header_filename.c_str());
	      exit(1);
	    }
          
	  for (size_t i = 0; i < (size_t)fh->header->n_targets; ++i)
            {
                const char* name = fh->header->target_name[i];
                uint32_t len  = fh->header->target_len[i];
                get_id(name, NULL, len);
                // fprintf(stderr, "SQ: %s - %u\n", name, len);
            }
	  // order_recs_lexicographically();
        }
    }
  
  // This function should NEVER return zero
  uint32_t get_id(const string& name,
		  Sequence* seq = NULL,
		  uint32_t len = 0)
  {
    pair<IDTable::iterator, bool> ret = 
      _by_name.insert(make_pair(name, SequenceInfo(_next_id, NULL, NULL, 0)));
    if (ret.second == true)
      {			
	char* _name = NULL;
	if (_keep_names)
	  _name = strdup(name.c_str());
	
	ret.first->second.name = _name;
	ret.first->second.seq = seq;
	ret.first->second.len = len;
	ret.first->second.observation_order = _next_id;
    //assert (_refid_to_hash.size() + 1 == _next_id);
	_refid_to_name.push_back (name);
			
	++_next_id;
      }
    else
      {
	if (seq)
	  {
	    ret.first->second.seq = seq;
	    ret.first->second.len = len;
	  }
      }

    return ret.first->second.observation_order;
  }
  
  const char* get_name(uint32_t ID) const
  {
    const string& name = _refid_to_name[ID-1];
    IDTable::const_iterator itr = _by_name.find(name);
    if (itr != _by_name.end())
      return itr->second.name;
    else
      return NULL;
  }
  
  uint32_t get_len(uint32_t ID) const
  {
    const string& name = _refid_to_name[ID-1];
    IDTable::const_iterator itr = _by_name.find(name);
    if (itr != _by_name.end())
      return itr->second.len;
    else
      return 0;
  }
  
  Sequence* get_seq(uint32_t ID) const
  {
    assert (ID > 0 && ID <= _refid_to_name.size());
    const string& name = _refid_to_name[ID-1];
    IDTable::const_iterator itr = _by_name.find(name);
    if (itr != _by_name.end())
      return itr->second.seq;
    else
      return NULL;
  }
  
  const SequenceInfo* get_info(uint32_t ID) const
  {
    assert (ID > 0 && ID <= _refid_to_name.size());
    const string& name = _refid_to_name[ID-1];
    IDTable::const_iterator itr = _by_name.find(name);
    if (itr != _by_name.end())
      {
	return &(itr->second);
      }
    else
      return NULL;
  }
  
  uint32_t observation_order(uint32_t ID) const
  {
    return ID;
  }
  
  iterator begin() { return _by_name.begin(); }
  iterator end() { return _by_name.end(); }
  
  const_iterator begin() const { return _by_name.begin(); }
  const_iterator end() const { return _by_name.end(); }
  
  size_t size() const { return _by_name.size(); }
	
  void clear()
  {
    _by_name.clear();
  }

  // strnum_cmp is taken from samtools.
  static inline int strnum_cmp(const string &a, const string &b)
  {
    char *pa = (char*)a.c_str(), *pb = (char*)b.c_str();
    while (*pa && *pb)
      {
      if (isdigit(*pa) && isdigit(*pb))
	{
	  long ai, bi;
	  ai = strtol(pa, &pa, 10);
	  bi = strtol(pb, &pb, 10);
	  if (ai != bi) return ai<bi? true : false;
	}
      else
	{
	  if (*pa != *pb) break;
	  ++pa; ++pb;
	}
      }

    if (*pa == *pb)
      return (pa-a.c_str()) < (pb-b.c_str())? true : false;
    
    return *pa<*pb? true : false;
  }

  void order_recs_lexicographically()
  {
    vector<string> vStr;
    for (IDTable::iterator i = _by_name.begin(); i != _by_name.end(); ++i)
      {
	vStr.push_back(i->first);
      }

    ::sort(vStr.begin(), vStr.end(), RefSequenceTable::strnum_cmp);
    
    _refid_to_name.clear();
    size_t new_order = 1;
    for (vector<string>::iterator i = vStr.begin(); i != vStr.end(); ++i, ++new_order)
      {
	_by_name.find(*i)->second.observation_order = new_order;
	_refid_to_name.push_back(*i);
      }
  }
  
private:
  uint32_t _next_id;
  bool _keep_names;
  IDTable _by_name;
  vector<string> _refid_to_name;
};
#endif

bool hit_insert_id_lt(const BowtieHit& h1, const BowtieHit& h2);

typedef vector<BowtieHit> HitList;

/* This class stores all the hits from a Bowtie map */

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

class HitStream;
class HitFactory {
   friend class HitStream;
  public:
    HitFactory(ReadTable& insert_table, 
               RefSequenceTable& reference_table) : 
    _insert_table(insert_table), _ref_table(reference_table)
                 {}
    virtual ~HitFactory() {}
    virtual void openStream(HitStream& hs)=0;
    virtual void rewind(HitStream& hs)=0;
    virtual void seek(HitStream& hs, int64_t offset)=0;
    virtual void closeStream(HitStream& hs)=0;
    BowtieHit create_hit(const string& insert_name, 
			 const string& ref_name,
			 const string& ref_name2,
			 int left,
			 const vector<CigarOp>& cigar,
			 bool antisense_aln,
			 bool antisense_splice,
			 unsigned char edit_dist,
			 unsigned char splice_mms,
			 bool end);
    
    BowtieHit create_hit(const string& insert_name, 
			 const string& ref_name,
			 uint32_t left,
			 uint32_t read_len,
			 bool antisense_aln,
			 unsigned char edit_dist,
			 bool end);
  
   virtual string hitfile_rec(HitStream& hs, const char* hit_buf)=0;
   virtual bool next_record(HitStream& hs, const char*& buf, size_t& buf_size) = 0;
   virtual bool get_hit_from_buf(const char* bwt_buf, 
                      BowtieHit& bh,
                      bool strip_slash,
                      char* name_out = NULL,
                      char* name_tags = NULL,
                      char* seq = NULL,
                      char* qual = NULL) = 0;
    
  protected:
    ReadTable& _insert_table;
    RefSequenceTable& _ref_table;
    HitStream* _hit_stream;
};//class HitFactory


class LineHitFactory  : public HitFactory {
//for text line-based formats like Bowtie and SAM
  public:
	LineHitFactory(ReadTable& insert_table,
			RefSequenceTable& reference_table) :
		HitFactory(insert_table, reference_table) {}

    string hitfile_rec(HitStream& hs, const char* hit_buf) {
      string r(hit_buf);
      return r;
      }
    void openStream(HitStream& hs);
    void rewind(HitStream& hs);
    void seek(HitStream&hs, int64_t offset);
    void closeStream(HitStream& hs);
    bool next_record(HitStream& hs, const char*& buf, size_t& buf_size);
   protected:
	static const size_t _hit_buf_max_sz = 10 * 1024;
	char _hit_buf[_hit_buf_max_sz];
	int _line_num;
  };

class BowtieHitFactory : public LineHitFactory {
  public:
	BowtieHitFactory(ReadTable& insert_table, 
			RefSequenceTable& reference_table) : 
		LineHitFactory(insert_table, reference_table) {}

	bool get_hit_from_buf(const char* bwt_buf, 
			      BowtieHit& bh,
			      bool strip_slash,
			      char* name_out = NULL,
			      char* name_tags = NULL,
			      char* seq = NULL,
			      char* qual = NULL);
};

class SplicedBowtieHitFactory : public LineHitFactory {
  public:
	SplicedBowtieHitFactory(ReadTable& insert_table, 
							RefSequenceTable& reference_table,
							int anchor_length) : 
	LineHitFactory(insert_table, reference_table),
	_anchor_length(anchor_length){}
	
	bool get_hit_from_buf(const char* bwt_buf, 
			      BowtieHit& bh,
			      bool strip_slash,
			      char* name_out = NULL,
			      char* name_tags = NULL,
			      char* seq = NULL,
			      char* qual = NULL);
  private:
	int _anchor_length;
	int _seg_offset;
	int _size_buf;
};

class SplicedSAMHitFactory : public LineHitFactory {
  public:
	SplicedSAMHitFactory(ReadTable& insert_table,
							RefSequenceTable& reference_table,
							int anchor_length) :
	LineHitFactory(insert_table, reference_table),
	_anchor_length(anchor_length){}

	bool get_hit_from_buf(const char* bwt_buf,
			      BowtieHit& bh,
			      bool strip_slash,
			      char* name_out = NULL,
			      char* name_tags = NULL,
			      char* seq = NULL,
			      char* qual = NULL);
  private:
	int _anchor_length;
	int _seg_offset;
	int _size_buf;
};


class SAMHitFactory : public LineHitFactory {
  public:
	SAMHitFactory(ReadTable& insert_table, 
				  RefSequenceTable& reference_table) : 
	   LineHitFactory(insert_table, reference_table) {}
	
	bool get_hit_from_buf(const char* bwt_buf, 
			      BowtieHit& bh,
			      bool strip_slash,
			      char* name_out = NULL,
			      char* name_tags = NULL,
			      char* seq = NULL,
			      char* qual = NULL);
};


/******************************************************************************
 BAMHitFactory turns SAM alignments into BowtieHits
 *******************************************************************************/
class BAMHitFactory : public HitFactory {
  public:

    BAMHitFactory(ReadTable& insert_table,
                  RefSequenceTable& reference_table) :
        HitFactory(insert_table, reference_table)
        {
         _sam_header=NULL;
        }
    void openStream(HitStream& hs);
    void rewind(HitStream& hs);
    void seek(HitStream& hs, int64_t offset);
    void closeStream(HitStream& hs);

    bool next_record(HitStream& hs, const char*& buf, size_t& buf_size);
	
	/*void mark_curr_pos()
	{ 
		_curr_pos = bgzf_tell(((samfile_t*)_hit_file)->x.bam);
	}

	
	void undo_hit() 
	{ 
		bgzf_seek(((samfile_t*)_hit_file)->x.bam, _curr_pos, SEEK_SET);
		_eof = false;
	}
    */
	string hitfile_rec(HitStream& hs, const char* hit_buf);

	bool get_hit_from_buf(const char* bwt_buf, 
			      BowtieHit& bh,
			      bool strip_slash,
			      char* name_out = NULL,
			      char* name_tags = NULL,
			      char* seq = NULL,
			      char* qual = NULL);
protected:
	//int64_t _curr_pos;
	//int64_t _beginning;
	bam1_t _next_hit; 
	bam_header_t* _sam_header;
    bool inspect_header(HitStream& hs);
};

class SplicedBAMHitFactory : public BAMHitFactory {
 public:
 SplicedBAMHitFactory(ReadTable& insert_table,
		      RefSequenceTable& reference_table,
		      int anchor_length) :
  BAMHitFactory(insert_table, reference_table),
    _anchor_length(anchor_length)
    {
    }

    bool get_hit_from_buf(const char* bwt_buf, 
			  BowtieHit& bh,
			  bool strip_slash,
			  char* name_out = NULL,
			  char* name_tags = NULL,
			  char* seq = NULL,
			  char* qual = NULL);
    
 private:
  int _anchor_length;
  int _seg_offset;
  int _size_buf;
};


struct HitsForRead
{
	HitsForRead() : insert_id(0) {}
	uint64_t insert_id;
	vector<BowtieHit> hits;
};

class HitStream
{
  friend class HitFactory;
  friend class LineHitFactory;
  friend class BAMHitFactory;
 //private:
  HitFactory* _factory;
  bool _spliced;
  bool _strip_slash;
  BowtieHit buffered_hit;
  bool _keep_bufs;
  bool _keep_seqs;
  bool _keep_quals;
  bool _from_bowtie;
  void* _hit_file;
  string _hit_file_name;
  FZPipe* _fzpipe;
  bool _eof;

public:
  HitStream(void* hit_file, //could be FILE* or samfile_t*
		  HitFactory* hit_factory, 
		  bool spliced, 
		  bool strip_slash, 
		  bool keep_bufs,
		  bool keep_seqs = false,
		  bool keep_quals = false,
		  bool from_bowtie = false) :
	_factory(hit_factory),
	_spliced(spliced),
	_strip_slash(strip_slash),
	buffered_hit(BowtieHit()),
	  _keep_bufs(keep_bufs),
	  _keep_seqs(keep_seqs),
	  _keep_quals(keep_quals),
	  _from_bowtie(from_bowtie),
      _hit_file(hit_file),
      _hit_file_name(),
      _fzpipe(NULL),
      _eof(false)
	{
		primeStream();
	}

	HitStream(const string& hit_filename,
		  HitFactory* hit_factory,
		  bool spliced,
		  bool strip_slash,
		  bool keep_bufs,
		  bool keep_seqs = false,
		  bool keep_quals = false,
		  bool from_bowtie = false) :
	_factory(hit_factory),
	_spliced(spliced),
	_strip_slash(strip_slash),
	buffered_hit(BowtieHit()),
	  _keep_bufs(keep_bufs),
	  _keep_seqs(keep_seqs),
	  _keep_quals(keep_quals),
	  _from_bowtie(from_bowtie),
      _hit_file(NULL),
      _hit_file_name(hit_filename),
      _fzpipe(NULL),
      _eof(false)
	{
	    _factory->openStream(*this);
	    primeStream();
	}

    HitStream(FZPipe& hit_filepipe,
          HitFactory* hit_factory,
          bool spliced,
          bool strip_slash,
          bool keep_bufs,
          bool keep_seqs = false,
          bool keep_quals = false,
          bool from_bowtie = false) :
    _factory(hit_factory),
    _spliced(spliced),
    _strip_slash(strip_slash),
    buffered_hit(BowtieHit()),
      _keep_bufs(keep_bufs),
      _keep_seqs(keep_seqs),
      _keep_quals(keep_quals),
      _from_bowtie(from_bowtie),
      _hit_file(NULL),
      _hit_file_name(),
      _fzpipe(&hit_filepipe),
      _eof(false)
    {
        _hit_file=_fzpipe->file;
        primeStream();
    }

    void primeStream() { //why?
      // Prime the stream by reading a single hit into the buffered_hit
      HitsForRead dummy = HitsForRead();
      next_read_hits(dummy);
    }
    bool eof() { return _eof; }
    bool ready() { return (_hit_file!=NULL); }
    void reset() {
		_factory->rewind(*this);
		_eof=false;
		// re-prime the stream;
		buffered_hit = BowtieHit();
		primeStream();
	}
    void seek(int64_t offset) {
      _factory->seek(*this, offset);
      _eof = false;
      buffered_hit = BowtieHit();
      primeStream();
    }
	
	bool next_read_hits(HitsForRead& hits_for_read)
	{
	  hits_for_read.hits.clear();
	  hits_for_read.insert_id = 0; 
	  
	  //if (!_hit_file || (feof(_hit_file) && buffered_hit.insert_id() == 0))
	  //  return false;
	  if (!this->ready())
	      //err_die("Error: next_read_hits() called on HitFactory with no file handle\n");
	      return false;
	  if (this->eof() && buffered_hit.insert_id() == 0) {
	          return false;
	          }

	  //char bwt_buf[2048]; bwt_buf[0] = 0;
	  char bwt_seq[2048]; bwt_seq[0] = 0;
	  char bwt_qual[2048]; bwt_qual[0] = 0;
	  
	  char* seq = _keep_seqs ? bwt_seq : NULL;
	  char* qual = _keep_quals ? bwt_qual : NULL;
	  
	  hits_for_read.insert_id = buffered_hit.insert_id();
	  if (hits_for_read.insert_id)
	    hits_for_read.hits.push_back(buffered_hit);
	  const char* hit_buf;
      size_t hit_buf_size = 0;
	  while (true) {

		if (!_factory->next_record(*this, hit_buf, hit_buf_size)) {
		              buffered_hit = BowtieHit();
		              break; }

		  //string clean_buf = bwt_buf;
		  // Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		if (_factory->get_hit_from_buf(hit_buf, bh, _strip_slash, 
		                         NULL, NULL, seq, qual)) {
		  if (_keep_bufs)
		    bh.hitfile_rec(_factory->hitfile_rec(*this, hit_buf));
		  
		  if (_keep_seqs)
		      bh.seq(seq);
		  
		  if (_keep_quals) {
		    // when it comes to convert from qual in color to qual in bp,
		    // we need to fill in the two extream qual values using the adjacent qual values.
		    size_t qual_len = strlen(qual);
		    if (color && qual_len > 2) {
			  qual[0] = qual[1];
			  qual[qual_len-1] = qual[qual_len-2];
			  }
			bh.qual(qual);
		    }
		  
		  if (bh.insert_id() == hits_for_read.insert_id) {
		      hits_for_read.hits.push_back(bh);
		      }
		  else {
		      buffered_hit = bh;
		      break;
		      }
		  } //hit parsed
	  } //while reading hits
	  
	  return (!hits_for_read.hits.empty() && hits_for_read.insert_id != 0);
	}
	
	uint64_t next_group_id() const 
	{
		return buffered_hit.insert_id();
	}
	bool fromBowtie() { return _from_bowtie; }
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

void add_hits_to_coverage(const HitList& hits, vector<unsigned short>& DoC);
void add_hit_to_coverage(const BowtieHit& bh, vector<unsigned int>& DoC);

void accept_all_hits(HitTable& hits);

//print BowtieHit as BAM record
void print_bamhit(GBamWriter& wbam,
		  const char* read_name,
		  const BowtieHit& bh,
		  const char* ref_name,
		  const char* ref_name2,
		  const char* sequence,
		  const char* qualities,
		  bool from_bowtie = false,
		  const vector<string>* extra_fields = NULL);


void extract_partial_hits(const BowtieHit& bh, const string& seq, const string& qual,
			  char* cigar1, char* cigar2, string& seq1, string& seq2,
			  string& qual1, string& qual2, int& left1, int& left2);

/**
 * Convert a vector of CigarOps to a string representation 
 */
std::string print_cigar(const vector<CigarOp>& bh_cigar);

/**
 * Calculate bowtie (1 or 2) related extra SAM fields such as
 * AS:i (alignment score)
 * MD:Z 
 * NM:i
 * etc
 */
void bowtie_sam_extra(const BowtieHit& bh, const RefSequenceTable& rt, vector<string>& fields);

#endif
