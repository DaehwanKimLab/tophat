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
#include <seqan/sequence.h>

#include <bam/sam.h>
using namespace std;

#include "common.h"
#define _FBUF_SIZE 10*1024
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
	
	bool operator==(const CigarOp& rhs) const 
	{ 
		return opcode == rhs.opcode && length == rhs.length; 
	}
	
};

typedef uint32_t ReadID;

class RefSequenceTable;

/*  Stores the information from a single record of the bowtie map. A given read
    may have many of these.  Reads up to 255bp are supported. 
*/
struct BowtieHit
{
	BowtieHit() : 
		_ref_id(0),
		_insert_id(0),
		_left(0),
		_antisense_splice(false),
		_antisense_aln(false),
		_edit_dist(0),
		_splice_mms(0),
		_end(false){}
	
	BowtieHit(uint32_t ref_id,
		  ReadID insert_id, 
		  int left, 
		  int read_len, 
		  bool antisense,
		  unsigned char edit_dist,
		  bool end) :
		_ref_id(ref_id),
		_insert_id(insert_id), 
		_left(left), 
		_cigar(vector<CigarOp>(1,CigarOp(MATCH,read_len))),
		_antisense_splice(false),
		_antisense_aln(antisense),
		_edit_dist(edit_dist),
		_splice_mms(0),
		_end(end)
	{
		assert(_cigar.capacity() == _cigar.size());
	}
	
	BowtieHit(uint32_t ref_id,
		  ReadID insert_id, 
		  int left,  
		  const vector<CigarOp>& cigar,
		  bool antisense_aln,
		  bool antisense_splice,
		  unsigned char edit_dist,
		  unsigned char splice_mms,
		  bool end) : 
		_ref_id(ref_id),
		_insert_id(insert_id), 	
		_left(left),
		_cigar(cigar),
		_antisense_splice(antisense_splice),
		_antisense_aln(antisense_aln),
		_edit_dist(edit_dist),
		_splice_mms(splice_mms),
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

  bool is_spliced()
  {
    for (size_t i = 0; i < _cigar.size(); ++i)
      {
	const CigarOp& op = _cigar[i];
	
	if (op.opcode == REF_SKIP)
	  return true;
      }
    return false;
  }
	
	bool antisense_splice()	const		{ return _antisense_splice; }
	bool antisense_align() const		{ return _antisense_aln;	}

	unsigned char edit_dist() const		{ return _edit_dist;		}
	unsigned char splice_mms() const	{ return _splice_mms;		}
	
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
	
	const string& hitfile_rec() const { return _hitfile_rec; }
	void hitfile_rec(const string& rec) { _hitfile_rec = rec; }

  	const string& seq() const { return _seq; }
	void seq(const string& seq) { _seq = seq; }

  	const string& qual() const { return _qual; }
	void qual(const string& qual) { _qual = qual; }
  
        bool end() const { return _end; }
  void end(bool end) { _end = end; }

  // this is for debugging purpose
  bool check_editdist_consistency(const RefSequenceTable& rt);

private:
	
  uint32_t _ref_id;
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
			return 0xFFFFFFFF;
		return ID;
	}
	
	
	size_t size() const { return _next_id; }
	
private:
	size_t _next_id;
};

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
					Sequence* seq,
                    uint32_t len)
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
    virtual void closeStream(HitStream& hs)=0;
    BowtieHit create_hit(const string& insert_name, 
                 const string& ref_name,
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
private:
	//int64_t _curr_pos;
	//int64_t _beginning;
	bam1_t _next_hit; 
	bam_header_t* _sam_header;
    bool inspect_header(HitStream& hs);
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

	HitStream(string& hit_filename,
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
		    if (color && !color_out && qual_len > 2) {
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

//print BowtieHit as SAM line
void print_hit(FILE* fout, 
	       const char* read_name,
	       const BowtieHit& bh,
	       const char* ref_name,
	       const char* sequence,
	       const char* qualities,
	       bool from_bowtie = false);

//print BowtieHit as BAM record
void print_bamhit(GBamWriter& wbam,
           const char* read_name,
           const BowtieHit& bh,
           const char* ref_name,
           const char* sequence,
           const char* qualities,
           bool from_bowtie = false);

/**
 * Convert a vector of CigarOps to a string representation 
 */
std::string print_cigar(vector<CigarOp>& bh_cigar);

#endif
