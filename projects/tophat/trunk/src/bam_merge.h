#ifndef BAM_MERGE_H
#define BAM_MERGE_H

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <queue>
#include <algorithm>

#include "common.h"

using namespace std;

extern bool raw_merge;

struct CBamLine {
  int filenum;
  uint64_t read_id;
  bam1_t* b;
  
CBamLine(int fno=-1, bam1_t* br=NULL, bam_header_t* header = NULL) :
  filenum(fno),
    read_id(0), b(br) {
    b_init(header);
  }
  
  void b_init(bam_header_t* header = NULL);
  void b_free();
};

struct equal_bam {
  bool operator() (const CBamLine& first, const CBamLine& second) const {
    if (raw_merge) return false;
    if (first.read_id != second.read_id)
      return false;

    if (first.b->core.tid != second.b->core.tid)
      return false;

    if (first.b->core.pos != second.b->core.pos)
      return false;

    if (first.b->core.n_cigar != second.b->core.n_cigar)
      return false;

    for (int i = 0; i < first.b->core.n_cigar; ++i){
      if (bam1_cigar(first.b)[i] != bam1_cigar(second.b)[i])
	return false;
    }

    // for fusion alignments, two alignments are always not equal
    if (bam_aux_get(first.b, "XF") || bam_aux_get(second.b, "XF"))
      return false;

    return true;
  }
};

struct less_bam {
  bool rev_cmp; //reverse the comparison
  
  less_bam(bool reverse_cmp = false) {
    rev_cmp = reverse_cmp;
  }
  
  bool operator() (const CBamLine& f, const CBamLine& s) const {
    if (raw_merge) return false;
    const CBamLine* first = &f;
    const CBamLine* second = &s;
    if (rev_cmp) {
      first = &s;
      second = &f;
    }
    
    if (first->read_id != second->read_id)
      return first->read_id < second->read_id;
    
    if (first->b->core.tid != second->b->core.tid)
      return first->b->core.tid < second->b->core.tid;

    if (first->b->core.pos != second->b->core.pos)
      return first->b->core.pos < second->b->core.pos;

    if (first->b->core.n_cigar != second->b->core.n_cigar)
      return first->b->core.n_cigar < second->b->core.n_cigar;

    for (int i = 0; i < first->b->core.n_cigar; ++i){
      if (bam1_cigar(first->b)[i] != bam1_cigar(second->b)[i])
	return bam1_cigar(first->b)[i] < bam1_cigar(second->b)[i];
    }

    // prefer a record with XS attribute
    char strand1 = 0, strand2 = 0;
    uint8_t* ptr = bam_aux_get(first->b, "XS");
    if (ptr) strand1 = bam_aux2A(ptr);

    ptr = bam_aux_get(second->b, "XS");
    if (ptr) strand2 = bam_aux2A(ptr);

    if (strand1 != strand2)
      {
	if (strand1 == '+' || strand2 == 0)
	  return true;
	else
	  return false;
      }

    // prefer more aux fields.
    if (bowtie2)
      {
	if (first->b->data_len != second->b->data_len)
	  return first->b->data_len > second->b->data_len;
      }
    
    return false;
  }
};

class BamMerge
{
 public:
  BamMerge(const vector<string>& bam_fnames,
	   vector<int64_t> file_offsets = vector<int64_t>());
  
  ~BamMerge();

 public:
  bool next_bam_lines(vector<CBamLine>& bam_lines);
  bam_header_t* get_sam_header()
  {
    if (_src_files.size() > 0)
      return _src_files[0]->header;
    else
      return NULL;
  }


 private:
  vector<string> _bam_fnames;
  priority_queue<CBamLine, vector<CBamLine>, less_bam> _lines;
  vector <samfile_t*> _src_files; //array of SAM file handles
  uint64_t _last_id;
};

#endif
