#include "bam_merge.h"

#define ERR_BAM_OPEN "Error: bam_merge failed to open BAM file %s\n"

bool raw_merge = false;

void CBamLine::b_init(bam_header_t* header) {
  if (b) {
    char *name  = bam1_qname(b);
    if (raw_merge) {
      read_id=0;
      return;
    }
    read_id=(uint64_t)atol(name);
    if (read_id<1 && header) {
      char* samline=bam_format1(header, b);
      err_die("Error: invalid read Id (must be numeric) for BAM record:\n%s\n",
	      samline);
    }
  }
}

void CBamLine::b_free() {
  if (b!=NULL) {
    bam_destroy1(b);
    b=NULL;
  }
}

BamMerge::BamMerge(const vector<string>& bam_fnames,
		   vector<int64_t> file_offsets) :
  _bam_fnames(bam_fnames),
  _lines(less_bam(true)),
  _last_id(0)
{
  if (bam_fnames.size() <= 0)
    return;
  
  for (size_t i = 0; i < _bam_fnames.size(); ++i) {
    const char* fname = _bam_fnames[i].c_str();
    samfile_t* fp = samopen(fname, "rb", 0);
    if (fp==0) {
      warn_msg(ERR_BAM_OPEN, fname);
      exit(1);
    }

    if (bam_fnames.size() == file_offsets.size() &&
	file_offsets[i] > 0)
      bgzf_seek(fp->x.bam, file_offsets[i], SEEK_SET);

    bam1_t* b = bam_init1();
    if (samread(fp, b) > 0) {
      _src_files.push_back(fp);
      CBamLine brec(_lines.size(), b, fp->header);
      _lines.push(brec);
    }
    else { bam_destroy1(b); }
  }

  if (_lines.size() == 0) {
    warn_msg("Warning: no input BAM records found.\n");
    exit(1);
  }
}

BamMerge::~BamMerge()
{
  while (_lines.size() > 0)
    {
      CBamLine brec(_lines.top());
      brec.b_free();
      _lines.pop();
    };
  
  for (size_t i = 0; i < _src_files.size(); ++i)
    samclose(_src_files[i]);
  
  _src_files.clear();
}

bool BamMerge::next_bam_lines(vector<CBamLine>& bam_lines)
{
  if (_lines.size() <= 0)
    return false;
  
  bam_lines.clear();
  vector<CBamLine> temp_bam_lines;
  while (_lines.size() > 0) {
    CBamLine brec(_lines.top()); //should have the smallest read_id
    assert (brec.fileno>=0 && brec.b!=NULL);

    if ((raw_merge || _last_id != brec.read_id) && temp_bam_lines.size() > 0)
      {
	break;
      }

    _lines.pop();
    _last_id = brec.read_id;
    temp_bam_lines.push_back(brec);
    //reuse brec
    brec.b = bam_init1();
    if (samread(_src_files[brec.filenum], brec.b)>0) {
      brec.b_init(_src_files[brec.filenum]->header);
      _lines.push(brec);
    }
    else { //no more BAM records
      brec.b_free();
    }
  }

  if (temp_bam_lines.size() <= 0)
    return false;

  // we need to eliminate duplicate alignments, which can happen when using Bowtie2
  // as we may often map the same read against transcriptome, genome, and novel/known splice junctions.
  std::sort (temp_bam_lines.begin(), temp_bam_lines.end(), less_bam());
  bool sense_strand = false, antisense_strand = false;
  vector<bool> free_indexes(temp_bam_lines.size(), false);
  for (size_t i = 0; i < temp_bam_lines.size(); ++i) {
    bool do_write = true;
    CBamLine& bam_line = temp_bam_lines[i];

    uint8_t* ptr = bam_aux_get(bam_line.b, "XS");
    char strand = 0;
    if (ptr) strand = bam_aux2A(ptr);

    if (i > 0) {
      if (equal_bam()(temp_bam_lines[i-1], bam_line)) {
	if (strand == 0) {
	  do_write = false;
	} else {
	  if (strand == '+' && sense_strand) do_write = false;
	  else sense_strand = true;

	  if (strand == '-' && antisense_strand) do_write = false;
	  else antisense_strand = true;
	}
      } else {
	sense_strand = false;
	antisense_strand = false;
      }
    }

    if (strand == '+')
      sense_strand = true;
    else if (strand == '-')
      antisense_strand = true;

    if (do_write)
      bam_lines.push_back(bam_line);
    else
      free_indexes[i] = true;
  }

  for (size_t i = 0; i < free_indexes.size(); ++i)
    {
      if (free_indexes[i])
	temp_bam_lines[i].b_free();
    }
  
  return bam_lines.size() > 0;
}
