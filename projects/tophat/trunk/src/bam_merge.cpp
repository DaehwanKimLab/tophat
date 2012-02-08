#include <vector>
#include <algorithm>
using namespace std;

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <bam/bam.h>
#include <bam/sam.h>

#include "common.h"
#include "GBase.h"
#include "GList.hh"

#define USAGE "Usage: bam_merge <out.bam> <in1.bam> <in2.bam> [...]\n"

#define ERR_BAM_OPEN "Error: bam_merge failed to open BAM file %s\n"

void print_usage()
{
  fprintf(stderr, USAGE);
}

samfile_t **srcfiles; //array of SAM file handles
samfile_t *fw; //output SAM handle

class CBamLine {
 public:
  int fileno;
  long read_id;
  bam1_t* b;
  bool operator==(CBamLine& l){
    return (read_id==l.read_id);
  }
  bool operator>(CBamLine& l){
    return (read_id<l.read_id); //(last has lowest #)                                                                                                                                                                                                                           
  }
  bool operator<(CBamLine& l){
    return (read_id>l.read_id);
  }

  CBamLine(int fno=-1, bam1_t* br=NULL) {
    fileno=fno;
    read_id=0;
    b=br;
    b_init();
  }
  void b_init() {
    if (b) {
      char *name  = bam1_qname(b);
      read_id=atol(name);
      if (read_id<1) {
	char* samline=bam_format1(srcfiles[0]->header, b);
	GError("Error: invalid read Id (must be numeric) for BAM record:\n%s\n",
	       samline);
      }
    }
  }
  
  void b_free() {
    if (b!=NULL) {  bam_destroy1(b); }
  }
};

struct equal_bam {
  bool operator() (const CBamLine& first, const CBamLine& second) const {
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

    return true;
  }
};

struct less_bam {
  bool operator() (const CBamLine& first, const CBamLine& second) const {
    if (first.read_id != second.read_id)
      return first.read_id < second.read_id;

    if (first.b->core.tid != second.b->core.tid)
      return first.b->core.tid < second.b->core.tid;

    if (first.b->core.pos != second.b->core.pos)
      return first.b->core.pos < second.b->core.pos;

    if (first.b->core.n_cigar != second.b->core.n_cigar)
      return first.b->core.n_cigar < second.b->core.n_cigar; 

    for (int i = 0; i < first.b->core.n_cigar; ++i){
      if (bam1_cigar(first.b)[i] != bam1_cigar(second.b)[i])
	return bam1_cigar(first.b)[i] < bam1_cigar(second.b)[i];
    }

    // prefer a record with XS attribute
    char strand1 = 0, strand2 = 0;
    uint8_t* ptr = bam_aux_get(first.b, "XS");
    if (ptr) strand1 = bam_aux2A(ptr);

    ptr = bam_aux_get(second.b, "XS");
    if (ptr) strand2 = bam_aux2A(ptr);

    if (strand1 == strand2)
      return false;

    if (strand1 == '+' || strand2 == 0)
      return true;
    else
      return false;
  }
};

void write_bam_lines(samfile_t* fw, vector<CBamLine>& bam_lines, FILE* findex = NULL) {
  static size_t count = 0;
  std::sort (bam_lines.begin(), bam_lines.end(), less_bam());

  bool sense_strand = false, antisense_strand = false;
  vector<CBamLine>::size_type i = 0;
  for (; i < bam_lines.size(); ++i) {
    bool do_write = true;
    const CBamLine& bam_line = bam_lines[i];
    if (i == 0 && findex != NULL) {
	if (count >= 1000) {
	  int64_t offset = bam_tell(fw->x.bam);
	  int block_offset = offset & 0xFFFF;
	  
	  // daehwan - this is a bug in bgzf.c in samtools
	  // I'll report this bug with a temporary solution, soon.
	  if (block_offset <= 0xF000)
	    {
	      fprintf(findex, "%lu\t%ld\n", bam_line.read_id, offset);
	      count = 0;
	    }
	}
      }

    uint8_t* ptr = bam_aux_get(bam_line.b, "XS");
    char strand = 0;
    if (ptr) strand = bam_aux2A(ptr);

    if (i > 0) {
      if (equal_bam()(bam_lines[i-1], bam_line)) {
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

    if (do_write) {
      samwrite(fw, bam_line.b);
      ++count;
    }
  }

  for (i = 0; i < bam_lines.size(); ++i)
    bam_lines[i].b_free();
		
  bam_lines.clear();
}

GList<CBamLine> lines(true,true,false);

int main(int argc, char *argv[])
{
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  
  char* outfname=NULL;
  if (argc-optind<3) {
    print_usage();
    if (argc>1)
      fprintf(stderr, "Error: only %d arguments given.\n", argc-1);
    return -1;
  }
  outfname=argv[optind];
  int num_src=argc-optind-1;
  GMALLOC(srcfiles, (num_src*sizeof(samfile_t*)));
  for (int i=optind+1;i<argc;i++) {
    int fno=i-optind-1;
    samfile_t* fp=samopen(argv[i], "rb", 0);
    if (fp==0) {
      fprintf(stderr, ERR_BAM_OPEN, argv[i]);
      return 1;
    }
    bam1_t *b = bam_init1();
    if (samread(fp, b) > 0) {
      srcfiles[fno]=fp;
      lines.Add(new CBamLine(fno, b));
    }
  }
  if (lines.Count()==0) {
    GMessage("Warning: no input BAM records found.\n");
  }
  fw=samopen(outfname, "wb", srcfiles[lines[0]->fileno]->header);
  if (fw==NULL)
    GError("Error creating output file %s\n", outfname);
  
  FILE* findex = NULL;
  if (!index_outfile.empty()) {
    findex = fopen(index_outfile.c_str(), "w");
    if (findex == NULL)
      err_die("Error: cannot create file %s\n", index_outfile.c_str());
  }
  
  vector<CBamLine> bam_lines;
  
  int last;
  long last_id = 0;
  while ((last=lines.Count()-1)>=0) {
    CBamLine* from=lines[last]; //should have the smallest read_id
    if (from->fileno<0 || from->b==NULL)
      GError("Invalid processTopLine() call with uninitialized value.");
    
    // we need to eliminate duplicate records, which can happen when using Bowtie2
    // as we may map the same read against transcriptome and genome.
    if (last_id != from->read_id && bam_lines.size() > 0)
      write_bam_lines(fw, bam_lines, findex);
    
    last_id = from->read_id;
    bam_lines.push_back(*from);

    from->b = bam_init1();
    if (samread(srcfiles[from->fileno], from->b)>0) {
      from->b_init();
      //adjust the position in the sorted lines list
      if (last<7) {//                                                                                                                                                                                                                                                       
	int i=last;
	while (i>0 && lines[i-1]->read_id<lines[i]->read_id) {
	  CBamLine* tmp=lines[i-1];
	  lines.Put(i-1, lines[i], false);
	  lines.Put(i,tmp, false);
	  i--;
	}
      }
      else { //use quick-insert                                                                                                                                                                                                                                          
	lines.Pop();
	lines.Add(from);
      }
    }
    else { //no more BAM records
      lines.Pop();
      from->b_free();
    }
  }
  
  if (bam_lines.size() > 0)
    write_bam_lines(fw, bam_lines, findex);
  
  samclose(fw);
  for (int i=0;i<num_src;i++)
    samclose(srcfiles[i]);
  GFREE(srcfiles);
  
  if (findex != NULL)
    fclose(findex);
  
  return 0;
}
