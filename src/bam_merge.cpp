#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "common.h"
#include <queue>
#include <algorithm>

using namespace std;

#define USAGE "Usage: bam_merge [-Q] <out.bam> <in1.bam> <in2.bam> [...]\n"

#define ERR_BAM_OPEN "Error: bam_merge failed to open BAM file %s\n"

void print_usage()
{
  fprintf(stderr, USAGE);
}


bool raw_merge = false;

vector <samfile_t*> srcfiles; //array of SAM file handles

struct CBamLine {
  int fileno;
  uint64_t read_id;
  bam1_t* b;

  CBamLine(int fno=-1, bam1_t* br=NULL): fileno(fno),
    read_id(0), b(br) {
    b_init();
  }

  void b_init() {
    if (b) {
      char *name  = bam1_qname(b);
      if (raw_merge) {
    	read_id=0;
    	return;
      }
      read_id=(uint64_t)atol(name);
      if (read_id<1) {
        char* samline=bam_format1(srcfiles[0]->header, b);
        err_die("Error: invalid read Id (must be numeric) for BAM record:\n%s\n",
                  samline);
      }
    }
  }

  void b_free() {
    if (b!=NULL) {
      bam_destroy1(b);
      b=NULL;
      }
  }
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
  less_bam(bool reverse_cmp=false) {
     rev_cmp=reverse_cmp;
     }
  bool operator() (const CBamLine& f, const CBamLine& s) const {
	 if (raw_merge) return false;
	 const CBamLine* first = &f;
	 const CBamLine* second = &s;
     if (rev_cmp) {
        first=&s;
        second=&f;
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

    if (strand1 == strand2)
      return false;

    if (strand1 == '+' || strand2 == 0)
      return true;
    else
      return false;
  }
};

void write_bam_lines(GBamWriter& bw, vector<CBamLine>& bam_lines) {
  std::sort (bam_lines.begin(), bam_lines.end(), less_bam());

  bool sense_strand = false, antisense_strand = false;
  vector<CBamLine>::size_type i = 0;
  for (; i < bam_lines.size(); ++i) {
    bool do_write = true;
    const CBamLine& bam_line = bam_lines[i];

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
      //samwrite(fw, bam_line.b);
      bw.write(bam_line.b, bam_line.read_id);
    }
  }

  for (i = 0; i < bam_lines.size(); ++i)
    bam_lines[i].b_free();

  bam_lines.clear();
}


//GList<CBamLine> lines(true,true,true); //sorted, free items, unique items

priority_queue< CBamLine , vector<CBamLine>, less_bam > lines(less_bam(true));

int main(int argc, char *argv[])
{
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  if (quals) raw_merge=true; //hijack the -Q flag for this specific merging option
  char* outfname=NULL;
  if (argc-optind<3) {
    print_usage();
    if (argc>1)
      warn_msg("Error: only %d arguments given.\n", argc-1);
    return -1;
  }
  outfname=argv[optind];
  int num_src=argc-optind-1;
  srcfiles.resize(num_src, NULL);
  for (int i=optind+1;i<argc;i++) {
    int fno=i-optind-1;
    samfile_t* fp=samopen(argv[i], "rb", 0);
    if (fp==0) {
      warn_msg(ERR_BAM_OPEN, argv[i]);
      return 1;
    }
    bam1_t *b = bam_init1();
    if (samread(fp, b) > 0) {
      srcfiles[fno]=fp;
      CBamLine brec(fno, b);
      lines.push(brec);
      }
    else { bam_destroy1(b); }
  }
  if (lines.size()==0) {
    warn_msg("Warning: no input BAM records found.\n");
    return 0;
  }

  vector<CBamLine> bam_lines;
  GBamWriter bamwriter(outfname, srcfiles[0]->header, index_outfile);

  uint64_t last_id = 0;
  while (lines.size()>0) {
    CBamLine brec(lines.top()); //should have the smallest read_id
    lines.pop();
    assert (brec.fileno>=0 && brec.b!=NULL);
    // we need to eliminate duplicate records, which can happen when using Bowtie2
    // as we may map the same read against transcriptome, genome, and novel/known splice junctions.
    if (last_id != brec.read_id && bam_lines.size() > 0)
       write_bam_lines(bamwriter, bam_lines);

    last_id = brec.read_id;
    bam_lines.push_back(brec);
    //reuse brec
    brec.b = bam_init1();
    if (samread(srcfiles[brec.fileno], brec.b)>0) {
      brec.b_init();
      lines.push(brec);
    }
    else { //no more BAM records
      brec.b_free();
    }
  }

  if (bam_lines.size() > 0)
    write_bam_lines(bamwriter, bam_lines);
  for (int i=0;i<num_src;i++)
    samclose(srcfiles[i]);
  srcfiles.clear();
  return 0;
}
