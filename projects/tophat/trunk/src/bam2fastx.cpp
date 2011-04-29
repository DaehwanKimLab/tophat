#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "bam/bam.h"
#include "bam/sam.h"

static bool is_fastq=true; //default is fastq
static char qseq[2048];

const char *short_options = "aq";

enum {
   OPT_FASTA = 127,
   OPT_FASTQ
   };
   
static struct option long_options[] = {
  {"fasta", no_argument, 0, OPT_FASTA},
  {"fastq", no_argument, 0, OPT_FASTQ},
  {0, 0, 0, 0} // terminator
  };

int parse_options(int argc, char** argv)
{
  int option_index = 0;
  int next_option;
  do {
     next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
     switch (next_option) {
       case -1:    
         break;
       case 'a':
       case OPT_FASTA:
         is_fastq = false;
         break;
       case 'q':
       case OPT_FASTQ:
         is_fastq = true;
         break;
       default:
         return 1;
       }
   } while(next_option != -1);
  
  return 0;
}


void showfastq(const bam1_t *b, samfile_t* fp) {
  //uint32_t *cigar = bam1_cigar(b);
  //const bam1_core_t *c = &b->core;
  char *name  = bam1_qname(b);
  char *qual  = (char*)bam1_qual(b);
  char *s    = (char*)bam1_seq(b);
  int i;
  for(i=0;i<(b->core.l_qseq);i++) {
    int8_t v = bam1_seqi(s,i);
    qseq[i] = bam_nt16_rev_table[v];
  }
  qseq[i] = 0;
  
  printf("@%s\n%s\n",name, qseq);
  //print SAM text format of this BAM record:
  //char * samtext=bam_format1(fp_in->header, b);
  //printf("%s\n",samtext);

  //printf("s cigar: %s\n",cigar);

  //printf("qual : ");
  for(i=0;i<(b->core.l_qseq);i++) {
    //printf(" %d",qual[n]);
    qseq[i]=qual[i]+33;        
    }
  qseq[i]=0;
  printf("+\n%s\n",qseq);
}

void showfasta(const bam1_t *b, samfile_t* fp) {
  char *name  = bam1_qname(b);
  char *s    = (char*)bam1_seq(b);
  int i;
  for(i=0;i<(b->core.l_qseq);i++) {
    int8_t v = bam1_seqi(s,i);
    qseq[i] = bam_nt16_rev_table[v];
  }
  qseq[i] = 0;
  printf(">%s\n%s\n",name, qseq);
}


#define USAGE "Usage: bam2fastx [--fasta|-a|--fastq|-q] <in.bam>\n"

int main(int argc, char *argv[])
{
    samfile_t *fp;
    char* fname=NULL;
    
    if (parse_options(argc, argv) || optind>=argc) {
       fprintf(stderr, USAGE);
       return -1;
       }
    fname=argv[optind++];

    if (fname==NULL || fname[0]==0) {
        fprintf(stderr, USAGE);
        return 1;
        }
    
    if ((fp = samopen(fname, "rb", 0)) == 0) {
        fprintf(stderr, "Error: bam2fastx failed to open BAM file %s\n", fname);
        return 1;
        }
    
    bam1_t *b = bam_init1();
    if (is_fastq) {
        while (samread(fp, b) >= 0) showfastq(b, fp);
        }
      else {
        while (samread(fp, b) >= 0) showfasta(b, fp);
        }
    
    bam_destroy1(b);
    
    samclose(fp);
    return 0;
}
