#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <getopt.h>

#include "bam/bam.h"
#include "bam/sam.h"

bool is_fastq=true; //default is fastq
bool sam_input=false; //default is BAM
bool all_reads=false;
bool mapped_only=false;
char outfname[1024];

#define USAGE "Usage: bam2fastx [--fasta|-a|--fastq|-q] [--sam|-s|-t] \n\
        [-M|--mapped-only|-A|--all] [-o <outfile>] <in.bam>|<in.sam> \n"



char qseq[2048];

const char *short_options = "o:aqstMUA";

enum {
   OPT_FASTA = 127,
   OPT_FASTQ,
   OPT_SAM,
   OPT_MAPPED_ONLY,
   OPT_ALL
   };
   
struct option long_options[] = {
  {"fasta", no_argument, 0, OPT_FASTA},
  {"fastq", no_argument, 0, OPT_FASTQ},
  {"sam", no_argument, 0, OPT_SAM},
  {"mapped-only", no_argument, 0, OPT_MAPPED_ONLY},
  {"all", no_argument, 0, OPT_ALL},
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
       case 's':
       case 't':
       case OPT_SAM: //sam (text) input
         sam_input = true;
         break;
       case 'M':
       case OPT_MAPPED_ONLY:
         mapped_only = true;
         break;
       case 'A':
       case OPT_ALL:
         all_reads = true;
         break;
       case 'o':
         strcpy(outfname, optarg);
         break;
       default:
         return 1;
       }
   } while(next_option != -1);
  if (all_reads && mapped_only) {
    fprintf(stderr, "Error: incompatible options !\n");
    exit(2);
    }
  return 0;
}

#define bam_unmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)

void showfastq(const bam1_t *b, samfile_t* fp, FILE* fout) {
  char *name  = bam1_qname(b);
  char *qual  = (char*)bam1_qual(b);
  char *s    = (char*)bam1_seq(b);
  int i;
  bool ismapped=((b->core.flag & BAM_FUNMAP) == 0);
  if (ismapped && !all_reads) return;
  if (mapped_only && !ismapped) return;
   
  for(i=0;i<(b->core.l_qseq);i++) {
    int8_t v = bam1_seqi(s,i);
    qseq[i] = bam_nt16_rev_table[v];
    }
  qseq[i] = 0;
  
  fprintf(fout, "@%s\n%s\n",name, qseq);
  for(i=0;i<(b->core.l_qseq);i++) {
    qseq[i]=qual[i]+33;
    }
  qseq[i]=0;
  fprintf(fout, "+\n%s\n",qseq);
}

void showfasta(const bam1_t *b, samfile_t* fp, FILE* fout) {
  char *name  = bam1_qname(b);
  char *s    = (char*)bam1_seq(b);
  int i;
  for(i=0;i<(b->core.l_qseq);i++) {
    int8_t v = bam1_seqi(s,i);
    qseq[i] = bam_nt16_rev_table[v];
  }
  qseq[i] = 0;
  fprintf(fout, ">%s\n%s\n",name, qseq);
}



int main(int argc, char *argv[])
{
    samfile_t *fp;
    outfname[0]=0;
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
    if (sam_input)
        fp = samopen(fname, "r", 0);
      else 
        fp = samopen(fname, "rb", 0);
    if (fp == 0) {
        fprintf(stderr, "Error: bam2fastx failed to open BAM file %s\n", fname);
        return 1;
        }
    FILE* fout=stdout;
    if (outfname[0]) {
       fout=fopen(outfname, "w");
       if (fout==NULL) {
           fprintf(stderr, "Error creating output file %s\n", outfname);
           return 2;
           }
       }
    bam1_t *b = bam_init1();
    if (is_fastq) {
        while (samread(fp, b) >= 0) showfastq(b, fp, fout);
        }
      else {
        while (samread(fp, b) >= 0) showfasta(b, fp, fout);
        }
    
    bam_destroy1(b);
    
    samclose(fp);
    return 0;
}
