#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <string>

#include "bam/bam.h"
#include "bam/sam.h"

using namespace std;

bool is_fastq=true; //default is fastq
bool sam_input=false; //default is BAM
bool all_reads=false;
bool mapped_only=false;
bool add_matenum=false;
bool pairs=false;
bool ignoreQC; // ignore qc fail flag 0x400

string outfname;

#define USAGE "Usage: bam2fastx [--fasta|-a|--fastq|-q] [-Q] [--sam|-s|-t]\n\
   [-M|--mapped-only|-A|--all] [-o <outfile>] [-P|--paired] [-N] <in.bam>\n\
   \nNote: By default, reads flagged as not passing quality controls are\n\
   discarded; the -Q option can be used to ignore the QC flag.\n\
   \nUse the -N option if the /1 and /2 suffixes should be applied to\n\
   read names according to the SAM flags (otherwise these suffixes\n\
   are not added unless the -P option is used)\n"

const char *short_options = "o:ac:qstQMAPN";

enum {
   OPT_FASTA = 127,
   OPT_FASTQ,
   OPT_SAM,
   OPT_PAIRED,
   OPT_MAPPED_ONLY,
   OPT_ALL
   };
   
struct Read {
	string name;
	int mate;
	string seq;
	string qual;
	void clear() {
	  name.clear();
	  mate=0;
	  seq.clear();
	  qual.clear();
	  }
};

struct option long_options[] = {
  {"fasta", no_argument, 0, OPT_FASTA},
  {"fastq", no_argument, 0, OPT_FASTQ},
  {"sam", no_argument, 0, OPT_SAM},
  {"paired", no_argument, 0, OPT_PAIRED},
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
       case 'P':
       case OPT_PAIRED:
         pairs = true;
         add_matenum=true;
         break;
       case 'Q':
    	 ignoreQC = true;
    	 break;
       case 'o':
         outfname=optarg;
         break;
       case 'N':
    	 add_matenum=true;
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

void getRead(const bam1_t *b, samfile_t* fp, Read& rd) {
  static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
  rd.clear();
  char *name  = bam1_qname(b);
  rd.name=name;
  char *qual  = (char*)bam1_qual(b);
  char *s    = (char*)bam1_seq(b);
  int i;

  if ((b->core.flag & BAM_FQCFAIL) && !ignoreQC) return;
  bool ismapped=((b->core.flag & BAM_FUNMAP) == 0);
  if (ismapped && !all_reads) return;
  if (mapped_only && !ismapped) return;

  bool isreversed=((b->core.flag & BAM_FREVERSE) != 0);
  bool is_paired = ((b->core.flag & BAM_FPAIRED) != 0);
  if (is_paired || add_matenum) {
     if (b->core.flag & BAM_FREAD1)
         rd.mate=1;
     else if (b->core.flag & BAM_FREAD2)
         rd.mate=2;
     }
  int seqlen = b->core.l_qseq;
  if (seqlen>0) {
	rd.seq.resize(seqlen);
	for(i=0;i<seqlen;i++) {
	  rd.seq[i] = bam1_seqi(s,i);
	}
	// copied from sam_view.c in samtools.
	if (isreversed) {
	   int l=0;
	   int r=seqlen-1;
	   while (l<r) {
		  char c=rd.seq[l];
		  rd.seq[l]=rd.seq[r];
		  rd.seq[r]=c;
		  l++;r--;
		  }
	   for (i=0;i<seqlen;i++) {
		 rd.seq[i]=seq_comp_table[(int)rd.seq[i]];
	   }
	}

	for(i=0;i<seqlen;i++) {
	  rd.seq[i] = bam_nt16_rev_table[(int)rd.seq[i]];
	}

	if (!is_fastq) return;

	rd.qual.resize(seqlen);
	for(i=0;i<seqlen;i++) {
	  if (qual[i]==0xFF)
		rd.qual[i]='I';
	  else rd.qual[i]=qual[i]+33;
	  }

	if (isreversed) {
	  int l=0;
	  int r=seqlen-1;
	  while (l<r) {
		int8_t t = rd.qual[l];
		rd.qual[l] = rd.qual[r];
		rd.qual[r] = t;
		l++;r--;
	  }
	}
  }
}

void writeRead(Read& rd, int& wpair, FILE* fout) {
  // shouldn't get an empty sequence here
  //if (rd.seq.empty()) {
  //	    return;
  //      }
  if (is_fastq) {
    if (rd.mate>0) {
      fprintf(fout, "@%s/%d\n%s\n",rd.name.c_str(), rd.mate, rd.seq.c_str());
      wpair|=rd.mate;
    }
     else {
       fprintf(fout, "@%s\n%s\n",rd.name.c_str(), rd.seq.c_str());
	   wpair|=4;
     }
    fprintf(fout, "+\n%s\n",rd.qual.c_str());
  }
  else {
    if (rd.mate>0) {
      fprintf(fout, ">%s/%d\n%s\n",rd.name.c_str(), rd.mate, rd.seq.c_str());
      wpair|=rd.mate;
    }
	 else {
	   fprintf(fout, ">%s\n%s\n",rd.name.c_str(), rd.seq.c_str());
	   wpair|=4;
	 }
  }
}

void writePaired(Read& rd, int& wpair, FILE* fout, FILE* fout2) {
  if (rd.mate==1) {
	 writeRead(rd, wpair, fout);
	  //if (write_mapped && last1>rd.name)
	  //err_order(last1, rd.name);
  }
  else if (rd.mate==2) {
	writeRead(rd, wpair, fout2);
	//if (write_mapped && last2>rd.name)
	  //      err_order(last2, rd.name);
  }
  else {
	fprintf(stderr, "Error: unpaired read encountered (%s)\n", rd.name.c_str());
	exit(1);
  }
}


void err_order(string& last) {
  fprintf(stderr, "Error: couldn't retrieve both reads for pair %s. "
	  "Perhaps the input file is not sorted by name?\n"
	  "(using 'samtools sort -n' might fix this)\n", last.c_str());
  exit(1);
}


string getFBase(const string& s, string& ext, string& pocmd) {
   string fbase(s);
   ext="";
   pocmd="";
   if (s.empty() || s=="-") return fbase;
   //size_t slen=s.length();
   size_t p=s.rfind('.');
   size_t d=s.rfind('/');
   if (p==string::npos || (d!=string::npos && p<d)) return fbase;
   fbase=s.substr(0,p);
   ext=s.substr(p+1);
   string fext(ext);
   if (fext.length()>2) fext=fext.substr(0,2);
   for(size_t i=0; i!=fext.length(); i++)
        fext[i] = std::tolower(fext[i]);
   if (fext=="gz" || fext=="z") pocmd="gzip -c ";
           else if (fext=="bz") pocmd="bzip2 -c ";
   if (!pocmd.empty()) {
	 p=fbase.rfind('.');
	 d=fbase.rfind('/');
     if (p==string::npos || (d!=string::npos && p<d)) return fbase;
     ext=fbase.substr(p+1)+"."+ext;
     return fbase.substr(0,p);
     }
   return fbase;
   }


int main(int argc, char *argv[])
{
    samfile_t *fp;
    char* fname=NULL;
    bool use_pclose=false;
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
    FILE* fout2=NULL;
    if (pairs && outfname.empty()) {
      fprintf(stderr, "Error: paired output (-P) requires the -o option.\n");
      return 1;
      }
    if (!outfname.empty()) {
 	   string fext;
 	   string pocmd;
 	   string fbase=getFBase(outfname, fext, pocmd);
       if (pairs) {
    	 outfname=fbase+".1."+fext;
    	 string out2=fbase+".2."+fext;
         if (!pocmd.empty()) {
               out2=pocmd+">"+out2;
               fout2=popen(out2.c_str(),"w");
               }
            else {
               fout2=fopen(out2.c_str(),"w");
               }
         if (fout2==NULL) {
      	    fprintf(stderr, "Error opening file stream: %s\n", out2.c_str());
      	    return 1;
            }
         }
  	   string out1(outfname);
       if (!pocmd.empty()) {
             out1=pocmd+">"+out1;
             fout=popen(out1.c_str(),"w");
             use_pclose=true;
             }
          else {
             fout=fopen(out1.c_str(),"w");
             }
       if (fout==NULL) {
    	    fprintf(stderr, "Error opening file stream: %s\n", out1.c_str());
    	    return 1;
          }
       }

    bam1_t *b = bam_init1();
    Read rd;
    //bool write_mapped=(all_reads || mapped_only);
    string last;
    int wpair=0; //writing pair status bitmask (bit 1 set mate 1 was written,
                //      bit 2 set if mate 2 was written, bit 3 set if unpaired read was written)
    while (samread(fp, b) >= 0) {
       getRead(b, fp, rd);
       if (rd.seq.empty())
    	 continue; //skip secondary alignments with no sequence
       int pstatus=rd.mate;
       if (rd.mate==0) pstatus=4;
       if (last!=rd.name) {
   	      if (pairs && !last.empty() && wpair!=3)
   	    	 err_order(last);
    	  wpair=0;
    	  last=rd.name;
          }
       if ( (pstatus & wpair)==0) {
    	    if (pairs) {
    	       writePaired(rd, wpair, fout, fout2);
    	     } //paired
    	     else  { //single reads
    		   writeRead(rd, wpair, fout);
    	     }
		   } //new pair
      }
    if (fout!=stdout) {
      if (use_pclose) pclose(fout);
      else fclose(fout);
    }
    if (fout2) {
      if (use_pclose) pclose(fout2);
      else fclose(fout2);
    }
    bam_destroy1(b);
    
    samclose(fp);
    return 0;
}
