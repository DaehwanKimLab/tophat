#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "local"
#define SVN_REVISION "unknown"
#endif


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <assert.h>
#include <string>

#include "bam.h"
#include "sam.h"

using namespace std;

bool is_fastq=true; //default is fastq
bool sam_input=false; //default is BAM
bool all_reads=false;
bool mapped_only=false;
bool add_matenum=false;
bool pairs=false;
bool color=false;
bool ignoreQC=false; // ignore qc fail flag 0x400
bool ignoreOQ=false; // ignore OQ tag
long s_w=0; //in paired mode, number of singletons written

string outfname;

#define USAGE "bam2fastx v%s (%s) usage:\n\
 bam2fastx [--fasta|-a] [-C|--color] [-P|--paired] [-N]\n\
 [-A|--all|-M|--mapped-only] [-Q] [--sam|-s|-t] [-o <outfname>] <in.bam>\n\
 \nBy default, bam2fastx only converts the unmapped reads from the input file,\n\
  discarding those unmapped reads flagged as QC failed.\n\
  The input BAM/SAM file MUST be sorted by read name (-n option for samtools\n\
  sort). If the input file name is \"-\", stdin will be used instead.\n\
 \nOptions:\n\
 -A,--all         convert all reads (mapped and unmapped)\n\
                  (but discarding those flagged as QC failed, unless -Q)\n\
 -P               paired reads are expected and converted into two output\n\
                  files (see <outfname> comments below)\n\
 -Q               convert unmapped reads even when flagged as QC failed\n\
 -M,--maped-only  convert only mapped reads\n\
 -N               for -P, append  /1 and /2 suffixes to read names\n\
 -O               ignore the original quality values (OQ tag) and write the\n\
                  current quality values (default is to use OQ data if found)\n\
 -C,--color       reads are in ABI SOLiD color format\n\
 -s,-t,--sam      input is a SAM text file (default: BAM input expected)\n\
 -a,--fasta       output FASTA records, not FASTQ (discard quality values)\n\
 -o <outfname>    output file name or template (see below)\n\
\n\
 <outfname> serves as a name template when -P option is provided, as suffixes\n\
 .1 and .2 (and .s) will be automatically inserted before the file extension\n\
 in <outfname>, such that multiple output files will be created.\n\
 If <outfname> ends in .gz or .bz2 then bam2fastx will write the\n\
 output compressed by gzip or bzip2 respectively.\n\n\
 Example of converting all paired reads from a BAM file to FASTQ format:\n\
    bam2fastx -PANQ -o sample.fq.gz sample.sortedbyname.bam\n\
 In this example the output will be written in two files: \n\
   sample.1.fq.gz and sample.2.fq.gz (and optionally sample.s.fq.gz\n\
   for unpaired reads)\n\
"

const char *short_options = "o:ac:qvhstOQCMAPN";

enum {
   OPT_HELP = 127,
   OPT_VERSION,
   OPT_FASTA,
   OPT_FASTQ,
   OPT_SAM,
   OPT_PAIRED,
   OPT_MAPPED_ONLY,
   OPT_ALL,
   OPT_COLOR
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
  {"help", no_argument, 0, OPT_HELP},
  {"version", no_argument, 0, OPT_VERSION},
  {"fasta", no_argument, 0, OPT_FASTA},
  {"fastq", no_argument, 0, OPT_FASTQ},
  {"sam", no_argument, 0, OPT_SAM},
  {"paired", no_argument, 0, OPT_PAIRED},
  {"mapped-only", no_argument, 0, OPT_MAPPED_ONLY},
  {"all", no_argument, 0, OPT_ALL},
  {"color", no_argument, 0, OPT_COLOR},
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
      case 'h':
      case OPT_HELP:
         fprintf(stdout, USAGE, PACKAGE_VERSION, SVN_REVISION);
         exit(0);
      case 'v':
      case OPT_VERSION:
         fprintf(stdout, "%s\n", PACKAGE_VERSION);
         exit(0);
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
      case 'C':
      case OPT_COLOR:
        color = true;
        break;
     case 'P':
        case OPT_PAIRED:
        pairs = true;
        break;
      case 'Q':
        ignoreQC = true;
        break;
      case 'O':
        ignoreOQ = true;
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
    fprintf(stderr, "Error: incompatible options, use either -A/--all or -M/--mapped-only!\n");
    exit(2);
    }
  return 0;
}

void getRead(const bam1_t *b, samfile_t* fp, Read& rd) {
  static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
  rd.clear();
  char *name  = bam1_qname(b);
  rd.name=name;
  unsigned char *qual  = NULL;
  unsigned char *s    = (unsigned char*)bam1_seq(b);
  int i;

  if ((b->core.flag & BAM_FQCFAIL) && !ignoreQC) return;
  bool ismapped=((b->core.flag & BAM_FUNMAP) == 0);
  if (ismapped && !(all_reads || mapped_only)) return;
  if (mapped_only && !ismapped) return;

  bool isreversed=((b->core.flag & BAM_FREVERSE) != 0);
  // bool is_paired = ((b->core.flag & BAM_FPAIRED) != 0);
  //if (add_matenum) {
     if (b->core.flag & BAM_FREAD1)
         rd.mate=1;
     else if (b->core.flag & BAM_FREAD2)
         rd.mate=2;
  //   }
  int seqlen = b->core.l_qseq;
  if (seqlen>0) {
	rd.seq.resize(seqlen);
	for(i=0;i<seqlen;i++) {
	  rd.seq[i] = bam1_seqi(s,i);
	}
	// copied from sam_view.c in samtools.
	if (isreversed) {

	  // didn't implement this for colorspace reads
	  assert (!color);
	  
	   int l=0;
	   int r=seqlen-1;
	   while (l<r) {
		  char c=rd.seq[l];
		  rd.seq[l]=rd.seq[r];
		  rd.seq[r]=c;
		  l++;r--;
		  }

	   if (!color) {
	     for (i=0;i<seqlen;i++) {
	       rd.seq[i]=seq_comp_table[(int)rd.seq[i]];
	     }
	   }
	}

  	if (color) {
	  const static char *color_bam_nt16_rev_table = "4014244434444444";
	  rd.seq[0] = bam_nt16_rev_table[(int)rd.seq[0]];
	  for(i=1;i<seqlen;i++) {
	    rd.seq[i] = color_bam_nt16_rev_table[(int)rd.seq[i]];
	  }
	}
	else {
	  for(i=0;i<seqlen;i++) {
	    rd.seq[i] = bam_nt16_rev_table[(int)rd.seq[i]];
	  }
	}

	if (!is_fastq) return;
	if (color) {
      qual = (unsigned char*)bam1_qual(b);
	  rd.qual.resize(seqlen-1);
	  for(i=1;i<seqlen;i++) {
	    if (qual[i]==0xFF)
	      rd.qual[i]='I';
	    else rd.qual[i-1]=qual[i]+33;
	  }
	}
	else {
	  bool fromOQ=false;
	  if (!ignoreOQ) {
	   uint8_t* ptr = bam_aux_get(b, "OQ");
	   if (ptr) {
	     fromOQ=true;
	     qual = (unsigned char*)bam_aux2Z(ptr);
	   }
	  }
	  if (fromOQ) {
	    rd.qual = (char*)qual;
	    //for(i=0;i<seqlen;i++) {
	    //  rd.qual[i]=qual[i];
	    //}
	  }
	  else {
		qual = (unsigned char*)bam1_qual(b);
	    rd.qual.resize(seqlen);
	    for(i=0;i<seqlen;i++) {
	      if (qual[i]==0xFF)
	        rd.qual[i]='I';
	      else rd.qual[i]=qual[i]+33;
	    }
	  }
	}

	if (isreversed) {
	  int l=0;
	  int r=seqlen-1;
	  if (color)
	    r -= 1;
	  
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
  wpair |= (rd.mate>0 ? rd.mate : 4);
  if (is_fastq) {
    if (rd.mate && add_matenum) {
      fprintf(fout, "@%s/%d\n%s\n",rd.name.c_str(), rd.mate, rd.seq.c_str());
    }
    else {
      fprintf(fout, "@%s\n%s\n",rd.name.c_str(), rd.seq.c_str());
    }
    fprintf(fout, "+\n%s\n",rd.qual.c_str());
  }
  else {
    if (rd.mate && add_matenum) {
      fprintf(fout, ">%s/%d\n%s\n",rd.name.c_str(), rd.mate, rd.seq.c_str());
    }
	 else {
	  fprintf(fout, ">%s\n%s\n",rd.name.c_str(), rd.seq.c_str());
	 }
  }
}

void flushPair(Read* pbuf, int& wpair, FILE* fout, FILE* fout2, FILE* fout_s) {
	if (pbuf[0].seq.length()>0) {
		if (pbuf[1].seq.length()>0) { //both mates present
		  if (pbuf[0].mate==1) writeRead(pbuf[0], wpair, fout);
		  else if (pbuf[0].mate==2)
			writeRead(pbuf[0], wpair, fout2);
		  else { writeRead(pbuf[0], wpair, fout_s); ++s_w; }
		  if (pbuf[1].mate==1) writeRead(pbuf[1], wpair, fout);
		  else if (pbuf[1].mate==2)
			writeRead(pbuf[1], wpair, fout2);
		  else { writeRead(pbuf[1], wpair, fout_s); ++s_w; }
		} //paired reads
		else { writeRead(pbuf[0], wpair, fout_s); ++s_w; }
	}
	else if (pbuf[1].seq.length()>0)
		{ writeRead(pbuf[1], wpair, fout_s); ++s_w; }
	pbuf[0].clear();
	pbuf[1].clear();
}


void writePaired(Read& rd, int& wpair, FILE* fout, FILE* fout2) {
  if (rd.mate==1) {
	 writeRead(rd, wpair, fout);
  }
  else if (rd.mate==2) {
	writeRead(rd, wpair, fout2);
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
       fprintf(stderr, USAGE, PACKAGE_VERSION, SVN_REVISION);
       return -1;
       }
    fname=argv[optind++];

    if (fname==NULL || fname[0]==0) {
        fprintf(stderr, USAGE, PACKAGE_VERSION, SVN_REVISION);
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
    FILE* fout_s=NULL;
    if (pairs && outfname.empty()) {
      fprintf(stderr, "Error: paired output (-P) requires the -o option.\n");
      return 1;
      }
    string f_single;
    if (!outfname.empty()) {
    	string fext;
    	string pocmd;
    	string fbase=getFBase(outfname, fext, pocmd);
    	if (pairs) {
    		outfname=fbase+".1."+fext;
    		string out2=fbase+".2."+fext;
    		string out_s=fbase+".s."+fext;
    		f_single=out_s;
    		if (!pocmd.empty()) {
    			out2=pocmd+">"+out2;
    			fout2=popen(out2.c_str(),"w");
    			out_s=pocmd+">"+out_s;
    			fout_s=popen(out_s.c_str(),"w");
    		}
    		else {
    			fout2=fopen(out2.c_str(),"w");
    			fout_s=fopen(out_s.c_str(),"w");
    		}
    		if (fout2==NULL) {
    			fprintf(stderr, "Error opening file stream: %s\n", out2.c_str());
    			return 1;
    		}
    		if (fout_s==NULL) {
    			fprintf(stderr, "Error opening file stream: %s\n", out_s.c_str());
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
    Read pbuf[2];
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
       if (last!=rd.name) { //new read/pair
   	      //if (pairs && !last.empty() && wpair!=3)
   	    	 //err_order(last); //previous pair was incomplete
    	  //flushPair(Read* pbuf, int& wpair, FILE* fout, FILE* fout2, FILE* fout_s)
   	      flushPair(pbuf, wpair, fout, fout2, fout_s);
    	  wpair=0;
    	  last=rd.name;
          }
       if ( (pstatus & wpair)==0) {
    	    if (pairs) {
    	       //writePaired(rd, wpair, fout, fout2);
    	       //keep pair in buffer pbuf
    	    	if (rd.mate>0) {
    	    		wpair|=rd.mate;
    	    		if (rd.mate>2) {
    	    			fprintf(stderr, "Error invalid mate index (%d) for read %s\n", rd.mate, rd.name.c_str());
    	    			return 1;
    	    		}
    	    		pbuf[rd.mate-1]=rd; //copy read in the buffer
    	    	}
    	    	else {
    	    	  wpair |= 4;
    	    	  pbuf[0]=rd;
    	    	}
    	     } //paired I/O
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
    if (fout_s) {
       if (use_pclose) pclose(fout_s);
       else fclose(fout_s);
       if (s_w==0) unlink(f_single.c_str());
    }
    bam_destroy1(b);
    samclose(fp);
    return 0;
}
