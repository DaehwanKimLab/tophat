/*
 *  common.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/26/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <cstdarg>
#include <getopt.h>

#include "common.h"

using namespace std;

#ifdef MEM_DEBUG 
//function for debugging memory usage of current program in Linux

#include <unistd.h>
#include <ios>
#include <fstream>

//////////////////////////////////////////////////////////////////////////////
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;
   vm_usage     = 0.0;
   resident_set = 0.0;
   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);
   // dummy vars for leading entries in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

void print_mem_usage() {
  double vs, rs;
  process_mem_usage(vs,rs);
  vs/=1024;
  rs/=1024;
  fprintf(stderr, "VMSize: %6.1fMB\tRSize: %6.1fMB\n", vs, rs);
  }
#endif


unsigned int max_insertion_length = 3;
unsigned int max_deletion_length = 3;


int inner_dist_mean = 200;
int inner_dist_std_dev = 20;
int max_mate_inner_dist = -1; 

int min_anchor_len = 8;
int min_report_intron_length = 50;
int max_report_intron_length = 500000;

int min_closure_intron_length = 50;
int max_closure_intron_length = 5000;

int min_coverage_intron_length = 50;
int max_coverage_intron_length = 20000;

int min_segment_intron_length = 50;
int max_segment_intron_length = 500000;

uint32_t min_closure_exon_length = 100; 

int island_extension = 25;
int segment_length = 25;
int segment_mismatches = 2;
int max_read_mismatches = 2;
int max_splice_mismatches = 1;

ReadFormat reads_format = FASTQ;

bool verbose = false;

unsigned int max_multihits = 40;
bool no_closure_search = false;
bool no_coverage_search = false;
bool no_microexon_search = false;
bool butterfly_search = false;
int num_cpus = 1;
float min_isoform_fraction = 0.15f;

string output_dir = "tophat_out";
string aux_outfile = ""; //auxiliary output file name (e.g. prep_reads read stats)
string gene_filter = "";
string gff_file = "";
string ium_reads = "";
string sam_header = "";
string sam_readgroup_id = "";
string zpacker = "";
string samtools_path = "samtools";

bool solexa_quals = false;
bool phred64_quals = false;
bool quals = false;
bool integer_quals = false;
bool color = false;
bool color_out = false;

string gtf_juncs = "";

string flt_reads = "";
string flt_mappings = "";

eLIBRARY_TYPE library_type = LIBRARY_TYPE_NONE;

extern void print_usage();

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

int parseIntOpt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloatOpt(float lower, float upper, const char *errmsg, void (*print_usage)()) {
    float l;
    l = (float)atof(optarg);
	
    if (l < lower) {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    if (l > upper)
    {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    return l;
	
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/*
  this is from
  http://www.winehq.org/pipermail/wine-patches/2001-November/001322.html
 */
char* get_token(char** str, const char* delims)
{
  char* token;
  if (*str == NULL)
      return NULL;
  
  token = *str;
  while (**str != '\0')
    {
      if (strchr(delims, **str) != NULL)
	{
	  **str = '\0';
	  ++(*str);
	  return token;
	}
    
      ++(*str);
    }
  
  *str = NULL;
  return token;
}


const char *short_options = "QCp:z:N:";

enum
  {
    OPT_FASTA = 127,
    OPT_FASTQ,
    OPT_MIN_ANCHOR,
    OPT_SPLICE_MISMATCHES,
    OPT_VERBOSE,
    OPT_INSERT_LENGTH_MEAN,
    OPT_INSERT_LENGTH_STD_DEV,
    OPT_MIN_ISOFORM_FRACTION,
    OPT_OUTPUT_DIR,
    OPT_GENE_FILTER,
    OPT_GFF_ANNOTATIONS,
    OPT_MAX_MULTIHITS,
    OPT_NO_CLOSURE_SEARCH,
    OPT_NO_COVERAGE_SEARCH,
    OPT_NO_MICROEXON_SEARCH,
    OPT_SEGMENT_LENGTH,
    OPT_READ_MISMATCHES,
    OPT_SEGMENT_MISMATCHES,
    OPT_MIN_CLOSURE_EXON,
    OPT_MAX_CLOSURE_INTRON,
    OPT_MIN_CLOSURE_INTRON,
    OPT_MAX_COVERAGE_INTRON,
    OPT_MIN_COVERAGE_INTRON,
    OPT_MIN_SEGMENT_INTRON,
    OPT_MAX_SEGMENT_INTRON,
    OPT_MIN_REPORT_INTRON,
    OPT_MAX_REPORT_INTRON,
    OPT_IUM_READS,
    OPT_BUTTERFLY_SEARCH,
    OPT_SOLEXA_QUALS,
    OPT_PHRED64_QUALS,
    OPT_SAM_HEADER,
    OPT_SAM_READGROUP_ID,
    OPT_QUALS,
    OPT_INTEGER_QUALS,
    OPT_COLOR,
    OPT_COLOR_OUT,
    OPT_LIBRARY_TYPE,
    OPT_MAX_DELETION_LENGTH,
    OPT_MAX_INSERTION_LENGTH,
    OPT_NUM_CPUS,
    OPT_ZPACKER,
    OPT_SAMTOOLS,
    OPT_AUX_OUT,
    OPT_GTF_JUNCS,
    OPT_FILTER_READS,
    OPT_FILTER_HITS
  };

static struct option long_options[] = {
{"fasta",		no_argument,		0,	OPT_FASTA},
{"fastq",		no_argument,		0,	OPT_FASTQ},
{"min-anchor",		required_argument,	0,	OPT_MIN_ANCHOR},
{"sam-header",		required_argument,	0,	OPT_SAM_HEADER},
{"rg-id",		required_argument,	0,	OPT_SAM_READGROUP_ID},
{"splice-mismatches",	required_argument,	0,	OPT_SPLICE_MISMATCHES},
{"verbose",		no_argument,		0,	OPT_VERBOSE},
{"inner-dist-mean",	required_argument,	0,	OPT_INSERT_LENGTH_MEAN},
{"inner-dist-std-dev",	required_argument,	0,	OPT_INSERT_LENGTH_STD_DEV},
{"output-dir",		required_argument,	0,	OPT_OUTPUT_DIR},
{"gene-filter",		required_argument,	0,	OPT_GENE_FILTER},
{"gtf-annotations",	required_argument,	0,	OPT_GFF_ANNOTATIONS},
{"max-multihits",	required_argument,	0,  OPT_MAX_MULTIHITS},
{"no-closure-search",	no_argument,		0,  OPT_NO_CLOSURE_SEARCH},
{"no-coverage-search",	no_argument,		0,  OPT_NO_COVERAGE_SEARCH},
{"no-microexon-search",	no_argument,		0,  OPT_NO_MICROEXON_SEARCH},
{"segment-length",	required_argument,	0,  OPT_SEGMENT_LENGTH},
{"segment-mismatches",	required_argument,	0,  OPT_SEGMENT_MISMATCHES},
{"max-mismatches",  required_argument,  0,  OPT_READ_MISMATCHES},
{"min-closure-exon",	required_argument,	0,  OPT_MIN_CLOSURE_EXON},
{"min-closure-intron",	required_argument,	0,  OPT_MIN_CLOSURE_INTRON},
{"max-closure-intron",	required_argument,	0,  OPT_MAX_CLOSURE_INTRON},
{"min-coverage-intron",	required_argument,	0,  OPT_MIN_COVERAGE_INTRON},
{"max-coverage-intron",	required_argument,	0,  OPT_MAX_COVERAGE_INTRON},
{"min-segment-intron",	required_argument,	0,  OPT_MIN_SEGMENT_INTRON},
{"max-segment-intron",	required_argument,	0,  OPT_MAX_SEGMENT_INTRON},
{"min-report-intron",	required_argument,	0,  OPT_MIN_REPORT_INTRON},
{"max-report-intron",	required_argument,	0,  OPT_MAX_REPORT_INTRON},
{"min-isoform-fraction",required_argument,	0,  OPT_MIN_ISOFORM_FRACTION},
{"ium-reads",		required_argument,	0,  OPT_IUM_READS},
{"butterfly-search",	no_argument,		0,	OPT_BUTTERFLY_SEARCH},
{"solexa-quals",	no_argument,		0,	OPT_SOLEXA_QUALS},
{"phred64-quals",	no_argument,		0,	OPT_PHRED64_QUALS},
{"quals",		no_argument,		0,	OPT_QUALS},
{"integer-quals",	no_argument,		0,	OPT_INTEGER_QUALS},
{"color",		no_argument,		0,	OPT_COLOR},
{"color-out",		no_argument,		0,	OPT_COLOR_OUT},
{"library-type",	required_argument,	0,	OPT_LIBRARY_TYPE},
{"max-deletion-length", required_argument, 0, OPT_MAX_DELETION_LENGTH},
{"max-insertion-length", required_argument, 0, OPT_MAX_INSERTION_LENGTH},
{"num-threads", required_argument, 0, OPT_NUM_CPUS},
{"zpacker", required_argument, 0, OPT_ZPACKER},
{"samtools", required_argument, 0, OPT_SAMTOOLS},
{"aux-outfile", required_argument, 0, OPT_AUX_OUT},
{"gtf-juncs", required_argument, 0, OPT_GTF_JUNCS},
{"flt-reads",required_argument, 0, OPT_FILTER_READS},
{"flt-hits",required_argument, 0, OPT_FILTER_HITS},
{0, 0, 0, 0} // terminator
};


void str_appendInt(string& str, int v) {
 stringstream ss;
 ss << v;
 str.append(ss.str());
}

bool str_endsWith(string& str, const char* suffix) {
 if (str.empty() || str.length()<3) return false;
 size_t l=strlen(suffix);
 if (str.length()<=l) return false;
 if (str.rfind(suffix, str.length()-l-1)!=string::npos) return true;
 return false;
}

int parse_options(int argc, char** argv, void (*print_usage)())
{
  int option_index = 0;
  int next_option;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
    case -1:    
      break;
    case OPT_FASTA:
      reads_format = FASTA;
      break;
    case OPT_FASTQ:
      reads_format = FASTQ;
      break;
    case OPT_MIN_ANCHOR:
      min_anchor_len = (uint32_t)parseIntOpt(3, "--min-anchor arg must be at least 3", print_usage);
      break;
    case OPT_SPLICE_MISMATCHES:
      max_splice_mismatches = parseIntOpt(0, "--splice-mismatches arg must be at least 0", print_usage);
      break; 
    case OPT_VERBOSE:
      verbose = true;
      break;
    case OPT_INSERT_LENGTH_MEAN:
      inner_dist_mean = parseIntOpt(-1024, "--inner-dist-mean arg must be at least -1024", print_usage);
      break;
    case OPT_INSERT_LENGTH_STD_DEV:
      inner_dist_std_dev = parseIntOpt(0, "--inner-dist-std-dev arg must be at least 0", print_usage);
      break;
    case OPT_OUTPUT_DIR:
      output_dir = optarg;
      break;
    case OPT_GENE_FILTER:
      gene_filter = optarg;
      break;
    case OPT_GFF_ANNOTATIONS:
      gff_file = optarg;
      break;
    case OPT_MAX_MULTIHITS:
      max_multihits = parseIntOpt(1, "--max-multihits arg must be at least 1", print_usage);
      break;
    case OPT_NO_CLOSURE_SEARCH:
      no_closure_search = true;
      break;
    case OPT_NO_COVERAGE_SEARCH:
      no_coverage_search = true;
      break;
    case OPT_NO_MICROEXON_SEARCH:
      no_microexon_search = true;
      break;
    case OPT_SEGMENT_LENGTH:
      segment_length = parseIntOpt(4, "--segment-length arg must be at least 4", print_usage);
      break;
    case OPT_SEGMENT_MISMATCHES:
      segment_mismatches = parseIntOpt(0, "--segment-mismatches arg must be at least 0", print_usage);
      break;
    case 'N':
    case OPT_READ_MISMATCHES:
      max_read_mismatches = parseIntOpt(0, "--max-mismatches arg must be at least 0", print_usage);
      break;
    case OPT_MIN_CLOSURE_EXON:
      min_closure_exon_length = parseIntOpt(1, "--min-closure-exon arg must be at least 1", print_usage);
      break;
    case OPT_MIN_CLOSURE_INTRON:
      min_closure_intron_length = parseIntOpt(1, "--min-closure-intron arg must be at least 1", print_usage);
      break;
    case OPT_MAX_CLOSURE_INTRON:
      max_closure_intron_length = parseIntOpt(1, "--max-closure-intron arg must be at least 1", print_usage);
      break;
    case OPT_MIN_COVERAGE_INTRON:
      min_coverage_intron_length = parseIntOpt(1, "--min-coverage-intron arg must be at least 1", print_usage);
      break;
    case OPT_MAX_COVERAGE_INTRON:
      max_coverage_intron_length = parseIntOpt(1, "--max-coverage-intron arg must be at least 1", print_usage);
      break;
    case OPT_MIN_SEGMENT_INTRON:
      min_segment_intron_length = parseIntOpt(1, "--min-segment-intron arg must be at least 1", print_usage);
      break;
    case OPT_MAX_SEGMENT_INTRON:
      max_segment_intron_length = parseIntOpt(1, "--max-segment-intron arg must be at least 1", print_usage);
      break;
    case OPT_MIN_REPORT_INTRON:
      min_report_intron_length = parseIntOpt(1, "--min-report-intron arg must be at least 1", print_usage);
      break;
    case OPT_MAX_REPORT_INTRON:
      max_report_intron_length = parseIntOpt(1, "--max-report-intron arg must be at least 1", print_usage);
      break;
    case OPT_MIN_ISOFORM_FRACTION:
      min_isoform_fraction = parseFloatOpt(0.0f, 1.0f, "--min-isoform-fraction arg must be [0.0,1.0]", print_usage);
      break;
    case OPT_IUM_READS:
      ium_reads = optarg;
      break;
    case OPT_SAM_HEADER:
      sam_header = optarg;
      break;
    case OPT_SAM_READGROUP_ID:
        sam_readgroup_id = optarg;
        break;
    case OPT_BUTTERFLY_SEARCH:
      butterfly_search = true;
      break;
    case OPT_SOLEXA_QUALS:
      solexa_quals = true;
      break;
    case OPT_PHRED64_QUALS:
      phred64_quals = true;
      break;
    case 'Q':
    case OPT_QUALS:
      quals = true;
      break;
    case OPT_INTEGER_QUALS:
      integer_quals = true;
      break;
    case 'C':
    case OPT_COLOR:
      color = true;
      break;
    case OPT_COLOR_OUT:
      color_out = true;
      break;
    case OPT_LIBRARY_TYPE:
      if (strcmp(optarg, "fr-unstranded") == 0)
	library_type = FR_UNSTRANDED;
      else if (strcmp(optarg, "fr-firststrand") == 0)
	library_type = FR_FIRSTSTRAND;
      else if (strcmp(optarg, "fr-secondstrand") == 0)
	library_type = FR_SECONDSTRAND;
      else if (strcmp(optarg, "ff-unstranded") == 0)
	library_type = FF_UNSTRANDED;
      else if (strcmp(optarg, "ff-firststrand") == 0)
	library_type = FF_FIRSTSTRAND;
      else if (strcmp(optarg, "ff-secondstrand") == 0)
	library_type = FF_SECONDSTRAND;
      break;
    case OPT_MAX_DELETION_LENGTH:
      max_deletion_length = parseIntOpt(0, "--max-deletion-length must be at least 0", print_usage);
      break;
    case OPT_MAX_INSERTION_LENGTH:
      max_insertion_length = parseIntOpt(0, "--max-insertion-length must be at least 0", print_usage);
      break;
    case 'z':
    case OPT_ZPACKER:
      zpacker =  optarg;
      break;
    case OPT_SAMTOOLS:
      samtools_path =  optarg;
      break;
    case OPT_AUX_OUT:
      aux_outfile =  optarg;
      break;
    case 'p':
    case OPT_NUM_CPUS:
      num_cpus=parseIntOpt(1,"-p/--num-threads must be at least 1",print_usage);
      break;
    case OPT_GTF_JUNCS:
      gtf_juncs = optarg;
      break;
    case OPT_FILTER_READS:
      flt_reads = optarg;
      break;
    case OPT_FILTER_HITS:
      flt_mappings = optarg;
      break;
    default:
      print_usage();
      return 1;
    }
  } while(next_option != -1);
  
  return 0;
}


// Error routine (prints error message and exits!)
void err_exit(const char* format,...){
  va_list arguments;
  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
#ifdef DEBUG
  // trigger a core dump for later inspection
  abort();
#endif
  exit(1);
}

FILE* FZPipe::openRead(const char* fname, string& popencmd) {
  pipecmd=popencmd;
  filename=fname;
  if (pipecmd.empty()) {
       file=fopen(filename.c_str(), "r");
       }
     else {
       string pcmd(pipecmd);
       pcmd.append(" '");
       pcmd.append(filename);
       pcmd.append("'");
       file=popen(pcmd.c_str(), "r");
       }
  return file;
}

FILE* FZPipe::openRead(const char* fname) {
  string pcmd;
  return this->openRead(fname,pcmd);
}

FILE* FZPipe::openWrite(const char* fname, string& popencmd) {
  pipecmd=popencmd;
  filename=fname;
  if (pipecmd.empty()) {
       file=fopen(filename.c_str(), "w");
       }
     else {
       string pcmd(pipecmd);
       pcmd.append(" - > '");
       pcmd.append(filename.c_str());
       pcmd.append("'");
       file=popen(pcmd.c_str(), "w");
       }
  return file;
}


FILE* FZPipe::openWrite(const char* fname) {
  string pcmd;
  return this->openWrite(fname,pcmd);
  }

void FZPipe::rewind() {
  if (pipecmd.empty()) {
      if (file!=NULL) {
           ::rewind(file);
           return;
           }
      if (!filename.empty()) {
          file=fopen(filename.c_str(),"r");
          return;
          }
      }
  if (filename.empty())
      err_die("Error: FZStream::rewind() failed (missing filename)!\n");
  this->close();
  string pcmd(pipecmd);
  pcmd.append(" '");
  pcmd.append(filename);
  pcmd.append("'");
  file=popen(pcmd.c_str(), "r");
  if (file==NULL) {
    err_die("Error: FZStream::rewind() popen(%s) failed!\n",pcmd.c_str());
    }
 }


string getFext(const string& s) {
   string r("");
   //if (xpos!=NULL) *xpos=0;
   if (s.empty() || s=="-") return r;
   int slen=(int)s.length();
   int p=s.rfind('.');
   int d=s.rfind('/');
   if (p<=0 || p>slen-2 || p<slen-7 || p<d) return r;
   r=s.substr(p+1);
   //if (xpos!=NULL) *xpos=p+1;
   for(size_t i=0; i!=r.length(); i++)
        r[i] = std::tolower(r[i]);
   return r;
   }

string guess_packer(const string& fname, bool use_all_cpus) {
   //only needed for the primary input files (given by user)
   string picmd("");
   string fext=getFext(fname);
   if (fext=="bam") {
     picmd="bam2fastx";
     return picmd;
     }
   if (fext=="gz" || fext=="gzip" || fext=="z") {
      if (use_all_cpus && str_endsWith(zpacker,"pigz")) {
           picmd=zpacker;
           if (num_cpus<2) picmd.append(" -p1");
                else {
                  picmd.append(" -p");
                  str_appendInt(picmd, num_cpus);
                  //picmd.append(" -cd");
                  }
           }
        else picmd="gzip";
      }
     else if (fext=="bz2" || fext=="bzip2" || fext=="bz" || fext=="bzip") {
      if (use_all_cpus && str_endsWith(zpacker,"pbzip2")) {
            picmd=zpacker;
            if (num_cpus<2) picmd.append(" -p1");
                 else {
                   picmd.append(" -p");
                   str_appendInt(picmd, num_cpus);
                   //picmd.append(" -cd");
                   }
            }
         else picmd="bzip2";
      }
  return picmd;
}

/*
string getBam2SamCmd(const string& fname) {
   string pipecmd("");
   string fext=getFext(fname);
   if (fext=="bam") {
      pipecmd=samtools_path;
      pipecmd.append(" view");
      }
  return pipecmd;
}
*/

void err_die(const char* format,...) { // Error exit
  va_list arguments;
  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
  exit(1);
}

string getUnpackCmd(const string& fname, bool use_all_cpus) {
 //prep_reads should use guess_packer() instead
  //otherwise compressed files MUST have the .z extension,
  //as they are all internally generated
 string pipecmd("");
 string fext=getFext(fname);
 if (fext=="bam") {
    pipecmd="bam2fastx";
    return pipecmd;
    }
 if (zpacker.empty() || fext!="z") { 
      return pipecmd; //no packer used
      }
 pipecmd=zpacker;
 if (str_endsWith(pipecmd, "pigz") ||str_endsWith(pipecmd, "pbzip2")) {
          if (use_all_cpus==false) pipecmd.append(" -p1");
              else if (num_cpus>1) {
                    pipecmd.append(" -p");
                    str_appendInt(pipecmd,num_cpus);
                    }
          }
 if (!pipecmd.empty()) pipecmd.append(" -cd");
 return pipecmd;
}

void checkSamHeader() {
  if (sam_header.empty())
    err_die("Error: writeSamHeader() with empty sam_header string\n");
  //copy the SAM header
  FILE* fh=fopen(sam_header.c_str(), "r");
  if (fh==NULL)
       err_die("Error: cannot open SAM header file %s\n",sam_header.c_str());
  fclose(fh);
}

void writeSamHeader(FILE* fout) {
  if (fout==NULL)
    err_die("Error: writeSamHeader(NULL)\n");
  checkSamHeader();
  //copy the SAM header
  FILE* fh=fopen(sam_header.c_str(), "r");
  int ch=-1;
  while ((ch=fgetc(fh))!=EOF) {
    if (fputc(ch, fout)==EOF)
          err_die("Error copying SAM header\n");
    }
  fclose(fh);
}

//auxiliary functions for BAM record handling
uint8_t* realloc_bdata(bam1_t *b, int size) {
  if (b->m_data < size) {
        b->m_data = size;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
  if (b->data_len<size) b->data_len=size;
  return b->data;
}

uint8_t* dupalloc_bdata(bam1_t *b, int size) {
  //same as realloc_bdata, but does not free previous data
  //but returns it instead
  //it ALWAYS duplicates data
  b->m_data = size;
  kroundup32(b->m_data);
  uint8_t* odata=b->data;
  b->data = (uint8_t*)malloc(b->m_data);
  memcpy((void*)b->data, (void*)odata, b->data_len);
  b->data_len=size;
  return odata; //user must FREE this after
}

extern unsigned short bam_char2flag_table[];

GBamRecord::GBamRecord(const char* qname, int32_t gseq_tid,
                 int pos, bool reverse, const char* qseq, const char* cigar, const char* quals) {
   novel=true;
   b=bam_init1();
   b->core.tid=gseq_tid;
   if (pos<=0) {
               b->core.pos=-1; //unmapped
               //if (gseq_tid<0)
               b->core.flag |= BAM_FUNMAP;
               }
          else b->core.pos=pos-1; //BAM is 0-based
   b->core.qual=255;
   int l_qseq=strlen(qseq);
   //this may not be accurate, setting CIGAR is the correct way
   //b->core.bin = bam_reg2bin(b->core.pos, b->core.pos+l_qseq-1);
   b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
   memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
   set_cigar(cigar); //this will also set core.bin
   add_sequence(qseq, l_qseq);
   add_quals(quals); //quals must be given as Phred33
   if (reverse) { b->core.flag |= BAM_FREVERSE ; }
   }

GBamRecord::GBamRecord(const char* qname, int32_t flags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals,
             const vector<string>* aux_strings) {
  novel=true;
  b=bam_init1();
  b->core.tid=g_tid;
  b->core.pos = (pos<=0) ? -1 : pos-1; //BAM is 0-based
  b->core.qual=map_qual;
  int l_qseq=strlen(qseq);
  b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
  memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
  set_cigar(cigar); //this will also set core.bin
  add_sequence(qseq, l_qseq);
  add_quals(quals); //quals must be given as Phred33
  set_flags(flags);
  set_mdata(mg_tid, (int32_t)(mate_pos-1), (int32_t)insert_size);
  if (aux_strings!=NULL) {
    for (vector<string>::const_iterator itr=aux_strings->begin();
              itr!=aux_strings->end(); ++itr) {
       add_aux(itr->c_str());
       }
    }
}
 void GBamRecord::set_cigar(const char* cigar) {
   //requires b->core.pos and b->core.flag to have been set properly PRIOR to this call
   int doff=b->core.l_qname;
   uint8_t* after_cigar=NULL;
   int after_cigar_len=0;
   uint8_t* prev_bdata=NULL;
   if (b->data_len>doff) {
      //cigar string already allocated, replace it
      int d=b->core.l_qname + b->core.n_cigar * 4;//offset of after-cigar data
      after_cigar=b->data+d;
      after_cigar_len=b->data_len-d;
      }
   const char *s;
   char *t;
   int i, op;
   long x;
   b->core.n_cigar = 0;
   if (cigar != NULL && strcmp(cigar, "*") != 0) {
        for (s = cigar; *s; ++s) {
            if (isalpha(*s)) b->core.n_cigar++;
            else if (!isdigit(*s)) {
                 err_die("Error: invalid CIGAR character (%s)\n",cigar);
                 }
            }
        if (after_cigar_len>0) { //replace/insert into existing full data
             prev_bdata=dupalloc_bdata(b, doff + b->core.n_cigar * 4 + after_cigar_len);
             memcpy((void*)(b->data+doff+b->core.n_cigar*4),(void*)after_cigar, after_cigar_len);
             free(prev_bdata);
             }
           else {
             realloc_bdata(b, doff + b->core.n_cigar * 4);
             }
        for (i = 0, s = cigar; i != b->core.n_cigar; ++i) {
            x = strtol(s, &t, 10);
            op = toupper(*t);
            if (op == 'M' || op == '=' || op == 'X') op = BAM_CMATCH;
            else if (op == 'I') op = BAM_CINS;
            else if (op == 'D') op = BAM_CDEL;
            else if (op == 'N') op = BAM_CREF_SKIP;
            else if (op == 'S') op = BAM_CSOFT_CLIP;
            else if (op == 'H') op = BAM_CHARD_CLIP;
            else if (op == 'P') op = BAM_CPAD;
            else err_die("Error: invalid CIGAR operation (%s)\n",cigar);
            s = t + 1;
            bam1_cigar(b)[i] = x << BAM_CIGAR_SHIFT | op;
        }
        if (*s) err_die("Error: unmatched CIGAR operation (%s)\n",cigar);
        b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
    } else {//no CIGAR string given
        if (!(b->core.flag&BAM_FUNMAP)) {
            fprintf(stderr, "Warning: mapped sequence without CIGAR (%s)\n", (char*)b->data);
            b->core.flag |= BAM_FUNMAP;
        }
        b->core.bin = bam_reg2bin(b->core.pos, b->core.pos + 1);
    }
   } //set_cigar()

 void GBamRecord::add_sequence(const char* qseq, int slen) {
   //must be called AFTER set_cigar (cannot replace existing sequence for now)
   if (qseq==NULL) return; //should we ever care about this?
   if (slen<0) slen=strlen(qseq);
   int doff = b->core.l_qname + b->core.n_cigar * 4;
   if (strcmp(qseq, "*")!=0) {
       b->core.l_qseq=slen;
       if (b->core.n_cigar && b->core.l_qseq != (int32_t)bam_cigar2qlen(&b->core, bam1_cigar(b)))
           err_die("Error: CIGAR and sequence length are inconsistent!(%s)\n",
                  qseq);
       uint8_t* p = (uint8_t*)realloc_bdata(b, doff + (b->core.l_qseq+1)/2 + b->core.l_qseq) + doff;
       //also allocated quals memory
       memset(p, 0, (b->core.l_qseq+1)/2);
       for (int i = 0; i < b->core.l_qseq; ++i)
           p[i/2] |= bam_nt16_table[(int)qseq[i]] << 4*(1-i%2);
       } else b->core.l_qseq = 0;
   }

 void GBamRecord::add_quals(const char* quals) {
   //requires core.l_qseq already set
   //and must be called AFTER add_sequence(), which also allocates the memory for quals
   uint8_t* p = b->data+(b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq+1)/2);
   if (quals==NULL || strcmp(quals, "*") == 0) {
      for (int i=0;i < b->core.l_qseq; i++) p[i] = 0xff;
      return;
      }
   for (int i=0;i < b->core.l_qseq; i++) p[i] = quals[i]-33;
   }

 void GBamRecord::add_aux(const char* str) {
     //requires: being called AFTER add_quals()
     static char tag[2];
     static uint8_t abuf[512];
     //requires: being called AFTER add_quals()
     int strl=strlen(str);
     //int doff = b->core.l_qname + b->core.n_cigar*4 + (b->core.l_qseq+1)/2 + b->core.l_qseq + b->l_aux;
     //int doff0=doff;
     if (strl < 6 || str[2] != ':' || str[4] != ':')
         parse_error("missing colon in auxiliary data");
     tag[0] = str[0]; tag[1] = str[1];
     uint8_t atype = str[3];
     uint8_t* adata=abuf;
     int alen=0;
     if (atype == 'A' || atype == 'a' || atype == 'c' || atype == 'C') { // c and C for backward compatibility
         atype='A';
         alen=1;
         adata=(uint8_t*)&str[5];
         }
      else if (atype == 'I' || atype == 'i') {
         long long x=(long long)atoll(str + 5);
         if (x < 0) {
             if (x >= -127) {
                 atype='c';
                 abuf[0] =  (int8_t)x;
                 alen=1;
                 }
             else if (x >= -32767) {
                 atype = 's';
                 *(int16_t*)abuf = (int16_t)x;
                 alen=2;
                 }
             else {
                 atype='i';
                 *(int32_t*)abuf = (int32_t)x;
                 alen=4;
                 if (x < -2147483648ll)
                     fprintf(stderr, "Parse warning: integer %lld is out of range.",
                             x);
                 }
             } else { //x >=0
             if (x <= 255) {
                 atype = 'C';
                 abuf[0] = (uint8_t)x;
                 alen=1;
                 }
             else if (x <= 65535) {
                 atype='S';
                 *(uint16_t*)abuf = (uint16_t)x;
                 alen=2;
                 }
             else {
                 atype='I';
                 *(uint32_t*)abuf = (uint32_t)x;
                 alen=4;
                 if (x > 4294967295ll)
                     fprintf(stderr, "Parse warning: integer %lld is out of range.",
                             x);
                 }
             }
         } //integer type
         else if (atype == 'f') {
             *(float*)abuf = (float)atof(str + 5);
             alen = sizeof(float);
             }
         else if (atype == 'd') { //?
             *(float*)abuf = (float)atof(str + 9);
             alen=8;
             }
         else if (atype == 'Z' || atype == 'H') {
             if (atype == 'H') { // check whether the hex string is valid
                 if ((strl - 5) % 2 == 1) parse_error("length of the hex string not even");
                 for (int i = 0; i < strl - 5; ++i) {
                     int c = toupper(str[5 + i]);
                     if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')))
                         parse_error("invalid hex character");
                     }
                 }
             memcpy(abuf, str + 5, strl - 5);
             abuf[strl-5] = 0;
             alen=strl-4;
             } else parse_error("unrecognized aux type");
  this->add_aux(tag, atype, alen, adata);
  }//add_aux()
