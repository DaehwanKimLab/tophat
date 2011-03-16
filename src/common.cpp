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
#include <stdarg.h>
#include <getopt.h>

#include "common.h"

using namespace std;

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

int max_splice_mismatches = 1;

ReadFormat reads_format = FASTQ;

bool verbose = false;

int max_multihits = 40;
bool no_closure_search = false;
bool no_coverage_search = false;
bool no_microexon_search = false;
bool butterfly_search = false;
int num_cpus = 1;
float min_isoform_fraction = 0.15f;

string output_dir = "tophat_out";
string gene_filter = "";
string gff_file = "";
string ium_reads = "";
string sam_header = "";
string zpacker = "";

bool solexa_quals = false;
bool phred64_quals = false;
bool quals = false;
bool integer_quals = false;
bool color = false;
bool color_out = false;

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


const char *short_options = "QCp:z:";

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
    OPT_QUALS,
    OPT_INTEGER_QUALS,
    OPT_COLOR,
    OPT_COLOR_OUT,
    OPT_LIBRARY_TYPE,
    OPT_MAX_DELETION_LENGTH,
    OPT_MAX_INSERTION_LENGTH,
    OPT_NUM_CPUS,
    OPT_ZPACKER
  };

static struct option long_options[] = {
{"fasta",		no_argument,		0,	OPT_FASTA},
{"fastq",		no_argument,		0,	OPT_FASTQ},
{"min-anchor",		required_argument,	0,	OPT_MIN_ANCHOR},
{"sam-header",		required_argument,	0,	OPT_SAM_HEADER},
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
    case 'p':
    case OPT_NUM_CPUS:
      num_cpus=parseIntOpt(1,"-p/--num-threads must be at least 1",print_usage);
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
  if (filename.empty()) {
        fprintf(stderr, "Error: FZStream::rewind() failed (missing filename)!\n");
        exit(1);
        }
  this->close();
  string pcmd(pipecmd);
  pcmd.append(" '");
  pcmd.append(filename);
  pcmd.append("'");
  file=popen(pcmd.c_str(), "r");
  if (file==NULL) {
    fprintf(stderr, "Error: FZStream::rewind() popen(%s) failed!\n",pcmd.c_str());
    exit(1);
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

string getUnpackCmd(const string& fname, bool use_all_cpus) {
//prep_reads should use guess_packer() instead
 string pipecmd("");
 if (zpacker.empty() || getFext(fname)!="z") {
      return pipecmd;
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

void writeSamHeader(FILE* fout) {
  if (fout==NULL) {
    fprintf(stderr, "Error: writeSamHeader(NULL)\n");
    exit(1);
    }
  if (sam_header.empty()) {
    fprintf(stderr, 
        "Error: writeSamHeader() with empty sam_header string\n");
    exit(1);
    }
  //copy the SAM header
  FILE* fh=fopen(sam_header.c_str(), "r");
  if (fh==NULL) {
       fprintf(stderr, "Error: cannot open SAM header file %s\n",
       sam_header.c_str());
       exit(1);
       }
  int ch=-1;
  while ((ch=fgetc(fh))!=EOF) {
    if (fputc(ch, fout)==EOF) {
          fprintf(stderr, "Error copying SAM header\n");
          exit(1);
          }
    }
  fclose(fh);
}
