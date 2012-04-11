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
#include <limits>
#include <getopt.h>

#include <unistd.h>
#include <ios>
#include <fstream>

using namespace std;

#include "common.h"
#include "tokenize.h"

#ifdef MEM_DEBUG
//function for debugging memory usage of current program in Linux

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

bool bowtie2 = true;
int bowtie2_min_score = -10;
int bowtie2_max_penalty = 6;
int bowtie2_min_penalty = 2;
int bowtie2_penalty_for_N = 1;
int bowtie2_read_gap_open = 5;
int bowtie2_read_gap_cont = 3;
int bowtie2_ref_gap_open = 5;
int bowtie2_ref_gap_cont = 3;

// daehwan - temporary
bool parallel = true;

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

unsigned int max_multihits = 20;
unsigned int max_seg_multihits = 40;
bool no_closure_search = false;
bool no_coverage_search = false;
bool no_microexon_search = false;
bool butterfly_search = false;
int num_threads = 1;

float min_isoform_fraction = 0.15f;

string output_dir = "tophat_out";
string std_outfile = "";
string aux_outfile = ""; //auxiliary output file name (e.g. prep_reads read stats)
string index_outfile = "";
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
string gtf_juncs = "";

bool report_secondary_alignments = false;
bool report_discordant_pair_alignments = false;

string flt_reads = "";
string flt_mappings = "";

bool fusion_search = false;
size_t fusion_anchor_length = 20;
size_t fusion_min_dist = 10000000;
size_t fusion_read_mismatches = 2;
size_t fusion_multireads = 2;
size_t fusion_multipairs = 2;
std::vector<std::string> fusion_ignore_chromosomes;
bool fusion_do_not_resolve_conflicts = false;

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

const char *short_options = "QCp:z:N:w:";

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
    OPT_MAX_SEG_MULTIHITS,
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
    OPT_LIBRARY_TYPE,
    OPT_MAX_DELETION_LENGTH,
    OPT_MAX_INSERTION_LENGTH,
    OPT_NUM_THREADS,
    OPT_ZPACKER,
    OPT_SAMTOOLS,
    OPT_AUX_OUT,
    OPT_STD_OUT,
    OPT_INDEX_OUT,
    OPT_GTF_JUNCS,
    OPT_FILTER_READS,
    OPT_FILTER_HITS,
    OPT_REPORT_SECONDARY_ALIGNMENTS,
    OPT_REPORT_DISCORDANT_PAIR_ALIGNMENTS,
    OPT_FUSION_SEARCH,
    OPT_FUSION_ANCHOR_LENGTH,
    OPT_FUSION_MIN_DIST,
    OPT_FUSION_READ_MISMATCHES,
    OPT_FUSION_MULTIREADS,
    OPT_FUSION_MULTIPAIRS,
    OPT_FUSION_IGNORE_CHROMOSOMES,
    OPT_FUSION_DO_NOT_RESOLVE_CONFLICTS,
    OPT_BOWTIE1,
    OPT_BOWTIE2_MIN_SCORE,
    OPT_BOWTIE2_MAX_PENALTY,
    OPT_BOWTIE2_MIN_PENALTY,
    OPT_BOWTIE2_PENALTY_FOR_N,
    OPT_BOWTIE2_READ_GAP_OPEN,
    OPT_BOWTIE2_READ_GAP_CONT,
    OPT_BOWTIE2_REF_GAP_OPEN,
    OPT_BOWTIE2_REF_GAP_CONT,
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
{"max-seg-multihits",	required_argument,	0,  OPT_MAX_SEG_MULTIHITS},
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
{"library-type",	required_argument,	0,	OPT_LIBRARY_TYPE},
{"max-deletion-length", required_argument, 0, OPT_MAX_DELETION_LENGTH},
{"max-insertion-length", required_argument, 0, OPT_MAX_INSERTION_LENGTH},
{"num-threads", required_argument, 0, OPT_NUM_THREADS},
{"zpacker", required_argument, 0, OPT_ZPACKER},
{"samtools", required_argument, 0, OPT_SAMTOOLS},
{"aux-outfile", required_argument, 0, OPT_AUX_OUT},
{"outfile", required_argument, 0, OPT_STD_OUT},
{"index-outfile", required_argument, 0, OPT_INDEX_OUT},
{"gtf-juncs", required_argument, 0, OPT_GTF_JUNCS},
{"flt-reads",required_argument, 0, OPT_FILTER_READS},
{"flt-hits",required_argument, 0, OPT_FILTER_HITS},
{"report-secondary-alignments", no_argument, 0, OPT_REPORT_SECONDARY_ALIGNMENTS},
{"report-discordant-pair-alignments", no_argument, 0, OPT_REPORT_DISCORDANT_PAIR_ALIGNMENTS},
{"fusion-search", no_argument, 0, OPT_FUSION_SEARCH},
{"fusion-anchor-length", required_argument, 0, OPT_FUSION_ANCHOR_LENGTH},
{"fusion-min-dist", required_argument, 0, OPT_FUSION_MIN_DIST},
{"fusion-read-mismatches", required_argument, 0, OPT_FUSION_READ_MISMATCHES},
{"fusion-multireads", required_argument, 0, OPT_FUSION_MULTIREADS},
{"fusion-multipairs", required_argument, 0, OPT_FUSION_MULTIPAIRS},
{"fusion-ignore-chromosomes", required_argument, 0, OPT_FUSION_IGNORE_CHROMOSOMES},
{"fusion-do-not-resolve-conflicts", no_argument, 0, OPT_FUSION_DO_NOT_RESOLVE_CONFLICTS},
{"bowtie1", no_argument, 0, OPT_BOWTIE1},
{"bowtie2-min-score", required_argument, 0, OPT_BOWTIE2_MIN_SCORE},
{"bowtie2-max-penalty", required_argument, 0, OPT_BOWTIE2_MAX_PENALTY},
{"bowtie2-min-penalty", required_argument, 0, OPT_BOWTIE2_MIN_PENALTY},
{"bowtie2-penalty-for-N", required_argument, 0, OPT_BOWTIE2_PENALTY_FOR_N},
{"bowtie2-read-gap-open", required_argument, 0, OPT_BOWTIE2_READ_GAP_OPEN},
{"bowtie2-read-gap-cont", required_argument, 0, OPT_BOWTIE2_READ_GAP_CONT},
{"bowtie2-ref-gap-open", required_argument, 0, OPT_BOWTIE2_REF_GAP_OPEN},
{"bowtie2-ref-gap-cont", required_argument, 0, OPT_BOWTIE2_REF_GAP_CONT},
{0, 0, 0, 0} // terminator
};


void str_appendInt(string& str, int64_t v) {
 stringstream ss;
 ss << v;
 str.append(ss.str());
}

void str_appendUInt(string& str, uint64_t v) {
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
    case OPT_MAX_SEG_MULTIHITS:
      max_seg_multihits = parseIntOpt(1, "--max-seg-multihits arg must be at least 1", print_usage);
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
    case 'w':
    case OPT_STD_OUT:
      std_outfile =  optarg;
      break;
    case OPT_INDEX_OUT:
      index_outfile =  optarg;
      break;
    case 'p':
    case OPT_NUM_THREADS:
      num_threads=parseIntOpt(1,"-p/--num-threads must be at least 1",print_usage);
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
    case OPT_REPORT_SECONDARY_ALIGNMENTS:
      report_secondary_alignments = true;
      break;
    case OPT_REPORT_DISCORDANT_PAIR_ALIGNMENTS:
      report_discordant_pair_alignments = true;
      break;
    case OPT_FUSION_SEARCH:
      fusion_search = true;
      break;
    case OPT_FUSION_ANCHOR_LENGTH:
      fusion_anchor_length = parseIntOpt(10, "--fusion-anchor-length must be at least 10", print_usage);
      break;
    case OPT_FUSION_MIN_DIST:
      fusion_min_dist = parseIntOpt(0, "--fusion-min-dist must be at least 0", print_usage);
      break;
    case OPT_FUSION_READ_MISMATCHES:
      fusion_read_mismatches = parseIntOpt(0, "--fusion-read-mismatches must be at least 0", print_usage);
      break;
    case OPT_FUSION_MULTIREADS:
      fusion_multireads = parseIntOpt(1, "--fusion-multireads must be at least 1", print_usage);
      break;
    case OPT_FUSION_MULTIPAIRS:
      fusion_multipairs = parseIntOpt(1, "--fusion-multipairs must be at least 0", print_usage);
      break;
    case OPT_FUSION_IGNORE_CHROMOSOMES:
      tokenize(optarg, ",", fusion_ignore_chromosomes);
      break;
    case OPT_FUSION_DO_NOT_RESOLVE_CONFLICTS:
      fusion_do_not_resolve_conflicts = true;
      break;
    case OPT_BOWTIE1:
      bowtie2 = false;
      break;
    case OPT_BOWTIE2_MIN_SCORE:
      bowtie2_min_score = -1 * parseIntOpt(0, "--bowtie2-min-score must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_MAX_PENALTY:
      bowtie2_max_penalty = parseIntOpt(0, "--bowtie2-max-penalty must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_MIN_PENALTY:
      bowtie2_min_penalty = parseIntOpt(0, "--bowtie2-min-penalty must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_PENALTY_FOR_N:
      bowtie2_penalty_for_N = parseIntOpt(0, "--bowtie2-penalty-for-N must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_READ_GAP_OPEN:
      bowtie2_read_gap_open = parseIntOpt(0, "--bowtie2-read-gap-open must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_READ_GAP_CONT:
      bowtie2_read_gap_cont = parseIntOpt(0, "--bowtie2-read-gap-cont must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_REF_GAP_OPEN:
      bowtie2_ref_gap_open = parseIntOpt(0, "--bowtie2-ref-gap-open must be at least 0", print_usage);
      break;
    case OPT_BOWTIE2_REF_GAP_CONT:
      bowtie2_ref_gap_cont = parseIntOpt(0, "--bowtie2-ref-gap-cont must be at least 0", print_usage);
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
	   return file;
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
  if (is_bam && !filename.empty()) {
	if (bam_file) {
	        samclose(bam_file);
	        bam_file=NULL;
	        }
	bam_file=samopen(filename.c_str(), "rb", 0);
	return;
  }
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
      err_die("Error: FZPipe::rewind() failed (missing filename)!\n");
  this->close();
  string pcmd(pipecmd);
  pcmd.append(" '");
  pcmd.append(filename);
  pcmd.append("'");
  file=popen(pcmd.c_str(), "r");
  if (file==NULL) {
    err_die("Error: FZPipe::rewind() popen(%s) failed!\n",pcmd.c_str());
    }
 }


string getFext(const string& s) {
   string r("");
   //if (xpos!=NULL) *xpos=0;
   if (s.empty() || s=="-") return r;
   int slen=(int)s.length();
   size_t pos=s.rfind('.');
   if (pos==string::npos) return r;
   int p=(int)pos;
   int d=s.rfind('/');
   if (p<=0 || p>slen-2 || p<slen-7 || p<d) return r;
   r=s.substr(p+1);
   //if (xpos!=NULL) *xpos=p+1;
   for(size_t i=0; i!=r.length(); i++)
        r[i] = std::tolower(r[i]);
   return r;
   }

string getFdir(const string& s) {
   string r("");
   //if (xpos!=NULL) *xpos=0;
   if (s.empty() || s=="-") return r;
   size_t p=s.rfind('/');
   if (p==string::npos) return r;
   r=s.substr(0,p+1);
   return r;
   }

void err_die(const char* format,...) { // Error exit
  va_list arguments;
  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
  exit(1);
}

void warn_msg(const char* format,...) {
  // print message to stderr
  va_list arguments;
  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
}

string guess_packer(const string& fname, bool use_all_cpus) {
   //only needed for the primary input files (given by user)
   string picmd("");
   string fext=getFext(fname);
   //if (fext=="bam") {
   //  picmd="bam2fastx";
   //  return picmd;
   //  }
   if (fext=="gz" || fext=="gzip" || fext=="z") {
      if (use_all_cpus && str_endsWith(zpacker,"pigz")) {
           picmd=zpacker;
           if (num_threads<2) picmd.append(" -p1");
                else {
                  picmd.append(" -p");
                  str_appendInt(picmd, num_threads);
                  //picmd.append(" -cd");
                  }
           }
        else picmd="gzip";
      }
     else if (fext=="bz2" || fext=="bzip2" || fext=="bz" || fext=="bzip") {
      if (use_all_cpus && str_endsWith(zpacker,"pbzip2")) {
            picmd=zpacker;
            if (num_threads<2) picmd.append(" -p1");
                 else {
                   picmd.append(" -p");
                   str_appendInt(picmd, num_threads);
                   //picmd.append(" -cd");
                   }
            }
         else picmd="bzip2";
      }
  return picmd;
}

string getUnpackCmd(const string& fname, bool use_all_cpus) {
 //prep_reads should use guess_packer() instead
  //otherwise compressed files MUST have the .z extension,
  //as they are all internally generated
 string pipecmd("");
 string fext=getFext(fname);
 //if (fext=="bam") {
 //   pipecmd="bam2fastx --all";
 //   return pipecmd;
 //   }
 if (zpacker.empty() || fext!="z") {
      return pipecmd; //no packer used
      }
 pipecmd=zpacker;
 if (str_endsWith(pipecmd, "pigz") ||str_endsWith(pipecmd, "pbzip2")) {
          if (use_all_cpus==false) pipecmd.append(" -p1");
              else if (num_threads>1) {
                    pipecmd.append(" -p");
                    str_appendInt(pipecmd,num_threads);
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
   if (pos<=0 || gseq_tid<0) {
               b->core.pos=-1; //unmapped
               b->core.flag |= BAM_FUNMAP;
               gseq_tid=-1;
               }
          else b->core.pos=pos-1; //BAM is 0-based
   b->core.tid=gseq_tid;
   b->core.mtid=-1;
   b->core.mpos=-1;
   b->core.qual=255;
   int l_qseq=strlen(qseq);
   //this may not be accurate, setting CIGAR is the correct way
   //b->core.bin = bam_reg2bin(b->core.pos, b->core.pos+l_qseq-1);
   b->core.l_qname=strlen(qname)+1; //includes the \0 at the end
   memcpy(realloc_bdata(b, b->core.l_qname), qname, b->core.l_qname);
   
   set_cigar(cigar); //this will also set core.bin
   add_sequence(qseq, l_qseq);
   add_quals(quals); //quals must be given as Phred33
   if (reverse) { b->core.flag |= BAM_FREVERSE; }
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
   if (cigar && strcmp(cigar, "*")) {
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
        for (i = 0, s = cigar; i != (int)b->core.n_cigar; ++i) {
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
        if (*s) err_die("Error: unmatched CIGAR operation (%s)\n", cigar);
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
     // static char tag[2];
     // static uint8_t abuf[512];
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
                     fprintf(stderr, "Parse warning: integer %lld is out of range.\n",
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
                     fprintf(stderr, "Parse warning: integer %lld is out of range.\n",
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
             alen=strl-4; //making sure the len includes the terminal \0
             } else parse_error("unrecognized aux type");
  this->add_aux(tag, atype, alen, adata);
  }//add_aux()


 uint8_t* GBamRecord::find_tag(const char tag[2]) {
   return bam_aux_get(this->b, tag);
 }

 char GBamRecord::tag_char(const char tag[2]) { //retrieve tag data as single char
   uint8_t* s=find_tag(tag);
   if (s) return ( bam_aux2A(s) );
   return 0;
  }

 int GBamRecord::tag_int(const char tag[2]) { //get the numeric value of tag
   uint8_t *s=find_tag(tag);
   if (s) return ( bam_aux2i(s) );
   return 0;
   }

 string GBamRecord::tag_str(const char tag[2]) { //return string value for a tag
   string r("");
   uint8_t *sz=find_tag(tag);
   if (sz) {
	 r = bam_aux2Z(sz);
   }
   return r;
   }

 char GBamRecord::spliceStrand() { // '+', '-' from the XS tag, or 0 if no XS tag
   char c=tag_char("XS");
   if (c) return c;
     else return '.';
   }

 string GBamRecord::sequence() {
   char *s = (char*)bam1_seq(b);
   string qseq;
   qseq.resize(b->core.l_qseq);
   for (int i=0;i<(b->core.l_qseq);i++) {
     int8_t v = bam1_seqi(s,i);
     qseq[i] = bam_nt16_rev_table[v];
     }
   return qseq;
   }

 string GBamRecord::qualities() {
   char *qual  = (char*)bam1_qual(b);
   string qv;
   qv.resize(b->core.l_qseq);
   for(int i=0;i<(b->core.l_qseq);++i) {
     qv[i]=qual[i]+33;
     }
   return qv;
   }

 string GBamRecord::seqData(string* readquals) {
   static const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
   string seq;
   string squal;
   char *qual  = (char*)bam1_qual(b);
   char *s    = (char*)bam1_seq(b);
   int i;
   //bool ismapped=((b->core.flag & BAM_FUNMAP) == 0);
   bool isreversed=((b->core.flag & BAM_FREVERSE) != 0);
   bool is_paired = ((b->core.flag & BAM_FPAIRED) != 0);
   int mate_num=0;
   if (is_paired) {
      if (b->core.flag & BAM_FREAD1)
          mate_num=1;
      else if (b->core.flag & BAM_FREAD2)
          mate_num=2;
      }
   int seqlen = b->core.l_qseq;
   if (seqlen>0) {
 	seq.resize(seqlen);
 	for(i=0;i<seqlen;i++) {
 	  seq[i] = bam1_seqi(s,i);
 	}
 	// copied from sam_view.c in samtools.
 	if (isreversed) {
 	   int l=0;
 	   int r=seqlen-1;
 	   while (l<r) {
 		  char c=seq[l];
 		  seq[l]=seq[r];
 		  seq[r]=c;
 		  l++;r--;
 		  }
 	   for (i=0;i<seqlen;i++) {
 		 seq[i]=seq_comp_table[(int)seq[i]];
 	   }
 	}

 	for(i=0;i<seqlen;i++) {
 	  seq[i] = bam_nt16_rev_table[(int)seq[i]];
 	}
  if (readquals) {
  	squal.resize(seqlen);
 	  for(i=0;i<seqlen;i++) {
 	    if (qual[i]==0xFF)
 		  squal[i]='I';
 	    else squal[i]=qual[i]+33;
 	    }

 	  if (isreversed) {
 	    int l=0;
 	    int r=seqlen-1;
 	    while (l<r) {
 		   int8_t t = squal[l];
 		   squal[l] = squal[r];
 		   squal[r] = t;
 		   l++;r--;
 	    }
 	  }
 	 *readquals = squal;
   }
  }//seqlen>0
 return seq;
 }

