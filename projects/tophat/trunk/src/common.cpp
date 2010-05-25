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
#include <stdarg.h>
#include <getopt.h>

#include "common.h"


using namespace std;

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

//uint32_t min_intron_length = 50;
//uint32_t max_intron_length = 20000;

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

float min_isoform_fraction = 0.15f;

string output_dir = "tophat_out";
string gene_filter = "";
string gff_file = "";
string ium_reads = "";

bool solexa_quals = false;
bool phred64_quals = false;

extern void print_usage();

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
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
static float parseFloat(float lower, float upper, const char *errmsg, void (*print_usage)()) {
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



const char *short_options = "";

#define OPT_FASTA						33
#define OPT_FASTQ						34
#define OPT_MIN_ANCHOR					35
#define OPT_SPLICE_MISMATCHES			36
#define OPT_VERBOSE						37
#define OPT_INSERT_LENGTH_MEAN			38
#define OPT_INSERT_LENGTH_STD_DEV		39
//#define OPT_COV_ISLAND_EXTENSION		40
#define OPT_MIN_ISOFORM_FRACTION		41
#define OPT_OUTPUT_DIR					42
#define OPT_GENE_FILTER					43
#define OPT_GFF_ANNOTATIONS				44
#define OPT_MAX_MULTIHITS				46
#define OPT_NO_CLOSURE_SEARCH			47
#define OPT_NO_COVERAGE_SEARCH			48
#define OPT_NO_MICROEXON_SEARCH			49
#define OPT_SEGMENT_LENGTH				50
#define OPT_SEGMENT_MISMATCHES			51
#define OPT_MIN_CLOSURE_EXON			52
#define OPT_MAX_CLOSURE_INTRON			53
#define OPT_MIN_CLOSURE_INTRON			54
#define OPT_MAX_COVERAGE_INTRON			55
#define OPT_MIN_COVERAGE_INTRON			56
#define OPT_MIN_SEGMENT_INTRON			57
#define OPT_MAX_SEGMENT_INTRON			58
#define OPT_MIN_REPORT_INTRON			59
#define OPT_MAX_REPORT_INTRON			60
#define OPT_IUM_READS					61
#define OPT_BUTTERFLY_SEARCH			62
#define OPT_SOLEXA_QUALS			    63
#define OPT_PHRED64_QUALS				64

static struct option long_options[] = {
{"fasta",				no_argument,		0,	OPT_FASTA},
{"fastq",				no_argument,		0,	OPT_FASTQ},
{"min-anchor",			required_argument,	0,	OPT_MIN_ANCHOR},
{"splice-mismatches",	required_argument,	0,	OPT_SPLICE_MISMATCHES},
{"verbose",				no_argument,		0,	OPT_VERBOSE},
{"inner-dist-mean",		required_argument,	0,	OPT_INSERT_LENGTH_MEAN},
{"inner-dist-std-dev",	required_argument,	0,	OPT_INSERT_LENGTH_STD_DEV},
{"output-dir",			required_argument,	0,	OPT_OUTPUT_DIR},
{"gene-filter",			required_argument,	0,	OPT_GENE_FILTER},
{"gff-annotations",		required_argument,	0,	OPT_GFF_ANNOTATIONS},
{"max-multihits",		required_argument,	0,  OPT_MAX_MULTIHITS},
{"no-closure-search",	no_argument,		0,  OPT_NO_CLOSURE_SEARCH},
{"no-coverage-search",	no_argument,		0,  OPT_NO_COVERAGE_SEARCH},
{"no-microexon-search",	no_argument,		0,  OPT_NO_MICROEXON_SEARCH},
{"segment-length",		required_argument,	0,  OPT_SEGMENT_LENGTH},
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
{"ium-reads",			required_argument,	0,  OPT_IUM_READS},
{"butterfly-search",	no_argument,		0,	OPT_BUTTERFLY_SEARCH},
{"solexa-quals",		no_argument,		0,	OPT_SOLEXA_QUALS},
{"phred64-quals",		no_argument,		0,	OPT_PHRED64_QUALS},

{0, 0, 0, 0} // terminator
};

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
				min_anchor_len = (uint32_t)parseInt(3, "--min-anchor arg must be at least 3", print_usage);
				break;
			case OPT_SPLICE_MISMATCHES:
				max_splice_mismatches = parseInt(0, "--splice-mismatches arg must be at least 0", print_usage);
				break; 
			case OPT_VERBOSE:
				verbose = true;
				break;
			case OPT_INSERT_LENGTH_MEAN:
				inner_dist_mean = parseInt(-1024, "--inner-dist-mean arg must be at least -1024", print_usage);
				break;
			case OPT_INSERT_LENGTH_STD_DEV:
				inner_dist_std_dev = parseInt(0, "--inner-dist-std-dev arg must be at least 0", print_usage);
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
				max_multihits = parseInt(1, "--max-multihits arg must be at least 1", print_usage);
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
				segment_length = parseInt(4, "--segment-length arg must be at least 4", print_usage);
				break;
			case OPT_SEGMENT_MISMATCHES:
				segment_mismatches = parseInt(0, "--segment-mismatches arg must be at least 0", print_usage);
				break;
			case OPT_MIN_CLOSURE_EXON:
				min_closure_exon_length = parseInt(1, "--min-closure-exon arg must be at least 1", print_usage);
				break;
			case OPT_MIN_CLOSURE_INTRON:
				min_closure_intron_length = parseInt(1, "--min-closure-intron arg must be at least 1", print_usage);
				break;
			case OPT_MAX_CLOSURE_INTRON:
				max_closure_intron_length = parseInt(1, "--max-closure-intron arg must be at least 1", print_usage);
				break;
			case OPT_MIN_COVERAGE_INTRON:
				min_coverage_intron_length = parseInt(1, "--min-coverage-intron arg must be at least 1", print_usage);
				break;
			case OPT_MAX_COVERAGE_INTRON:
				max_coverage_intron_length = parseInt(1, "--max-coverage-intron arg must be at least 1", print_usage);
				break;
			case OPT_MIN_SEGMENT_INTRON:
				min_segment_intron_length = parseInt(1, "--min-segment-intron arg must be at least 1", print_usage);
				break;
			case OPT_MAX_SEGMENT_INTRON:
				max_segment_intron_length = parseInt(1, "--max-segment-intron arg must be at least 1", print_usage);
				break;
			case OPT_MIN_REPORT_INTRON:
				min_report_intron_length = parseInt(1, "--min-report-intron arg must be at least 1", print_usage);
				break;
			case OPT_MAX_REPORT_INTRON:
				max_segment_intron_length = parseInt(1, "--max-segment-intron arg must be at least 1", print_usage);
				break;
			case OPT_MIN_ISOFORM_FRACTION:
				min_isoform_fraction = parseFloat(0.0f, 1.0f, "--min-isoform-fraction arg must be [0.0,1.0]", print_usage);
				break;
			case OPT_IUM_READS:
				ium_reads = optarg;
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
