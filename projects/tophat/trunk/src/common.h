#ifndef COMMON_H
#define COMMON_H
/*
 *  common.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/26/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
#include <stdint.h>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>

/*
 * Maximum allowable length of an
 * an insertion. Used mainly in
 * segment_juncs.cpp
 */
extern unsigned int max_insertion_length;

/*
 * Maximum allowable length of a
 * deletion. Used mainly in segment_juncs.cpp
 * and long_spanning_reads.cpp
 */
extern unsigned int max_deletion_length;

extern int inner_dist_mean;
extern int inner_dist_std_dev;
extern int max_mate_inner_dist; 

extern int min_anchor_len;
extern int min_report_intron_length;
extern int max_report_intron_length;

extern int min_closure_intron_length;
extern int max_closure_intron_length;

extern int min_coverage_intron_length;
extern int max_coverage_intron_length;

extern int min_segment_intron_length;
extern int max_segment_intron_length;

extern uint32_t min_closure_exon_length;
extern int island_extension;
extern int num_cpus;
extern int segment_length; // the read segment length used by the pipeline
extern int segment_mismatches;

extern int max_splice_mismatches;

enum ReadFormat {FASTA, FASTQ};
extern ReadFormat reads_format;

extern bool verbose;
extern int max_multihits;
extern bool no_closure_search;
extern bool no_coverage_search;
extern bool no_microexon_search;
extern bool butterfly_search;

extern float min_isoform_fraction;

extern std::string output_dir;
extern std::string gff_file;
extern std::string gene_filter;

extern std::string ium_reads;
extern std::string sam_header;
extern std::string zpacker; //program to use to read/write compressed files (gzip, pigz, bzip2, pbzip2)
extern bool solexa_quals;
extern bool phred64_quals;
extern bool quals;
extern bool integer_quals;
extern bool color;
extern bool color_out;

enum eLIBRARY_TYPE
  {
    LIBRARY_TYPE_NONE = 0,

    FR_UNSTRANDED,
    FR_FIRSTSTRAND,
    FR_SECONDSTRAND,

    FF_UNSTRANDED,
    FF_FIRSTSTRAND,
    FF_SECONDSTRAND,

    NUM_LIBRARY_TYPE
  };

extern eLIBRARY_TYPE library_type;

bool str_endsWith(std::string& str, const char* suffix);
void str_appendInt(std::string& str, int v);
FILE* fzOpen(std::string& fname, const char* mode);

int parseIntOpt(int lower, const char *errmsg, void (*print_usage)());
int parse_options(int argc, char** argv, void (*print_usage)());

void err_exit(const char* format,...); // exit with an error

char* get_token(char** str, const char* delims);
std::string guess_packer(const std::string& fname, bool use_all_cpus);
std::string getUnpackCmd(const std::string& fname, bool use_all_cpus=false);

void writeSamHeader(FILE* fout);

struct FZPipe {
 FILE* file;
 std::string filename;
 std::string pipecmd;
 FZPipe():filename(),pipecmd() {
   file=NULL;
   }
 FZPipe(std::string& fname, std::string& pcmd):filename(fname),pipecmd(pcmd) {
   //open as a compressed file reader
   file=NULL;
   this->openRead(fname.c_str(), pipecmd);
   }
 void close() {
   if (file!=NULL) {
     if (pipecmd.empty()) fclose(file);
                     else pclose(file);
     file=NULL;
     }
   }
 FILE* openWrite(const char* fname, std::string& popencmd);
 FILE* openWrite(const char* fname);
 FILE* openRead(const char* fname, std::string& popencmd);

 FILE* openRead(const char* fname);
 FILE* openRead(const std::string fname, std::string& popencmd) {
   return this->openRead(fname.c_str(),popencmd);
   }
 FILE* openRead(const std::string fname) {
   return this->openRead(fname.c_str());
   }
 void rewind();
};



#endif
