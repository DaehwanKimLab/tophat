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
#include <vector>
#include "bam/bam.h"
#include "bam/sam.h"


#ifdef MEM_DEBUG
 void process_mem_usage(double& vm_usage, double& resident_set);
 void print_mem_usage();
#endif

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
extern std::string sam_readgroup_id;
extern std::string zpacker; //path to program to use for de/compression (gzip, pigz, bzip2, pbzip2)
extern std::string samtools_path; //path to samtools executable
extern std::string aux_outfile; //auxiliary output file name
extern bool solexa_quals;
extern bool phred64_quals;
extern bool quals;
extern bool integer_quals;
extern bool color;
extern bool color_out;

extern std::string gtf_juncs;

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
std::string getFext(const std::string& s); //returns file extension converted to lowercase
bool str_endsWith(std::string& str, const char* suffix);
void str_appendInt(std::string& str, int v);
FILE* fzOpen(std::string& fname, const char* mode);

int parseIntOpt(int lower, const char *errmsg, void (*print_usage)());
int parse_options(int argc, char** argv, void (*print_usage)());

void err_exit(const char* format,...); // exit with an error

char* get_token(char** str, const char* delims);
std::string guess_packer(const std::string& fname, bool use_all_cpus);
std::string getUnpackCmd(const std::string& fname, bool use_all_cpus=false);

void checkSamHeader();
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

void err_die(const char* format,...);

//uint8_t* realloc_bdata(bam1_t *b, int size);
//uint8_t* dupalloc_bdata(bam1_t *b, int size);

class GBamRecord {
   bam1_t* b;
   // b->data has the following strings concatenated:
   //  qname (including the terminal \0)
   //  +cigar (each event encoded on 32 bits)
   //   +seq  (4bit-encoded)
   //    +qual
   //     +aux
   bool novel;
 public:
   GBamRecord(bam1_t* from_b=NULL) {
      if (from_b==NULL) {
           b=bam_init1();
           novel=true;
           }
        else {
           b=from_b;
           novel=false;
           }
      }

     void clear() {
        if (novel) {
            bam_destroy1(b);
            }
        novel=true;
        b=bam_init1();
        }

    ~GBamRecord() {
       if (novel) {  bam_destroy1(b); }
       }

    void parse_error(const char* s) {
      err_die("BAM parsing error: %s\n", s);
      }

    bam1_t* get_b() { return b; }

    void set_mdata(int32_t mtid, int32_t m0pos, //0-based coordinate, -1 if not available
                     int32_t isize=0) { //mate info for current record
      b->core.mtid=mtid;
      b->core.mpos=m0pos; // should be -1 if '*'
      b->core.isize=isize; //should be 0 if not available
      }

    void set_flags(uint16_t flags) {
      b->core.flag=flags;
      }

    //creates a new record from 1-based alignment coordinate
    //quals should be given as Phred33
    //Warning: pos and mate_pos must be given 1-based!
    GBamRecord(const char* qname, int32_t gseq_tid,
                    int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* quals=NULL);
    GBamRecord(const char* qname, int32_t flags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals=NULL,
             const std::vector<std::string>* aux_strings=NULL);
    void set_cigar(const char* cigar); //converts and adds CIGAR string given in plain SAM text format
    void add_sequence(const char* qseq, int slen=-1); //adds the DNA sequence given in plain text format
    void add_quals(const char* quals); //quality values string in Phred33 format
    void add_aux(const char* str); //adds one aux field in plain SAM text format (e.g. "NM:i:1")
    void add_aux(const char tag[2], char atype, int len, uint8_t *data) {
      int ori_len = b->data_len;
      b->data_len += 3 + len;
      b->l_aux += 3 + len;
      if (b->m_data < b->data_len) {
        b->m_data = b->data_len;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
      }
      b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
      b->data[ori_len + 2] = atype;
      memcpy(b->data + ori_len + 3, data, len);
      }
};

class GBamWriter {
   samfile_t* bam_file;
   bam_header_t* bam_header;
 public:
   void create(const char* fname, bool uncompressed=false) {
      if (bam_header==NULL)
         err_die("Error: no bam_header for GBamWriter::create()!\n");
      if (uncompressed) {
         bam_file=samopen(fname, "wbu", bam_header);
         }
        else {
         bam_file=samopen(fname, "wb", bam_header);
         }
      if (bam_file==NULL)
         err_die("Error: could not create BAM file %s!\n",fname);
      //do we need to call bam_header_write() ?
      }

   GBamWriter(const char* fname, bam_header_t* bh, bool uncompressed=false) {
      bam_header=bh;
      create(fname, uncompressed);
      }

   GBamWriter(const char* fname, const char* samfname, bool uncompressed=false) {
      tamFile samf_in=sam_open(samfname);
      if (samf_in==NULL)
         err_die("Error: could not open SAM file %s\n", samfname);
      bam_header=sam_header_read(samf_in);
      if (bam_header==NULL)
         err_die("Error: could not read SAM header from %s!\n",samfname);
      sam_close(samf_in);
      create(fname, uncompressed);
      }

    ~GBamWriter() {
      samclose(bam_file);
      bam_header_destroy(bam_header);
      }
   bam_header_t* get_header() { return bam_header; }
   int32_t get_tid(const char *seq_name) {
      if (bam_header==NULL)
         err_die("Error: missing SAM header (get_tid())\n");
      return bam_get_tid(bam_header, seq_name);
      }

   //just a convenience function for creating a new record, but it's NOT written
   //given pos must be 1-based (so it'll be stored as pos-1 because BAM is 0-based)
   GBamRecord* new_record(const char* qname, const char* gseqname,
            int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* qual=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (bam_header->n_targets == 0) {
               err_die("Error: missing/invalid SAM header\n");
               } else
                   fprintf(stderr, "Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }

      return (new GBamRecord(qname, gseq_tid, pos, reverse, qseq, cigar, qual));
      }

   GBamRecord* new_record(const char* qname, int32_t flags, const char* gseqname,
         int pos, int map_qual, const char* cigar, const char* mgseqname, int mate_pos,
         int insert_size, const char* qseq, const char* quals=NULL,
                          const std::vector<std::string>* aux_strings=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (bam_header->n_targets == 0) {
               err_die("Error: missing/invalid SAM header\n");
               } else
                   fprintf(stderr, "Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }
      int32_t mgseq_tid=-1;
      if (mgseqname!=NULL) {
         if (strcmp(mgseqname, "=")==0) {
            mgseq_tid=gseq_tid;
            }
          else {
            mgseq_tid=get_tid(mgseqname);
            if (mgseq_tid < 0 && strcmp(mgseqname, "*")) {
                fprintf(stderr, "Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   mgseqname);
                }
            }
          }
      return (new GBamRecord(qname, flags, gseq_tid, pos, map_qual, cigar,
              mgseq_tid, mate_pos, insert_size, qseq, quals, aux_strings));
      }

   void write(GBamRecord* brec) {
      if (brec!=NULL)
          samwrite(this->bam_file,brec->get_b());
      }
};


#endif
