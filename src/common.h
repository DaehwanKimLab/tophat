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


#define MAX_READ_LEN 1024

//for fastq and bam indexing by read# (multi-threading)
#define INDEX_REC_COUNT 1000

#define VMAXINT32 0xFFFFFFFF

#ifdef MEM_DEBUG
 void process_mem_usage(double& vm_usage, double& resident_set);
 void print_mem_usage();
#endif

extern bool bowtie2;
extern int bowtie2_min_score;
extern int bowtie2_max_penalty;
extern int bowtie2_min_penalty;
extern int bowtie2_penalty_for_N;
extern int bowtie2_read_gap_open;
extern int bowtie2_read_gap_cont;
extern int bowtie2_ref_gap_open;
extern int bowtie2_ref_gap_cont;

// daehwan - temporary for parallelization
extern bool parallel;

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
extern int num_threads;
extern int segment_length; // the read segment length used by the pipeline
extern int segment_mismatches;
extern int max_read_mismatches;

extern int max_splice_mismatches;

enum ReadFormat {FASTA, FASTQ};
extern ReadFormat reads_format;

extern bool verbose;
extern unsigned int max_multihits;
extern unsigned int max_seg_multihits;
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
extern std::string std_outfile; //main output file that some modules can use instead of stdout
extern std::string aux_outfile; //auxiliary output file name
extern std::string index_outfile; //index output file name
extern bool solexa_quals;
extern bool phred64_quals;
extern bool quals;
extern bool integer_quals;
extern bool color;
extern std::string gtf_juncs;

extern bool report_secondary_alignments;
extern bool report_discordant_pair_alignments;

//prep_reads only: --flt-reads <bowtie-fastq_for--max>
//  filter out reads if their numeric ID is in this fastq file
// OR if flt_mappings was given too, filter out reads if their ID
// is NOT in this fastq file
extern std::string flt_reads;

//prep_reads special usage: filter out mappings whose read ID
//is NOT found in the flt_reads file, and write them into
// aux_outfile; also reverses the flt_reads filter itself
extern std::string flt_mappings;

extern bool fusion_search;
extern size_t fusion_anchor_length;
extern size_t fusion_min_dist;
extern size_t fusion_read_mismatches;
extern size_t fusion_multireads;
extern size_t fusion_multipairs;
extern std::vector<std::string> fusion_ignore_chromosomes;
extern bool fusion_do_not_resolve_conflicts;

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
std::string getFdir(const std::string& s); //returns file extension converted to lowercase
bool str_endsWith(std::string& str, const char* suffix);
void str_appendInt(std::string& str, int64_t v);
void str_appendUInt(std::string& str, uint64_t v);
FILE* fzOpen(std::string& fname, const char* mode);

int parseIntOpt(int lower, const char *errmsg, void (*print_usage)());
int parse_options(int argc, char** argv, void (*print_usage)());

void err_exit(const char* format,...); // exit with an error

char* get_token(char** str, const char* delims);
std::string guess_packer(const std::string& fname, bool use_all_cpus=false);
std::string getUnpackCmd(const std::string& fname, bool use_all_cpus=false);

void checkSamHeader();
void writeSamHeader(FILE* fout);

class FZPipe {
 public:
	 union {
	   FILE* file;
	   samfile_t* bam_file;
	 };
	 std::string filename;
	 std::string pipecmd;
	 bool is_bam;
	 FZPipe(const std::string& fname, bool guess=false):filename(fname),pipecmd() {
	   //this constructor is only to use FZPipe as a READER
       //also accepts/recognizes BAM files (without needing pipes)
	   //for which it only stores the filename
	   openRead(fname, guess);
	   }

   FILE* openRead(const std::string& fname, bool guess=false) {
     filename=fname;
     pipecmd="";
     is_bam=false;
     if (getFext(fname) == "bam") {
           //file=(FILE*)this; //just to be non-NULL;
           is_bam=true;
           bam_file=samopen(filename.c_str(), "rb", 0);
           return file;
           }
     if (guess) {
       pipecmd=guess_packer(fname);
       if (!pipecmd.empty()) pipecmd.append(" -cd");
     }
     else {
       pipecmd=getUnpackCmd(fname);
     }
     if (pipecmd.empty()) {
           //this->openRead(fname.c_str());
           file=fopen(filename.c_str(), "r");
           return file;
           }
     this->openRead(fname.c_str(), pipecmd);
     return file;
     }

     FILE* openRead(const char* fname, std::string& popencmd);

	 FILE* openRead(const char* fname) {
	       std::string s(fname);
	       return openRead(s);
	       }
	 FILE* openRead(const std::string fname, std::string& popencmd) {
	   return this->openRead(fname.c_str(),popencmd);
	   }


	 FZPipe():filename(),pipecmd() {
	   is_bam=false;
	   file=NULL;
	   }
	 FZPipe(std::string& fname, std::string& pcmd):filename(fname),pipecmd(pcmd) {
	   //open as a compressed file reader
	   if (pipecmd.empty()) {
		 this->openRead(fname);
		 return;
	   }
	   is_bam=false;
	   file=NULL;
	   this->openRead(fname.c_str(), pipecmd);
	   }

	 void close() {
	   if (file!=NULL) {
		 if (is_bam) {
	        samclose(bam_file);
	        bam_file=NULL;
		 }
		 else {
		   if (pipecmd.empty()) fclose(file);
						 else pclose(file);
		   }
		 file=NULL;
		 }
	   }
	 FILE* openWrite(const char* fname, std::string& popencmd);
	 FILE* openWrite(const char* fname);

	 void rewind();

	 void seek(int64_t foffset)
	 {
	   if (is_bam) {
		 bgzf_seek(bam_file->x.bam, foffset, SEEK_SET);
	   }
	   else {
		 fseek(this->file, foffset, SEEK_SET);
	   }
	 }
};

void err_die(const char* format,...);
void warn_msg(const char* format,...);

class GBamRecord {
   bam1_t* b;
   // b->data has the following strings concatenated:
   //  qname (including the terminal \0)
   //  +cigar (each event encoded on 32 bits)
   //   +seq  (4bit-encoded)
   //    +qual
   //     +aux
   bool novel;

   char tag[2];
   uint8_t abuf[1024];
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

    void set_flag(uint16_t flag) { //use BAM_F* constants
      b->core.flag |= flag;
    }

    void unset_flag(uint16_t flag) { //use BAM_F* constants
      b->core.flag &= ~flag;
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
      //IMPORTANT:  strings (Z,H) should include the terminal \0
      int addz=0;
      if ((atype=='Z' || atype=='H') && data[len-1]!=0) {
    	addz=1;
        }
      int ori_len = b->data_len;
      b->data_len += 3 + len + addz;
      b->l_aux += 3 + len + addz;
      if (b->m_data < b->data_len) {
        b->m_data = b->data_len;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
      }
      b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
      b->data[ori_len + 2] = atype;
      if (addz) {
    	b->data[ori_len+len+3]=0;
        }
      memcpy(b->data + ori_len + 3, data, len);
      }
   //--reading back aux tags:
	 uint8_t* find_tag(const char tag[2]);
	 //returns pointer at the beginning of tag data, or NULL if tag not found
	 //the returned pointer can then be used by bam_aux2*() functions

	 std::string tag_str(const char tag[2]); //return tag value for tag type 'Z'
	 int tag_int(const char tag[2]); //return numeric value of tag (for numeric types)
	 char tag_char(const char tag[2]); //return char value of tag (for type 'A')
	 char spliceStrand(); // '+', '-' from the XS tag, or '.' if no XS tag
	//--
	 std::string qualities(); //return quality string, as is (ignores BAM_FREVERSE)
	 std::string sequence();  //return read sequence as is (ignores BAM_FREVERSE)
	 std::string seqData(std::string* readquals=NULL); //return seq and qv, reversed if BAM_FREVERSE
};

class GBamWriter {
   samfile_t* bam_file;
   bam_header_t* bam_header;
   FILE* findex;
   uint64_t wcount;
   uint64_t idxcount;
   int64_t idx_last_id;
   bool external_header;
 public:
   void create(const char* fname, bool uncompressed=false) {
	  findex=NULL;
	  wcount=0;
	  idxcount=0;
	  idx_last_id=0;
	  external_header=false;
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

   void create(const char* fname, std::string& idxfile) {
	  findex=NULL;
	  wcount=0;
	  idxcount=0;
	  idx_last_id=0;
	  external_header=false;
      if (bam_header==NULL)
         err_die("Error: no bam_header for GBamWriter::create()!\n");
      bam_file=samopen(fname, "wb", bam_header);
      if (bam_file==NULL)
         err_die("Error: could not create BAM file %s!\n",fname);
      if (!idxfile.empty()) {
    	   findex = fopen(idxfile.c_str(), "w");
    	   if (findex == NULL)
    	     err_die("Error: cannot create file %s\n", idxfile.c_str());
           }
      }

   void create(const char* fname, bam_header_t* bh, bool uncompressed=false) {
     findex=NULL;
	 wcount=0;
	 idxcount=0;
	 idx_last_id=0;
     external_header=false;
	 bam_header=bh;
	 create(fname, uncompressed);
     }

   GBamWriter(const char* fname, bam_header_t* bh, bool uncompressed=false) {
      create(fname, bh, uncompressed);
      external_header=true;
      }

   GBamWriter(const char* fname, bam_header_t* bh, std::string& idxfile) {
	  bam_header=bh;
      create(fname, idxfile);
      external_header=true;
      }

   GBamWriter(std::string& fname, std::string& idxfile) {
	 //create BAM with empty header
     external_header=false;
	 bam_header=bam_header_init();
	 create(fname.c_str());
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

   GBamWriter(const char* fname, const char* samfname, std::string idxfile) {
      tamFile samf_in=sam_open(samfname);
      if (samf_in==NULL)
         err_die("Error: could not open SAM file %s\n", samfname);
      bam_header=sam_header_read(samf_in);
      if (bam_header==NULL)
         err_die("Error: could not read SAM header from %s!\n",samfname);
      sam_close(samf_in);
      create(fname, idxfile);
      }

    ~GBamWriter() {
      samclose(bam_file);
      if (bam_header && !external_header)
        bam_header_destroy(bam_header);
      if (findex != NULL)
        fclose(findex);
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
	 if (gseqname==NULL || strcmp(gseqname, "*")==0) {
	   //probably an unmapped read
	   //if (pos>0) err_die("Error: genomic position given for unmapped read!\n");
	   return (new GBamRecord(qname, -1, 0, false, qseq, cigar, qual));
	  }
	 else {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0) {
            if (bam_header->n_targets == 0) {
               err_die("Error: missing/invalid SAM header\n");
               } else
                   fprintf(stderr, "Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }
      return (new GBamRecord(qname, gseq_tid, pos, reverse, qseq, cigar, qual));
	  }
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
      if (brec!=NULL) {
        samwrite(this->bam_file,brec->get_b());
        wcount++;
      }

      }
   void write(bam1_t* b, int64_t read_id=0) {
     int64_t pre_block_addr=0; //offsets after last write()
     int     pre_block_offs=0; //but before this write()
     int64_t pre_pos=0;
     bool write_index=false;
     if (findex && read_id) {
       if (idxcount >= INDEX_REC_COUNT && read_id != idx_last_id) {
	 pre_pos = this->tell();
	 pre_block_offs = pre_pos & 0xFFFF;
	 pre_block_addr = (pre_pos >> 16) & 0xFFFFFFFFFFFFLL;
	 write_index=true;
       }
       idx_last_id=read_id;
       idxcount++;
     }
     
     samwrite(this->bam_file, b);
     wcount++;
     if (write_index) {
       int64_t offset = this->tell();
       int     post_block_offs = offset & 0xFFFF; //offsets after this write()
       int64_t post_block_addr = (offset >> 16) & 0xFFFFFFFFFFFFLL;
       int data_len = b->data_len+BAM_CORE_SIZE;
       if (post_block_addr != pre_block_addr &&
	   post_block_offs>=data_len)
	 //all data written in this block
	 //WARNING: this check fails for very large BAM records (> 64K)
	 {
	   pre_pos = post_block_addr << 16;
	 }

       fprintf(findex, "%lu\t%ld\n", read_id, pre_pos);
       idxcount = 0;
     }
   }
   

   int64_t tell() {
     return bam_tell(this->bam_file->x.bam);
   }


   int64_t writtenCount() { return wcount; }

   void flush() {
     bgzf_flush(this->bam_file->x.bam);
   }

   void seek(int64_t offset) {
     bam_seek(this->bam_file->x.bam, offset, SEEK_SET);
   }
};

#endif
