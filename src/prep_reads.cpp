/*
 *  polyA_reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/2/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *	Derived from maq "catfilter", by Heng Li at Sanger
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <cassert>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "common.h"
#include "reads.h"
#include "tokenize.h"
#include "qual.h"

//bool fastq_db = true;

using namespace std;

void format_qual_string(string& qual_str)
{
  for (size_t i = 0; i < qual_str.size(); ++i)
    {
      qual_str[i] = charToPhred33(qual_str[i], 
				      solexa_quals, 
				      phred64_quals);
    }
}

bool readmap_loaded=false;
vector<bool> readmap; //for filter_multihits


void load_readmap() {
  readmap_loaded=false;
  if (flt_reads.empty()) return;
  //FZPipe(filter_multihits, false);
  FILE* rfile=fopen(flt_reads.c_str(),"r");
  if (rfile==NULL)
     err_die("Error: cannot read file %s\n", flt_reads.c_str());
  FLineReader fr(rfile);
  skip_lines(fr); //should not be needed, bowtie doesn't add extra header lines for this
  Read read;
  while (next_fastx_read(fr, read, reads_format)) {
    uint32_t rnum=(uint32_t)atol(read.name.c_str());
    if (rnum>=(uint32_t) readmap.size())
         readmap.resize(rnum+1, false);
    readmap[rnum] = true;
    }
  readmap_loaded=true;
  fclose(rfile);
}

bool check_readmap(uint32_t read_id) {
 if (read_id>=readmap.size()) return false;
    else return readmap[read_id];
}
void flt_reads_and_hits(vector<FZPipe>& reads_files) {
  if (!readmap_loaded)
      err_die("Error: filtering reads not enabled, aborting.");
  if (aux_outfile.empty())
       err_die("Error: auxiliary output file not provided.");
  // -- filter mappings:
  string fext=getFext(flt_mappings);
  if (fext=="bam") {
    //TODO
    //WRITEME: use samtools C API to write a quick bam filter
    err_die("Error: cannot handle BAM files yet![WRITEME]");
    return;
    }
  bool is_sam=false;
  string s(flt_mappings);
  for(size_t i=0; i!=s.length(); ++i)
          s[i] = std::tolower(s[i]);
  if (fext=="sam" || s.rfind(".sam.")+fext.length()+5==s.length()) is_sam=true;
  string unzcmd=getUnpackCmd(flt_mappings);
  FZPipe hitsfile(flt_mappings, unzcmd);
  FLineReader fh(hitsfile);
  FZPipe outfile;
  outfile.openWrite(aux_outfile.c_str(), zpacker);
  if (outfile.file==NULL) err_die("Error: cannot create file %s", aux_outfile.c_str());
  const char* line;
  while ((line=fh.nextLine())) {
    if (is_sam && line[0]=='@') {
        fprintf(outfile.file, "%s\n", line);
        continue;
        }
    char* tab=strchr((char*)line, '\t');
    if (tab==NULL)
        err_die("Error: cannot find tab character in %s mappings line:\n%s",
            flt_mappings.c_str(),line);
    *tab=0;
    uint32_t rid = (uint32_t) atol(line);
    if (rid==0)
      err_die("Error: invalid read ID (%s) parsed from mapping file %s",
            line, flt_mappings.c_str());
    *tab='\t';
    if (check_readmap(rid)) {
       fprintf(outfile.file, "%s\n", line);
       }
    }
  outfile.close();
  hitsfile.close();
  // -- now filter reads
  for (size_t fi = 0; fi < reads_files.size(); ++fi) {
      //only one file expected here, this is not the initial prep_reads
      Read read;
      FLineReader fr(reads_files[fi]);
      //skip_lines(fr);
      while (next_fastx_read(fr, read)) {
        uint32_t rnum=(uint32_t)atol(read.name.c_str());
        if (check_readmap(rnum))
          printf("@%s\n%s\n+%s\n%s\n",
             read.name.c_str(),
             read.seq.c_str(),
             read.alt_name.c_str(),
             read.qual.c_str());
        }
    }
}


void process_reads(vector<FZPipe>& reads_files, vector<FZPipe>& quals_files)
{	
   //TODO: add the option to write the garbage reads into separate file(s)
  int num_reads_chucked = 0;
  int multimap_chucked=0;
  int min_read_len = 20000000;
  int max_read_len = 0;
  uint32_t next_id = 0;
  FILE* fw=NULL;
  if (!aux_outfile.empty()) {
    fw=fopen(aux_outfile.c_str(), "w");
    if (fw==NULL)
       err_die("Error: cannot create file %s\n",aux_outfile.c_str());
    }

  for (size_t fi = 0; fi < reads_files.size(); ++fi)
    {
      Read read;
      FLineReader fr(reads_files[fi]);
      skip_lines(fr);
      FZPipe fq;
      if (quals)
	      fq = quals_files[fi];
      FLineReader frq(fq);
      skip_lines(frq);
      
      while (!fr.isEof()) {
	    //read.clear();
	    // Get the next read from the file
        if (!next_fastx_read(fr, read, reads_format, ((quals) ? &frq : NULL)))
            break;

        ++next_id;
        //IMPORTANT: to keep paired reads in sync, this must be
        //incremented BEFORE any reads are chucked !

        if (read.seq.length()<12) {
            ++num_reads_chucked;
            continue;
            }
        if ((int)read.seq.length()<min_read_len)
             min_read_len=read.seq.length();
        if ((int)read.seq.length()>max_read_len)
             max_read_len=read.seq.length();

        // daehwan - check this later, it's due to bowtie
        if (color && read.seq[1] == '4') {
          ++num_reads_chucked;
          continue;
          }

        if (readmap_loaded && check_readmap(next_id)) {
          ++num_reads_chucked;
          ++multimap_chucked;
          continue;
          }
	      format_qual_string(read.qual);
        std::transform(read.seq.begin(), read.seq.end(), read.seq.begin(), ::toupper);
        char counts[256];
        memset(counts, 0, sizeof(counts));
        // Count up the bad characters
        for (unsigned int i = 0; i != read.seq.length(); ++i)
          {
            char c = (char)toupper(read.seq[i]);
            counts[(size_t)c]++;
          }

        double percent_A = (double)(counts[(size_t)'A']) / read.seq.length();
        double percent_C = (double)(counts[(size_t)'C']) / read.seq.length();
        double percent_G = (double)(counts[(size_t)'G']) / read.seq.length();
        double percent_T = (double)(counts[(size_t)'T']) / read.seq.length();
        double percent_N = (double)(counts[(size_t)'N']) / read.seq.length();
        double percent_4 = (double)(counts[(size_t)'4']) / read.seq.length();

        // Chuck the read if there are at least 5 'N's or if it's mostly
        // (>90%) 'N's and 'A's

        if (percent_A > 0.9 ||
            percent_C > 0.9 ||
            percent_G > 0.9 ||
            percent_T > 0.9 ||
            percent_N >= 0.1 ||
            percent_4 >= 0.1)
          {
            ++num_reads_chucked;
          }
        else
          {
            if (reads_format == FASTQ or (reads_format == FASTA && quals))
              {

              printf("@%u\n%s\n+%s\n%s\n",
               next_id,
               read.seq.c_str(),
               read.name.c_str(),
               read.qual.c_str());
              }
            else if (reads_format == FASTA)
              {
                string qual;
                if (color)
                  qual = string(read.seq.length()-1, 'I').c_str();
                else
                  qual = string(read.seq.length(), 'I').c_str();

                printf("@%u\n%s\n+%s\n%s\n",
                   next_id,
                   read.seq.c_str(),
                   read.name.c_str(),
                   qual.c_str());
              }
          }
      } //while !fr.isEof()
    fr.close();
    frq.close();
    } //for each input file
  fprintf(stderr, "%u out of %u reads have been filtered out\n",
	  num_reads_chucked, next_id);
  if (readmap_loaded)
    fprintf(stderr, "\t(%u filtered out due to %s)\n",
        multimap_chucked, flt_reads.c_str());
  if (fw!=NULL) {
    fprintf(fw, "min_read_len=%d\n",min_read_len - (color ? 1 : 0));
    fprintf(fw, "max_read_len=%d\n",max_read_len - (color ? 1 : 0));
    fprintf(fw, "reads_in =%d\n",next_id);
    fprintf(fw, "reads_out=%d\n",next_id-num_reads_chucked);
    fclose(fw);
    }
}

void print_usage()
{
  fprintf(stderr, "Usage:\n prep_reads [--filter-multi <multi.fq>] <reads1.fa/fq,...,readsN.fa/fq>\n");
  //alternate usage (--ctv_to_num) : doesn't filter out any reads,
  //but simply convert read names to numeric ids
}


int main(int argc, char *argv[])
{
  fprintf(stderr, "prep_reads v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
  fprintf(stderr, "---------------------------\n");
  
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string reads_file_list = argv[optind++];
  
  vector<string> reads_file_names;
  vector<FZPipe> reads_files;
  tokenize(reads_file_list, ",",reads_file_names);
  for (size_t i = 0; i < reads_file_names.size(); ++i)
    {
      string fname=reads_file_names[i];
      string pipecmd=guess_packer(fname, true);
      if (!pipecmd.empty()) pipecmd.append(" -cd");
      FZPipe seg_file(fname,pipecmd);
      if (seg_file.file == NULL) {
	  fprintf(stderr, "Error: cannot open reads file %s for reading\n",
	      reads_file_names[i].c_str());
	      exit(1);
          }
      reads_files.push_back(seg_file);
    }
  
  vector<string> quals_file_names;
  vector<FZPipe> quals_files;
  if (quals)
    {
      string quals_file_list = argv[optind++];
      tokenize(quals_file_list, ",", quals_file_names);
      for (size_t i = 0; i < quals_file_names.size(); ++i)
	{
      string fname(quals_file_names[i]);
      string pipecmd=guess_packer(fname, true);
      FZPipe seg_file(fname, pipecmd);
	  if (seg_file.file == NULL)
	    {
	      fprintf(stderr, "Error: cannot open reads file %s for reading\n",
		      quals_file_names[i].c_str());
	      exit(1);
	      }
	  quals_files.push_back(seg_file);
	}
    }
  load_readmap();
  // Only print to standard out the good reads
  //TODO: a better, more generic read filtering protocol
  if (flt_mappings.empty())
     process_reads(reads_files, quals_files);
  else //special use case: filter previous bowtie results
     flt_reads_and_hits(reads_files);
    
  return 0;
}

