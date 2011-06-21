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


void filter_garbage_reads(vector<FZPipe>& reads_files, vector<FZPipe>& quals_files)
{	

  int num_reads_chucked = 0, num_reads = 0;
  int min_read_len = 2000000;
  int max_read_len = 0;
  int next_id = 0;
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
        if (read.seq.length()<12) {
            ++num_reads_chucked;
            continue;
            }
        if ((int)read.seq.length()<min_read_len) min_read_len=read.seq.length();
        if ((int)read.seq.length()>max_read_len) max_read_len=read.seq.length();
	    format_qual_string(read.qual);
        std::transform(read.seq.begin(), read.seq.end(), read.seq.begin(), ::toupper);
        ++num_reads;
        ++next_id;
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

        if (reads_format == FASTQ &&
            ((!color && read.seq.length() != read.qual.length()) ||
             (color && read.seq.length() != read.qual.length()+1)) )
          {
            ++num_reads_chucked;
            continue;
          }

        // daehwan - check this later, it's due to bowtie
        if (color && read.seq[1] == '4') {
          continue;
          }
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
         /*   if (!fastq_db)
          {
            if (reads_format == FASTQ  or (reads_format == FASTA && quals))
              printf("@%s\n%s\n+\n%s\n",
                 read.name.c_str(), read.seq.c_str(),read.qual.c_str());
            else if (reads_format == FASTA)
              printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
          }
            else
          { */
            if (reads_format == FASTQ or (reads_format == FASTA && quals))
              {
                printf("@%d\n%s\n+%s\n%s\n",
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

                printf("@%d\n%s\n+%s\n%s\n",
                   next_id,
                   read.seq.c_str(),
                   read.name.c_str(),
                   qual.c_str());
              }
          //} only used if fastq_db is false (?)
          }
      } //while !fr.isEof()
    fr.close();
    frq.close();
    } //for each input file
  
  fprintf(stderr, "%d out of %d reads have been filtered out\n", 
	  num_reads_chucked, num_reads);
  if (fw!=NULL) {
    fprintf(fw, "min_read_len=%d\n",min_read_len - (color ? 1 : 0));
    fprintf(fw, "max_read_len=%d\n",max_read_len - (color ? 1 : 0));
    fprintf(fw, "reads_in =%d\n",num_reads);
    fprintf(fw, "reads_out=%d\n",num_reads-num_reads_chucked);
    fclose(fw);
    }
}

void print_usage()
{
  fprintf(stderr, "Usage:   prep_reads <reads1.fa/fq,...,readsN.fa/fq>\n");
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
  
  // Only print to standard out the good reads
  //TODO: a better, more generic read filtering protocol
  filter_garbage_reads(reads_files, quals_files);
    
  return 0;
}

