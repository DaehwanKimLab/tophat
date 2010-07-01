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

bool fastq_db = true;
using namespace std;

void format_qual_string(const string& orig_qual_str,
						string& out_qual_str)
{
	out_qual_str = orig_qual_str;
	for (size_t i = 0; i < orig_qual_str.size(); ++i)
	{
		out_qual_str[i] = charToPhred33(orig_qual_str[i], 
										solexa_quals, 
										phred64_quals);
	}
}


void filter_garbage_reads(vector<FILE*> reads_files)
{	
	int num_reads_chucked = 0, num_reads = 0;
	int next_id = 0;
	for (size_t fi = 0; fi < reads_files.size(); ++fi)
	{
		Read read;
		FILE* fa = reads_files[fi];
		FLineReader fr(fa);
		//while(!feof(fa))
		while (!fr.isEof())
		{
			read.clear();
			
			// Get the next read from the file
			if (reads_format == FASTA)
			{
				if (!next_fasta_record(fr, read.name, read.seq))
					break;
			}
			else if (reads_format == FASTQ)
			{
				string orig_qual;
				if (!next_fastq_record(fr, read.name, read.seq, read.alt_name, orig_qual))
					break;
				format_qual_string(orig_qual, read.qual);
			}
			
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
			
			if (reads_format == FASTQ &&
				read.seq.length() != read.qual.length())
			{
				++num_reads_chucked;
				continue;
			}
			
			// Chuck the read if there are at least 5 'N's or if it's mostly
			// (>90%) 'N's and 'A's
			
			if (percent_A > 0.9 ||
				percent_C > 0.9 ||
				percent_G > 0.9 ||
				percent_T > 0.9 ||
				percent_N >= 0.1)
			{
				++num_reads_chucked;
			} 
			else
			{
				if (!fastq_db)
				{
					if (reads_format == FASTA)
						printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
					else if (reads_format == FASTQ)
						printf("@%s\n%s\n+\n%s\n", 
							   read.name.c_str(), read.seq.c_str(),read.qual.c_str());
				}
				else
				{
					if (reads_format == FASTA)
					{
						printf("@%d\n%s\n+%s\n%s\n",
							   next_id,
							   read.seq.c_str(),
							   read.name.c_str(),
							   string(read.seq.length(), 'I').c_str());
					}
					else if (reads_format == FASTQ)
					{
						printf("@%d\n%s\n+%s\n%s\n",
							   next_id,
							   read.seq.c_str(),
							   read.name.c_str(),
							   read.qual.c_str());
					}
					
				}
			}
		}
	}
	
	fprintf(stderr, "%d out of %d reads have been filtered out\n", 
			num_reads_chucked, num_reads);
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
    vector<FILE*> reads_files;
    tokenize(reads_file_list, ",",reads_file_names);
    for (size_t i = 0; i < reads_file_names.size(); ++i)
    {
        FILE* seg_file = fopen(reads_file_names[i].c_str(), "r");
        if (seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                    reads_file_names[i].c_str());
            exit(1);
        }
        reads_files.push_back(seg_file);
    }
	
	// Only print to standard out the good reads
    //TODO: a better, more generic read filtering protocol
	filter_garbage_reads(reads_files);
	
	return 0;
}

