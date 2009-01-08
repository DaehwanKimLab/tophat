/*
 *  polyA_reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/2/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *	Derived from maq "catfilter", by Heng Li at Sanger
 */

#include <stdio.h>
#include <cassert>
#include "reads.h"

enum Format {FASTA, FASTQ};
Format format = FASTA;

void filter_garbage_reads(FILE *fa)
{
	int num_reads_chucked = 0, num_reads = 0;
	 
	Read read;
	while(!feof(fa))
	{
		read.clear();
		
		// Get the next read from the file
		if (format == FASTA)
		{
			if (!next_fasta_record(fa, read.name, read.seq))
				break;
		}
		else if (format == FASTQ)
		{
			if (!next_fastq_record(fa, read.name, read.seq, read.qual))
				break;
		}
		
		
		
		
		++num_reads;
		
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
			if (format == FASTA)
				printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
			else if (format == FASTQ)
				printf("@%s\n%s\n+\n%s\n", 
					   read.name.c_str(), read.seq.c_str(),read.qual.c_str());
					   
		}
	}
	
	fprintf(stderr, "%d out of %d reads have been filtered out\n", 
			num_reads_chucked, num_reads);
}



int main(int argc, char *argv[])
{
	FILE *fp;
	int c;
	
	// Parse command line options
	while ((c = getopt(argc, argv, "hqf")) >= 0) {
		switch (c) {
			case 'h': 
			{
				fprintf(stderr, "This program filters bad reads from stdin\n");
				return 1;
			}
			case 'q':
			{
				format = FASTQ;
				break;
			}
			case 'f':
			{
				format = FASTA;
				break;
			}
			default: break;
		}
	}

	fp = stdin;
	//fp = (strcmp(argv[optind], "-") == 0)? stdin : fopen(argv[optind], "r");
	assert(fp);
	
	// Only print to standard out the good reads
	filter_garbage_reads(fp);
	fclose(fp);
	return 0;
}

