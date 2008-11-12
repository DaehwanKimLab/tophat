/*
 *  read_ids.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/25/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <cassert>
#include <stdio.h>
#include "reads.h"

enum Format {FASTA, FASTQ};
Format format = FASTA;

void print_read_ids(FILE *fa)
{
	
	Read read;
	int num_reads = 0;
	
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
		printf("%s\n", read.name.c_str());
	}
	
	fprintf(stderr, "Reported %d read ids\n",num_reads);
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
				fprintf(stderr, "This program extracts the FAST(A/Q) read identifiers from from stdin\n");
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
	print_read_ids(fp);
	fclose(fp);
	return 0;
}
