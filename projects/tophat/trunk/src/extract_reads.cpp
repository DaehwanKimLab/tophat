/*
 *  extract_reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/25/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <stdio.h>
#include <iostream>
#include <set>
#include <algorithm>
#include <string>
#include <cassert>
#include "reads.h"

using namespace std;

enum Format {FASTA, FASTQ};
Format format = FASTA;

bool invert_selection = false;

struct r_to_l_lexcompare
{
	bool operator()(const string& s1, const string& s2) const
	{
		return lexicographical_compare(s1.rbegin(), s1.rend(),
									   s2.rbegin(), s2.rend());
	}
};


typedef set<string, r_to_l_lexcompare> READSET;

void extract_reads(FILE *fa, const READSET& selected_reads)
{
	int num_reads = 0;
	int reads_examined = 0;
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
		reads_examined++;
		bool read_in_set = selected_reads.find(read.name) != selected_reads.end();
		if (read_in_set || (!read_in_set && invert_selection))
		{
			++num_reads;
			if (format == FASTA)
				printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
			else if (format == FASTQ)
				printf("@%s\n%s\n+\n%s\n", 
					   read.name.c_str(), read.seq.c_str(),read.qual.c_str());

		}
	}
}


void get_next_chunk(FILE* fp, 
					READSET& selected_reads, 
					unsigned int chunk_size)
{
	char buf[2048];
	
	unsigned int id_count = 0;
	while (!feof(fp) && fgets(buf, sizeof(buf), fp)) 
	{
		// Chomp trailing newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		// Chomp leading whitespace
		char* p = buf;
		while(*p && isspace(*p)) ++p;
		if (!*p)
			continue;
		 
		selected_reads.insert(p);
		if (++id_count >= chunk_size)
			break;
	}
	fprintf(stderr, "chunk contains %d read ids\n", id_count);
}


static void print_usage() 
{
	cerr << "Usage: extract_reads [options] read_ids < reads.f*" << endl;
	cerr << "    -f       reads are in FASTA format" << endl;
	cerr << "    -q       reads are in FASTQ format" << endl;
	cerr << "    -r <int> select reads in ARG sized chunks to reduce peak memory." << endl;
	cerr << "    -v       invert the selection (i.e. select reads NOT in <read_ids>." << endl;
	cerr << "\nSelects reads listed in read_ids from the standard input" << endl;
}



/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			fprintf(stderr,"%s\n",errmsg);
			print_usage();
			exit(1);
		}
		return (int32_t)l;
	}
	fprintf(stderr,"%s\n",errmsg);
	print_usage();
	exit(1);
	return -1;
}


int main(int argc, char *argv[])
{
	FILE *fp;
	int c;
	unsigned int chunk_size = 0xFFFFFFFF;
	// Parse command line options
	while ((c = getopt(argc, argv, "hqfr:")) >= 0) {
		switch (c) {
			case 'h': 
			{
				print_usage();
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
			case 'r':
			{
				chunk_size = parseInt(1,"-r arg must be at least 1");
				break;
			}
			case 'v':
			{
				invert_selection = true;
				break;
			}
			default: break;
		}
	}
	
	if(optind >= argc) {
		print_usage();
		return 1;
	} 
	
	FILE* f_selected_reads = fopen(argv[optind],"r");
	
	fp = stdin;
	//fp = (strcmp(argv[optind], "-") == 0)? stdin : fopen(argv[optind], "r");
	assert(fp);
	
	READSET selected_reads;
	//unsigned int chunk_number = 1;
	do {
		selected_reads.clear();
		rewind(fp);
		
		//fprintf(stderr,"Starting chunk # %d\n", chunk_number++);
		get_next_chunk(f_selected_reads, selected_reads, chunk_size);
		// Only print to standard out the good reads
		extract_reads(fp, selected_reads);
		
	}while(!feof(f_selected_reads));

	fclose(fp);
	return 0;
}



