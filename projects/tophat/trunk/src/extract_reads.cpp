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
#include <cstring>
#include <cassert>
#include "reads.h"
#include "common.h"

using namespace std;

ReadFormat format = FASTA;

bool invert_selection = false;
bool use_alt_name = false;

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
	FLineReader fr(fa);
	//while(!feof(fa))
	while(!fr.isEof())
	{
		read.clear();
		
		// Get the next read from the file
		if (!next_fasta_record(fr, read.name, read.seq, format))
		  break;
		if (format == FASTQ)
		{
		  if (!next_fastq_record(fr, read.seq, read.alt_name, read.qual, format))
		    break;
		}
		reads_examined++;
		const string& name_selected = use_alt_name ? read.alt_name : read.name;
		bool read_in_set = selected_reads.find(name_selected) != selected_reads.end();
		if ((read_in_set && !invert_selection) || (!read_in_set && invert_selection))
		{
			++num_reads;
			if (format == FASTA)
				printf(">%s\n%s\n", read.name.c_str(), read.seq.c_str());
			else if (format == FASTQ)
				printf("@%s\n%s\n+%s\n%s\n", 
					   read.name.c_str(), read.seq.c_str(),read.alt_name.c_str(),read.qual.c_str());

		}
	}
}


void get_next_chunk(FILE* fp, 
					READSET& selected_reads, 
					unsigned int chunk_size)
{
	static int buf_size = 2048;
	char buf[buf_size];
	
	unsigned int id_count = 0;
	while (!feof(fp) && fgets(buf, buf_size, fp)) 
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
	cerr << "    -a       Select reads based on alternate name field (FASTQ only)" << endl;
	cerr << "    -r <int> select reads in ARG sized chunks to reduce peak memory." << endl;
	cerr << "    -v       invert the selection (i.e. select reads NOT in <read_ids>." << endl;
	cerr << "\nSelects reads listed in read_ids from the standard input" << endl;
}

int main(int argc, char *argv[])
{
	FILE *fp;
	int c;
	unsigned int chunk_size = 0xFFFFFFFF;
	// Parse command line options
	while ((c = getopt(argc, argv, "hqfr:va")) >= 0) {
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
			case 'a':
			{
				use_alt_name = true;
				break;
			}
			case 'f':
			{
				format = FASTA;
				break;
			}
			case 'r':
			{
				chunk_size = parseIntOpt(1,"-r arg must be at least 1", print_usage);
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
	
	string selected_reads_name = argv[optind];
	FILE* f_selected_reads = fopen(selected_reads_name.c_str(),"r");
	if (!f_selected_reads)
	{
		fprintf(stderr, "Error: could not open %s\n", selected_reads_name.c_str());
		exit(1);
	}
	
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



