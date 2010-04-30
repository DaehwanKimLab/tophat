#ifndef READS_H
#define READS_H
/*
 *  reads.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/2/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <string>
#include "common.h"

using std::string;

static const int max_read_bp = 64;

// Note: qualities are not currently used by TopHat
struct Read
{
	Read() 
	{
		seq.reserve(max_read_bp);
		qual.reserve(max_read_bp);
	}
	
	string name;
	string seq;
	string alt_name;
	string qual;
	
	bool lengths_equal() { return seq.length() == qual.length(); }
	void clear() 
	{ 
		name.clear(); 
		seq.clear(); 
		qual.clear(); 
		alt_name.clear();
	}
};

bool next_fasta_record(FILE* fp, string& defline, string& seq);
bool next_fastq_record(FILE* fp, string& defline, string& seq, string& alt_name, string& qual);
void reverse_complement(string& seq);

class ReadTable;

bool get_read_from_stream(uint64_t insert_id,
						  FILE* reads_file,
						  ReadFormat reads_format,
						  bool strip_slash,
						  char read_name [], 
						  char read_seq  [],
						  char read_alt_name [],
						  char read_qual []);

#endif
