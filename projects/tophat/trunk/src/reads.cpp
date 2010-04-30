/*
 *  reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/2/48.
 *  Copyright 2448 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <string>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include "reads.h"
#include "bwt_map.h"

using namespace std;

bool next_fasta_record(FILE* fp, 
					   string& defline, 
					   string& seq)
{
	static int buf_size = 2048;
	char buf[buf_size];
	
	while (!feof(fp) && fgets(buf, buf_size, fp)) 
	{
		// Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		if (buf[0] == '>')
		{
			defline = buf + 1;
			string::size_type space_pos = defline.find_first_of(" \t\r");
			if (space_pos != string::npos)
			{
				defline.resize(space_pos);
			}
			break;
		}
	}
	
	if (feof(fp))
		return false;
	
	int c;
	while ((c = fgetc(fp)) != EOF)
	{
		if (c == '>')
		{
			ungetc(c, fp);
			break;
		}
		else if (isalpha(c))
		{
			c = toupper(c);
			seq.push_back(c);
		}
	}
	
	return !(seq.empty());
}

bool next_fastq_record(FILE* fp, 
					   string& defline, 
					   string& seq, 
					   string& alt_name,
					   string& qual)
{
	static int buf_size = 2048;
	char buf[buf_size];
	defline.clear();
	seq.clear();
	alt_name.clear();
	qual.clear();
	// Put the name of the record into defline
	while (!feof(fp) && fgets(buf, buf_size, fp)) 
	{
		// Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		// Is this line the start of a new record?
		if (buf[0] == '@')
		{
			defline = buf + 1;
			string::size_type space_pos = defline.find_first_of(" \t\r");
			if (space_pos != string::npos)
			{
				defline.resize(space_pos);
			}
			break;
		}
	}
	if (feof(fp))
		return false;
	
	// Put the sequece for the record in seq
	int c;
	while ((c = fgetc(fp)) != EOF)
	{
		// We might have multiple newlines after the sequence, but the "+" line
		// means we've reached the end of the sequence part of the record 
		if (c == '+')
		{
			ungetc(c, fp);
			break;
		}
		else if (isalpha(c))
		{
			seq.push_back(c);
		}
	}
	
	if (c == EOF)
		return false;
	
	while (!feof(fp) && fgets(buf, buf_size, fp))
	{
		// Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		if (buf[0] == '+')
		{
			alt_name = buf + 1;
			break;
		}	
	}
	
	if (feof(fp))
		return false;
	
	// Put the qualities of the record into qual
	while ((c = fgetc(fp)) != EOF)
	{
		// If we hit a "@" char, and we already have enough qualities for the 
		// sequence string we read, then this is the start of a new record.
		if (c == '@' && qual.length() >= seq.length())
		{
			ungetc(c, fp);
			break;
		}
		// Only accept quality chars in the Sanger-scaled (but printable) range
		//else if (c >= '!' && c <= 'I')
		//{
		if (qual.length() < seq.length())
			qual.push_back(c);
		//}
	}
	
	return true;
}

// This could be faster.
void reverse_complement(string& seq)
{
	//fprintf(stderr,"fwd: %s\n", seq.c_str());
	for (string::size_type i = 0; i < seq.length(); ++i)
	{
		switch(seq[i])
		{
			case 'A' : seq[i] = 'T'; break;
			case 'T' : seq[i] = 'A'; break;
			case 'C' : seq[i] = 'G'; break;
			case 'G' : seq[i] = 'C'; break;
			default: seq[i]   = 'N'; break;
		}
	}
	reverse(seq.begin(), seq.end());
	//fprintf(stderr, "rev: %s\n", seq.c_str());
}

bool get_read_from_stream(uint64_t insert_id,
						  FILE* reads_file,
						  ReadFormat reads_format,
						  bool strip_slash,
						  char read_name [], 
						  char read_seq  [],
						  char read_alt_name [], 
						  char read_qual [])
{
	
	Read read;
	while(!feof(reads_file))
	{
		read.clear();
		
		// Get the next read from the file
		if (reads_format == FASTA)
		{
			if (!next_fasta_record(reads_file, read.name, read.seq))
				break;
		}
		else if (reads_format == FASTQ)
		{
			if (!next_fastq_record(reads_file, read.name, read.seq, read.alt_name, read.qual))
				break;
		}
		
		
		if (strip_slash)
		{
			string::size_type slash = read.name.rfind("/");
			if (slash != string::npos)
				read.name.resize(slash);
		}
		
		if ((uint64_t)atoi(read.name.c_str()) == insert_id)
		{
			if (read_name) strcpy(read_name, read.name.c_str());
			if (read_seq) strcpy(read_seq, read.seq.c_str());
			if (read_alt_name) strcpy(read_alt_name, read.alt_name.c_str());
			if (read_qual) strcpy(read_qual, read.qual.c_str());
			return true;
		}
		
        //rt.get_id(read.name, ref_str);
    }	
	
	return false;
}


