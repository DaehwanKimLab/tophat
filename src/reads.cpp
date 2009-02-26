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
#include <cstring>
#include "reads.h"

using namespace std;

bool next_fasta_record(FILE* fp, 
					   string& defline, 
					   string& seq)
{
	char buf[2048];
	
	while (!feof(fp) && fgets(buf, sizeof(buf), fp)) 
	{
		// Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		if (buf[0] == '>')
		{
			defline = buf + 1;
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

bool next_fastq_record(FILE* fp, string& defline, string& seq, string& qual)
{
	char buf[2048];
	
	// Put the name of the record into defline
	while (!feof(fp) && fgets(buf, sizeof(buf), fp)) 
	{
		// Chomp the newline
		char* nl = strrchr(buf, '\n');
		if (nl) *nl = 0;
		
		// Is this line the start of a new record?
		if (buf[0] == '@')
		{
			defline = buf + 1;
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
	
	// Discard the optional secondary name (the one on the "+" line)
	while (!feof(fp) && fgets(buf, sizeof(buf), fp))
	{
		if (buf[0] == '+')
		{
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


