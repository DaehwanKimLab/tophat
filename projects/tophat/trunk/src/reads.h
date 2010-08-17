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

void reverse_complement(string& seq);
string convert_color_to_bp(const string& color);
string convert_bp_to_color(const string& bp, bool remove_primer = false);

class ReadTable;

bool get_read_from_stream(uint64_t insert_id,
						  FILE* reads_file,
						  ReadFormat reads_format,
						  bool strip_slash,
						  char read_name [], 
						  char read_seq  [],
						  char read_alt_name [],
						  char read_qual []);

class FLineReader { //simple text line reader class, buffering last line read
  int len;
  int allocated;
  char* buf;
  bool isEOF;
  FILE* file;
  bool pushed; //pushed back
  int lcount; //counting all lines read by the object
public:
  char* chars() { return buf; }
  char* line() { return buf; }
  int readcount() { return lcount; } //number of lines read
  int length() { return len; } //length of the last line read
  bool isEof() {return isEOF; }
  char* nextLine();
  void pushBack() { if (lcount>0) pushed=true; } // "undo" the last getLine request
           // so the next call will in fact return the same line
  FLineReader(FILE* stream=NULL) {
    len=0;
    isEOF=false;
    allocated=512;
    buf=(char*)malloc(allocated);
    lcount=0;
    buf[0]=0;
    file=stream;
    pushed=false;
    }
  ~FLineReader() {
    free(buf);
    }
};

bool next_fasta_record(FLineReader& fr, string& defline, string& seq, ReadFormat reads_format);
bool next_fastq_record(FLineReader& fr, const string& seq, string& alt_name, string& qual, ReadFormat reads_format);

#endif
