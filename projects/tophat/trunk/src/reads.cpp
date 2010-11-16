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
#include "tokenize.h"

using namespace std;

char* FLineReader::nextLine() {
   if(!file) return NULL;
   if (pushed) { pushed=false; return buf; }
   //reads a char at a time until \n and/or \r are encountered
   len=0;
   int c=0;
   while ((c=getc(file))!=EOF) {
     if (len>=allocated-1) {
        allocated+=512;
        buf=(char*)realloc(buf,allocated);
     }
     if (c=='\n' || c=='\r') {
       buf[len]='\0';
       if (c=='\r') { //DOS file: double-char line terminator, skip the second one
          if ((c=getc(file))!='\n')
              ungetc(c,file); //this will always happen on Mac
          }
       lcount++;
       return buf;
       }
     buf[len]=(char)c;
     len++;
     }
   if (c==EOF) {
     isEOF=true;
     if (len==0) return NULL;
     }
   buf[len]='\0';
   lcount++;
   return buf;
}

void skip_lines(FLineReader& fr)
{
  char* buf = NULL;
  while ((buf = fr.nextLine()) != NULL) {
    if (buf[0] == '\0') continue;
    if (buf[0] == '>' || buf[0] == '@')
      {
	fr.pushBack();
	break;
      }
  }
}

bool next_fasta_record(FLineReader& fr,
		       string& defline, 
		       string& seq,
		       ReadFormat reads_format)

{
  seq.clear();
  defline.clear();
  char* buf=NULL;
  while ((buf=fr.nextLine())!=NULL) {
    if (buf[0]==0) continue; //skip empty lines
    if ((reads_format == FASTA && buf[0] == '>') || (reads_format == FASTQ && (buf[0] == '+' || buf[0] == '@'))) { //next record
        if (seq.length()>0) { //current record ending
           fr.pushBack();
           return true;
           }
        defline=buf+1;
        string::size_type space_pos = defline.find_first_of(" \t");
        if (space_pos != string::npos) {
            defline.resize(space_pos);
            }
        continue;
        } //defline
    // sequence line
    seq.append(buf);
    } //line reading loop

    replace(seq.begin(),seq.end(),'.','N'); //shouldn't really be needed for FASTA files
	return !(seq.empty());
}

bool next_fastq_record(FLineReader& fr,
		       const string& seq,
		       string& alt_name,
		       string& qual,
		       ReadFormat reads_format)
{
  alt_name.clear();
  qual.clear();
  char* fline=fr.nextLine();
  if (fline==NULL) return false;
  while (fline[0]==0) { //skip empty lines
    fline=fr.nextLine();
    if (fline==NULL) return false;
    }
  //must be on '+' line here
  if (fline==NULL || (reads_format == FASTQ && fline[0] != '+') || (reads_format == FASTA && quals && fline[0] != '>')) {
     err_exit("Error: '+' not found for fastq record %s\n",fline);
     return false;
     }
  alt_name=fline+1;
  string::size_type space_pos = alt_name.find_first_of(" \t");
  if (space_pos != string::npos) alt_name.resize(space_pos);
   //read qv line(s) now:
  while ((fline=fr.nextLine())!=NULL) {
    if (integer_quals)
      {
	vector<string> integer_qual_values;
	tokenize(string(fline), " ", integer_qual_values);

	string temp_qual;
	for (size_t i = 0; i < integer_qual_values.size(); ++i)
	  {
	    int qual_value = atoi(integer_qual_values[i].c_str());
	    if (qual_value < 0) qual_value = 0;
	    temp_qual.push_back((char)(qual_value + 33));
	  }

	qual.append(temp_qual);
      }
    else
      qual.append(fline);
    
    if ((!color && qual.length()>=seq.length()) || (color && qual.length()+1>=seq.length())) break;
     }
  // final check
  if ((!color && seq.length()!=qual.length()) || (color && seq.length()!=qual.length()+1)) {
           err_exit("Error: qual length (%d) differs from seq length (%d) for fastq record %s!\n",
               qual.length(), seq.length(), alt_name.c_str());
           return false;
           }
  //

  return !(qual.empty());
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

string convert_color_to_bp(const string& color)
{
  if (color.length() <= 0)
    return "";

  char base = color[0];
  string bp;
  for (string::size_type i = 1; i < color.length(); ++i)
    {
      char next = color[i];
      switch(base)
	{
	  // 'A0':'A', 'A1':'C', 'A2':'G', 'A3':'T', 'A4':'N', 'A.':'N',
	case 'A':
	  {
	    switch(next)
	      {
	      case '0': next = 'A'; break;
	      case '1': next = 'C'; break;
	      case '2': next = 'G'; break;
	      case '3': next = 'T'; break;
	      default: next = 'N'; break;
	      }
	  }
	  break;
	case 'C':
	  {
	    // 'C0':'C', 'C1':'A', 'C2':'T', 'C3':'G', 'C4':'N', 'C.':'N',
	     switch(next)
	      {
	      case '0': next = 'C'; break;
	      case '1': next = 'A'; break;
	      case '2': next = 'T'; break;
	      case '3': next = 'G'; break;
	      default: next = 'N'; break;
	      }
	  }
	  break;
	case 'G':
	  {
	    // 'G0':'G', 'G1':'T', 'G2':'A', 'G3':'C', 'G4':'N', 'G.':'N',
	     switch(next)
	      {
	      case '0': next = 'G'; break;
	      case '1': next = 'T'; break;
	      case '2': next = 'A'; break;
	      case '3': next = 'C'; break;
	      default: next = 'N'; break;
	      }
	  }
	  break;
	case 'T':
	  {
	    // 'T0':'T', 'T1':'G', 'T2':'C', 'T3':'A', 'T4':'N', 'T.':'N',
	     switch(next)
	      {
	      case '0': next = 'T'; break;
	      case '1': next = 'G'; break;
	      case '2': next = 'C'; break;
	      case '3': next = 'A'; break;
	      default: next = 'N'; break;
	      }
	  }
	  break;
	default: next = 'N'; break;
	}

      bp.push_back(next);
      base = next;
    }

  return bp;
}

#define check_color(b1, b2, c1, c2) ((b1 == c1 && b2 == c2) || (b1 == c2 && b2 == c1))

string convert_bp_to_color(const string& bp, bool remove_primer)
{
  if (bp.length() <= 1)
    return "";

  char base = toupper(bp[0]);
  string color;
  if (!remove_primer)
    color.push_back(base);
  
  for (string::size_type i = 1; i < bp.length(); ++i)
    {
      char next = toupper(bp[i]);
      
      if ((base == 'A' || base == 'G' || base == 'C' || base == 'T') && base == next)
	color.push_back('0');
      else if (check_color(base, next, 'A', 'C') || check_color(base, next, 'G', 'T'))
	color.push_back('1');
      else if (check_color(base, next, 'A', 'G') || check_color(base, next, 'C', 'T'))
	color.push_back('2');
      else if (check_color(base, next, 'A', 'T') || check_color(base, next, 'C', 'G'))
	color.push_back('3');
      else
	color.push_back('4');

      base = next;
    }

  return color;
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
	FLineReader fr(reads_file);
	//while(!feof(reads_file))
    while(!fr.isEof())
	{
		read.clear();
		
		// Get the next read from the file
		if (!next_fasta_record(fr, read.name, read.seq, reads_format))
		  break;

		if (reads_format == FASTQ)
		{
		  if (!next_fastq_record(fr, read.seq, read.alt_name, read.qual, reads_format))
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


