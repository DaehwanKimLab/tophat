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

#include <iostream>
#include <cassert>
#include <string>
#include <algorithm>
#include <cstring>
#include <cstdlib>

#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

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
              ungetc(c,file);
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
  if (fr.fhandle() == NULL) return;
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

  replace(seq.begin(), seq.end(), '.', color ? '4' : 'N'); //shouldn't really be needed for FASTA files
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
    
     //if ((!color && qual.length()>=seq.length()) || (color && qual.length()+1>=seq.length())) break;
     if (qual.length()+1>=seq.length()) break;
     }
  // final check
  if (color && qual.length()==seq.length()-1) {
	    string tmpq("!");
	    tmpq+=qual;
	    qual=tmpq;
        }
  if (seq.length()!=qual.length()) {
           err_exit("Error: qual string length (%d) differs from seq length (%d) for read %s!\n",
               qual.length(), seq.length(), alt_name.c_str());
           return false;
           }
  //
  return !(qual.empty());
}

bool next_fastx_read(FLineReader& fr, Read& read, ReadFormat reads_format, 
                        FLineReader* frq) {
  read.clear();
  char* buf=NULL;
  while ((buf=fr.nextLine())!=NULL) {
    if (buf[0]==0) continue; //skip empty lines
    if ((reads_format == FASTA && buf[0] == '>') ||
          (reads_format == FASTQ && (buf[0] == '+' || buf[0] == '@'))) { //next record
        if (read.seq.length()>0) { //current record ending
           fr.pushBack();
           break;
           }
        read.name=buf+1;
        string::size_type space_pos = read.name.find_first_of(" \t");
        if (space_pos != string::npos) {
            read.name.resize(space_pos);
            }
        continue;
        } //defline
    // sequence line
    read.seq.append(buf);
    } //line reading loop

  replace(read.seq.begin(), read.seq.end(), '.', color ? '4' : 'N'); //shouldn't really be needed for FASTA files
  if (reads_format != FASTQ && frq==NULL)
      return (!read.seq.empty());
  if (frq==NULL) frq=&fr; //FASTQ
  //FASTQ or quals in a separate file -- now read quality values
  buf=frq->nextLine();
  if (buf==NULL) return false;
  while (buf[0]==0) { //skip empty lines
    buf=frq->nextLine();
    if (buf==NULL) return false;
    }
  //must be on '+' line here
  if (buf==NULL || (reads_format == FASTQ && buf[0] != '+') ||
           (reads_format == FASTA && buf[0] != '>')) {
     err_exit("Error: beginning of quality values record not found! (%s)\n",buf);
     return false;
     }
  read.alt_name=buf+1;
  string::size_type space_pos = read.alt_name.find_first_of(" \t");
  if (space_pos != string::npos) read.alt_name.resize(space_pos);
   //read qv line(s) now:
  while ((buf=frq->nextLine())!=NULL) {
    if (integer_quals)
      {
        vector<string> integer_qual_values;
        tokenize(string(buf), " ", integer_qual_values);
        string temp_qual;
        for (size_t i = 0; i < integer_qual_values.size(); ++i)
          {
            int qual_value = atoi(integer_qual_values[i].c_str());
            if (qual_value < 0) qual_value = 0;
            temp_qual.push_back((char)(qual_value + 33));
          }
        read.qual.append(temp_qual);
      }
    else {
      // if (color && read.qual.length()==0 && buf[0]=='!')
      //   read.qual.append(&(buf[1])); //some color qual strings start with '!' for the adaptor
      //  else
         read.qual.append(buf);
      }
    if ((!color && read.qual.length()>=read.seq.length())
          || (color && read.qual.length()+1>=read.seq.length())) break;
    } //while qv lines
  
  // final check
  // final check
  if (color && read.qual.length()==read.seq.length()-1) {
	    string tmpq("!");
	    tmpq+=read.qual;
	    read.qual=tmpq;
        }
  if (read.seq.length()!=read.qual.length()) {
  //if ((!color && read.seq.length()!=read.qual.length()) || (color && read.seq.length()!=read.qual.length()+1)) {
           err_exit("Error: qual length (%d) differs from seq length (%d) for fastq record %s!\n",
               read.qual.length(), read.seq.length(), read.alt_name.c_str());
           return false;
           }
  //
  return !(read.qual.empty());
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

// daehwan - reduce code redundancy!
seqan::String<char> convert_color_to_bp(char base, const seqan::String<char>& color)
{
  if (seqan::length(color) <= 0)
    return "";

  string bp;
  for (string::size_type i = 0; i < seqan::length(color); ++i)
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

#define two_bps_to_color(b1, b2, c)				      \
  if (((b1) == 'A' || (b1) == 'G' || (b1) == 'C' || (b1) == 'T') && (b1) == (b2)) \
  c = '0'; \
  else if (check_color((b1), (b2), 'A', 'C') || check_color((b1), (b2), 'G', 'T')) \
  c = '1'; \
  else if (check_color((b1), (b2), 'A', 'G') || check_color((b1), (b2), 'C', 'T')) \
  c = '2'; \
  else if (check_color((b1), (b2), 'A', 'T') || check_color((b1), (b2), 'C', 'G')) \
  c = '3'; \
  else \
  c = '4';
 

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

      char c = '0';
      two_bps_to_color(base, next, c);
      color.push_back(c);

      base = next;
    }

  return color;
}

// daehwan - check this - 
seqan::String<char> convert_bp_to_color(const seqan::String<char>& bp, bool remove_primer)
{
  if (seqan::length(bp) <= 1)
    return "";

  char base = toupper(bp[0]);
  string color;
  if (!remove_primer)
    color.push_back(base);
  
  for (string::size_type i = 1; i < seqan::length(bp); ++i)
    {
      char next = toupper(bp[i]);

      char c = '0';
      two_bps_to_color(base, next, c);
      color.push_back(c);

      base = next;
    }

  return color;
}

/*
 */
void BWA_decode(const string& color, const string& qual, const string& ref, string& decode)
{
  assert(color.length() == ref.length() - 1);
  
  const size_t max_length = 256;
  const unsigned int max_value = max_length * 0xff;
  size_t length = color.length();
  if (length < 1 || length + 1 > max_length)
    {
      return;
    }

  unsigned int f[max_length * 4];
  char ptr[max_length * 4];

  unsigned int q_prev = 0;
  for (unsigned int i = 0; i < length + 1; ++i)
    {
      unsigned int q = (unsigned int) (qual.length() <= i ? 'I' : qual[i]) - 33;
      for (unsigned int j = 0; j < 4; ++j)
	{
	  size_t i_j = i * 4 + j;
	  if (i == 0)
	    {
	      f[i_j] = "ACGT"[j] == ref[i] ? 0 : q;
	      ptr[i_j] = 4;
	      continue;
	    }

	  f[i_j] = max_value;
	  char base = "ACGT"[j];
	  for (unsigned int k = 0; k < 4; ++k)
	    {
	      char base_prev = "ACGT"[k];
	      char ref_color;
	      two_bps_to_color(base_prev, base, ref_color);

	      char base_prev_prev = "ACGTN"[ptr[(i-1)*4 + k]];
	      char ref_color_prev;
	      two_bps_to_color(base_prev_prev, base_prev, ref_color_prev);

	      char color_curr = color[i-1];
	      char color_prev = i >= 2 ? color[i-2] : '4';	      

	      int q_hat = 0;
	      if (color_prev == ref_color_prev && color_prev != '4')
		{
		  if (color_curr == ref_color)
		    q_hat = q + q_prev;
		  else
		    q_hat = q_prev - q;
		}
	      else if (color_curr == ref_color)
		{
		  q_hat = q - q_prev;
		}

	      unsigned int f_k = f[(i-1) * 4 + k] +
		(base == ref[i] ? 0 : q_hat) +
		(color_curr == ref_color ? 0 : q);

	      if (f_k < f[i_j])
		{
		  f[i_j] = f_k;
		  ptr[i_j] = k;
		}
	    }
	}

      q_prev = q;
    }

  unsigned int min_index = 0;
  unsigned int min_f = f[length * 4];
  for (unsigned int i = 1; i < 4; ++i)
    {
      unsigned int temp_f = f[length * 4 + i];
      if (temp_f < min_f)
	{
	  min_f = temp_f;
	  min_index = i;
	}
    }

  decode.resize(length + 1);
  decode[length] = "ACGT"[min_index];
  for (unsigned int i = length; i > 0; --i)
    {
      min_index = ptr[i * 4 + min_index];
      decode[i-1] = "ACGT"[min_index];
    }
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
  while(!fr.isEof())
	{
	read.clear();
	  
	  // Get the next read from the file
	if (!next_fastx_read(fr, read, reads_format))
	    break;

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
