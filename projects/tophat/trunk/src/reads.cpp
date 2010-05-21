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

char* FLineReader::nextLine() {
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




bool next_fasta_record(FLineReader& fr,
					   string& defline, 
					   string& seq)

{
  seq.clear();
  defline.clear();
  char* buf=NULL;
  while ((buf=fr.nextLine())!=NULL) {
    if (buf[0]==0) continue; //skip empty lines
    if (buf[0] == '>') { //next record
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
  /*
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
	*/
    replace(seq.begin(),seq.end(),'.','N'); //shouldn't really be needed for FASTA files
	return !(seq.empty());
}

bool next_fastq_record(FLineReader& fr,
					   string& defline, 
					   string& seq, 
					   string& alt_name,
					   string& qual)
{
  seq.clear();
  defline.clear();
  alt_name.clear();
  qual.clear();
  char* fline=fr.nextLine();
  if (fline==NULL) return false;
  while (fline[0]==0) { //skip empty lines
    fline=fr.nextLine();
    if (fline==NULL) return false;
    }
  if (fline[0] != '@') { //first non-empty line must be a defline
    err_exit("Error: fastq record must start with '@'!\n",fline);
    return false;
    }
  defline=fline+1;
  string::size_type space_pos = defline.find_first_of(" \t\r");
  if (space_pos != string::npos) defline.resize(space_pos);
  //parse sequence:
  while ((fline=fr.nextLine())!=NULL) {
    if (fline[0]=='+') {
       break;
       }
    seq.append(fline);
    }
  //must be on '+' line here
  if (fline==NULL || fline[0]!='+') {
     err_exit("Error: '+' not found for fastq record %s\n",defline.c_str());
     return false;
     }
  alt_name=fline+1;
  space_pos = alt_name.find_first_of(" \t");
  if (space_pos != string::npos) alt_name.resize(space_pos);
   //read qv line(s) now:
  while ((fline=fr.nextLine())!=NULL) {
     qual.append(fline);
     if (qual.length()>=seq.length()) break;
     }
  // final check
  if (seq.length()!=qual.length()) {
           err_exit("Error: qual length (%d) differs from seq length (%d) for fastq record %s!\n",
               qual.length(), seq.length(), defline.c_str());
           return false;
           }
  //

  replace(seq.begin(),seq.end(),'.','N'); //replace dots with Ns in the sequence
  return !(seq.empty());

  /*
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
	*/
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
	FLineReader fr(reads_file);
	//while(!feof(reads_file))
    while(!fr.isEof())
	{
		read.clear();
		
		// Get the next read from the file
		if (reads_format == FASTA)
		{
			if (!next_fasta_record(fr, read.name, read.seq))
				break;
		}
		else if (reads_format == FASTQ)
		{
			if (!next_fastq_record(fr, read.name, read.seq, read.alt_name, read.qual))
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


