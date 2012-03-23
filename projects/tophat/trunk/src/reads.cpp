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
  if (fline==NULL || (reads_format == FASTQ && fline[0] != '+') ||
      (reads_format == FASTA && quals && fline[0] != '>')) {
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
      if (qual.length()>=seq.length()-1) break;
     }
  // final check
  if (color) {
     if (seq.length()==qual.length()) {
        //discard first qv
        qual=qual.substr(1);
        }
     if (seq.length()!=qual.length()+1) {
        err_exit("Error: length of quality string does not match seq length (%d) for color read %s!\n",
           seq.length(), alt_name.c_str());
        }
     }
  else {
    if (seq.length()!=qual.length()) {
           err_exit("Error: qual string length (%d) differs from seq length (%d) for read %s!\n",
               qual.length(), seq.length(), alt_name.c_str());
           //return false;
           }
    }
  //
  return !(qual.empty());
}

bool next_fastx_read(FLineReader& fr, Read& read, ReadFormat reads_format, 
                        FLineReader* frq) {
  /*
  if (fr.pushed_read)
    {
      read = fr.last_read;
      fr.pushed_read = false;
      return true;
    }
  */
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
        read.qual.append(buf);
      }
    if (read.qual.length()>=read.seq.length()-1)
          break;
    } //while qv lines
  
  // final check
  if (color) {
     if (read.seq.length()==read.qual.length()) {
        //discard first qv
        read.qual=read.qual.substr(1);
        }
     if (read.seq.length()!=read.qual.length()+1) {
        err_exit("Error: length of quality string does not match sequence length (%d) for color read %s!\n",
            read.seq.length(), read.alt_name.c_str());
        }
     }
  else {
   if (read.seq.length()!=read.qual.length()) {
           err_exit("Error: qual length (%d) differs from seq length (%d) for fastq record %s!\n",
               read.qual.length(), read.seq.length(), read.alt_name.c_str());
           return false;
           }
    }

  //fr.last_read = read;
  return !(read.seq.empty());
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
  
  static const size_t max_length = MAX_READ_LEN;
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


void bam2Read(bam1_t *b, Read& rd, bool alt_name=false) {
  GBamRecord bamrec(b);
  rd.clear();
  rd.seq=bamrec.seqData(&rd.qual);
  rd.name=bam1_qname(b);
  if (alt_name)
    rd.alt_name=bamrec.tag_str("ZN");
}


bool ReadStream::next_read(QReadData& rdata, ReadFormat read_format) {
  while (read_pq.size()<ReadBufSize && !r_eof) {
    //keep the queue topped off
    Read rf;
    if (get_direct(rf, read_format)) {
      uint64_t id = (uint64_t)atol(rf.name.c_str());
      QReadData rdata(id, rf, last_b());
      read_pq.push(rdata);
      }
    }
  if (read_pq.size()==0)
     return false;
  const QReadData& t = read_pq.top();
  rdata=t; //copy strings and duplicate b pointer!
  read_pq.pop();
  return true;
}

bool ReadStream::get_direct(Read& r, ReadFormat read_format) {
  if (fstream.file==NULL) return false;
  if (fstream.is_bam) {
	 bool got_read=false;
	 while (!got_read) {
       if (samread(fstream.bam_file, b) < 0) {
          r_eof=true;
          return false;
       }
       else {
    	if (bam_ignoreQC || (b->core.flag & BAM_FQCFAIL)==0)
            got_read=true;
       }
	 }
     bam2Read(b, r, bam_alt_name);
     return true;
   }
   if (!next_fastx_read(*flseqs, r, read_format, flquals)) {
        r_eof=true;
        return false;
        }
  return true;
}

// reads must ALWAYS be requested in increasing order of their ID
bool ReadStream::getRead(uint64_t r_id,
			 Read& read,
			 ReadFormat read_format,
			 bool strip_slash,
			 uint64_t begin_id,
			 uint64_t end_id,
			 GBamWriter* um_out, //unmapped reads output
			 bool um_write_found //write the found ones
			 ) {
  if (!fstream.file)
       err_die("Error: calling ReadStream::getRead() with no file handle!");
  if (r_id<last_id)
      err_die("Error: ReadStream::getRead() called with out-of-order id#!");
  last_id=r_id;
  bool found=false;
  read.clear();
  while (!found) {
    QReadData rdata;
    if (!next_read(rdata, read_format))
        break;
    /*
    if (strip_slash) {
       string::size_type slash = rdata.read.name.rfind("/");
       if (slash != string::npos)
          rdata.read.name.resize(slash);
       }
    uint64_t id = (uint64_t)atol(read.name.c_str());
    */
    if (rdata.id >= end_id)
      return false;

    if (rdata.id < begin_id)
      continue;

    if (rdata.id == r_id)
      {
      read=rdata.read;
	  found=true;
      }
    else if (rdata.id > r_id)
      { //can't find it, went too far
      //only happens when reads [mates] were removed for some reason
      //read_pq.push(make_pair(id, read));
      read_pq.push(rdata);
      break;
      }
    if (um_out && ((um_write_found && found) ||
    	           (!um_write_found && !found))) {
       //write unmapped reads
       //fprintf(um_out, "@%s\n%s\n+\n%s\n", read.alt_name.c_str(),
       //                        read.seq.c_str(), read.qual.c_str());
    	string rname(rdata.read.alt_name);

    	GBamRecord bamrec(rname.c_str(), -1, 0, false, rdata.read.seq.c_str(),
    		NULL, rdata.read.qual.c_str());
    	if (rdata.matenum) {
    	  bamrec.set_flag(BAM_FPAIRED);
    	  if (rdata.matenum==1) bamrec.set_flag(BAM_FREAD1);
    	  else bamrec.set_flag(BAM_FREAD2);
    	}
    	if (rdata.trashCode) {
    	  if (rdata.trashCode!='M') {
    		   //multi-mapped reads did not really QC-fail
  		       bamrec.set_flag(BAM_FQCFAIL);
    	  }
    	  bamrec.add_aux("ZT", 'A', 1, (uint8_t*)&rdata.trashCode);
    	}
        um_out->write(&bamrec);
      }
    } //while reads
  return found;
}
