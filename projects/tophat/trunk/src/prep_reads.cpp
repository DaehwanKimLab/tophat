/*
 *  prep_reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/2/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *	Derived from maq "catfilter", by Heng Li at Sanger
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <cassert>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "common.h"
#include "reads.h"
#include "tokenize.h"
#include "qual.h"

using namespace std;

void format_qual_string(string& qual_str)
{
  for (size_t i = 0; i < qual_str.size(); ++i)
    {
      qual_str[i] = charToPhred33(qual_str[i], 
				      solexa_quals, 
				      phred64_quals);
    }
}

vector<string> flt_reads_fnames;
bool readmap_loaded=false;
vector<bool> readmap; //for filter_multihits
vector<bool> mate_readmap; //for filter_multihits

void load_readmap(string& flt_fname, vector<bool>& rmap) {
  //readmap_loaded=false;
  if (flt_fname.empty()) return;
  ReadStream rdstream(flt_fname, NULL, true);
  Read read;
  while (rdstream.get_direct(read, reads_format)) {
    uint32_t rnum=(uint32_t)atol(read.name.c_str());
    if (rnum>=(uint32_t) rmap.size())
         rmap.resize(rnum+1, false);
    rmap[rnum] = true;
    }
  readmap_loaded=true;
}

bool check_readmap(vector<bool>& rmap, uint32_t read_id) {
 if (read_id>=rmap.size()) return false;
    else return rmap[read_id];
}

void flt_reads_and_hits(vector<string>& reads_files) {
  if (!readmap_loaded)
      err_die("Error: filtering reads not enabled, aborting.");
  if (aux_outfile.empty())
       err_die("Error: auxiliary output file not provided.");
  // -- filter mappings:
  string fext=getFext(flt_mappings);
  if (fext=="bam") {
	samfile_t* fbam=samopen(flt_mappings.c_str(), "rb", 0);
    if (fbam==NULL)
        err_die("Error opening BAM file %s!\n", flt_mappings.c_str());
    bam1_t *b = bam_init1();
    string aux_index(aux_outfile);
    aux_index+=".index";
    GBamWriter wbam(aux_outfile.c_str(), fbam->header, aux_index);
    while (samread(fbam, b) > 0) {
      char* rname  = bam1_qname(b);
      uint32_t rid=(uint32_t)atol(rname);
      if (check_readmap(readmap, rid)) {
    	 //write this mapping into the output file
    	 wbam.write(b, rid);
         }
      }
    bam_destroy1(b);
    samclose(fbam);
    }
  else {
	bool is_sam=false;
	string s(flt_mappings);
	for(size_t i=0; i!=s.length(); ++i)
			s[i] = std::tolower(s[i]);
	if (fext=="sam" || s.rfind(".sam.")+fext.length()+5==s.length()) is_sam=true;
	string unzcmd=getUnpackCmd(flt_mappings);
	FZPipe hitsfile(flt_mappings, unzcmd);
	FLineReader fh(hitsfile);
	FZPipe outfile;
	outfile.openWrite(aux_outfile.c_str(), zpacker);
	if (outfile.file==NULL) err_die("Error: cannot create file %s", aux_outfile.c_str());
	const char* line;
	while ((line=fh.nextLine())) {
	  if (is_sam && line[0]=='@') {
		  //copy the header
		  fprintf(outfile.file, "%s\n", line);
		  continue;
		  }
	  char* tab=strchr((char*)line, '\t');
	  if (tab==NULL)
		  err_die("Error: cannot find tab character in %s mappings line:\n%s",
			  flt_mappings.c_str(),line);
	  *tab=0;
	  uint32_t rid = (uint32_t) atol(line);
	  if (rid==0)
		err_die("Error: invalid read ID (%s) parsed from mapping file %s",
			  line, flt_mappings.c_str());
	  *tab='\t';
	  if (check_readmap(readmap, rid)) {
		 fprintf(outfile.file, "%s\n", line);
		 }
	  }
	outfile.close();
	hitsfile.close();
	}
  // -- now filter reads
  //FILE* findex = NULL;
  GBamWriter* wbam=NULL;
  FILE* fout=NULL;
  if (std_outfile.empty()) {
    fout=stdout;
  }
  else { //output file name explicitely given
	if (getFext(std_outfile)=="bam") {
	  if (sam_header.empty()) err_die("Error: sam header file not provided.\n");
	  wbam = new GBamWriter(std_outfile.c_str(), sam_header.c_str(), index_outfile);
	  //wbam = new GBamWriter(std_outfile, index_outfile);
	}
	else {
	  fout = fopen(std_outfile.c_str(), "w");
	  if (fout==NULL)
	       err_die("Error: cannot create file %s\n", std_outfile.c_str());
	}
  }
  /*
  if (wbam==NULL && !index_outfile.empty()) {
    findex = fopen(index_outfile.c_str(), "w");
    if (findex == NULL)
       err_die("Error: cannot create file %s\n", index_outfile.c_str());
    }
    */
  for (size_t fi = 0; fi < reads_files.size(); ++fi) {
      //only one file expected here, this is not the initial prep_reads
      Read read;
      ReadStream readstream(reads_files[fi], NULL, true);
      //skip_lines(fr);
      while (readstream.get_direct(read)) {
        uint32_t rnum=(uint32_t)atol(read.name.c_str());
        if (check_readmap(readmap, rnum)) {
          if (wbam) {
        	GBamRecord bamrec(read.name.c_str(), -1, 0, false,
        		read.seq.c_str(), NULL, read.qual.c_str());
        	wbam->write(bamrec.get_b(), rnum);
          }
          else {
        	fprintf(fout, "@%s\n%s\n+%s\n%s\n",
                read.name.c_str(),
                read.seq.c_str(),
                read.alt_name.c_str(),
                read.qual.c_str());
          }
        }
      }
  }
 if (wbam) delete wbam;
 if (fout && fout!=stdout) fclose(fout);
}


void writePrepBam(GBamWriter* wbam, Read& read, uint32_t rid, char trashcode=0, int matenum=0) {
  if (wbam==NULL) return;
  string rnum;
  str_appendUInt(rnum, rid);
  string rname(read.name);

  // attach a primer tag and the following color to the end of the read name for colorspace reads
  // also, attach a quality value
  if (color)
    {
      rnum.push_back(read.seq[0]);
      rnum.push_back(read.seq[1]);

      if (!read.qual.empty())
	rnum.push_back(read.qual[1]);
      else
	rnum.push_back('!');
    }
    
  GBamRecord bamrec(rnum.c_str(), -1, 0, false, read.seq.c_str(), NULL, read.qual.c_str());
  if (matenum) {
	bamrec.set_flag(BAM_FPAIRED);
	if (matenum==1) bamrec.set_flag(BAM_FREAD1);
	else if (matenum==2) bamrec.set_flag(BAM_FREAD2);
  }

  bamrec.add_aux("ZN", 'Z', rname.length(), (uint8_t*)rname.c_str());

  if (trashcode) {
	 bamrec.set_flag(BAM_FQCFAIL);
     bamrec.add_aux("ZT", 'A', 1, (uint8_t*)&trashcode);
     }
  wbam->write(bamrec.get_b(), rid);
}

bool processRead(int matenum, Read& read, uint32_t next_id,  int& num_reads_chucked,
	int& multimap_chucked, GBamWriter* wbam, FILE* fout, FILE* fqindex,
	int& min_read_len, int& max_read_len, uint64_t& fout_offset, vector<bool>& rmap) {
  if (read.seq.length()<12) {
            ++num_reads_chucked;
	  writePrepBam(wbam, read, next_id, 'S', matenum);
	  return false;
	}
        if ((int)read.seq.length()<min_read_len)
             min_read_len=read.seq.length();
        if ((int)read.seq.length()>max_read_len)
             max_read_len=read.seq.length();

        if (color && read.seq[1] == '4') {
          ++num_reads_chucked;
	writePrepBam(wbam, read, next_id, 'c', matenum);
	return false;
          }

  if (readmap_loaded && check_readmap(rmap, next_id)) {
          ++num_reads_chucked;
          ++multimap_chucked;
	writePrepBam(wbam, read, next_id, 'M', matenum);
	return false;
          }
        format_qual_string(read.qual);
        std::transform(read.seq.begin(), read.seq.end(), read.seq.begin(), ::toupper);
        char counts[256];
        memset(counts, 0, sizeof(counts));
        // Count up the bad characters
        for (unsigned int i = 0; i != read.seq.length(); ++i)
          {
            char c = (char)toupper(read.seq[i]);
            counts[(size_t)c]++;
          }

        double percent_A = (double)(counts[(size_t)'A']) / read.seq.length();
        double percent_C = (double)(counts[(size_t)'C']) / read.seq.length();
        double percent_G = (double)(counts[(size_t)'G']) / read.seq.length();
        double percent_T = (double)(counts[(size_t)'T']) / read.seq.length();
        double percent_N = (double)(counts[(size_t)'N']) / read.seq.length();
        double percent_4 = (double)(counts[(size_t)'4']) / read.seq.length();

        // Chuck the read if there are at least 5 'N's or if it's mostly
        // (>90%) 'N's and 'A's
        char trash_code=0;
        if (percent_A > 0.9 ||
            percent_C > 0.9 ||
            percent_G > 0.9 ||
            percent_T > 0.9)
              trash_code='L';
        else if (percent_N >= 0.1 || percent_4 >=0.1)
              trash_code='N';
        if (trash_code) {
          ++num_reads_chucked;
	writePrepBam(wbam, read, next_id, trash_code, matenum);
	return false;
        }

        if (wbam) {
	  if (reads_format == FASTA && !quals)
	    read.qual = string(read.seq.length(), 'I').c_str();
	  else if (color && quals)
	    read.qual = "!" + read.qual;
	  
          writePrepBam(wbam, read, next_id, 0, matenum);
        }
        else {
		  // daehwan - we should not use buf in printf function
		  // because it may contain some control characters such as "\" from quality values.
		  // Here, buf is only used for calculating the file offset
		  char buf[2048] = {0};
		  if (reads_format == FASTQ or (reads_format == FASTA && quals))
				{
		  sprintf(buf, "@%u\n%s\n+%s\n%s\n",
			  next_id,
			  read.seq.c_str(),
			  read.name.c_str(),
			  read.qual.c_str());

		  fprintf(fout, "@%u\n%s\n+%s\n%s\n",
				 next_id,
				 read.seq.c_str(),
				 read.name.c_str(),
				 read.qual.c_str());
				}
			  else if (reads_format == FASTA)
				{
				  string qual;
				  if (color)
					qual = string(read.seq.length()-1, 'I').c_str();
				  else
					qual = string(read.seq.length(), 'I').c_str();

				  sprintf(buf, "@%u\n%s\n+%s\n%s\n",
			  next_id,
			  read.seq.c_str(),
			  read.name.c_str(),
			  qual.c_str());

		  fprintf(fout, "@%u\n%s\n+%s\n%s\n",
			  next_id,
			  read.seq.c_str(),
			  read.name.c_str(),
			  qual.c_str());
				}
		  else
			{
		  assert(0);
			}

	if (fqindex != NULL)
			{
		  if ((next_id - num_reads_chucked) % INDEX_REC_COUNT == 0)
		fprintf(fqindex, "%d\t%lu\n", next_id, (long unsigned)fout_offset);
	}

	fout_offset += strlen(buf);
  }
  return true;
} //validate read


const char* ERR_FILE_CREATE="Error: cannot create file %s\n";

void process_reads(vector<string>& reads_fnames, vector<FZPipe>& quals_files,
	    vector<string>& mate_fnames, vector<FZPipe>& mate_quals_files)
{
   //TODO: add the option to write the garbage reads into separate file(s)
  int num_reads_chucked = 0;
  int multimap_chucked = 0;
  int mate_num_reads_chucked = 0;
  int mate_multimap_chucked = 0;
  int min_read_len = 20000000;
  int max_read_len = 0;
  int mate_min_read_len = 20000000;
  int mate_max_read_len = 0;
  uint32_t next_id = 0;
  uint32_t num_left = 0;
  uint32_t num_mates = 0;

  FILE* fw=NULL; //aux output file
  string outfname; //std_outfile after instancing template
  string mate_outfname;
  string idxfname; //index_outfile after instancing template
  string mate_idxfname;
  bool have_mates = (mate_fnames.size() > 0);

  if (!aux_outfile.empty()) {
    fw=fopen(aux_outfile.c_str(), "w");
    if (fw==NULL)
       err_die(ERR_FILE_CREATE,aux_outfile.c_str());
    }

  FILE* fqindex = NULL; //fastq index
  FILE* mate_fqindex = NULL;
  GBamWriter* wbam=NULL;
  GBamWriter* mate_wbam=NULL;
  FILE* fout=NULL;
  FILE* mate_fout=NULL;
  uint64_t fout_offset = 0;
  uint64_t mate_fout_offset = 0;
  if (std_outfile.empty()) {
    fout=stdout;
    //for PE reads, flt_side will decide which side is printed (can't be both)
    if (have_mates && flt_side==2)
      err_die("Error: --flt-side option required for PE reads directed to stdout!\n");
    mate_fout=stdout;
  }
  else { //output file name explicitely given
	//could be a template
	if (std_outfile.find("%side%") != string::npos) {
	  outfname=str_replace(std_outfile, "%side%", "left");
	  if (have_mates)
	    mate_outfname=str_replace(std_outfile, "%side%", "right");
	}
	else {
	  outfname=std_outfile;
	}
	if (index_outfile.find("%side%") != string::npos) {
	  idxfname=str_replace(index_outfile, "%side%", "left");
	  if (have_mates)
	    mate_idxfname=str_replace(index_outfile, "%side%", "right");
	}
	else {
	  idxfname=index_outfile;
	}
	if (getFext(outfname)=="bam") {
	  if (sam_header.empty()) err_die("Error: sam header file not provided.\n");
	  wbam = new GBamWriter(outfname.c_str(), sam_header.c_str(), idxfname);
	  if (!mate_outfname.empty()) {
		mate_wbam = new GBamWriter(mate_outfname.c_str(), sam_header.c_str(), mate_idxfname);
	  }
	}
	else { //fastq output
	  fout = fopen(outfname.c_str(), "w");
	  if (fout==NULL)
	       err_die(ERR_FILE_CREATE, outfname.c_str());
	  mate_fout = fopen(mate_outfname.c_str(), "w");
	  if (mate_fout==NULL)
	       err_die(ERR_FILE_CREATE, mate_outfname.c_str());
	}
  }
  if (wbam==NULL && !idxfname.empty()) {
	//fastq file output, indexed
	fqindex = fopen(idxfname.c_str(), "w");
	if (fqindex == NULL)
	  err_die(ERR_FILE_CREATE, idxfname.c_str());
	if (!mate_idxfname.empty()) {
	  mate_fqindex = fopen(mate_idxfname.c_str(), "w");
	  if (mate_fqindex == NULL)
		err_die(ERR_FILE_CREATE, mate_idxfname.c_str());
			}

  }
  bool possible_mate_mismatch=false;
  size_t max_files=max(reads_fnames.size(), mate_fnames.size());
  for (size_t fi = 0; fi < max_files; ++fi)
  {
	Read read;
	Read mate_read;
    ReadStream* reads=NULL;
	ReadStream* mate_reads=NULL;
	FZPipe* fq=NULL;
	FZPipe* mate_fq=NULL;
	bool have_l_reads=(fi<reads_fnames.size());
	bool have_r_reads=(have_mates && fi<mate_fnames.size());
	if (have_l_reads) {
	  if (quals)
		fq = & quals_files[fi];
	  reads=new ReadStream(reads_fnames[fi], fq, true);
	}
	if (have_r_reads) {
	  if (quals)
		 mate_fq = & mate_quals_files[fi];
	  mate_reads=new ReadStream(mate_fnames[fi], mate_fq, true);
	}

	while (have_l_reads || have_r_reads) {
	  if (have_l_reads && (have_l_reads=reads->get_direct(read, reads_format)))
		 num_left++;
	  // Get the next read from the file
	  int matenum=0; // 0 = unpaired, 1 = left, 2 = right
	  if (have_r_reads && (have_r_reads=mate_reads->get_direct(mate_read, reads_format)) ) {
		num_mates++;
	  }
	  if (have_l_reads && have_r_reads) {
		matenum = 1; //read is first in a pair
		if (have_l_reads && have_r_reads && !possible_mate_mismatch) {
  		  //check if reads are paired correctly
		  int nl=read.name.length();
		  bool mate_match=(nl==(int)mate_read.name.length());
		  int m_len=0, c=0;
		  while (c<nl) {
			if (read.name[c]==mate_read.name[c]) {
			  m_len++; //number of matching chars in the same position
			}
			c++;
		  }
		  if (mate_match && nl-m_len > 2) //more than 2 chars differ
			mate_match=false;
		  if (!mate_match) {
			fprintf(stderr, "WARNING: read pairing issues detected (check prep_reads.log) !\n"
				" Pair #%d name mismatch: %s vs %s\n",
				next_id+1, read.name.c_str(), mate_read.name.c_str());
			possible_mate_mismatch=true;
		  }
		} //mate check
	  } //paired reads
	  if (have_l_reads || have_r_reads) {
		  //IMPORTANT: to keep paired reads in sync, this must be
		  //incremented BEFORE any reads are chucked !
		  ++next_id;
	  }
	  if ((flt_side & 1)==0 && have_l_reads)
		  //for unpaired reads or left read in a pair
	      processRead(matenum, read, next_id,  num_reads_chucked, multimap_chucked,
		     wbam, fout, fqindex, min_read_len,  max_read_len,  fout_offset, readmap);
	  if (flt_side>0 && have_r_reads) {
		  matenum = have_l_reads ? 2 : 0;
		  processRead(matenum, mate_read, next_id,  mate_num_reads_chucked, mate_multimap_chucked,
			  mate_wbam, mate_fout, mate_fqindex, mate_min_read_len,  mate_max_read_len,
			  mate_fout_offset, mate_readmap);
      }
    } //while !fr.isEof()
	if (reads)
	  delete reads;
	if (mate_reads)
	  delete mate_reads;
    } //for each input file
  if (fout!=stdout || (flt_side & 1) == 0) {
	fprintf(stderr, "%u out of %u reads have been filtered out\n",
		num_reads_chucked, num_left);
	if (readmap_loaded)
	  fprintf(stderr, "\t(%u filtered out due to %s)\n",
		  multimap_chucked, flt_reads_fnames[0].c_str());
  }

  if (have_mates && (fout!=stdout || flt_side>0)) {
	fprintf(stderr, "%u out of %u read mates have been filtered out\n",
		mate_num_reads_chucked, num_mates);
	if (readmap_loaded && mate_multimap_chucked)
	  fprintf(stderr, "\t(%u mates filtered out due to %s)\n",
		  mate_multimap_chucked, flt_reads_fnames[1].c_str());
  }

  if (wbam) { delete wbam; }
  if (mate_wbam) { delete mate_wbam; }
  if (fout && fout!=stdout) fclose(fout);
  if (mate_fout) fclose(mate_fout);
  if (fw!=NULL) {
	string side("");
	if (have_mates)
	   side="left_";
    fprintf(fw, "%smin_read_len=%d\n", side.c_str(), min_read_len - (color ? 1 : 0));
    fprintf(fw, "%smax_read_len=%d\n", side.c_str(), max_read_len - (color ? 1 : 0));
    fprintf(fw, "%sreads_in =%d\n", side.c_str(), num_left);
    fprintf(fw, "%sreads_out=%d\n", side.c_str(), num_left-num_reads_chucked);
    if (have_mates) {
      side="right_";
      fprintf(fw, "%smin_read_len=%d\n", side.c_str(), mate_min_read_len - (color ? 1 : 0));
      fprintf(fw, "%smax_read_len=%d\n", side.c_str(), mate_max_read_len - (color ? 1 : 0));
      fprintf(fw, "%sreads_in =%d\n", side.c_str(), num_mates);
      fprintf(fw, "%sreads_out=%d\n", side.c_str(), num_mates-mate_num_reads_chucked);
    }
    fclose(fw);
    }

  if (fqindex) fclose(fqindex);
  if (mate_fqindex) fclose(mate_fqindex);
}

void print_usage()
{
  fprintf(stderr, "Usage:\n prep_reads [--filter-multi <multi.fq>] <reads1.fa/fq,...,readsN.fa/fq> [<read_qual_files>,..] \\"
	  "[<mate_reads1.fa/fq,...,mate_readsN.fa/fq> [<mate_qual_files>,..]\n");
}

void open_qual_files(vector<FZPipe>& quals_files, string& quals_file_list) {
  vector<string> quals_file_names;
  tokenize(quals_file_list, ",", quals_file_names);
  for (size_t i = 0; i < quals_file_names.size(); ++i)
  {
	FZPipe seg_file(quals_file_names[i], true);
	if (seg_file.file == NULL)
	{
	  fprintf(stderr, "Error: cannot open qual. file %s\n",
		  quals_file_names[i].c_str());
	  exit(1);
	}
	quals_files.push_back(seg_file);
  }
}

int main(int argc, char *argv[])
{
  fprintf(stderr, "prep_reads v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
  fprintf(stderr, "---------------------------\n");

  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string reads_file_list(argv[optind++]);
  vector<string> reads_filenames;
  tokenize(reads_file_list, ",",reads_filenames);
  vector<FZPipe> quals_files;
  if (quals)
    {
    if (optind>=argc) {
      err_die("Error: quality value file(s) not provided !\n");
    }
      string quals_file_list = argv[optind++];
    open_qual_files(quals_files, quals_file_list);
    if (quals_files.size()!=reads_filenames.size())
        err_die("Error: number of quality value files must much the number of read files!\n");
  }
  string mate_file_list;
  vector<string> mate_filenames;
  vector<FZPipe> mate_quals_files;
  if (optind<argc)
	mate_file_list = argv[optind++];

  if (!mate_file_list.empty()) {
	tokenize(mate_file_list, ",",mate_filenames);
	/*
	if (reads_filenames.size()!=mate_filenames.size()) {
	    err_die("Error: same number of input files required for paired reads!\n");
	  }
	*/
	if (quals)
	{
	  if (optind>=argc) {
	      err_die("Error: mate quality value file(s) not provided !\n");
	  }
	  string mate_quals_file_list = argv[optind++];
	  open_qual_files(mate_quals_files, mate_quals_file_list);
	  if (mate_quals_files.size()!=mate_filenames.size())
	      err_die("Error: number of quality value files must much the number of read files!\n");
	      }
	}
  if (!flt_reads.empty()) {
	//for multi-mapped prefiltering usage
	readmap_loaded = false;
	tokenize(flt_reads, ",", flt_reads_fnames);
	load_readmap(flt_reads_fnames[0], readmap);
	if (flt_reads_fnames.size()==2)
	  load_readmap(flt_reads_fnames[1], mate_readmap);
    }
  if (flt_mappings.empty())
     process_reads(reads_filenames, quals_files, mate_filenames, mate_quals_files);
  else //special use case: filter previous mappings (when prefiltering)
     flt_reads_and_hits(reads_filenames);
    
  return 0;
}

