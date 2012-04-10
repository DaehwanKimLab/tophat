/*
 *  juncs_db.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "fusions.h"


using namespace std;
using namespace seqan;
using std::set;


// Length of the outer dimension of a single insert from the paired-end library

static int read_length = -1;
void print_usage()
{
    fprintf(stderr, "Usage:   juncs_db <min_anchor> <read_length> <splice_coords1,...,splice_coordsN> <insertion_coords1,...,insertion_coordsN> <deletion_coords1,...,deletion_coordsN> <ref.fa>\n");
}

typedef vector<string> Mapped;

bool splice_junc_lt(const pair<size_t, size_t>& lhs, 
					const pair<size_t, size_t>& rhs)
{
	if (lhs.first < rhs.first)
		return true;
	else
		return lhs.second < rhs.second;
}

/**
 * Given an insertion set, this code will print FASTA entries
 * for the surrounding sequence. The names of the FASTA entries
 * contain information about the exact location and nature of the
 * insertion. The entry is generally of the form
 * <contig>|<left end of fasta sequence>|<position of insertion>-<sequence of insertion>|<right end of fasta sequence>|ins|<[fwd|rev]>
 */ 
template<typename TStr>
void print_insertion(const Insertion& insertion,
		     int read_len,
		     TStr& ref_str,
		     const string& ref_name,
		     ostream& splice_db)
{
  int half_splice_len = read_len - min_anchor_len;
  
  size_t left_start, right_start;
  size_t left_end, right_end;

  size_t ref_len = length(ref_str);
  if (insertion.left >= 0 && insertion.left <= ref_len)
    {
      left_start = (int)insertion.left - half_splice_len + 1 >= 0 ? (int)insertion.left - half_splice_len + 1 : 0;
      left_end = left_start + half_splice_len;
      
      right_start = (int)left_end; 
      right_end = right_start + half_splice_len  < ref_len ? right_start + half_splice_len : ref_len;

      Infix<RefSequenceTable::Sequence>::Type left_splice = infix(ref_str, left_start, left_end);
      Infix<RefSequenceTable::Sequence>::Type right_splice = infix(ref_str, right_start, right_end); 
      
      splice_db << ">" << ref_name << "|" << left_start << "|" << insertion.left << "-" << insertion.sequence 
		<< "|" << right_end << "|ins|" << ("fwd")  << endl;
      
      splice_db << left_splice << insertion.sequence << right_splice << endl;
    }
}


template<typename TStr>
void print_splice(const Junction& junction,
		  int read_len,
		  const string& tag,
		  TStr& ref_str,
		  const string& ref_name,
		  ostream& splice_db)
{
  // daehwan - this is tentative, let's think about this more :)
  // int half_splice_len = read_len - min_anchor_len;
  int half_splice_len = read_len;
  size_t left_start, right_start;
  size_t left_end, right_end;

  size_t ref_len = length(ref_str);
  
  if (junction.left >= 0 && junction.left <= ref_len &&
      junction.right >= 0 && junction.right <= ref_len)
    {
      left_start = (int)junction.left - half_splice_len + 1 >= 0 ? (int)junction.left - half_splice_len + 1 : 0;
      left_end = left_start + half_splice_len;
      
      right_start = junction.right;
      right_end = right_start + half_splice_len < ref_len ? right_start + half_splice_len : ref_len;

      Infix<RefSequenceTable::Sequence>::Type left_splice = infix(ref_str,
								  left_start, 
								  left_end);
      Infix<RefSequenceTable::Sequence>::Type right_splice = infix(ref_str, 
								   right_start, 
								   right_end);
      
      splice_db << ">" << ref_name << "|" << left_start << "|" << junction.left <<
	"-" << junction.right << "|" << right_end << "|" << tag << endl;
      
      splice_db << left_splice << right_splice << endl;
    }
}

template<typename TStr>
void print_fusion(const Fusion& fusion,
		  int read_len,
		  TStr& left_ref_str,
		  TStr& right_ref_str,
		  const char* left_ref_name,
		  const char* right_ref_name,
		  ostream& fusion_db)
{  
  int half_splice_len = read_len - min_anchor_len;
  
  size_t left_start, right_start;
  size_t left_end, right_end;
  
  if (fusion.left >= 0 && fusion.left < length(left_ref_str) &&
      fusion.right >= 0 && fusion.right < length(right_ref_str))
    {
      if (fusion.dir == FUSION_FF || fusion.dir == FUSION_FR)
	{
	  left_start = fusion.left + 1 >= half_splice_len ? fusion.left - half_splice_len + 1 : 0;
	  left_end = left_start + half_splice_len;
	}
      else
	{
	  left_start = fusion.left;
	  left_end = left_start + half_splice_len < length(left_ref_str) ? left_start + half_splice_len : length(left_ref_str);
	}

      if (fusion.dir == FUSION_FF || fusion.dir == FUSION_RF)
	{
	  right_start = fusion.right;
	  right_end = right_start + half_splice_len < length(right_ref_str) ? right_start + half_splice_len : length(right_ref_str);
	}
      else
	{
	  right_end = fusion.right + 1;
	  right_start = right_end >= half_splice_len ? right_end - half_splice_len : 0;
	}

      seqan::Dna5String left_splice = infix(left_ref_str, left_start, left_end);
      seqan::Dna5String right_splice = infix(right_ref_str, right_start, right_end);

      if (fusion.dir == FUSION_RF || fusion.dir == FUSION_RR)
	{
	  seqan::reverseComplement(left_splice);
	  left_start = left_end - 1;
	}

      if (fusion.dir == FUSION_FR || fusion.dir == FUSION_RR)
	{
	  seqan::reverseComplement(right_splice);
	  right_end = right_start - 1;
	}

      const char* dir = "ff";
      if (fusion.dir == FUSION_FR)
	dir = "fr";
      else if (fusion.dir == FUSION_RF)
	dir = "rf";
      else if (fusion.dir == FUSION_RR)
	dir = "rr";
      
      fusion_db << ">" << left_ref_name << "-" << right_ref_name << "|"
		<< left_start << "|"
		<< fusion.left << "-" << fusion.right << "|"
		<< right_end << "|fus|" << dir <<  endl;
      
      fusion_db << left_splice << right_splice << endl;
    }
}


/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parse_oInt(int lower, char* arg, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			print_usage();
			exit(1);
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	print_usage();
	exit(1);
	return -1;
}
//
//int parse_options(int argc, char** argv)
//{
//	int option_index = 0;
//	int next_option; 
//	do { 
//		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);		
//		switch (next_option) {
//			case 'v':
//				verbose = true;
//				break;
//			case -1: /* Done with options. */
//				break;
//			default: 
//				print_usage();
//				return 1;
//		}
//	} while(next_option != -1);
//	
//	return 0;
//}

void get_seqs(istream& ref_stream,
	      RefSequenceTable& rt,
	      bool keep_seqs = true)
{    
    while(ref_stream.good() &&
          !ref_stream.eof())
    {
      RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence();
        string name;
        readMeta(ref_stream, name, Fasta());
	string::size_type space_pos = name.find_first_of(" \t\r");
	if (space_pos != string::npos)
	  {
	    name.resize(space_pos);
	  }
	fprintf(stderr, "\tLoading %s...", name.c_str());
	seqan::read(ref_stream, *ref_str, Fasta());
	fprintf(stderr, "done\n");
        rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
	if (!keep_seqs)
	  delete ref_str;
    }
}

void driver(const vector<FILE*>& splice_coords_files,
	    const vector<FILE*>& insertion_coords_files,
	    const vector<FILE*>& deletion_coords_files,
	    const vector<FILE*>& fusion_coords_files, 
	    ifstream& ref_stream)
{	
	char splice_buf[2048];
	RefSequenceTable rt(sam_header, true);
	get_seqs(ref_stream, rt, true);

	JunctionSet junctions;
	for (size_t i = 0; i < splice_coords_files.size(); ++i)
	{
		FILE* splice_coords = splice_coords_files[i];
		if (!splice_coords)
			continue;
		while (fgets(splice_buf, 2048, splice_coords))
		{
			char* nl = strrchr(splice_buf, '\n');
			char* buf = splice_buf;
			if (nl) *nl = 0;
			
			/**
			 Fields are:
			 1) reference name
			 2) left coord of splice (last char of the left exon)
			 3) right coord of splice (first char of the right exon)
			 */
			
			char* ref_name                   = get_token((char**)&buf, "\t");
			char* scan_left_coord            = get_token((char**)&buf, "\t");
			char* scan_right_coord           = get_token((char**)&buf, "\t");
			char* orientation				 = get_token((char**)&buf, "\t");
			
			if (!scan_left_coord || !scan_right_coord || !orientation)
			{
				fprintf(stderr,"Error: malformed splice coordinate record\n");
				exit(1);
			}
			uint32_t ref_id = rt.get_id(ref_name, NULL, 0);
			uint32_t left_coord = atoi(scan_left_coord);
			uint32_t right_coord = atoi(scan_right_coord);
			bool antisense = *orientation == '-';
			junctions.insert(make_pair<Junction, JunctionStats>(Junction(ref_id, left_coord, right_coord, antisense), JunctionStats()));
		}
	}


	/*
	 * Read in the deletion coordinates
	 * and store in a set
	 */	
	std::set<Deletion> deletions;
	for(size_t i=0; i < deletion_coords_files.size(); ++i){
		FILE* deletion_coords = deletion_coords_files[i];
		if(!deletion_coords){
			continue;
		} 
		while (fgets(splice_buf, 2048, deletion_coords))
		{
			char* nl = strrchr(splice_buf, '\n');
			char* buf = splice_buf;
			if (nl) *nl = 0;
			
			/**
			 Fields are:
			 1) reference name
			 2) left coord of splice (last char of the left exon)
			 3) right coord of splice (first char of the right exon)
			 */
			
			char* ref_name                   = get_token((char**)&buf, "\t");
			char* scan_left_coord            = get_token((char**)&buf, "\t");
			char* scan_right_coord           = get_token((char**)&buf, "\t");
			
			if (!scan_left_coord || !scan_right_coord)
			{
				fprintf(stderr,"Error: malformed deletion coordinate record\n");
				exit(1);
			}

			/*
			 * Note that when reading in a deletion, the left co-ord is the position of the 
			 * first deleted based. Since we are co-opting the junction data structure, need
			 * to fix up this location
			 */
			uint32_t ref_id = rt.get_id(ref_name, NULL, 0);
			uint32_t left_coord = atoi(scan_left_coord);
			uint32_t right_coord = atoi(scan_right_coord);
			deletions.insert(Deletion(ref_id, left_coord - 1, right_coord, false));
		}
	}

	/*
	 * Read in the insertion coordinates
	 * and store in a set
	 */
	std::set<Insertion> insertions;
	for(size_t i=0; i < insertion_coords_files.size(); ++i){
		FILE* insertion_coords = insertion_coords_files[i];
		if(!insertion_coords){
			continue;
		} 
		while(fgets(splice_buf, 2048, insertion_coords)){
			char* nl = strrchr(splice_buf, '\n');
			char* buf = splice_buf;
			if (nl) *nl = 0;
			
			char* ref_name = get_token((char**)&buf, "\t");
			char* scan_left_coord = get_token((char**)&buf, "\t");
			char* scan_right_coord = get_token((char**)&buf, "\t");
			char* scan_sequence = get_token((char**)&buf, "\t");

			if (!scan_left_coord || !scan_sequence || !scan_right_coord)
			{
				fprintf(stderr,"Error: malformed insertion coordinate record\n");
				exit(1);
			}
			
			seqan::Dna5String sequence = seqan::Dna5String(scan_sequence);
			bool containsN = false;
			for(size_t index = 0; index < seqan::length(sequence); index += 1){
				/*
				 * Don't allow any ambiguities in the insertion
				 */
				if(sequence[index] == 'N'){
					containsN = true;
					break;	
				}
			}
			if(containsN){
				continue;
			}
			seqan::CharString charSequence = sequence;
			uint32_t ref_id = rt.get_id(ref_name,NULL,0);
			uint32_t left_coord = atoi(scan_left_coord);
			insertions.insert(Insertion(ref_id, left_coord, seqan::toCString(charSequence)));
		}
	}

	std::set<Fusion> fusions;
	for(size_t i=0; i < fusion_coords_files.size(); ++i){
		FILE* fusion_coords = fusion_coords_files[i];
		if(!fusion_coords){
			continue;
		} 
		while(fgets(splice_buf, 2048, fusion_coords)){
			char* nl = strrchr(splice_buf, '\n');
			char* buf = splice_buf;
			if (nl) *nl = 0;
			
			char* ref_name1 = strsep((char**)&buf, "\t");
			char* scan_left_coord = strsep((char**)&buf, "\t");
			char* ref_name2 = strsep((char**)&buf, "\t");
			char* scan_right_coord = strsep((char**)&buf, "\t");
			char* scan_dir = strsep((char**)&buf, "\t");

			if (!ref_name1 || !scan_left_coord || !ref_name2 || !scan_right_coord || !scan_dir)
			{
			  fprintf(stderr,"Error: malformed insertion coordinate record\n");
			  exit(1);
			}
			
			uint32_t ref_id1 = rt.get_id(ref_name1, NULL, 0);
			uint32_t ref_id2 = rt.get_id(ref_name2, NULL, 0);
			uint32_t left_coord = atoi(scan_left_coord);
			uint32_t right_coord = atoi(scan_right_coord);
			uint32_t dir = FUSION_FF;
			if (strcmp(scan_dir, "fr") == 0)
			  dir = FUSION_FR;
			else if(strcmp(scan_dir, "rf") == 0)
			  dir = FUSION_RF;
			else if(strcmp(scan_dir, "rr") == 0)
			  dir = FUSION_RR;
		  
			fusions.insert(Fusion(ref_id1, ref_id2, left_coord, right_coord, dir));
		}
	}

	{
	  JunctionSet::iterator itr = junctions.begin();
	  for (; itr != junctions.end(); ++itr)
	    {
	      RefSequenceTable::Sequence* ref_str = rt.get_seq(itr->first.refid);
	      if (ref_str == NULL) continue;
	      
	      const char* name = rt.get_name(itr->first.refid);
	      print_splice(itr->first, read_length, itr->first.antisense ? "GTAG|rev" : "GTAG|fwd", *ref_str, name, cout);
	    }
	}

	{
	  std::set<Deletion>::iterator itr = deletions.begin();
	  for (; itr != deletions.end(); ++itr)
	    {
	      RefSequenceTable::Sequence* ref_str = rt.get_seq(itr->refid);
	      if (ref_str == NULL) continue;
	      
	      const char* name = rt.get_name(itr->refid);
	      print_splice((Junction)*itr, read_length, itr->antisense ? "del|rev" : "del|fwd", *ref_str, name, cout);
	    }
	}

	{
	  std::set<Insertion>::iterator itr  = insertions.begin();
	  for (; itr != insertions.end(); ++itr){
	    RefSequenceTable::Sequence* ref_str = rt.get_seq(itr->refid);
	    if (ref_str == NULL) continue;
	    
	    const char* name = rt.get_name(itr->refid);
	    print_insertion(*itr, read_length, *ref_str, name, cout);	
	  }
	}

	{
	  std::set<Fusion>::iterator itr = fusions.begin();
	  for (; itr != fusions.end(); ++itr){
	    RefSequenceTable::Sequence* left_ref_str = rt.get_seq(itr->refid1);
	    RefSequenceTable::Sequence* right_ref_str = rt.get_seq(itr->refid2);

	    if (left_ref_str == NULL || right_ref_str == NULL) continue;
	    
	    const char* left_ref_name = rt.get_name(itr->refid1);
	    const char* right_ref_name = rt.get_name(itr->refid2);
	    print_fusion(*itr, read_length, *left_ref_str, *right_ref_str, left_ref_name, right_ref_name, cout);	
	  }
	}
}

int main(int argc, char** argv)
{
	fprintf(stderr, "juncs_db v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "---------------------------\n");
	
	int parse_ret = parse_options(argc, argv, print_usage);
	if (parse_ret)
		return parse_ret;
	
    if(optind >= argc) 
    {
		print_usage();
		return 1;
	}
	
	min_anchor_len = parse_oInt(3, argv[optind++], "anchor length must be at least 3");
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}

	read_length = parse_oInt(4, argv[optind++], "read length must be at least 4");
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}
	
	string splice_coords_file_list = argv[optind++];
	vector<string> splice_coords_file_names;
	vector<FILE*> coords_files;
	tokenize(splice_coords_file_list, ",", splice_coords_file_names);
	for (size_t s = 0; s < splice_coords_file_names.size(); ++s)
	{
	
		FILE* coords_file = fopen(splice_coords_file_names[s].c_str(), "r");
		
		if (!coords_file)
		{
			fprintf(stderr, "Warning: cannot open %s for reading\n",
					splice_coords_file_names[s].c_str());
			continue;
		}
		coords_files.push_back(coords_file);
	}
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}

	/*
	 * Read in the insertion co-ordinates
	 */
	string insertion_coords_file_list = argv[optind++];
	vector<string> insertion_coords_file_names;
	vector<FILE*> insertion_coords_files;
	tokenize(insertion_coords_file_list, ",", insertion_coords_file_names);
	for(size_t s = 0; s < insertion_coords_file_names.size(); ++s)
	{
		FILE* insertion_coords_file = fopen(insertion_coords_file_names[s].c_str(),"r");
		if(!insertion_coords_file)
		{
			fprintf(stderr, "Warning: cannot open %s for reading\n",
					insertion_coords_file_names[s].c_str());
			continue;
		}
		insertion_coords_files.push_back(insertion_coords_file);
	}
	if(optind >= argc)
	{
		print_usage();
		return 1;
	}


	/*
	 * Read in the deletion co-ordinates
	 */
	string deletion_coords_file_list = argv[optind++];
	vector<string> deletion_coords_file_names;
	vector<FILE*> deletion_coords_files;
	tokenize(deletion_coords_file_list, ",", deletion_coords_file_names);
	for(size_t s = 0; s < deletion_coords_file_names.size(); ++s)
	{
		FILE* deletion_coords_file = fopen(deletion_coords_file_names[s].c_str(),"r");
		if(!deletion_coords_file)
		{
			fprintf(stderr, "Warning: cannot open %s for reading\n",
					deletion_coords_file_names[s].c_str());
			continue;
		}
		deletion_coords_files.push_back(deletion_coords_file);
	}
	if(optind >= argc)
	{
		print_usage();
		return 1;
	}


	/*
	 */
	string fusion_coords_file_list = argv[optind++];
	vector<string> fusion_coords_file_names;
	vector<FILE*> fusion_coords_files;
	tokenize(fusion_coords_file_list, ",", fusion_coords_file_names);
	for(size_t s = 0; s < fusion_coords_file_names.size(); ++s)
	{
		FILE* fusion_coords_file = fopen(fusion_coords_file_names[s].c_str(),"r");
		if(!fusion_coords_file)
		{
			fprintf(stderr, "Warning: cannot open %s for reading\n",
					fusion_coords_file_names[s].c_str());
			continue;
		}
		fusion_coords_files.push_back(fusion_coords_file);
	}
	if(optind >= argc)
	{
		print_usage();
		return 1;
	}

	
	string ref_file_name = argv[optind++];
	ifstream ref_stream(ref_file_name.c_str());
	
	if (!ref_stream.good())
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				ref_file_name.c_str());
		exit(1);
	}
    
	driver(coords_files, insertion_coords_files, deletion_coords_files, fusion_coords_files, ref_stream);
    return 0;
}
