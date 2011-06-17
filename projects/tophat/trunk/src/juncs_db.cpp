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
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"

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
  if (insertion.left >= 0 && insertion.left <= length(ref_str))
    {
      left_start = (int)insertion.left - half_splice_len + 1 >= 0 ? (int)insertion.left - half_splice_len + 1 : 0;
      left_end = left_start + half_splice_len;
      
      right_start = (int)left_end; 
      right_end = right_start + half_splice_len  < length(ref_str) ? right_start + half_splice_len : length(ref_str);
      
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
  
  if (junction.left >= 0 && junction.left <= length(ref_str) &&
      junction.right >= 0 && junction.right <= length(ref_str))
    {
      left_start = (int)junction.left - half_splice_len + 1 >= 0 ? (int)junction.left - half_splice_len + 1 : 0;
      left_end = left_start + half_splice_len;
      
      right_start = junction.right;
      right_end = right_start + half_splice_len < length(ref_str) ? right_start + half_splice_len : length(ref_str) - right_start;
      
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

void driver(const vector<FILE*>& splice_coords_files,
			const vector<FILE*>& insertion_coords_files,
			const vector<FILE*>& deletion_coords_files, 
			ifstream& ref_stream)
{	
	char splice_buf[2048];
	RefSequenceTable rt(true);
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


	typedef RefSequenceTable::Sequence Reference;
	
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;

		readMeta(ref_stream, name, Fasta());
		string::size_type space_pos = name.find_first_of(" \t\r");
		if (space_pos != string::npos)
		{
			name.resize(space_pos);
		}
		
		read(ref_stream, ref_str, Fasta());
		
		uint32_t refid = rt.get_id(name, NULL, 0);
		Junction dummy_left(refid, 0, 0, true);
		Junction dummy_right(refid, 0xFFFFFFFF, 0xFFFFFFFF, true);
		
		pair<JunctionSet::iterator, JunctionSet::iterator> r;
		r.first = junctions.lower_bound(dummy_left);
		r.second = junctions.upper_bound(dummy_right);
		
		JunctionSet::iterator itr = r.first;
		
		while(itr != r.second && itr != junctions.end())
		{
			print_splice(itr->first, read_length, itr->first.antisense ? "GTAG|rev" : "GTAG|fwd", ref_str, name, cout);
			++itr;
		}
	}


	ref_stream.clear();
	ref_stream.seekg(0, ios::beg);


	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;

		readMeta(ref_stream, name, Fasta());
		string::size_type space_pos = name.find_first_of(" \t\r");
		if (space_pos != string::npos)
		{
			name.resize(space_pos);
		}
		
		read(ref_stream, ref_str, Fasta());
		
		uint32_t refid = rt.get_id(name, NULL,0);
		Deletion dummy_left(refid, 0, 0, true);
		Deletion dummy_right(refid, 0xFFFFFFFF, 0xFFFFFFFF, true);
		
		pair<std::set<Deletion>::iterator, std::set<Deletion>::iterator> r;
		r.first = deletions.lower_bound(dummy_left);
		r.second = deletions.upper_bound(dummy_right);
		
		std::set<Deletion>::iterator itr = r.first;
		
		while(itr != r.second && itr != deletions.end())
		{
			print_splice((Junction)*itr, read_length, itr->antisense ? "del|rev" : "del|fwd", ref_str, name, cout);
			++itr;
		}
	}

	ref_stream.clear();
	ref_stream.seekg(0, ios::beg);



	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;

		readMeta(ref_stream, name, Fasta());
		string::size_type space_pos = name.find_first_of(" \t\r");
		if (space_pos != string::npos)
		{
			name.resize(space_pos);
		}
		
		read(ref_stream, ref_str, Fasta());
		
		uint32_t refid = rt.get_id(name, NULL,0);
		Insertion dummy_left(refid, 0, "");
		Insertion dummy_right(refid, 0xFFFFFFFF, "");
	
		std::set<Insertion>::iterator itr = insertions.lower_bound(dummy_left);
		std::set<Insertion>::iterator upper   = insertions.upper_bound(dummy_right);

		while(itr != upper && itr != insertions.end()){
			print_insertion(*itr, read_length, ref_str, name, cout);	
			++itr;
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


	
	string ref_file_name = argv[optind++];
	ifstream ref_stream(ref_file_name.c_str());
	
	if (!ref_stream.good())
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				ref_file_name.c_str());
		exit(1);
	}
    
	driver(coords_files, insertion_coords_files, deletion_coords_files, ref_stream);
    return 0;
}
