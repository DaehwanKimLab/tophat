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
//#include "closures.h"
#include "tokenize.h"
#include "junctions.h"

#ifdef PAIRED_END
#include "inserts.h"
#endif

using namespace std;
using namespace seqan;
using std::set;


// Length of the outer dimension of a single insert from the paired-end library

static bool verbose = false;
static int read_length = -1;
void print_usage()
{
    fprintf(stderr, "Usage:   juncs_db <min_anchor> <read_length> <splice_coords1,...,splice_coordsN> <ref.fa>\n");
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

template<typename TStr>
void print_splice(const Junction& junction,
				  int read_len,
				  const string& tag,
				  TStr& ref_str,
				  const string& ref_name,
				  ostream& splice_db)
{
	int half_splice_len = read_len - min_anchor_len;

	size_t left_start, right_start;
	uint32_t left_end, right_end;
	
	left_start = (int)junction.left - half_splice_len + 1 >= 0 ? (int)junction.left - half_splice_len + 1 : 0;
	left_end = left_start + half_splice_len;
	
	right_start = junction.right;
	right_end = right_start + half_splice_len < length(ref_str) ? right_start + half_splice_len : length(ref_str) - right_start;
	
	
	
	Infix<String<Dna5, Alloc<> > >::Type left_splice = infix(ref_str,
															 left_start, 
															 left_end);
	Infix<String<Dna5, Alloc<> > >::Type right_splice = infix(ref_str, 
															  right_start, 
															  right_end);
	
	splice_db << ">" << ref_name << "|" << left_start << "|" << junction.left <<
		"-" << junction.right << "|" << right_end << "|" << tag << endl;

	splice_db << left_splice << right_splice << endl;
		
}


//const char *short_options = "r:I:d:s:va:e:i:";
const char *short_options = "v";

static struct option long_options[] = {
{"verbose",		no_argument,	0,							'v'},
//{"insert-len",      required_argument,       0,            'I'},
//{"insert-stddev",      required_argument,       0,            's'},
//{"max-dist",       required_argument,       0,            'd'},
//{"min-intron",       required_argument,       0,            'i'},
//{"min-exon",       required_argument,       0,            'e'},
{0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, char* arg, const char *errmsg) {
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

int parse_options(int argc, char** argv)
{
	int option_index = 0;
	int next_option; 
	do { 
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);		
		switch (next_option) {
			case 'v':
				verbose = true;
				break;
			case -1: /* Done with options. */
				break;
			default: 
				print_usage();
				return 1;
		}
	} while(next_option != -1);
	
	return 0;
}

void driver(const vector<FILE*>& splice_coords_files, 
			ifstream& ref_stream)
{
	char splice_buf[2048];
	SequenceTable rt(true);
	JunctionSet junctions;
	for (size_t i = 0; i < splice_coords_files.size(); ++i)
	{
		FILE* splice_coords = splice_coords_files[i];
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
			
			char* ref_name                   = strsep((char**)&buf, "\t");
			char* scan_left_coord            = strsep((char**)&buf, "\t");
			char* scan_right_coord           = strsep((char**)&buf, "\t");
			char* orientation				 = strsep((char**)&buf, "\t");
			
			if (!scan_left_coord || !scan_right_coord || !orientation)
			{
				fprintf(stderr,"Error: malformed splice coordinate record\n");
				exit(1);
			}
			uint32_t ref_id = rt.get_id(ref_name);
			uint32_t left_coord = atoi(scan_left_coord);
			uint32_t right_coord = atoi(scan_right_coord);
			bool antisense = *orientation == '-';
			junctions.insert(make_pair<Junction, JunctionStats>(Junction(ref_id, left_coord, right_coord, antisense), JunctionStats()));
		}
	}
	
	typedef String< Dna5, Alloc<> > Reference;
	
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;
		readMeta(ref_stream, name, Fasta());
		read(ref_stream, ref_str, Fasta());
		
		uint32_t refid = rt.get_id(name);
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
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
	if (parse_ret)
		return parse_ret;
	
    if(optind >= argc) 
    {
		print_usage();
		return 1;
	}
	
	min_anchor_len = parseInt(3, argv[optind++], "read length must be at least 3");
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}

	read_length = parseInt(20, argv[optind++], "read length must be at least 20");
	
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
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_coords_file_names[s].c_str());
			exit(1);
		}
		coords_files.push_back(coords_file);
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
    
	driver(coords_files, ref_stream);
    return 0;
}
