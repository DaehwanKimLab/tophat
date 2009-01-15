/*
 *  tophat_reports.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/20/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

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
#include "junctions.h"
#include "fragments.h"
#include "wiggles.h"

#ifdef PAIRED_END
#include "inserts.h"
#endif

using namespace std;
using namespace seqan;
using std::set;

const char *short_options = "r:I:d:s:va:AF:";

static int min_anchor_len = 5;
static int min_intron_len = 40;

static bool filter_junctions = true;
static float min_isoform_fraction = 0.15;
static bool accept_all = true;

bool verbose = false;

void print_usage()
{
#ifdef PAIRED_END
    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <map1.bwtout> [splice_map1.sbwtout] [map2.bwtout] [splice_map2.sbwtout]\n");
#else
	fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <map1.bwtout> [splice_map1.sbwtout]\n");
#endif
}

#ifdef PAIRED_END
void insert_best_pairings(SequenceTable& rt,
						  HitTable& hits1,
						  HitTable& hits2,
						  BestInsertAlignmentTable& best_pairings)						
{
	for(SequenceTable::const_iterator ci = rt.begin();
		ci != rt.end();
		++ci) 
	{		
		
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
		vector<size_t> map1_singletons;
		vector<size_t> map2_singletons;
		vector<pair<size_t, size_t> > happy_mates;
		
		string name = ci->first;
		uint32_t ref_id = rt.get_id(name);
		HitList* hits1_in_ref = hits1.get_hits(ref_id);
		HitList* hits2_in_ref = hits2.get_hits(ref_id);
		
		if (!hits1_in_ref || !hits2_in_ref)
			continue;
		
		if (verbose)
			fprintf(stderr, "Looking for best insert mappings in %s\n", name.c_str());
		
		best_insert_mappings(ref_id,
							 name, 
							 *hits1_in_ref, 
							 *hits2_in_ref, 
							 best_pairings);
	}	
}
#endif

void fragment_best_alignments(SequenceTable& rt,
							  HitTable& hits1,
							  BestFragmentAlignmentTable& best_alignments)						
{
	for(SequenceTable::const_iterator ci = rt.begin();
		ci != rt.end();
		++ci) 
	{		
		
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
				
		string name = ci->first;
		uint32_t ref_id = rt.get_id(name);
		HitList* hits_in_ref = hits1.get_hits(ref_id);
		
		if (!hits_in_ref)
			continue;
		
		if (verbose)
			fprintf(stderr, "Looking for best alignments in %s\n", name.c_str());
		
		best_fragment_mappings(ref_id,
							 name, 
							 *hits_in_ref, 
							 best_alignments);
	}	
}


void driver(FILE* map1, 
			FILE* splice_map1,
			FILE* map2, 
			FILE* splice_map2,
			FILE* coverage_out,
			FILE* junctions_out)
{
	
	bool paired_end = map2 && splice_map2;
	SequenceTable it(paired_end);
	SequenceTable rt(true);
    HitTable hits1(it,rt);
    HitTable hits2(it,rt);
    
    get_mapped_reads(map1, hits1, false);
	if (splice_map1)
		get_mapped_reads(splice_map1, hits1, true);
	
	JunctionSet junctions;
	
	print_wiggle_header(coverage_out);
	
	if (!paired_end)
	{
		fprintf(stderr, "Finished reading alignments\n");
		
		if (accept_all)
		{
			accept_all_hits(hits1);
		}
		else
		{
			BestFragmentAlignmentTable best_alignments(it.size());
			fragment_best_alignments(rt, hits1, best_alignments);
			
			accept_unique_hits(best_alignments);
		}
		junctions_from_alignments(hits1, junctions);
		
		for (SequenceTable::const_iterator ci = rt.begin();
			 ci != rt.end();
			 ++ci)
		{
			vector<short> DoC;
			const HitList* h1 = hits1.get_hits(ci->second);
			if (h1)
				add_hits_to_coverage(*h1, DoC);
			
			if (filter_junctions)
				accept_valid_junctions(junctions, ci->second, DoC, min_isoform_fraction);
			else
				accept_all_junctions(junctions, ci->second);
			
			print_wiggle_for_ref(coverage_out, ci->first, DoC);
		}
		
		print_junctions(junctions_out, junctions, rt);
	}
#ifdef PAIRED_END	
	else
	{
		get_mapped_reads(map2, hits2, false);
		get_mapped_reads(splice_map2, hits2, true);
		
		BestInsertAlignmentTable best_pairings(it.size());
		insert_best_pairings(rt, hits1, hits2, best_pairings);
		
		accept_valid_hits(best_pairings);
		
		junctions_from_alignments(hits1, junctions);
		junctions_from_alignments(hits2, junctions);
		
		print_junctions(junctions_out, junctions, rt);
		
		for (SequenceTable::const_iterator ci = rt.begin();
			 ci != rt.end();
			 ++ci)
		{
			vector<short> DoC;
			const HitList* h1 = hits1.get_hits(ci->second);
			if (h1)
				add_hits_to_coverage(*h1, DoC);
			
			const HitList* h2 = hits2.get_hits(ci->second);
			if (h2)
				add_hits_to_coverage(*h2, DoC);
			
			print_wiggle_for_ref(coverage_out, ci->first, DoC);
		}
    }
#endif
	
	uint32_t accepted_junctions = 0;
	for (JunctionSet::iterator itr = junctions.begin(); itr != junctions.end(); ++itr)
	{
		if(itr->second.accepted)
			accepted_junctions++;
	}
	
    fprintf(stderr, "Found %d junctions from happy spliced reads\n", accepted_junctions);
}

static struct option long_options[] = {
{"insert-len",      required_argument,       0,            'I'},
{"insert-stddev",      required_argument,       0,            's'},
{"read-len",       required_argument,       0,            'r'},
{"max-dist",       required_argument,       0,            'd'},
{"min-anchor",       required_argument,       0,            'a'},
{"min-intron",       required_argument,       0,            'i'},
{"min-isoform-fraction",       required_argument,       0,            'F'},
{"verbose",		no_argument,	0,							'v'},
{"accept-all-hits",      no_argument,       0,            'A'},
{0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
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

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloat(float lower, float upper, const char *errmsg) {
	float l;
	l = atof(optarg);
	
	if (l < lower) {
		cerr << errmsg << endl;
		print_usage();
		exit(1);
	}
	
	if (l > upper)
	{
		cerr << errmsg << endl;
		print_usage();
		exit(1);
	}
	
	return l;
	
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
			case 'd':
	   			max_mate_inner_dist = (uint32_t)parseInt(0, "-d/--max-dist arg must be at least 0");
	   			break;
			case 'I':
	   			insert_len = (uint32_t)parseInt(1, "-I/--insert-len arg must be at least 1");
	   			break;
			case 's':
	   			insert_len_std_dev = (uint32_t)parseInt(1, "-s/--insert-stddev arg must be at least 1");
	   			break;
			case 'v':
				verbose = true;
				break;
			case 'a':
				min_anchor_len = (uint32_t)parseInt(4, "-a/--min-anchor arg must be at least 4");
				break;
			case 'i':
				min_intron_len = (uint32_t)parseInt(1, "-a/--min-intron arg must be at least 1");
				break;
			case 'F':
				min_isoform_fraction = parseFloat(0.0, 1.0, "-a/--min-isoform-fraction arg must be [0.0,1.0]");
				if (min_isoform_fraction = 0.0)
				{
					filter_junctions = false;
				}
				break;
			case 'A':
				accept_all = true;
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
	
	string coverage_file_name = argv[optind++];
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}
	
	string junctions_file_name = argv[optind++];
	
    if(optind >= argc) 
    {
		print_usage();
		return 1;
	}
	
	string map1_file_name = argv[optind++];
	
	string splice_map1_file_name;
	string map2_file_name;
	string splice_map2_file_name;
	
	if (optind < argc) 
		splice_map1_file_name = argv[optind++];
	
	if (optind < argc)
		map2_file_name = argv[optind++];

	if (optind < argc) 
		splice_map2_file_name = argv[optind++];
		
	// Open the approppriate files

	FILE* coverage_file = fopen(coverage_file_name.c_str(), "w");
	if (coverage_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for writing\n",
				coverage_file_name.c_str());
		exit(1);
	}
	
	FILE* junctions_file = fopen(junctions_file_name.c_str(), "w");
	if (junctions_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for writing\n",
				junctions_file_name.c_str());
		exit(1);
	}
	
    FILE* map1_file = fopen(map1_file_name.c_str(), "r");
	if (map1_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map1_file_name.c_str());
		exit(1);
	}
	
	FILE* splice_map1_file = NULL;
	
	if (!splice_map1_file_name.empty())
	{
		splice_map1_file = fopen(splice_map1_file_name.c_str(), "r");
		if (splice_map1_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map1_file_name.c_str());
			exit(1);
		}
	}
	
	if ((splice_map2_file_name.empty() && !map2_file_name.empty()) ||
		(!splice_map2_file_name.empty() && map2_file_name.empty()))
	{
		fprintf(stderr, "Error: please specify both a contiguous and a spliced map file for paired-end reports\n");
		exit(1);
	}
	
	FILE* map2_file = NULL;
	FILE* splice_map2_file = NULL;
		
	if (!splice_map2_file_name.empty() && !map2_file_name.empty())
	{
		map2_file = fopen(map2_file_name.c_str(), "r");
		if (map2_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					map2_file_name.c_str());
			exit(1);
		}
		
		splice_map2_file = fopen(splice_map2_file_name.c_str(), "r");
		if (splice_map2_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map2_file_name.c_str());
			exit(1);
		}
	}
    
    driver(map1_file, 
		   splice_map1_file, 
		   map2_file,
		   splice_map2_file, 
		   coverage_file, 
		   junctions_file);
    
    return 0;
}


