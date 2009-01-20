/*
 *  library_stats.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
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
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"

using namespace std;

void print_usage()
{
    fprintf(stderr, "Usage:   library_stats <map1.bwtout> <map2.bwtout> [splice_map1.bwtout] [splice_map2.bwtout]\n");
}

struct LibraryStats
{
	LibraryStats() : 
		m1_aligned_reads(0), 
		m2_aligned_reads(0), 
		m1_singletons(0), 
		m2_singletons(0), 
		spliced_reads(0),
		too_far_inserts(0), 
		too_close_inserts(0), 
		same_strand_inserts(0), 
		happy_inserts(0) {}
	
	// A name for this statistics object
	string ref_region;
	
	// Read-level alignment counters
	uint32_t m1_aligned_reads;
	uint32_t m2_aligned_reads;
	uint32_t m1_singletons;
	uint32_t m2_singletons;
	uint32_t spliced_reads;
	
	// Insert category counts
	uint32_t too_far_inserts;
	uint32_t too_close_inserts;
	uint32_t same_strand_inserts;
	uint32_t happy_inserts;
	

	
	// Roll the stats for rhs into this object
	LibraryStats& operator+=(const LibraryStats& rhs)
	{
		m1_aligned_reads += rhs.m1_aligned_reads;
		m2_aligned_reads += rhs.m2_aligned_reads;
		m1_singletons += rhs.m1_singletons;
		m2_singletons += rhs.m2_singletons;
		
		too_far_inserts += rhs.too_far_inserts;
		too_close_inserts += rhs.too_close_inserts;
		same_strand_inserts += rhs.same_strand_inserts;
		happy_inserts += rhs.happy_inserts;
		
		spliced_reads += rhs.spliced_reads;
		
		return *this;
	}
		
	uint32_t contiguously_aligned_reads() const
	{
		return m2_aligned_reads + m1_aligned_reads - spliced_reads;
	}
};


ostream& operator<<(ostream& os, const LibraryStats& L)
{
	os << L.ref_region << ":" << endl;
	os << "\t" << "Aligned reads:\t" << L.m2_aligned_reads + L.m1_aligned_reads << endl;
	os << "\t" << "Singletons:\t" << L.m2_singletons + L.m1_singletons << endl;
	os << "\t" << "Happy reads:\t" << 2 * L.happy_inserts << endl;
	os << "\t" << "Unhappy inserts, too far:\t" << 2 * L.too_far_inserts << endl;
	os << "\t" << "Unhappy inserts, too close:\t" << 2 * L.too_close_inserts << endl;
	os << "\t" << "Unhappy inserts, same strand:\t" << 2 * L.same_strand_inserts << endl;
	os << "\t" << "Spliced read alignments:\t" << L.spliced_reads << endl;
	os << "\t" << "Contiguous read alignments:\t" << L.contiguously_aligned_reads() << endl;
	
	return os;
}


// Print a simple listing of inner distances on the best insert mappings.
// Useful if you want to build a histogram of insert sizes, or determine 
// the mean and std dev of your library insert size empirically.
void print_mate_distances(const string& dist_file_name,
						  const vector<InsertStatus>& best_status_for_inserts)
{
	ofstream dist_file(dist_file_name.c_str());
	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
	{
		if (best_status_for_inserts[i].mask < SINGLETON)
		{
			dist_file << best_status_for_inserts[i].inner_dist << endl;
		}
	}
}

// Compute some statistics about the mapping for this library, and return them
// in a LibraryStats object.

//TODO: eliminate this functions dependence on the AlignStatus enum, unify with
// status masks
LibraryStats compute_stats(const vector<InsertStatus>& best_status_for_inserts)
{
	LibraryStats best_pairing_stats;
	best_pairing_stats.ref_region = "Best alignment pairings";

	// Analyze the table of best pairings for each insert and report the same 
	// statistics
	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
	{
		switch (best_status_for_inserts[i].mask)
		{
			case HAPPY:
				best_pairing_stats.happy_inserts++; break;
			case TOO_FAR:
				best_pairing_stats.too_far_inserts++; break;
			case TOO_CLOSE:
				best_pairing_stats.too_close_inserts++; break;
			case WRONG_ORIENTATION:
				best_pairing_stats.same_strand_inserts++; break;
			default: break;
		}
		
		// Now integrate the read-level stastitics (spliced/contiguous, etc)
		// for this insert into the LibraryStats object.
		const InsertStatus& st = best_status_for_inserts[i];
		for (size_t j = 0; j < st.best_alignments.size(); ++j)
		{
			const InsertAlignment& a = st.best_alignments[j];
			AlignStatus left_status = status(a.left_alignment);
			AlignStatus right_status = status(a.right_alignment);
			
			// If just one read is aligned, count it as a singleton
			if (left_status == UNALIGNED && right_status != UNALIGNED)
			{
				best_pairing_stats.m1_aligned_reads++;
				best_pairing_stats.m1_singletons++;
			}
			else if (right_status == UNALIGNED && left_status != UNALIGNED)
			{
				best_pairing_stats.m2_aligned_reads++;
				best_pairing_stats.m2_singletons++;
			}
			else
			{
				best_pairing_stats.m1_aligned_reads++;
				best_pairing_stats.m2_aligned_reads++;
			}
			if (left_status == SPLICED)
				best_pairing_stats.spliced_reads++;
			if (right_status == SPLICED)
				best_pairing_stats.spliced_reads++;
			
		}
			 
	}
	return best_pairing_stats;
}

void driver(FILE* map1, FILE* map2, FILE* splice_map1 = NULL, FILE* splice_map2 = NULL)
{
	// A sequence table for the inserts
	SequenceTable it;
	
	// A sequence table for the reference contigs
	SequenceTable rt;
	
	// The table of hits to reference contigs by inserts from the left ends of 
	// the library
    HitTable hits1(it,rt);
	
	// The table of hits to reference contigs by inserts from the right ends of
	// the library
    HitTable hits2(it,rt);

	// Populate the left end hit table from a genomic bowtie map of the left 
	// ends
	get_mapped_reads(map1, hits1, false);
	
	// Populate the left end hit table from a genomic bowtie map of the left 
	// ends
    get_mapped_reads(map2, hits2, false);
	
	// If we have valid splices files for both ends of the library as well,
	// then the user is running the library stats program after having 
	// generated possible splice junctions, and mapped reads to them.
	if (splice_map1 && splice_map2)
	{
		get_mapped_reads(splice_map1, hits1, true);
		get_mapped_reads(splice_map2, hits2, true);
    }
	
	// A table to hold the best alignment(s) for each insert
	vector<InsertStatus> insert_best_pairings(it.size());
	
	// Fill the above best hit table by examining hits to each reference 
	// contig independently.
	for(SequenceTable::const_iterator ci = rt.begin();
		ci != rt.end();
		++ci) 
	{
		string name = ci->first;
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
		vector<size_t> map1_singletons;
		vector<size_t> map2_singletons;
		vector<pair<size_t, size_t> > happy_mates;
		
		// Get the unique internal ID of this reference contig
		uint32_t ref_id = ci->second;
		
		// Get the left and right end hits to this reference contig
		const HitList* p_hits1_in_ref = hits1.get_hits(ref_id);
		const HitList* p_hits2_in_ref = hits2.get_hits(ref_id);
		
		// If we have hits for both ends to this contig, we can classify
		// inserts by happiness, etc.  If not, skip this contig
		if (!p_hits1_in_ref || !p_hits2_in_ref)
			continue;
		const HitList& hits1_in_ref = *p_hits1_in_ref;
		const HitList& hits2_in_ref = *p_hits2_in_ref;
		
		// For each insert, determine a set of best alignments, if any, to this
		// contig
		best_insert_mappings(ref_id,
							 name, 
							 hits1_in_ref, 
							 hits2_in_ref, 
							 insert_best_pairings);
	}
	
	// Print a simple listing of inner distances on the best insert mappings.
	// Useful if you want to build a histogram of insert sizes, or determine 
	// the mean and std dev of your library insert size empirically.
	print_mate_distances("mate_dist.txt", insert_best_pairings);
	
	// Compute and report some statistics about the mapping for this library.
	LibraryStats best_pairing_stats = compute_stats(insert_best_pairings);
	cout << best_pairing_stats;
}

const char *short_options = "r:I:d:s:";
static struct option long_options[] = {
{"insert-len",      required_argument,       0,            'I'},
{"insert-stddev",      required_argument,       0,            's'},
{"read-len",       required_argument,       0,            'r'},
{"max-dist",       required_argument,       0,            'd'},
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
	
	string map1_file_name = argv[optind++];
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}
	
	string map2_file_name = argv[optind++];
	
    FILE* map1_file = fopen(map1_file_name.c_str(), "r");
	if (map1_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map1_file_name.c_str());
		exit(1);
	}
	
	
	FILE* map2_file = fopen(map2_file_name.c_str(), "r");
	if (map2_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map2_file_name.c_str());
		exit(1);
	}
	
	FILE* splice_map1_file = NULL;
	FILE* splice_map2_file = NULL;
	
	if (optind + 1< argc)
	{
		string splice_map1_file_name = argv[optind++];
		splice_map1_file = fopen(splice_map1_file_name.c_str(), "r");
		if (splice_map1_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map1_file_name.c_str());
			exit(1);
		}
		
		string splice_map2_file_name = argv[optind++];
		splice_map2_file = fopen(splice_map2_file_name.c_str(), "r");
		if (splice_map2_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map1_file_name.c_str());
			exit(1);
		}
	}
	    
    driver(map1_file, map2_file, splice_map1_file, splice_map2_file);
    
    return 0;
}
