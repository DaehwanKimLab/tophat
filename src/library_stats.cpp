/*
 *  library_stats.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
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
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "inserts.h"
#include "tokenize.h"

using namespace std;

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
//void print_mate_distances(const string& dist_file_name,
//						  const vector<InsertStatus>& best_status_for_inserts)
//{
//	ofstream dist_file(dist_file_name.c_str());
//	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
//	{
//		if (best_status_for_inserts[i].mask < SINGLETON)
//		{
//			dist_file << best_status_for_inserts[i].inner_dist << endl;
//		}
//	}
//}

// Compute some statistics about the mapping for this library, and return them
// in a LibraryStats object.



LibraryStats compute_stats(const BestInsertAlignmentTable best_status_for_inserts)
{
	LibraryStats best_pairing_stats;
	//best_pairing_stats.ref_region = "Best alignment pairings";

	// Analyze the table of best pairings for each insert and report the same 
	// statistics
	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
	{
		const pair<InsertAlignmentGrade, vector<InsertAlignment> >& st = best_status_for_inserts[i];
		
		for (size_t j = 0; j < st.second.size(); ++j)
		{
			
			InsertAlignment a = st.second[j];
			if (a.left_alignment && a.right_alignment)
			{
				pair<int, int> distances = pair_distances(*(a.left_alignment),
														  *(a.right_alignment));
				int inner_dist = distances.second;
				printf("%d\n", inner_dist);
			}
		}
	}
	return best_pairing_stats;
}

void load_hits(const vector<string>& filenames,
			   RefSequenceTable& rt,
			   ReadTable& it,
               HitTable& hits)
{
    for (size_t i = 0; i < filenames.size(); ++i)
    {
        const string& filename = filenames[i];
		HitFactory* hit_factory = NULL;
		
        if (filename.rfind(".sam") != string::npos)
		{
			hit_factory = new SAMHitFactory(it,rt);
		}
		else
		{
			hit_factory = new BowtieHitFactory(it,rt);
		}
		
		
        FILE* map = fopen(filename.c_str(), "r");
		
        if (map == NULL)
        {
            fprintf(stderr, "Error: could not open %s\n", filename.c_str());
            exit(1);
        }
        fprintf(stderr, "Loading hits from %s\n", filename.c_str());
        size_t num_hits_before_load = hits.total_hits();
        get_mapped_reads(map, hits, *hit_factory, true);
        fprintf(stderr, "Loaded %d hits from %s\n", 
                (int)hits.total_hits() - (int)num_hits_before_load, 
                filename.c_str());
		delete hit_factory;
    }
}

void load_hits(RefSequenceTable& rt,
			   ReadTable& it,
			   const string& left_read_maplist,
               HitTable& left_hits,
               const string& right_read_maplist,
               HitTable& right_hits)
{
    vector<string> left_filenames;
    tokenize(left_read_maplist, ",", left_filenames);
    load_hits(left_filenames, rt, it, left_hits);
	vector<string> right_filenames;
	tokenize(right_read_maplist, ",", right_filenames);
	load_hits(right_filenames, rt, it, right_hits);
}

void driver(const string& left_maps,
			const string& right_maps)
{
	// Load the set of left maps, and if provided, the set of right maps
    ReadTable it;
    RefSequenceTable rt(true);
	
	//HitFactory hit_factory(it,rt);
	
    HitTable left_hits;
    HitTable right_hits;
	
	load_hits(rt, it, left_maps, left_hits, right_maps, right_hits);
	
	BestInsertAlignmentTable best_pairings(it.size() + 1);
	
	insert_best_pairings(rt, it, left_hits, right_hits, best_pairings, true);
	
	// Print a simple listing of inner distances on the best insert mappings.
	// Useful if you want to build a histogram of insert sizes, or determine 
	// the mean and std dev of your library insert size empirically.
	//print_mate_distances("mate_dist.txt", insert_best_pairings);
	
	// Compute and report some statistics about the mapping for this library.
	LibraryStats best_pairing_stats = compute_stats(best_pairings);
	//cout << best_pairing_stats;
}


void print_usage()
{
    fprintf(stderr, "Usage:   library_stats <left_map1,...,left_mapN> <right_map1,...,right_mapN>\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

int main(int argc, char** argv)
{	
	int parse_ret = parse_options(argc, argv, print_usage);
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

    driver(map1_file_name, map2_file_name);
    
    return 0;
}
