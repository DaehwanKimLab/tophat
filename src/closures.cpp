/*
 *  closures.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/15/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bwt_map.h"
#include "inserts.h"
#include "closures.h"

bool possible_cotranscript(const BowtieHit& h1, const BowtieHit& h2, bool check_strand = true)
{
	if (h1.insert_id != h2.insert_id) 
		return false;
	int min_mate_inner_dist = insert_len - h1.read_len() - 
	h2.read_len() - insert_len_std_dev;
	if (max_mate_inner_dist == -1)
	{
		max_mate_inner_dist = insert_len - h1.read_len() - 
		h2.read_len() + insert_len_std_dev;
	}
	
	InsertAlignmentGrade grade(h1,h2, min_mate_inner_dist, max_mate_inner_dist);
	return (!grade.too_far && !grade.too_close && grade.opposite_strands);
}


void check_mates(const HitList& hits1_in_ref,
				 const HitList& hits2_in_ref,
				 vector<pair<size_t, size_t> >& happy_mates,
				 vector<size_t>& map1_singletons,
				 vector<size_t>& map2_singletons)
{
	std::set<size_t> marked;
	// TODO: if this shows up on the profile, replace it with a linear
	// time algorithm.  This one is 2*n*lg(n).
	HitList::const_iterator last_good = hits2_in_ref.begin();
	
	// Sanity checking to verify hits are sorted by insert id.
#if !NDEBUG		
	for (size_t i = 1; i < hits1_in_ref.size(); ++i)
	{
		const BowtieHit& h1 = hits1_in_ref[i];
		const BowtieHit& h2 = hits1_in_ref[i-1];
		assert(h1.insert_id >= h2.insert_id);
	}
	
	for (size_t i = 1; i < hits2_in_ref.size(); ++i)
	{
		const BowtieHit& h1 = hits2_in_ref[i];
		const BowtieHit& h2 = hits2_in_ref[i-1];
		assert(h1.insert_id >= h2.insert_id);
	}
#endif
	
	for (size_t i = 0; i < hits1_in_ref.size(); ++i)
	{
		pair<HitList::const_iterator, HitList::const_iterator> range_pair;
		range_pair = equal_range(last_good, hits2_in_ref.end(),
								 hits1_in_ref[i], hit_insert_id_lt);
		bool found_hit = false;
		if (range_pair.first != range_pair.second)
			last_good = range_pair.first;
		for (HitList::const_iterator f = range_pair.first;
			 f != range_pair.second;
			 ++f)
		{
			if (possible_cotranscript(hits1_in_ref[i], *f))
			{
				happy_mates.push_back(make_pair(i,f - hits2_in_ref.begin()));
				marked.insert(f - hits2_in_ref.begin());
				found_hit = true;
			}
		}
		if (!found_hit)
			map1_singletons.push_back(i);
	}
	
	for (size_t i = 0; i < hits2_in_ref.size(); ++i)
	{
		if (marked.find(i) == marked.end())
		{
			map2_singletons.push_back(i);
		}
	}	
}

void closure_driver(FILE* map1, FILE* map2, ifstream& ref_stream)
{
	typedef String< Dna5, Alloc<> > Reference;
	
	SequenceTable it(true);
	SequenceTable rt(true);
    HitTable hits1(it,rt);
    HitTable hits2(it,rt);
    
    get_mapped_reads(map1, hits1, false);
    get_mapped_reads(map2, hits2, false);
    
	uint32_t num_happy = 0;
	
	uint32_t num_potential_splices = 0;
	
	ofstream splice_db("possible_junc_db.fa");
	//ofstream splice_db = cout;
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;
		readMeta(ref_stream, name, Fasta());
		read(ref_stream, ref_str, Fasta());
		
		uint32_t ref_id = rt.get_id(name);
		
		if (verbose)
			fprintf(stderr, "Looking for happy mates in %s\n", name.c_str());
		
		CoordSet fwd_splices;
		CoordSet rev_splices;
		
		closure_search(ref_id, ref_str, hits1, hits2, read_length * 1.5, fwd_splices, rev_splices);
		
		print_splices(fwd_splices, read_length, "GTAG|fwd", ref_str, name, splice_db);
		print_splices(rev_splices, read_length, "GTAG|rev", ref_str, name, splice_db);
		
		num_potential_splices += rev_splices.size();
		num_potential_splices += fwd_splices.size();
		if (verbose)
		{
			fprintf(stderr, "\tFound %d possible splices\n", (int)rev_splices.size() + (int)fwd_splices.size());
		}
		//num_happy += happy_mates.size();
	}
	
    //fprintf(stderr, "Found %d happy mates\n", num_happy);
	fprintf(stderr, "Found %d total possible splices\n", num_potential_splices);
}
