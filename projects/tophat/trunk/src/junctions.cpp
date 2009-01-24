/*
 *  junctions.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
 
#include <cassert>
#include "junctions.h"

extern bool verbose;

// This routine DOES NOT set the real refid!  
pair<Junction, JunctionStats> junction_from_spliced_hit(const BowtieHit& h)
{
	
	assert (h.splice_pos_left != -1 && h.splice_pos_right != -1); 
	Junction junc;
	junc.refid = 0xFFFFFFFF;
	junc.left = h.left + h.splice_pos_left;
	junc.right = h.right - h.splice_pos_right;
	junc.antisense = h.antisense_splice;
	
	JunctionStats stats;
	stats.left_extent = h.splice_pos_left;
	stats.right_extent = h.splice_pos_right;
	stats.num_reads = 1;
	return make_pair(junc, stats);
}

void print_junction(FILE* junctions_out, 
					const string& name, 
					const Junction& j, 
					const JunctionStats& s, 
					uint32_t junc_id)
{
	fprintf(stdout,
			"%s\t%d\t%d\tJUNC%08d\t%d\t%c\t%d\t%d\t255,0,0\t2\t%d,%d\t0,%d\n",
			name.c_str(),
			j.left - s.left_extent,
			j.right + s.right_extent - 1,
			junc_id,
			s.num_reads,
			j.antisense ? '-' : '+',
			j.left - s.left_extent,
			j.right + s.right_extent - 1,
			s.left_extent,
			s.right_extent,
			j.right - (j.left - s.left_extent));
}

void junction_from_alignment(const BowtieHit& spliced_alignment,
							 uint32_t refid,
							 JunctionSet& junctions)
{
	pair<Junction, JunctionStats> junc;
	junc = junction_from_spliced_hit(spliced_alignment);
	junc.first.refid = refid;
	JunctionSet::iterator itr = junctions.find(junc.first);
	
	if (itr != junctions.end())
	{
		JunctionStats& j = itr->second;
		j.left_extent = max(j.left_extent, junc.second.left_extent);
		j.right_extent = max(j.right_extent, junc.second.right_extent);
		j.num_reads++;
	}
	else
	{
		assert(junc.first.refid != 0xFFFFFFFF);
		junctions[junc.first] = junc.second;
	}
	
}

#if !NDEBUG
void validate_junctions(const JunctionSet& junctions)
{
	uint32_t invalid_juncs = 0;
	for (JunctionSet::const_iterator i = junctions.begin();
		 i != junctions.end();
		 ++i)
	{
		if (!i->first.valid())
			invalid_juncs++;
	}
	fprintf(stderr, "Found %d invalid junctions\n", invalid_juncs);
}
#endif

int rejected = 0;
int rejected_spliced = 0;
int total_spliced = 0;
int total = 0;
void junctions_from_alignments(const HitTable& hits,
							   JunctionSet& junctions)
{

	std::set<pair< int, int > > splice_coords;
	//JunctionSet raw_junctions;
	for (HitTable::const_iterator ci = hits.begin();
		 ci != hits.end();
		 ++ci)
	{
		const HitList& rh = ci->second;
		if (rh.size() == 0)
			continue;
		for (size_t i = 0; i < rh.size(); ++i)
		{
			AlignStatus s = status(&rh[i]);
			total++;
			if (s == SPLICED)
				total_spliced++;
			if (!rh[i].accepted)
			{
				rejected++;
				if (s == SPLICED)
					rejected_spliced++;
				continue;
			}
			
			if (s == SPLICED)
			{
				junction_from_alignment(rh[i], ci->first, junctions); 
			}
		}
	}
	
	
}

bool accept_if_valid(const Junction& j, 
					 JunctionStats& s, 
					 const vector<short>& DoC,
					 float min_isoform_fraction)
{
	uint32_t left_junc_doc = 0;
	uint32_t right_junc_doc = 0;
	
	for (uint32_t i = j.left - s.left_extent + 1; i <= j.left; ++i)
	{
		left_junc_doc += DoC[i];
	}
	for (uint32_t i = j.right; i < j.right + s.right_extent; ++i)
	{
		right_junc_doc += DoC[i];
	}
	
	uint32_t junc_doc = 0;
	uint8_t extent = 0;
	if (left_junc_doc > right_junc_doc)
	{
		junc_doc = left_junc_doc;
		extent = s.left_extent;
	}
	else
	{
		junc_doc = right_junc_doc;
		extent = s.right_extent;
	}
	
	float avg_junc_doc = junc_doc / (float)(extent);
	
	//if (avg_junc_doc / (float) s.num_reads > 100.0)
	if (s.num_reads / avg_junc_doc < min_isoform_fraction)
	{
		s.accepted = false;
	}
	else
	{
		s.accepted = true;
	}
	return s.accepted;
}

void accept_valid_junctions(JunctionSet& junctions,
							const uint32_t refid,
							const vector<short>& DoC,
							float min_isoform_fraction)
{
	Junction dummy_left(refid, 0, 0, true);
	Junction dummy_right(refid, 0xFFFFFFFF, 0xFFFFFFFF, true);
	
	pair<JunctionSet::iterator, JunctionSet::iterator> r;
	r.first = junctions.lower_bound(dummy_left);
	r.second = junctions.upper_bound(dummy_right);
	
	JunctionSet::iterator itr = r.first;

	while(itr != r.second && itr != junctions.end())
	{
		accept_if_valid(itr->first, itr->second, DoC, min_isoform_fraction);
		++itr;
	}
}

void accept_all_junctions(JunctionSet& junctions,
						  const uint32_t refid)
{
	for (JunctionSet::iterator itr = junctions.begin(); itr != junctions.end(); ++itr)
	{
		itr->second.accepted = true;
	}
}

void print_junctions(FILE* junctions_out, 
					 const JunctionSet& junctions,
					 SequenceTable& ref_sequences)
{
	uint32_t junc_id = 1;
	ref_sequences.invert();
	fprintf(junctions_out, "track name=junctions description=\"TopHat junctions\"\n");
	for (JunctionSet::const_iterator i = junctions.begin();
		 i != junctions.end();
		 ++i)
	{
		const pair<Junction, JunctionStats>& j_itr = *i; 
		const Junction& j = j_itr.first;
		const JunctionStats& s = j_itr.second;			
		if (!s.accepted)
			continue;
		assert(ref_sequences.get_name(j.refid));
		//fprintf(stdout,"%d\t%d\t%d\t%c\n", j.refid, j.left, j.right, j.antisense ? '-' : '+');
		
		print_junction(junctions_out, 
					   *(ref_sequences.get_name(j.refid)),
					   j,
					   s, 
					   junc_id++);
	}
	fprintf(stderr, "Rejected %d / %d alignments, %d / %d spliced\n", rejected, total, rejected_spliced, total_spliced);
}
