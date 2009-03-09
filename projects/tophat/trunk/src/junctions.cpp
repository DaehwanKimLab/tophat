/*
 *  junctions.cpp
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
#include "junctions.h"

extern bool verbose;


void junctions_from_spliced_hit(BowtieHit& h, 
								vector<pair<Junction, JunctionStats> >& new_juncs)

{
	assert (!h.contiguous()); 
	const vector<CigarOp>& cigar = h.cigar();
	int j = h.left();
	
	for (size_t c = 0 ; c < cigar.size(); ++c)
	{
		Junction junc;
		JunctionStats stats;
		switch(cigar[c].opcode)
		{
			case REF_SKIP:
				
				junc.refid = h.ref_id();
				junc.left = j;
				junc.right = junc.left + cigar[c].length;
				junc.antisense = h.antisense_splice();
				
				
				assert (c > 0 && c < cigar.size() - 1); 
				
				stats.left_extent = cigar[c - 1].length;
				stats.right_extent = cigar[c + 1].length;
				
				stats.supporting_hits.insert(&h);
				new_juncs.push_back(make_pair(junc, stats));
				//fall through this case to MATCH is intentional
			case MATCH:
				j += cigar[c].length;
				break;
			default:
				break;
		}
	}

}

void print_junction(FILE* junctions_out, 
					const char* name, 
					const Junction& j, 
					const JunctionStats& s, 
					uint64_t junc_id)
{
	fprintf(junctions_out,
			"%s\t%d\t%d\tJUNC%08d\t%d\t%c\t%d\t%d\t255,0,0\t2\t%d,%d\t0,%d\n",
			name,
			j.left - s.left_extent,
			j.right + s.right_extent,

			(int)junc_id,
			(int)(s.supporting_hits.size()),

			j.antisense ? '-' : '+',
			j.left - s.left_extent,
			j.right + s.right_extent,
			s.left_extent,
			s.right_extent,
			j.right - (j.left - s.left_extent));
}


void junctions_from_alignment(BowtieHit& spliced_alignment,

							 JunctionSet& junctions)
{
	vector<pair<Junction, JunctionStats> > juncs;
	junctions_from_spliced_hit(spliced_alignment, juncs);
	
	for (size_t i = 0; i < juncs.size(); ++i)
	{
		pair<Junction, JunctionStats>& junc = juncs[i];
		JunctionSet::iterator itr = junctions.find(junc.first);
		
		if (itr != junctions.end())
		{
			JunctionStats& j = itr->second;
			j.left_extent = max(j.left_extent, junc.second.left_extent);
			j.right_extent = max(j.right_extent, junc.second.right_extent);
			j.supporting_hits.insert(&spliced_alignment);
		}
		else
		{
			assert(junc.first.refid != 0xFFFFFFFF);
			junctions[junc.first] = junc.second;
		}
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
void junctions_from_alignments(HitTable& hits,
							   JunctionSet& junctions)
{

	std::set<pair< int, int > > splice_coords;
	//JunctionSet raw_junctions;
	for (HitTable::iterator ci = hits.begin();
		 ci != hits.end();
		 ++ci)
	{
		HitList& rh = ci->second;
		if (rh.size() == 0)
			continue;
		for (size_t i = 0; i < rh.size(); ++i)
		{
			BowtieHit& bh = rh[i];
			AlignStatus s = status(&bh);
			total++;
			if (s == SPLICED)
				total_spliced++;
			if (!bh.accepted())
			{
				rejected++;
				if (s == SPLICED)
					rejected_spliced++;
				continue;
			}
			
			if (s == SPLICED)
			{
				junctions_from_alignment(bh, junctions); 
			}
		}
	}
	
	
}

bool accept_if_valid(const Junction& j, 
					 JunctionStats& s, 
					 const vector<short>& DoC,
					 double min_isoform_fraction)
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
	
	double avg_junc_doc = junc_doc / (double)(extent);
	
	//if (avg_junc_doc / (float) s.num_reads > 100.0)
	if (s.supporting_hits.size() / avg_junc_doc < min_isoform_fraction)
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
							double min_isoform_fraction)
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
					 RefSequenceTable& ref_sequences)
{
	uint64_t junc_id = 1;
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
					   ref_sequences.get_name(j.refid),
					   j,
					   s, 
					   junc_id++);
	}
	fprintf(stderr, "Rejected %d / %d alignments, %d / %d spliced\n", rejected, total, rejected_spliced, total_spliced);
}
