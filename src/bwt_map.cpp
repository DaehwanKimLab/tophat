/*
 *  bwt_map.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <iostream>
#include <set>
#include <vector>

#include "common.h"
#include "bwt_map.h"

using namespace std;



bool hit_insert_id_lt(const BowtieHit& h1, const BowtieHit& h2)
{
	return h1.insert_id < h2.insert_id;
}

void get_mapped_reads(FILE* bwtf, HitTable& hits, bool spliced, bool verbose)
{
    const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";
	static const int buf_size = 256;
	char orientation;
	char name[buf_size];
	int bwtf_ret = 0;
	//uint32_t seqid = 0;
	char text_name[buf_size];
	unsigned int text_offset;
	char sequence[buf_size];
	
	char qualities[buf_size];
	unsigned int other_occs;
	char mismatches[buf_size];
	
    char bwt_buf[2048];
	uint32_t reads_extracted = 0;
	
	while (fgets(bwt_buf, 2048, bwtf))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		memset(mismatches, 0, sizeof(mismatches));
		// Get a new record from the tab-delimited Bowtie map
		bwtf_ret = sscanf(bwt_buf,
						  bwt_fmt_str,
						  name,
						  &orientation,
						  text_name,   // name of reference sequence
						  &text_offset,
						  sequence,
						  qualities,
						  &other_occs,
						  mismatches);
		
		// If we didn't get enough fields, this record is bad, so skip it
		if (bwtf_ret > 0 && bwtf_ret < 6)
		{
			fprintf(stderr, "Warning: found malformed record, skipping\n");
			continue;
		}
        
		// Stripping the slash and number following it gives the insert name
        char* slash = strrchr(name, '/');
        if (slash)
        {
            *slash = 0;
		}
		int read_len = strlen(sequence);
		
		// Add this alignment to the table of hits for this half of the
		// Bowtie map
		if (spliced)
		{
			// Parse the text_name field to recover the splice coords
			vector<string> toks;
			char* pch = strtok (text_name,"|-");
			while (pch != NULL)
			{
				toks.push_back(pch);
				pch = strtok (NULL, "|-");
			}
			
			if (toks.size() == 7)
			{
				string contig = toks[0];
				
				uint32_t left = atoi(toks[1].c_str()) + text_offset;
				uint32_t spliced_read_len = strlen(sequence);
				int8_t left_splice_pos = atoi(toks[2].c_str()) - left;
				int8_t right_splice_pos = spliced_read_len - left_splice_pos;
				
				uint32_t right = atoi(toks[3].c_str()) + right_splice_pos;
				atoi(toks[4].c_str());
				
				assert (string(toks[6]) == "rev" || string(toks[6]) == "fwd");
				assert (orientation == '-' || orientation == '+');
				
				//vector<string> mismatch_toks;
				char* pch = strtok (mismatches,",");
				bool mismatch_in_anchor = false;
				while (pch != NULL)
				{
					char* colon = strchr(pch, ':');
					if (colon) 
					{
						*colon = 0;
						int mismatch_pos = atoi(pch);
						if ((orientation == '+' && abs(mismatch_pos - left_splice_pos) < 5) ||
							(orientation == '-' && abs(((int)spliced_read_len - left_splice_pos + 1) - mismatch_pos)) < 5)
							mismatch_in_anchor = true;
					}
					//mismatch_toks.push_back(pch);
					pch = strtok (NULL, ",");
				}
				
				if (!mismatch_in_anchor)
				{
					hits.add_spliced_hit(name, 
										 contig, 
										 left, 
										 right, 
										 left_splice_pos, 
										 right_splice_pos, 
										 spliced_read_len, 
										 orientation == '-', 
										 string(toks[6]) == "rev");
				}
			}
		}
		else
		{
			hits.add_hit(name, text_name, text_offset, read_len, orientation == '-');
		}
		reads_extracted++;
	}
	
	// This will sort the map by insert id.
	hits.finalize();
	if (verbose)
	{
        fprintf(stderr, "Extracted %d alignments from Bowtie map\n", reads_extracted);
    }
}

pair<int, int> pair_distances(const BowtieHit& h1, 
							  const BowtieHit& h2)
{
	int minor_hit_start, major_hit_start;
	int minor_hit_end, major_hit_end;
	if (h1.left < h2.left)
	{
		minor_hit_start = (int)h1.left;
		minor_hit_end = (int)h1.right;
		major_hit_start = (int)h2.left;
		major_hit_end = (int)h2.right;
	}
	else
	{
		minor_hit_start = (int)h2.left;
		minor_hit_end = (int)h2.right;
		major_hit_start = (int)h1.left;
		major_hit_end = (int)h1.right;
	}
	
	int inner_dist = major_hit_start - minor_hit_end;
	int outer_dist = major_hit_end - minor_hit_start;
	return make_pair(outer_dist, inner_dist);
}

AlignStatus status(const BowtieHit* align)
{
	if (!align)
		return UNALIGNED;
	if (align->splice_pos_left == -1)
		return CONTIGUOUS;
	return SPLICED;
}

//bool left_status_better(MateStatusMask left, MateStatusMask right)
//{
//	return left < right;
//}
//
//bool status_equivalent(MateStatusMask left, MateStatusMask right)
//{	
//	return left == right;
//}

void best_insert_mappings(uint32_t refid,
						  const string& name,
						  HitList& hits1_in_ref,
						  HitList& hits2_in_ref,
						  BestInsertAlignmentTable& best_status_for_inserts)
{	
	std::set<size_t> marked;
	HitList::iterator last_good = hits2_in_ref.begin();
	
	for (size_t i = 0; i < hits1_in_ref.size(); ++i)
	{
		BowtieHit& h1 = hits1_in_ref[i];
		pair<HitList::iterator, HitList::iterator> range_pair;
		range_pair = equal_range(last_good, hits2_in_ref.end(),
								 h1, hit_insert_id_lt);
		bool found_hit = false;
		if (range_pair.first != range_pair.second)
			last_good = range_pair.first;

		for (HitList::iterator f = range_pair.first;
			 f != range_pair.second;
			 ++f)
		{
			BowtieHit& h2 = *f;

			if (h1.insert_id == h2.insert_id)
			{
				uint32_t insert_id = h2.insert_id;
				
				int min_mate_inner_dist = insert_len - h1.read_len() - 
				h2.read_len() - insert_len_std_dev;
				if (max_mate_inner_dist == -1)
				{
					max_mate_inner_dist = insert_len - h1.read_len() - 
					h2.read_len() + insert_len_std_dev;
				}
				InsertAlignmentGrade s(h1, h2, max_mate_inner_dist, min_mate_inner_dist);
				
				pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
					= best_status_for_inserts[insert_id];
				InsertAlignmentGrade& current = insert_best.first;
				// Is the new status better than the current best one?
				if (current < s)
				{
					insert_best.second.clear();
					current = s;
					insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
				}
				else if (! (s < current))
				{
					insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
				}
				
				marked.insert(f - hits2_in_ref.begin());
				found_hit = true;
			}
			
		}
		if (!found_hit)
		{
			pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
				= best_status_for_inserts[h1.insert_id];
			InsertAlignmentGrade& current = insert_best.first;	
			
			InsertAlignmentGrade s(h1);
			
			if (current < s)
			{
				insert_best.second.clear();
				current = s;
				insert_best.second.push_back(InsertAlignment(refid, &h1, NULL));
			}
			else if (! (s < current))
			{
				insert_best.second.push_back(InsertAlignment(refid, &h1, NULL));
			}
			
		}
	}
	
	for (size_t i = 0; i < hits2_in_ref.size(); ++i)
	{
		BowtieHit& h2 = hits2_in_ref[i];

		pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
			= best_status_for_inserts[h2.insert_id];
		InsertAlignmentGrade& current = insert_best.first;	
		
		InsertAlignmentGrade s(h2);
		// Did we include h2 as part of a pairing already, or is this first time
		// we've seen it?  If so, it's a singleton.
		if (marked.find(i) == marked.end())
		{
			if (current < s)
			{
				insert_best.second.clear();
				current = s;
				insert_best.second.push_back(InsertAlignment(refid, NULL, &h2));
			}
			else if (! (s < current))
			{
				insert_best.second.push_back(InsertAlignment(refid, NULL, &h2));
			}
		}
	}	
	
}

void best_fragment_mappings(uint32_t refid,
							const string& name,
							HitList& hits_in_ref,
							BestFragmentAlignmentTable& best_status_for_fragments)
{	
	for (size_t i = 0; i < hits_in_ref.size(); ++i)
	{
		BowtieHit& h1 = hits_in_ref[i];
		
		uint32_t fragment_id = h1.insert_id;
				
		FragmentAlignmentGrade s(h1);
				
		pair<FragmentAlignmentGrade, vector<FragmentAlignment> >& fragment_best
				= best_status_for_fragments[fragment_id];
		FragmentAlignmentGrade& current = fragment_best.first;
		// Is the new status better than the current best one?
		if (current < s)
		{
			fragment_best.second.clear();
			current = s;
			fragment_best.second.push_back(FragmentAlignment(refid, &h1));
		}
		else if (! (s < current)) // is it just as good?
		{
			fragment_best.second.push_back(FragmentAlignment(refid, &h1));
		}
	}
}


int long_spliced = 0;
int short_spliced = 0;
int singleton_splices = 0;


bool valid_insert_alignment(const InsertAlignmentGrade& g, const InsertAlignment& a)
{
	if (!a.left_alignment || !a.right_alignment)
		return false;
	if (g.both_mapped)
	{ 
		// Take all the contiguously mapped pairs
		if (!g.one_spliced && !g.both_spliced)
			return true;
	
		// Take the pairs that include one or more spliced reads as long as 
		// the inner dist isn't too big
		if (g.one_spliced || g.both_spliced)
		{
			if (g.too_far || g.too_close)
				return false;
		}
		return true;
	}
	return false;
}

bool valid_fragment_alignment(const FragmentAlignmentGrade& g, const FragmentAlignment& a)
{
	// stub
	return true;
}

void accept_valid_hits(BestInsertAlignmentTable& best_status_for_inserts)
{
	for (size_t i = 0; i < best_status_for_inserts.size(); ++i)
	{
		const pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
			= best_status_for_inserts[i];
		for (size_t j = 0; j < insert_best.second.size(); ++j)
		{
			const InsertAlignment& a = insert_best.second[j];
			bool valid = valid_insert_alignment(insert_best.first, a);
			if (a.left_alignment)
			{
				a.left_alignment->accepted = valid;
			}
			if (a.right_alignment)
			{
				a.right_alignment->accepted = valid;
			}
		}
	}
}

void accept_unique_hits(BestFragmentAlignmentTable& best_status_for_fragments)
{
	for (size_t i = 0; i < best_status_for_fragments.size(); ++i)
	{
		const pair<FragmentAlignmentGrade, vector<FragmentAlignment> >& fragment_best
			= best_status_for_fragments[i];
		if (fragment_best.first.spliced)
		{
			if (fragment_best.second.size() == 1)
				fragment_best.second[0].alignment->accepted = true;
			else
				fragment_best.second[0].alignment->accepted = false;
		}
		else
		{
			for (size_t j = 0; j < fragment_best.second.size(); ++j)
			{
				const FragmentAlignment& a = fragment_best.second[j];
				a.alignment->accepted = true;
			}
		}
			
	}
}

void accept_valid_hits(BestFragmentAlignmentTable& best_status_for_fragments)
{
	for (size_t i = 0; i < best_status_for_fragments.size(); ++i)
	{
		const pair<FragmentAlignmentGrade, vector<FragmentAlignment> >& fragment_best
		= best_status_for_fragments[i];

		for (size_t j = 0; j < fragment_best.second.size(); ++j)
		{
			const FragmentAlignment& a = fragment_best.second[j];
			bool valid = valid_fragment_alignment(fragment_best.first, a);
			a.alignment->accepted = valid;
		}
	}
}

void accept_all_hits(HitTable& hits)
{
	for (HitTable::iterator i = hits.begin(); 
		 i != hits.end();
		 ++i)
	{
		for (size_t j = 0;
			 j != i->second.size();
			 ++j)
		{
			i->second[j].accepted = true;
		}
	}
}

void add_hits_to_coverage(const HitList& hits, vector<short>& DoC)
{
	int max_hit_pos = -1;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		max_hit_pos = max((int)hits[i].right,max_hit_pos);
	}
	
	if ((int)DoC.size() < max_hit_pos)
		DoC.resize(max_hit_pos);
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const BowtieHit& bh = hits[i];
		
		if (!bh.accepted)
			continue;
		
		if (bh.splice_pos_left != -1) //spliced alignment?
		{
			// split up the coverage contibution for this reads
			size_t j = bh.left;
			while (j <= bh.left + bh.splice_pos_left)
			{
				DoC[j++]++;
			}
			
			j = bh.right - bh.splice_pos_right;
			while (j < bh.right)
			{
				DoC[j++]++;
			}
		}
		else // contiguous alignment
		{
			size_t j = bh.left;
			while (j < bh.right)
			{
				DoC[j++]++;
			}
		}
	}
}