/*
 *  tophat_reports.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/20/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
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
#include "junctions.h"
#include "fragments.h"
#include "wiggles.h"
#include "tokenize.h"
#include "genes.h"
#include "reads.h"

#include "inserts.h"

using namespace std;
using namespace seqan;
using std::set;

void fragment_best_alignments(const HitsForRead& hits_for_read,
							  FragmentAlignmentGrade& best_grade,
                              HitsForRead& best_hits)
{
	const vector<BowtieHit>& hits = hits_for_read.hits;
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		FragmentAlignmentGrade g(hits[i]);
		// Is the new status better than the current best one?
		if (best_grade < g)
		{
			best_hits.hits.clear();
			best_grade = g;
			best_hits.hits.push_back(hits[i]);
		}
		else if (!(g < best_grade)) // is it just as good?
		{
			best_grade.num_alignments++;
			best_hits.hits.push_back(hits[i]);
		}
	}
}

void insert_best_alignments(const HitsForRead& left_hits,
							const HitsForRead& right_hits,
							InsertAlignmentGrade& best_grade,
							HitsForRead& left_best_hits,
							HitsForRead& right_best_hits)
{
	// max mate inner distance (genomic)
	int min_mate_inner_dist = inner_dist_mean - inner_dist_std_dev;
	if (max_mate_inner_dist == -1)
	{
		max_mate_inner_dist = inner_dist_mean + inner_dist_std_dev;
	}
	
	const vector<BowtieHit>& left = left_hits.hits;
	const vector<BowtieHit>& right = right_hits.hits;
	
	for (size_t i = 0; i < left.size(); ++i)
	{
		for (size_t j = 0; j < right.size(); ++j)
		{
		
			const BowtieHit& lh = left[i];
			const BowtieHit& rh = right[j];
			
			InsertAlignmentGrade g(lh, rh, min_mate_inner_dist, max_mate_inner_dist);

			// Is the new status better than the current best one?
			if (best_grade < g)
			{
				left_best_hits.hits.clear();
				right_best_hits.hits.clear();
				
				best_grade = g;
				left_best_hits.hits.push_back(lh);
				right_best_hits.hits.push_back(rh);
			}
			else if (!(g < best_grade))
			{
				left_best_hits.hits.push_back(lh);
				right_best_hits.hits.push_back(rh);
				best_grade.num_alignments++;
			}
		}
	}
}

enum FragmentType {FRAG_UNPAIRED, FRAG_LEFT, FRAG_RIGHT};

bool rewrite_sam_hit(const RefSequenceTable& rt,
                     const BowtieHit& bh,
					 const char* bwt_buf,
					 char* rebuf, 
					 char* read_alt_name,
					 const FragmentAlignmentGrade& grade,
					 FragmentType insert_side,
                     int num_hits,
                     const BowtieHit* next_hit)
{
	// Rewrite this hit, filling in the alt name, mate mapping
	// and setting the pair flag
	vector<string> sam_toks;
	tokenize(bwt_buf, "\t", sam_toks);
	char* slash_pos = strrchr(read_alt_name,'/');
	if (slash_pos)
	{
		*slash_pos = 0;
	}
	
	*rebuf = 0;
	
	for (size_t t = 0; t < sam_toks.size(); ++t)
	{
		switch(t)
		{
			case 0: //QNAME
			{
				sam_toks[t] = read_alt_name;
				break;
			}
			
			case 1:
			{
				// SAM FLAG
				if (insert_side != FRAG_UNPAIRED)
				{
					int flag = 0;
					// mark this guys as a singleton mate
					flag |= 0x0001;
					if (insert_side == FRAG_LEFT)
						flag |= 0x0040;
					else if (insert_side == FRAG_RIGHT)
						flag |= 0x0080;
					flag |= 0x0008;
					
					char flag_buf[64];
					sprintf(flag_buf, "%d", flag);
					sam_toks[t] = flag_buf;
				}
				break;
			}
				
			case 4: //MAPQ
			{
				int mapQ;
				if (grade.num_alignments > 1)
				{
					double err_prob = 1 - (1.0 / grade.num_alignments);
					mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
				}
				else
				{
					mapQ = 255;
				}
				char mapq_buf[64];
				sprintf(mapq_buf, "%d", mapQ);
				sam_toks[t] = mapq_buf;
				break;
			}

			default:
				break;
		}	
		strcat (rebuf, sam_toks[t].c_str());
		if (t != sam_toks.size() - 1)
			strcat(rebuf, "\t");
	}
    
    char nh_buf[2048];
    
    sprintf(nh_buf, 
            "\tNH:i:%d", 
            num_hits);
    
    strcat(rebuf, nh_buf);
    
    if (next_hit)
    {
        const char* nh_ref_name = rt.get_name(next_hit->ref_id());
        assert (nh_ref_name != NULL);
        const char* curr_ref_name = rt.get_name(bh.ref_id());
        assert (curr_ref_name != NULL);
        
        char mate_buf[2048];
        bool same_contig = !strcmp(curr_ref_name, nh_ref_name);
        
        sprintf(mate_buf, 
                "\tCC:Z:%s\tCP:i:%d", 
                same_contig ? "=" : nh_ref_name, 
                next_hit->left() + 1);
        strcat(rebuf, mate_buf);
    }
    strcat(rebuf, "\n");
    
	return true;
}

bool rewrite_sam_hit(const RefSequenceTable& rt,
                     const BowtieHit& bh,
					 const char* bwt_buf,
					 char* rebuf, 
					 char* read_alt_name,
					 const InsertAlignmentGrade& grade,
					 FragmentType insert_side,
					 const BowtieHit* partner,
                     int num_hits,
                     const BowtieHit* next_hit)
{
	// Rewrite this hit, filling in the alt name, mate mapping
	// and setting the pair flag
	vector<string> sam_toks;
	tokenize(bwt_buf, "\t", sam_toks);
	char* slash_pos = strrchr(read_alt_name,'/');
	if (slash_pos)
	{
		*slash_pos = 0;
	}
	
	*rebuf = 0;
	
	for (size_t t = 0; t < sam_toks.size(); ++t)
	{
		switch(t)
		{
			case 0: //QNAME
			{
				sam_toks[t] = read_alt_name;
				break;
			}
			case 1: //SAM FLAG
			{
				// 0x0010 (strand of query) is assumed to be set correctly
				// to begin with
				
				int flag = atoi(sam_toks[1].c_str());
				flag |= 0x0001;
				if (insert_side == FRAG_LEFT)
					flag |= 0x0040;
				else if (insert_side == FRAG_RIGHT)
					flag |= 0x0080;
				

				if (grade.happy() && partner)
					flag |= 0x0002;
				
				if (partner)
				{
					if (partner->antisense_align())
						flag |= 0x0020;
				}
				else
				{
					flag |= 0x0008;
				}
				
				char flag_buf[64];
				sprintf(flag_buf, "%d", flag);
				sam_toks[t] = flag_buf;
				break;
			}

			case 4: //MAPQ
			{
				int mapQ;
				if (grade.num_alignments > 1)
				{
					double err_prob = 1 - (1.0 / grade.num_alignments);
					mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
				}
				else
				{
					mapQ = 255;
				}
				char mapq_buf[64];
				sprintf(mapq_buf, "%d", mapQ);
				sam_toks[t] = mapq_buf;
				break;
			}
			case 6: //MRNM
			{
				if (partner)
				{
					//FIXME: this won't be true forever.  Someday, we will report 
					//alignments of pairs not on the same contig.  
					sam_toks[t] = "=";
				}
				else
				{
					sam_toks[t] = "*";
				}
				break;
			}
			case 7: //MPOS
			{
				if (partner)
				{
					char pos_buf[64];
					int partner_pos = partner->left() + 1;  // SAM is 1-indexed
					
					sprintf(pos_buf, "%d", partner_pos);
					sam_toks[t] = pos_buf;
					break;
				}
				else
				{
					sam_toks[t] = "0";
				}
			}
			default:
				break;
		}	
		strcat (rebuf, sam_toks[t].c_str());
		if (t != sam_toks.size() - 1)
			strcat(rebuf, "\t");
	}
    
    char nh_buf[2048];
    
    sprintf(nh_buf, 
            "\tNH:i:%d", 
            num_hits);
    
    strcat(rebuf, nh_buf);
    
    
    if (next_hit)
    {
        const char* nh_ref_name = rt.get_name(next_hit->ref_id());
        assert (nh_ref_name != NULL);
        const char* curr_ref_name = rt.get_name(bh.ref_id());
        assert (curr_ref_name != NULL);
        
        char mate_buf[2048];
        bool same_contig = !strcmp(curr_ref_name, nh_ref_name);
        
        assert (num_hits > 1);
        
        sprintf(mate_buf, 
                "\tCC:Z:%s\tCP:i:%d",
                same_contig ? "=" : nh_ref_name, 
                next_hit->left() + 1);
        
        strcat(rebuf, mate_buf);
    }
    strcat(rebuf, "\n");
    
	return true;
}

struct lex_hit_sort
{
    lex_hit_sort(const RefSequenceTable& rt) : _rt(rt) {}
    
    bool operator()(const BowtieHit& lhs, const BowtieHit& rhs) const
    {
        return (strcmp(_rt.get_name(lhs.ref_id()), _rt.get_name(rhs.ref_id())) < 0);
    }
    
    const RefSequenceTable& _rt;
};

void print_sam_for_hits(const RefSequenceTable& rt,
                        const HitsForRead& hits,
						const FragmentAlignmentGrade& grade,
						FragmentType frag_type,
						FILE* reads_file,
						FILE* fout)
{
	static const int buf_size = 2048;
	char read_name[buf_size];
	char read_seq[buf_size];
	char read_alt_name[buf_size];
	char read_quals[buf_size];
	
	char rebuf[buf_size];
    
	HitsForRead sorted_hits = hits;
    lex_hit_sort s(rt);
    sort(sorted_hits.hits.begin(), sorted_hits.hits.end(), s);
    
	bool got_read = get_read_from_stream(hits.insert_id, 
										 reads_file,
										 reads_format,
										 false,
										 read_name, 
										 read_seq,
										 read_alt_name,
										 read_quals);
	
	assert (got_read);
	
	for (size_t i = 0; i < sorted_hits.hits.size(); ++i)
	{
		const BowtieHit& bh = sorted_hits.hits[i];
		if (rewrite_sam_hit(rt, 
                            bh, 
                            bh.hitfile_rec().c_str(), 
                            rebuf, 
                            read_alt_name, 
                            grade, 
                            frag_type,
                            sorted_hits.hits.size(),
                            (i < sorted_hits.hits.size() - 1) ? &(sorted_hits.hits[i+1]) : NULL))
        {
			fprintf(fout, "%s", rebuf);
        }
	}
}

void print_sam_for_hits(const RefSequenceTable& rt,
                        const HitsForRead& left_hits,
						const HitsForRead& right_hits,
						const InsertAlignmentGrade& grade,
						FILE* left_reads_file,
						FILE* right_reads_file,
						FILE* fout)
{
	assert (left_hits.insert_id == right_hits.insert_id);
	
	static const int buf_size = 2048;
	char left_read_name[buf_size];
	char left_read_seq[buf_size];
	char left_read_alt_name[buf_size];
	char left_read_quals[buf_size];
	char left_rebuf[buf_size];
	
	char right_read_name[buf_size];
	char right_read_seq[buf_size];
	char right_read_alt_name[buf_size];
	char right_read_quals[buf_size];
	char right_rebuf[buf_size];
	
    lex_hit_sort s(rt);
    
    HitsForRead left_sorted_hits = left_hits;
    sort(left_sorted_hits.hits.begin(), left_sorted_hits.hits.end(), s);
    
    HitsForRead right_sorted_hits = right_hits;
    sort(right_sorted_hits.hits.begin(), right_sorted_hits.hits.end(), s);
    
	bool got_left_read = get_read_from_stream(left_sorted_hits.insert_id, 
											  left_reads_file,
											  reads_format,
											  false,
											  left_read_name, 
											  left_read_seq,
											  left_read_alt_name,
											  left_read_quals);
	
	bool got_right_read = get_read_from_stream(right_sorted_hits.insert_id, 
											   right_reads_file,
											   reads_format,
											   false,
											   right_read_name, 
											   right_read_seq,
											   right_read_alt_name,
											   right_read_quals);
	
	assert (left_sorted_hits.hits.size() == right_sorted_hits.hits.size() ||
			(left_sorted_hits.hits.empty() || right_sorted_hits.hits.empty()));
	
	if (left_sorted_hits.hits.size() == right_sorted_hits.hits.size())
	{
		assert (got_left_read && got_right_read);
		for (size_t i = 0; i < right_sorted_hits.hits.size(); ++i)
		{
			const BowtieHit& right_bh = right_sorted_hits.hits[i];
			if (rewrite_sam_hit(rt,
                                right_bh, 
								right_bh.hitfile_rec().c_str(), 
								right_rebuf, 
								right_read_alt_name, 
								grade, 
								FRAG_RIGHT, 
								&left_sorted_hits.hits[i],
                                right_sorted_hits.hits.size(),
                                (i < right_sorted_hits.hits.size() - 1) ? &(right_sorted_hits.hits[i+1]) : NULL))
            {
				fprintf(fout, "%s", right_rebuf);
			}
			
			const BowtieHit& left_bh = left_sorted_hits.hits[i];
			if (rewrite_sam_hit(rt,
                                left_bh, 
								left_bh.hitfile_rec().c_str(), 
								left_rebuf, 
								left_read_alt_name, 
								grade, 
								FRAG_LEFT, 
								&right_sorted_hits.hits[i],
                                left_sorted_hits.hits.size(),
                                (i < left_sorted_hits.hits.size() - 1) ? &(left_sorted_hits.hits[i+1]) : NULL))
			{
				fprintf(fout, "%s", left_rebuf);
			}
		}
	}
	else if (left_sorted_hits.hits.empty())
	{
		for (size_t i = 0; i < right_sorted_hits.hits.size(); ++i)
		{
			const BowtieHit& bh = right_sorted_hits.hits[i];
			if (rewrite_sam_hit(rt,
                                bh, 
								bh.hitfile_rec().c_str(), 
								right_rebuf, 
								right_read_alt_name, 
								grade, 
								FRAG_RIGHT, 
								NULL,
                                right_sorted_hits.hits.size(),
                                (i < right_sorted_hits.hits.size() - 1) ? &(right_sorted_hits.hits[i+1]) : NULL))
				fprintf(fout, "%s", right_rebuf);
		}
	}
	else if (right_sorted_hits.hits.empty())
	{
		for (size_t i = 0; i < left_sorted_hits.hits.size(); ++i)
		{
			const BowtieHit& bh = left_sorted_hits.hits[i];
			if (rewrite_sam_hit(rt,
                                bh, 
								bh.hitfile_rec().c_str(), 
								left_rebuf, 
								left_read_alt_name, 
								grade, 
								FRAG_LEFT, 
								NULL,
                                left_sorted_hits.hits.size(),
                                (i < left_sorted_hits.hits.size() - 1) ? &(left_sorted_hits.hits[i+1]) : NULL))
				fprintf(fout, "%s", left_rebuf);
		}
	}
	else
	{
		assert (false);
	}
}

// Extracts junctions from all the SAM hits (based on REF_SKIPs) in the hit file
// resets the stream when finished.
void get_junctions_from_hits(HitStream& hit_stream, 
							 ReadTable& it, 
							 JunctionSet& junctions)
{
	HitsForRead curr_hit_group;
	hit_stream.next_read_hits(curr_hit_group);

	uint32_t curr_obs_order = it.observation_order(curr_hit_group.insert_id);
	
	while(curr_obs_order != 0xFFFFFFFF)
	{
		for (size_t i = 0; i < curr_hit_group.hits.size(); ++i)
		{
			BowtieHit& bh = curr_hit_group.hits[i];
			if (!bh.contiguous())
			{
				junctions_from_alignment(bh, junctions);
			}
			hit_stream.next_read_hits(curr_hit_group);
			curr_obs_order = it.observation_order(curr_hit_group.insert_id);
		}
	}
	
	hit_stream.reset();
}

void update_junctions(const HitsForRead& hits,
					  JunctionSet& junctions)
{
	for (size_t i = 0; i < hits.hits.size(); ++i)
	{
		const BowtieHit& bh = hits.hits[i];
		junctions_from_alignment(bh, junctions);
	}
}

// Extracts junctions from all the SAM hits (based on REF_SKIPs) in the hit file
// resets the stream when finished.
void exclude_hits_on_filtered_junctions(const JunctionSet& junctions,
										HitsForRead& hits)
{
	HitsForRead remaining;
	remaining.insert_id = hits.insert_id;

	for (size_t i = 0; i < hits.hits.size(); ++i)
	{
		BowtieHit& bh = hits.hits[i];
		bool filter_hit = false;
		if (!bh.contiguous())
		{
			JunctionSet bh_juncs;
			junctions_from_alignment(bh, bh_juncs);
			for (JunctionSet::iterator itr = bh_juncs.begin();
				 itr != bh_juncs.end();
				 itr++)
			{
				const Junction& j = itr->first;
				JunctionSet::const_iterator target = junctions.find(j);
				if (target == junctions.end() || !target->second.accepted)
				{
					filter_hit = true;
					break;
				}
			}
		}
		if (!filter_hit)
		{
			remaining.hits.push_back(bh);
		}
		else
		{
			//fprintf(stderr, "Filtering hit\n");
		}
	}
	hits = remaining;
}

void get_junctions_from_best_hits(HitStream& left_hs,
								  HitStream& right_hs,
								  ReadTable& it, 
								  JunctionSet& junctions)
{
	HitsForRead curr_left_hit_group;
	HitsForRead curr_right_hit_group;
	
	left_hs.next_read_hits(curr_left_hit_group);
	right_hs.next_read_hits(curr_right_hit_group);
	
	uint32_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	uint32_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	
	// While we still have unreported hits...
	while(curr_left_obs_order != 0xFFFFFFFF || 
		  curr_right_obs_order != 0xFFFFFFFF)
	{		
		// Chew up left singletons
		while (curr_left_obs_order < curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_left_obs_order;
			FragmentAlignmentGrade grade;
			
			// Process hits for left singleton, select best alignments
			fragment_best_alignments(curr_left_hit_group, grade, best_hits);
			
			update_junctions(best_hits, junctions);
			
			// Get next hit group
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		}
		
		// Chew up right singletons
		while (curr_left_obs_order > curr_right_obs_order &&
			   curr_right_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_right_obs_order;
			FragmentAlignmentGrade grade;
			
			// Process hit for right singleton, select best alignments
			fragment_best_alignments(curr_right_hit_group,grade, best_hits);
			
			update_junctions(best_hits, junctions);
			
			// Get next hit group
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
		
		// Since we have both left hits and right hits for this insert,
		// Find the best pairing and print both
		while (curr_left_obs_order == curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{			 
			if (curr_left_hit_group.hits.empty())
			{
				HitsForRead right_best_hits;
				right_best_hits.insert_id = curr_right_obs_order;
				
				FragmentAlignmentGrade grade;
				fragment_best_alignments(curr_right_hit_group, grade, right_best_hits);
				
				update_junctions(right_best_hits, junctions);
			}
			else if (curr_right_hit_group.hits.empty())
			{
				HitsForRead left_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				
				FragmentAlignmentGrade grade;
				// Process hits for left singleton, select best alignments
				fragment_best_alignments(curr_left_hit_group, grade, left_best_hits);
				
				update_junctions(left_best_hits, junctions);
			}
			else
			{		
				HitsForRead left_best_hits;
				HitsForRead right_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				right_best_hits.insert_id = curr_right_obs_order;
				
				InsertAlignmentGrade grade;
				insert_best_alignments(curr_left_hit_group, 
									   curr_right_hit_group, 
									   grade,
									   left_best_hits,
									   right_best_hits);
				
				update_junctions(left_best_hits, junctions);
				update_junctions(right_best_hits, junctions);
			}
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
	}
	
	
	left_hs.reset();
	right_hs.reset();
}


void driver(FILE* left_map,
			FILE* left_reads,
            FILE* right_map,
			FILE* right_reads,
            FILE* junctions_out,
			FILE* accepted_hits_out)
{	
    ReadTable it;
    RefSequenceTable rt(true);
	
    SAMHitFactory hit_factory(it,rt);
	
	HitStream left_hs(left_map, &hit_factory, false, true, true);
	HitStream right_hs(right_map, &hit_factory, false, true, true);
	
    JunctionSet junctions;
		
	get_junctions_from_best_hits(left_hs, right_hs, it, junctions);
	
	size_t num_unfiltered_juncs = junctions.size();
	
	fprintf(stderr, "Loaded %lu junctions\n", num_unfiltered_juncs); 
	
	HitsForRead curr_left_hit_group;
	HitsForRead curr_right_hit_group;
	
	left_hs.next_read_hits(curr_left_hit_group);
	right_hs.next_read_hits(curr_right_hit_group);

	uint32_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	uint32_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	
	// Read hits, extract junctions, and toss the ones that arent strongly enough supported.
	
	filter_junctions(junctions);
	//size_t num_juncs_after_filter = junctions.size();
	//fprintf(stderr, "Filtered %lu junctions\n", num_unfiltered_juncs - num_juncs_after_filter);
	
	size_t small_overhangs = 0;
	for (JunctionSet::iterator i = junctions.begin(); i != junctions.end(); ++i)
	{
		if (i->second.accepted && 
			(i->second.left_extent < min_anchor_len || i->second.right_extent < min_anchor_len))
		{
			small_overhangs++;
		}
	}
	
	if (small_overhangs >0)
		fprintf(stderr, "Warning: %lu small overhang junctions!\n", small_overhangs);
	
	JunctionSet final_junctions; // the junctions formed from best hits
	
	fprintf (stderr, "Reporting final accepted alignments...");
	
	// While we still have unreported hits...
	while(curr_left_obs_order != 0xFFFFFFFF || 
		  curr_right_obs_order != 0xFFFFFFFF)
	{		
		// Chew up left singletons
		while (curr_left_obs_order < curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_left_obs_order;
			FragmentAlignmentGrade grade;
			
			exclude_hits_on_filtered_junctions(junctions, curr_left_hit_group);
			
			// Process hits for left singleton, select best alignments
			fragment_best_alignments(curr_left_hit_group, grade, best_hits);
			
			update_junctions(best_hits, final_junctions);
			
			print_sam_for_hits(rt,
                               best_hits, 
							   grade,
							   right_map ? FRAG_LEFT : FRAG_UNPAIRED,
							   left_reads, 
							   accepted_hits_out);
			
			// Get next hit group
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		}
		
		// Chew up right singletons
		while (curr_left_obs_order > curr_right_obs_order &&
			   curr_right_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_right_obs_order;
			FragmentAlignmentGrade grade;
			
			exclude_hits_on_filtered_junctions(junctions, curr_right_hit_group);
			
			// Process hit for right singleton, select best alignments
			fragment_best_alignments(curr_right_hit_group,grade, best_hits);
			
			update_junctions(best_hits, final_junctions);
			
			print_sam_for_hits(rt,
                               best_hits, 
							   grade, 
							   FRAG_RIGHT,
							   right_reads, 
							   accepted_hits_out);
			
			// Get next hit group
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
		
		// Since we have both left hits and right hits for this insert,
		// Find the best pairing and print both
		while (curr_left_obs_order == curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{			 
			exclude_hits_on_filtered_junctions(junctions, curr_left_hit_group);
			exclude_hits_on_filtered_junctions(junctions, curr_right_hit_group);
			
			if (curr_left_hit_group.hits.empty())
			{
				HitsForRead right_best_hits;
				right_best_hits.insert_id = curr_right_obs_order;
				
				FragmentAlignmentGrade grade;
				fragment_best_alignments(curr_right_hit_group, grade, right_best_hits);
				
				update_junctions(right_best_hits, final_junctions);
				
				print_sam_for_hits(rt,
                                   right_best_hits, 
								   grade, 
								   FRAG_RIGHT,
								   right_reads, 
								   accepted_hits_out);	
			}
			else if (curr_right_hit_group.hits.empty())
			{
				HitsForRead left_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				
				FragmentAlignmentGrade grade;
				// Process hits for left singleton, select best alignments
				fragment_best_alignments(curr_left_hit_group, grade, left_best_hits);
				
				update_junctions(left_best_hits, final_junctions);
				
				print_sam_for_hits(rt,
                                   left_best_hits, 
								   grade,
								   FRAG_LEFT,
								   left_reads, 
								   accepted_hits_out);
			}
			else
			{		
				HitsForRead left_best_hits;
				HitsForRead right_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				right_best_hits.insert_id = curr_right_obs_order;
				
				InsertAlignmentGrade grade;
				insert_best_alignments(curr_left_hit_group, 
									   curr_right_hit_group, 
									   grade,
									   left_best_hits,
									   right_best_hits);
				
				update_junctions(left_best_hits, final_junctions);
				update_junctions(right_best_hits, final_junctions);
				
				print_sam_for_hits(rt,
                                   left_best_hits,
								   right_best_hits,
								   grade,
								   left_reads, 
								   right_reads, 
								   accepted_hits_out);
			}
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
	}
	
	fprintf (stderr, "done\n");
	
	small_overhangs = 0;
	for (JunctionSet::iterator i = final_junctions.begin(); i != final_junctions.end();)
	{
		const JunctionStats& stats = i->second;
		if (i->second.supporting_hits == 0 || 
			i->second.left_extent < 8 || 
			i->second.right_extent < 8)
		{
			final_junctions.erase(i++);
		}
		else
		{
			++i;
		}
	}
	
	
	
//	if (small_overhangs > 0)
//		fprintf(stderr, "Warning: %lu small overhang junctions!\n", small_overhangs);
	
	fprintf (stderr, "Printing junction BED track...");
	print_junctions(junctions_out, final_junctions, rt);
	fprintf (stderr, "done\n");
	
    fprintf(stderr, "Found %lu junctions from happy spliced reads\n", final_junctions.size());
}

void print_usage()
{
    fprintf(stderr, "Usage:   tophat_reports <junctions.bed> <accepted_hits.sam> <left_map1,...,left_mapN> <left_reads.fq>  [right_map1,...,right_mapN] [right_reads.fq]\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

int main(int argc, char** argv)
{
	fprintf(stderr, "tophat_reports v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "---------------------------------------\n");
	
	reads_format = FASTQ;
	
    int parse_ret = parse_options(argc, argv, print_usage);
    if (parse_ret)
        return parse_ret;

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
	
	string accepted_hits_file_name = argv[optind++];
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }

    string left_map_filename = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
	FILE* left_map = fopen(left_map_filename.c_str(), "r");
	if (!left_map)
	{
		fprintf(stderr, "Error: cannot open map file %s for reading\n",
				left_map_filename.c_str());
		exit(1);
	}
	
    string left_reads_filename = argv[optind++];
	
    string* right_map_filename = NULL;
	FILE* right_map = NULL;
	string* right_reads_filename = NULL;
	FILE* right_reads_file = NULL;
	
    if (optind < argc)
	{
        right_map_filename = new string(argv[optind++]);
		
		if(optind >= argc)
		{
			print_usage();
			return 1;
		}
		
		right_map = fopen(right_map_filename->c_str(), "r");
		if (!right_map)
		{
			fprintf(stderr, "Error: cannot open map file %s for reading\n",
					right_map_filename->c_str());
			exit(1);
		}
		
		right_reads_filename = new string(argv[optind++]);
		right_reads_file = fopen(right_reads_filename->c_str(), "r");
		if (!right_reads_file)
		{
			fprintf(stderr, "Error: cannot open reads file %s for reading\n",
					right_reads_filename->c_str());
			exit(1);
		}
	}

    FILE* junctions_file = fopen((output_dir + "/" + junctions_file_name).c_str(), "w");
    if (junctions_file == NULL)
    {
        fprintf(stderr, "Error: cannot open BED file %s for writing\n",
                junctions_file_name.c_str());
        exit(1);
    }

    FILE* accepted_hits_file = fopen((output_dir + "/" + accepted_hits_file_name).c_str(), "w");
    if (accepted_hits_file == NULL)
    {
        fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
                accepted_hits_file_name.c_str());
        exit(1);
    }
	
	FILE* left_reads_file = fopen(left_reads_filename.c_str(), "r");
    if (!left_reads_file)
    {
        fprintf(stderr, "Error: cannot open reads file %s for reading\n",
                left_reads_filename.c_str());
        exit(1);
    }

    driver(left_map,
		   left_reads_file,
           right_map,
		   right_reads_file,
           junctions_file,
		   accepted_hits_file);

    return 0;
}


