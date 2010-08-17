#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  long_spanning_reads.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/5/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <cassert>
#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <cstring>
#include <bitset>
//#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "segments.h"
#include "reads.h"

#include "junctions.h"

//#include "fragments.h"
//#include "wiggles.h"
//#include "tokenize.h"
//#include "genes.h"
//#include "FSA/gff.h"
//#include "FSA/sequence.h"

using namespace seqan;
using namespace std;

void print_usage()
{
    fprintf(stderr, "Usage:   long_spanning_reads <reads.fa/.fq> <possible_juncs1,...,possible_juncsN> <seg1.bwtout,...,segN.bwtout> [spliced_seg1.bwtout,...,spliced_segN.bwtout]\n");
}

//const char *short_options = "a:i:I:fqm:";

// This is the maximum number of bowtie mismatches allower per segment hit
//static const int num_bowtie_mismatches = 2;

//int max_splice_mismatches = 0;

//static struct option long_options[] = {
//{"min-anchor",       required_argument,       0,            'a'},
//{"min-intron",       required_argument,       0,            'i'},
//{"splice-mismatches",       required_argument,       0,            'm'},
//{0, 0, 0, 0} // terminator
//};
//
//int parse_options(int argc, char** argv)
//{
//    int option_index = 0;
//    int next_option;
//    do {
//        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
//        switch (next_option) {
//			case 'a':
//				min_anchor_len = (uint32_t)parseInt(3, "-a/--min-anchor arg must be at least 3", print_usage);
//				break;
//			case 'i':
//				min_intron_length = (uint32_t)parseInt(1, "-a/--min-intron arg must be at least 1", print_usage);
//				break;
//			case 'I':
//				max_intron_length = parseInt(1,"-I arg must be at least 1", print_usage);
//				break;
//			case 'm':
//				max_splice_mismatches = parseInt(0,"-I arg must be at least 1", print_usage);
//				break;
//			case 'f':
//				reads_format = FASTA;
//				break;
//			case 'q':
//				reads_format = FASTQ;
//				break;    
//			case -1:     /* Done with options. */
//                break;
//            default:
//                print_usage();
//                return 1;
//        }
//    } while(next_option != -1);
//    
//    return 0;
//}

bool key_lt(const pair<uint32_t, HitsForRead>& lhs,
			const pair<uint32_t, HitsForRead>& rhs)
{
	return lhs.first < rhs.first;
}

void get_seqs(istream& ref_stream,
			  RefSequenceTable& rt,
			  bool keep_seqs = true,
			  bool strip_slash = false)
{    
    while(ref_stream.good() &&
          !ref_stream.eof())
    {
		RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence();
        string name;
        readMeta(ref_stream, name, Fasta());
		string::size_type space_pos = name.find_first_of(" \t\r");
		if (space_pos != string::npos)
		{
			name.resize(space_pos);
		}
		seqan::read(ref_stream, *ref_str, Fasta());
		
        rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
    }	
}


void look_left_for_hit_group(ReadTable& unmapped_reads,
							 vector<HitStream>& contig_hits,
							 size_t curr_file,
							 vector<HitStream>& spliced_hits,
							 const HitsForRead& targets,
							 vector<HitsForRead>& seg_hits_for_read)
{
	int left_file = curr_file - 1;
	
	HitStream& left = contig_hits[left_file];
	//HitStream& curr = seg_files[curr_file];
	
	uint64_t curr_next_group_id = targets.insert_id;

	int curr_order = unmapped_reads.observation_order(curr_next_group_id);
	
	assert (curr_order != -1);
	while(true)
	{
		HitsForRead hit_group;
		uint64_t left_next_group_id = left.next_group_id();
		int left_order = unmapped_reads.observation_order(left_next_group_id);
		
		// If we would have seen the hits by now, bail out.
		if (curr_order < left_order || left_order == -1)
		{
			break;
		}
		if (left.next_read_hits(hit_group))
		{
			if (hit_group.insert_id == targets.insert_id)
			{
				// Some of the targets may be missing, we need to 
				// process them individually
				seg_hits_for_read[left_file] = hit_group;
				break;
			}
		}
	}
	
	HitsForRead& curr_seg_hits = seg_hits_for_read[left_file];
	
	if (left_file < (int)spliced_hits.size() && left_file >= 0)
	{
//		vector<pair<uint32_t, HitsForRead> >::iterator splice_itr;
//		pair<uint32_t, HitsForRead> dummy = make_pair(curr_order, HitsForRead());
//		splice_itr = lower_bound(spliced_hits[left_file].begin(),
//								 spliced_hits[left_file].end(),
//								 dummy,
//								 key_lt);
//		
//		if (splice_itr != spliced_hits[left_file].end() &&
//			splice_itr->first == (uint32_t)curr_order)
//		{
//			HitsForRead& spliced_group = splice_itr->second;
//			curr_seg_hits.insert_id = splice_itr->second.insert_id;
//			
//			curr_seg_hits.hits.insert(curr_seg_hits.hits.end(),
//									  spliced_group.hits.begin(),
//									  spliced_group.hits.end());
//		}
		
		// Scan forward in the spliced hits file for this hit group
		HitsForRead spliced_group;
		HitsForRead curr_spliced_group;
		while (spliced_hits[left_file].next_group_id() > 0 && 
			   spliced_hits[left_file].next_group_id() <= (uint32_t)curr_order)
		{
			spliced_hits[left_file].next_read_hits(curr_spliced_group);
			
			if (curr_spliced_group.insert_id == (uint32_t)curr_order)
			{
				spliced_group = curr_spliced_group;
				break;
			}
		}
		
		if (!spliced_group.hits.empty())
		{
			curr_seg_hits.insert_id = spliced_group.insert_id;
			curr_seg_hits.hits.insert(curr_seg_hits.hits.end(),
									  spliced_group.hits.begin(),
									  spliced_group.hits.end());
		}
	}
	
	if (curr_seg_hits.hits.empty())
		return;
	else if (left_file > 0)
	{
		look_left_for_hit_group(unmapped_reads, 
								contig_hits, 
								curr_file - 1,
								spliced_hits,
								curr_seg_hits,
								seg_hits_for_read);
	}
}

BowtieHit merge_sense_chain(std::set<Junction>& possible_juncs,
							list<BowtieHit>& hit_chain)
{
	
	uint32_t reference_id = hit_chain.front().ref_id();
	uint64_t insert_id = hit_chain.front().insert_id();
	
	int	left = hit_chain.front().left();	
	
	list<BowtieHit>::iterator prev_hit = hit_chain.begin();
	list<BowtieHit>::iterator curr_hit = ++(hit_chain.begin());

	string seq;
	string qual;
	int old_read_length = 0;
	for(list<BowtieHit>::iterator i = hit_chain.begin();
		i != hit_chain.end(); 
		++i)
	{
	  seq += i->seq();
	  qual += i->qual();
	  old_read_length += i->read_len();
	}
	
	
	while(curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
	{
		int gap  = curr_hit->left() - prev_hit->right();
		if (gap > 0 && 
			(gap < min_report_intron_length || gap > max_report_intron_length))
		{
			return BowtieHit();
		}
		
		++prev_hit;
		++curr_hit;
	}
	
	prev_hit = hit_chain.begin();
	curr_hit = ++(hit_chain.begin());
	
	while(curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
	{
//		BowtieHit& lh = *left_hit;
//		BowtieHit& rh = *right_hit;
		//printf ("rh = [%d,%d], lh = [%d,%d]\n", rh.left(), rh.right(), lh.left(), lh.right());
		if (prev_hit->right() != curr_hit->left())
		{
			if (curr_hit->left() - prev_hit->right() < (int)min_report_intron_length)
				return BowtieHit();
			
			std::set<Junction>::iterator lb,ub;
			
			int left_boundary = prev_hit->right() - 4;
			int right_boundary = curr_hit->left() + 4;
			
			lb = possible_juncs.upper_bound(Junction(reference_id, left_boundary, right_boundary - 8, true));
			ub = possible_juncs.lower_bound(Junction(reference_id, left_boundary + 8, right_boundary, false));
			
			
			int new_left = -1;
			bool antisense_closure = false;
			vector<CigarOp> new_cigar;
			bool found_closure = false;
			
			while (lb != ub && lb != possible_juncs.end())
			{

//				int lbl = lb->left;
//				int lbr = lb->right;
//				int rhr = right_hit->right();
//				int rhl = right_hit->left();
//				int lhl = left_hit->left();
//				int lhr = left_hit->right();
//				fprintf(stderr, "[%d]\t%d\t%d\n", (int)insert_id, lb->left, lb->right);
				int dist_to_left = lb->left - prev_hit->right() + 1;
				int dist_to_right = lb->right - curr_hit->left();
				if ( abs(dist_to_left) <= 4 && abs(dist_to_right) <= 4 && dist_to_left == dist_to_right)
				{
					antisense_closure = lb->antisense;
					if (found_closure)
					{
						fprintf(stderr, "Warning: multiple closures found for read # %d\n", (int)insert_id);
						return BowtieHit();
					}
					
					found_closure = true;
					
					new_left = prev_hit->left();
					new_cigar = prev_hit->cigar();
					
//					int new_right_back_len = new_cigar.back().length;
//					new_right_back_len += dist_to_left;
//					
//					vector<CigarOp> new_left_cig = left_hit->cigar();
//					int new_left_front_len = new_left_cig.front().length;
//					new_left_front_len -= dist_to_right;
					
					int new_left_back_len = new_cigar.back().length;
					new_left_back_len += dist_to_left;
					
					vector<CigarOp> new_right_cig = curr_hit->cigar();
					int new_right_front_len = new_right_cig.front().length;
					new_right_front_len -= dist_to_right;
					
					if (new_right_front_len > 0 && new_left_back_len > 0)
					{
						//new_cigar.back().length += dist_to_left; //fix the left endpoint
						new_cigar.push_back(CigarOp(REF_SKIP, lb->right - lb->left - 1));
						
						//vector<CigarOp> new_right_cig = right_hit->cigar();
						new_right_cig.front().length = new_right_front_len;
						new_cigar.back().length = new_left_back_len;
						
						for (size_t c = 0; c < new_right_cig.size(); ++c)
							new_cigar.push_back(new_right_cig[c]);
					}
				}
				++lb;
			}
			
			if (found_closure)
			{
				// If we got here, it means there's exactly one closure.
				BowtieHit merged_hit(reference_id,
									 insert_id,
									 new_left,
									 new_cigar,
									 false,
									 antisense_closure,
									 prev_hit->edit_dist() + curr_hit->edit_dist(),
									 prev_hit->splice_mms() + curr_hit->splice_mms());
				
				prev_hit = hit_chain.erase(prev_hit,++curr_hit);
				hit_chain.insert(prev_hit, merged_hit);
				curr_hit = prev_hit;
				curr_hit++;
				continue;
			}
			else
			{
				//fprintf (stderr, "Couldn't get sense closure for read #%d\n", (int)insert_id);
				// If we get here, we couldn't even find a closure for the hits
				return BowtieHit();
			}
		}
		
		++prev_hit;
		++curr_hit;
	}
	
	
	bool saw_antisense_splice = false;
	bool saw_sense_splice = false;
	vector<CigarOp> long_cigar;
	int num_mismatches = 0;
	int num_splice_mms = 0;
	for (list<BowtieHit>::iterator s = hit_chain.begin(); s != hit_chain.end(); ++s)
	{
		num_mismatches += s->edit_dist();
		num_splice_mms += s->splice_mms();
		if (!s->contiguous())
		{
			if (s->antisense_splice())
			{
				if (saw_sense_splice)
					return BowtieHit();
				saw_antisense_splice = true;
			}
			else
			{
				if (saw_antisense_splice)
					return BowtieHit();
				saw_sense_splice = true;
			}
		}
		const vector<CigarOp>& cigar = s->cigar();
		if (long_cigar.empty())
		{
			long_cigar = cigar;
		}
		else
		{
			CigarOp& last = long_cigar.back();
			last.length += cigar[0].length;
			for (size_t b = 1; b < cigar.size(); ++b)
			{
				long_cigar.push_back(cigar[b]);
			}
		}
	}
	
	BowtieHit new_hit(reference_id,
					  insert_id, 
					  left, 
					  long_cigar, 
					  false,
					  saw_antisense_splice,
					  num_mismatches,
					  num_splice_mms);

	new_hit.seq(seq);
	new_hit.qual(qual);
	
	int new_read_len = new_hit.read_len();
	if (new_read_len != old_read_length)
	{
		//fprintf(stderr, "Warning: malformed closure\n");
		return BowtieHit();
	}
	return new_hit;
	
}

BowtieHit merge_antisense_chain(std::set<Junction>& possible_juncs,
								list<BowtieHit>& hit_chain)
{
	assert(hit_chain.size() > 1);
	uint32_t reference_id = hit_chain.front().ref_id();
	uint64_t insert_id = hit_chain.front().insert_id();
	
	int	left = hit_chain.back().left();	
	
	list<BowtieHit>::iterator prev_hit = hit_chain.begin();
	list<BowtieHit>::iterator curr_hit = ++(hit_chain.begin());

	string seq;
	string qual;
	int old_read_length = 0;
	for(list<BowtieHit>::iterator i = hit_chain.begin();
		i != hit_chain.end(); 
		++i)
	{
	    seq = i->seq() + seq;
	    qual = i->qual() + qual;
	    old_read_length += i->read_len();
	}
	
	while(curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
	{
		int gap  = prev_hit->left() - curr_hit->right();
		if (gap > 0 && 
			(gap < min_report_intron_length || gap > max_report_intron_length))
		{
			return BowtieHit();
		}
		
		++prev_hit;
		++curr_hit;
	}
	
	prev_hit = hit_chain.begin();
	curr_hit = ++(hit_chain.begin());
	
	while(curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
	{
		//BowtieHit& lh = *prev_hit;
		//BowtieHit& rh = *curr_hit;
		//printf ("rh = [%d,%d], lh = [%d,%d]\n", rh.left(), rh.right(), lh.left(), lh.right());
//		int r = right_hit->right();
//		int l = left_hit->left();
		if (curr_hit->right() < prev_hit->left())
		{
			if (prev_hit->left() - curr_hit->right() < min_report_intron_length)
				return BowtieHit();
			
			std::set<Junction>::iterator lb,ub;
			//printf("** %d\t%d\n", right_hit->right() - 2, left_hit->left() + 2);
			
			int left_boundary = curr_hit->right() - 4;
			int right_boundary = prev_hit->left() + 4;
			// 
			lb = possible_juncs.upper_bound(Junction(reference_id, left_boundary, right_boundary - 8, true));
			ub = possible_juncs.lower_bound(Junction(reference_id, left_boundary + 8, right_boundary, false));
			
			bool found_closure = false;
			bool antisense_closure = false;
			int new_left = -1;
			vector<CigarOp> new_cigar;
			
			while (lb != ub && lb != possible_juncs.end())
			{
				
				//<curr_hit> ---------- <prev_hit>
				//        left---------right
//				int lbl = lb->left;
//				int lbr = lb->right;
//				int curr_right = curr_hit->right();
//				int curr_left = curr_hit->left();
//				int prev_left = prev_hit->left();
//				int prev_right = prev_hit->right();
//				fprintf(stderr, "[%d]\t%d\t%d, %d\t%d\n", (int)insert_id, lb->left, lb->right, prev_hit->left(), curr_hit->right());
				int dist_to_left = (int)lb->left - (int)curr_hit->right() + 1;
				int dist_to_right = (int)(lb->right) - (int)(prev_hit->left());
				
				if (abs(dist_to_left) <= 4 && abs(dist_to_right) <= 4 && dist_to_left == dist_to_right)
				{
					//printf("%d\t%d\n", lb->left, lb->right);
					
					new_left = curr_hit->left();
					new_cigar = curr_hit->cigar();
					
					int new_curr_len = new_cigar.back().length;
					new_curr_len += dist_to_left;
					
					vector<CigarOp> new_left_cig = prev_hit->cigar();
					int new_prev_len = new_left_cig.front().length;
					new_prev_len -= dist_to_right;
					
					if (new_prev_len > 0 && new_curr_len > 0)
					{
						new_cigar.back().length = new_curr_len; //fix the left endpoint
						new_cigar.push_back(CigarOp(REF_SKIP, lb->right - lb->left - 1));
						
						new_left_cig.front().length = new_prev_len;
						
						for (size_t c = 0; c < new_left_cig.size(); ++c)
							new_cigar.push_back(new_left_cig[c]);
						
						antisense_closure = lb->antisense;

						if (found_closure)
						{
							fprintf(stderr, "Warning: multiple closures found for read # %d\n", (int)insert_id);
							return BowtieHit();
						}
						
						found_closure = true;
					}
				}
				++lb;
			}
			
			if (found_closure)
			{				
				BowtieHit merged_hit(reference_id,
									 insert_id,
									 new_left,
									 new_cigar,
									 true,
									 antisense_closure,
									 prev_hit->edit_dist() + curr_hit->edit_dist(),
									 prev_hit->splice_mms() + curr_hit->splice_mms());
				
				prev_hit = hit_chain.erase(prev_hit,++curr_hit);
				hit_chain.insert(prev_hit, merged_hit);
				curr_hit = prev_hit;
				curr_hit++;
				continue;
			}
			else
			{
				//fprintf (stderr, "Couldn't get antisense closure for read #%d\n", (int)insert_id);
				// If we get here, we couldn't even find a closure for the hits
				return BowtieHit();
			}
		}
		++prev_hit;
		++curr_hit;
	}
	
	bool saw_antisense_splice = false;
	bool saw_sense_splice = false;
	vector<CigarOp> long_cigar;
	int num_mismatches = 0;
	int num_splice_mms = 0;
	for (list<BowtieHit>::reverse_iterator s = hit_chain.rbegin(); s != hit_chain.rend(); ++s)
	{
		num_mismatches += s->edit_dist();
		num_splice_mms += s->splice_mms();
		if (!s->contiguous())
		{
			if (s->antisense_splice())
			{
				if (saw_sense_splice)
					return BowtieHit();
				saw_antisense_splice = true;
			}
			else
			{
				if (saw_antisense_splice)
					return BowtieHit();
				saw_sense_splice = true;
			}
		}
		const vector<CigarOp>& cigar = s->cigar();
		if (long_cigar.empty())
		{
			long_cigar = cigar;
		}
		else
		{
			CigarOp& last = long_cigar.back();
			last.length += cigar[0].length;
			for (size_t b = 1; b < cigar.size(); ++b)
			{
				long_cigar.push_back(cigar[b]);
			}
		}
	}
	
	BowtieHit new_hit(reference_id,
					  insert_id, 
					  left, 
					  long_cigar, 
					  true,
					  saw_antisense_splice,
					  num_mismatches,
					  num_splice_mms);
	new_hit.seq(seq);
	new_hit.qual(qual);
	
	int new_read_len = new_hit.read_len();
	if (new_read_len != old_read_length)
	{
		//fprintf(stderr, "Warning: malformed closure\n");
		return BowtieHit();
	}
	return new_hit;
}

int multi_closure = 0;
int anchor_too_short = 0;
int gap_too_short = 0;

bool valid_hit(const BowtieHit& bh)
{
	if (bh.insert_id())
	{		
		// validate the cigar chain - no gaps shorter than an intron, etc.
		for (size_t i = 0; i < bh.cigar().size(); ++i)
		{
			const CigarOp& cig = bh.cigar()[i];
			if (cig.opcode == REF_SKIP && cig.length < (uint64_t)min_report_intron_length)
			{
				gap_too_short++;
				return false;
			}
		}
		if (bh.cigar().front().opcode != MATCH || 
			bh.cigar().back().opcode != MATCH /* ||
			(int)bh.cigar().front().length < min_anchor_len||
			(int)bh.cigar().back().length < min_anchor_len*/ )
		{
			anchor_too_short++;
			return false;
		}	
	}
	else
	{
		multi_closure++;
		return false;
	}
	
	return true;
}

void merge_segment_chain(std::set<Junction>& possible_juncs,
						 vector<BowtieHit>& hits,
						 vector<BowtieHit>& merged_hits)
{
	if (hits.size() == 0)
		return;
	
	BowtieHit bh;
	if (hits.size() > 1)
	{
		list<BowtieHit> hit_chain;
	
		copy(hits.begin(), hits.end(), back_inserter(hit_chain));

		if (hit_chain.front().antisense_align())
		{
			bh = merge_antisense_chain(possible_juncs, hit_chain);
		}
		else
		{
			bh = merge_sense_chain(possible_juncs, hit_chain);
		}
	}
	else
	{
		bh = hits[0];
	}
	
//	if (valid_hit(bh))
//		print_hit(stdout, read_name, bh, ref_name, read_seq, read_quals);
	if (valid_hit(bh))
	{
		merged_hits.push_back(bh);
	}
}

bool dfs_seg_hits(std::set<Junction>& possible_juncs, 
				  vector<HitsForRead>& seg_hits_for_read,
				  size_t curr,
				  vector<BowtieHit>& seg_hit_stack,
				  vector<BowtieHit>& joined_hits)
{
	assert (!seg_hit_stack.empty());
	bool join_success = false;
	
	if (curr < seg_hits_for_read.size())
	{
		for (size_t i = 0; i < seg_hits_for_read[curr].hits.size(); ++i)
		{
			BowtieHit& bh = seg_hits_for_read[curr].hits[i];
			BowtieHit& back = seg_hit_stack.back();
			bool consistent_sense = bh.antisense_align() == back.antisense_align();
			bool same_contig = bh.ref_id() == back.ref_id();
			
			// FIXME: when we have stranded reads, we need to fix this condition
			//bool consistent_strand = (bh.contiguous() || back.contiguous() || 
			//						  (bh.antisense_splice() == back.antisense_splice()));
			if (consistent_sense && same_contig /*&& consistent_strand*/)
			{
				if (bh.antisense_align())
				{
					int bh_r = bh.right();
					int back_left = seg_hit_stack.back().left();
					if (bh_r + max_report_intron_length >= back_left &&
						back_left >= bh_r)
					{
						// these hits are compatible, so push bh onto the 
						// stack, recurse, and pop it when done.
						seg_hit_stack.push_back(bh);
						bool success = dfs_seg_hits(possible_juncs,
													seg_hits_for_read, 
													curr + 1,
													seg_hit_stack,
													joined_hits);
						if (success)
							join_success = true;
						
						seg_hit_stack.pop_back();
					}
				}
				else
				{
					int bh_l = bh.left();
					int back_right = seg_hit_stack.back().right();
					if (back_right + max_report_intron_length >= bh_l  &&
						bh_l >= back_right)
					{
						// these hits are compatible, so push bh onto the 
						// stack, recurse, and pop it when done.
						seg_hit_stack.push_back(bh);
						bool success = dfs_seg_hits(possible_juncs,
													seg_hits_for_read, 
													curr + 1,
													seg_hit_stack,
													joined_hits);
						if (success)
							join_success = true;
						
						seg_hit_stack.pop_back();
					}
				}
			}
		}
	}
	else
	{
		// We have a hit, so chain the hits at currently on the stack and 
		// print...
		//const char* ref_name = rt.get_name(seg_hit_stack[0].ref_id());
		merge_segment_chain(possible_juncs,
							seg_hit_stack,
							joined_hits);
		return join_success = true;
	}
	return join_success;
}

bool join_segments_for_read(std::set<Junction>& possible_juncs,
							vector<HitsForRead>& seg_hits_for_read,
							vector<BowtieHit>& joined_hits)
{	
	vector<BowtieHit> seg_hit_stack;
	bool join_success = false;
	
	for (size_t i = 0; i < seg_hits_for_read[0].hits.size(); ++i)
	{
		BowtieHit& bh = seg_hits_for_read[0].hits[i];
		seg_hit_stack.push_back(bh);
		bool success = dfs_seg_hits(possible_juncs,
									seg_hits_for_read, 
									1, 
									seg_hit_stack,
									joined_hits);
		if (success)
			join_success = true;
		seg_hit_stack.pop_back();
	}
	
	return join_success;
}

void join_segment_hits(std::set<Junction>& possible_juncs, 
					   ReadTable& unmapped_reads,
					   RefSequenceTable& rt,
					   FILE* reads_file,
					   vector<HitStream>& contig_hits,
					   vector<HitStream>& spliced_hits)
{
	uint32_t curr_contig_obs_order = 0xFFFFFFFF;
	HitStream* last_seg_contig_stream = NULL;
	uint64_t next_contig_id = 0;
	
	if (contig_hits.size())
	{
		last_seg_contig_stream = &(contig_hits.back());
		next_contig_id = last_seg_contig_stream->next_group_id();
		curr_contig_obs_order = unmapped_reads.observation_order(next_contig_id);
	}
	
	HitsForRead curr_hit_group;
	
	uint32_t curr_spliced_obs_order = 0xFFFFFFFF;
	HitStream* last_seg_spliced_stream = NULL;
	uint64_t next_spliced_id = 0;
	
	if (spliced_hits.size())
	{
		last_seg_spliced_stream = &(spliced_hits.back());
		next_spliced_id = last_seg_spliced_stream->next_group_id();
		curr_spliced_obs_order = unmapped_reads.observation_order(next_spliced_id);	
	}
	
	//size_t spliced_num = last_seg_spliced_stream->size();
	
	while(curr_contig_obs_order != 0xFFFFFFFF || 
		  curr_spliced_obs_order != 0xFFFFFFFF)
	{
		uint32_t read_in_process;
		vector<HitsForRead> seg_hits_for_read;
		seg_hits_for_read.resize(contig_hits.size());
		
		if (curr_contig_obs_order < curr_spliced_obs_order)
		{
			// Get hit group
			last_seg_contig_stream->next_read_hits(curr_hit_group);
			seg_hits_for_read.back() = curr_hit_group;
			
			//uint64_t last_id = next_contig_id;
			next_contig_id = last_seg_contig_stream->next_group_id();
			uint32_t next_order = unmapped_reads.observation_order(next_contig_id);
			
			read_in_process = curr_contig_obs_order;
			curr_contig_obs_order = next_order;
		}
		else if  (curr_spliced_obs_order < curr_contig_obs_order)
		{
			
			// Get hit group
//			curr_hit_group = splice_itr->second;
//			seg_hits_for_read.back() = curr_hit_group;
//			
//			read_in_process = curr_spliced_obs_order;
//			
//			++splice_itr;
//			if (splice_itr != last_seg_spliced_stream->end())
//			{
//				curr_spliced_obs_order = splice_itr->first;
//			}
//			else
//			{
//				curr_spliced_obs_order = 0xFFFFFFFF;
//			}
			
			last_seg_spliced_stream->next_read_hits(curr_hit_group);
			seg_hits_for_read.back() = curr_hit_group;
			
			//uint64_t last_id = next_contig_id;
			next_spliced_id = last_seg_spliced_stream->next_group_id();
			uint32_t next_order = unmapped_reads.observation_order(next_spliced_id);
			
			read_in_process = curr_spliced_obs_order;
			curr_spliced_obs_order = next_order;
		}
		else if (curr_contig_obs_order == curr_spliced_obs_order &&
				 curr_contig_obs_order != 0xFFFFFFFF && 
				 curr_spliced_obs_order != 0xFFFFFFFF)
		{
			last_seg_contig_stream->next_read_hits(curr_hit_group);
			
			HitsForRead curr_spliced_group;
			last_seg_spliced_stream->next_read_hits(curr_spliced_group);
			
			curr_hit_group.hits.insert(curr_hit_group.hits.end(),
									   curr_spliced_group.hits.begin(),
									   curr_spliced_group.hits.end());
			seg_hits_for_read.back() = curr_hit_group;
			
			
			read_in_process = curr_spliced_obs_order;
			
			next_contig_id = last_seg_contig_stream->next_group_id();
			uint32_t next_order = unmapped_reads.observation_order(next_contig_id);
			
			next_spliced_id = last_seg_spliced_stream->next_group_id();
			uint32_t next_spliced_order = unmapped_reads.observation_order(next_spliced_id);
			
			curr_spliced_obs_order = next_spliced_order;
			curr_contig_obs_order = next_order;
		}
		else
		{
			break;
		}
		
		if (contig_hits.size() > 1)
		{
			look_left_for_hit_group(unmapped_reads, 
									contig_hits, 
									contig_hits.size() - 1,
									spliced_hits,
									curr_hit_group,
									seg_hits_for_read);
		}
		
		if (!seg_hits_for_read.empty() && !seg_hits_for_read[0].hits.empty())
		{
			uint64_t insert_id = seg_hits_for_read[0].hits[0].insert_id();
			char read_name[256];
			char read_seq[256];
			char read_alt_name[256];
			char read_quals[256];
			
			if (get_read_from_stream(insert_id,  
									 reads_file,
									 reads_format,
									 false,
									 read_name, 
									 read_seq,
									 read_alt_name,
									 read_quals))
			{
				vector<BowtieHit> joined_hits;
				join_segments_for_read(possible_juncs,
									   seg_hits_for_read, 
									   joined_hits);
				
				sort(joined_hits.begin(), joined_hits.end());
				vector<BowtieHit>::iterator new_end = unique(joined_hits.begin(), joined_hits.end());
				joined_hits.erase(new_end, joined_hits.end());

				for (size_t i = 0; i < joined_hits.size(); i++)
				{
					const char* ref_name = rt.get_name(joined_hits[i].ref_id());
					if (color && !color_out)
					  print_hit(stdout, read_name, joined_hits[i], ref_name, joined_hits[i].seq().c_str(), joined_hits[i].qual().c_str(), true);
					else
					  print_hit(stdout, read_name, joined_hits[i], ref_name, read_seq, read_quals, false);
				}
			}
			else
			{
				fprintf(stderr, "Error: could not get read # %d from stream\n", read_in_process);
				exit(1);
			}
		}
		else
		{
			//fprintf(stderr, "Warning: couldn't join segments for read # %d\n", read_in_process);
		}	
	}
}					   

void driver(vector<FILE*> possible_juncs_files,
			vector<FILE*>& spliced_seg_files,
			vector<FILE*>& seg_files,
			FILE* reads_file)
{	
	if (seg_files.size() == 0)
	{
		fprintf(stderr, "No hits to process, exiting\n");
		exit(0);
	}

	RefSequenceTable rt(true, true);
	ReadTable it;

	bool need_seq = false;
	bool need_qual = false;
	if (color && !color_out)
	  {
	    need_seq = true;
	    need_qual = true;
	  }
	
	//HitFactory hit_factory(it,rt);
	
	rewind(reads_file);
	
	//vector<HitTable> seg_hits;
	vector<HitStream> contig_hits;
	vector<HitStream> spliced_hits;
	
	vector<HitFactory*> factories;
	for (size_t i = 0; i < seg_files.size(); ++i)
	{
		HitFactory* fac = new BowtieHitFactory(it, rt);
		factories.push_back(fac);
		HitStream hs(seg_files[i], fac, false, false, false, need_seq, need_qual);
		contig_hits.push_back(hs);
	}
	
	fprintf(stderr, "Loading spliced hits...");
	
	//vector<vector<pair<uint32_t, HitsForRead> > > spliced_hits;
	for (size_t i = 0; i < spliced_seg_files.size(); ++i)
	{
		int anchor_length = 0;
		if (i == 0 || i == spliced_seg_files.size() - 1)
			anchor_length = min_anchor_len;
		//HitFactory* fac = new SplicedBowtieHitFactory(it, rt, i == 0 || i == spliced_seg_files.size() - 1);
		HitFactory* fac = new SplicedBowtieHitFactory(it, 
													  rt, 
													  anchor_length);
		factories.push_back(fac);
		
		HitStream hs(spliced_seg_files[i], fac, true, false, false, need_seq, need_qual);
		spliced_hits.push_back(hs);
		
		//HitsForRead hit_group;

//		spliced_hits.push_back(vector<pair<uint32_t, HitsForRead> >());
//		while (hs.next_read_hits(hit_group))
//		{
//			if (hit_group.insert_id)
//			{
//				uint32_t key  = it.observation_order(hit_group.insert_id);
//				spliced_hits[i].push_back(make_pair(key,hit_group));
//			}
//		}
		
		//sort(spliced_hits[i].begin(), spliced_hits[i].end(), key_lt);
		//fprintf(stderr, "%d items\n", (int)spliced_hits[i].size());
			
	}
	fprintf(stderr, "done\n");
	
	fprintf(stderr, "Loading junctions...");
	std::set<Junction> possible_juncs;
	
	for (size_t i = 0; i < possible_juncs_files.size(); ++i)
	{
		char buf[2048];
		while(!feof(possible_juncs_files[i]) &&
			  fgets(buf, sizeof(buf), possible_juncs_files[i]))
		{
			char junc_ref_name[256];
			int left;
			int right; 
			char orientation;
			int ret = sscanf(buf, "%s %d %d %c", junc_ref_name, &left, &right, &orientation);
			if (ret != 4)
				continue;
			uint32_t ref_id = rt.get_id(junc_ref_name, NULL, 0);
			possible_juncs.insert(Junction(ref_id, left, right, orientation == '-'));
		}
	}
	fprintf(stderr, "done\n");
	
	join_segment_hits(possible_juncs, it, rt, reads_file, contig_hits, spliced_hits);

	for (size_t fac = 0; fac < factories.size(); ++fac)
	{
		delete factories[fac];
	}
}

int main(int argc, char** argv)
{
	fprintf(stderr, "long_spanning_reads v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "--------------------------------------------\n");
	
    int parse_ret = parse_options(argc, argv, print_usage);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string reads_file_name = argv[optind++];
    
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string juncs_file_list = argv[optind++];
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string segment_file_list = argv[optind++];
	
	
    string spliced_segment_file_list;
	if(optind < argc)
    {
		spliced_segment_file_list = argv[optind++];
    }
	
	fprintf(stderr, "Opening %s for reading\n",
					reads_file_name.c_str());
	FILE* reads_file = fopen(reads_file_name.c_str(), "r");
    if (!reads_file)
    {
        fprintf(stderr, "Error: cannot open %s for reading\n",
                reads_file_name.c_str());
        exit(1);
    }
    
	vector<string> juncs_file_names;
    vector<FILE*> juncs_files;
    tokenize(juncs_file_list, ",",juncs_file_names);
    for (size_t i = 0; i < juncs_file_names.size(); ++i)
    {
		fprintf(stderr, "Opening %s for reading\n",
						juncs_file_names[i].c_str());
        FILE* juncs_file = fopen(juncs_file_names[i].c_str(), "r");
        if (juncs_file == NULL)
        {
            fprintf(stderr, "Warning: cannot open %s for reading\n",
                    juncs_file_names[i].c_str());
            continue;
        }
        juncs_files.push_back(juncs_file);
    }
	
    vector<string> segment_file_names;
    vector<FILE*> segment_files;
    tokenize(segment_file_list, ",",segment_file_names);
    for (size_t i = 0; i < segment_file_names.size(); ++i)
    {
		fprintf(stderr, "Opening %s for reading\n",
						segment_file_names[i].c_str());
        FILE* seg_file = fopen(segment_file_names[i].c_str(), "r");
        if (seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open %s for reading\n",
                    segment_file_names[i].c_str());
            exit(1);
        }
        segment_files.push_back(seg_file);
    }
	
	vector<string> spliced_segment_file_names;
    vector<FILE*> spliced_segment_files;
    tokenize(spliced_segment_file_list, ",",spliced_segment_file_names);
    for (size_t i = 0; i < spliced_segment_file_names.size(); ++i)
    {
		fprintf(stderr, "Opening %s for reading\n",
						spliced_segment_file_names[i].c_str());
        FILE* spliced_seg_file = fopen(spliced_segment_file_names[i].c_str(), "r");
        if (spliced_seg_file == NULL)
        {
            fprintf(stderr, "Error: cannot open %s for reading\n",
                    spliced_segment_file_names[i].c_str());
            exit(1);
        }
        spliced_segment_files.push_back(spliced_seg_file);
    }
	
	if (spliced_segment_files.size())
	{
		if (spliced_segment_files.size() != segment_files.size())
		{
			fprintf(stderr, "Error: each segment file must have a corresponding spliced segment file\n");
			exit(1);
		}
	}
	
    driver(juncs_files, spliced_segment_files, segment_files, reads_file);
    
    return 0;
}
