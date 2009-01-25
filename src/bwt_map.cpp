/*
 *  bwt_map.cpp
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
				int8_t left_splice_pos = atoi(toks[2].c_str()) - left + 1;
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

AlignStatus status(const BowtieHit* align)
{
	if (!align)
		return UNALIGNED;
	if (align->splice_pos_left == -1)
		return CONTIGUOUS;
	return SPLICED;
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
