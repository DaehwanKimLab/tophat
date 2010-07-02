/*
 *  mask_sam.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 9/10/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"

void print_usage()
{
    fprintf(stderr, "Usage: mask_sam hits.sam repeatmask.out\n");
}

struct Repeat
{
	Repeat(uint32_t i, int l, int r) : ref_id(i), left(l), right(r) {}
	uint32_t ref_id;
	int left;
	int right;
	
	bool operator<(const Repeat& rhs) const
	{
		if (ref_id != rhs.ref_id)
			return ref_id < rhs.ref_id;
		if (left != rhs.left)
			return left < rhs.left;
		if (right != rhs.right)
			return right < rhs.right;
		return false;
	}	
};

void driver(FILE* map_file, FILE* repeat_file)
{
	RefSequenceTable rt(true, true);
	ReadTable it;
	
	SAMHitFactory factory(it, rt);
	
	vector<Repeat > repeats;
	
	char buf[1024 * 16];
	while(!feof(repeat_file) &&
		  fgets(buf, sizeof(buf), repeat_file))
	{
		char chrom[1024];
		
		int left;
		int right; 
		int score;
		float t1,t2,t3;
		int ret = sscanf(buf, "%d %f %f %f %s %d %d", &score,&t1,&t2,&t3,chrom,&left,&right);
		if (ret < 7)
			continue;
		
		int id = rt.get_id(chrom, NULL, 0);
		repeats.push_back(Repeat(id, left, right));
	}
	
	
	sort(repeats.begin(), repeats.end());
	
	vector<Repeat>::iterator curr_repeat = repeats.begin();
	
	static int buf_size = 4096;
	char bwt_buf[buf_size];
	
	uint32_t last_ref_id_seen = 0;
	int last_pos_seen = 0;
	
	while (fgets(bwt_buf, 4096, map_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		if (*bwt_buf == 0)
			continue;
		
		string clean_buf = bwt_buf;
		// Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		
		if (factory.get_hit_from_buf(bwt_buf, bh, true))
		{

			if (bh.left() < last_pos_seen)
			{
				fprintf(stderr, "Error: this SAM file doesn't appear to be correctly sorted!\n");
				fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
						rt.get_name(bh.ref_id()),
						bh.left(),
						rt.get_name(last_ref_id_seen),
						last_pos_seen);
				
				exit(1);
			}
			
			while (curr_repeat != repeats.end() &&
				   (curr_repeat->ref_id < bh.ref_id() ||
					curr_repeat->ref_id == bh.ref_id() && curr_repeat->right < bh.left()))
			{
				curr_repeat++;
			}
			
			vector<Repeat>::iterator next_repeat = curr_repeat;
			
			bool found_mask = false;
			while (next_repeat != repeats.end() &&
				   next_repeat->left < bh.right())
			{
				if (next_repeat->ref_id == bh.ref_id() &&
					next_repeat->left <= bh.left() && next_repeat->right >= bh.right())
				{
					found_mask = true;
					break;
				}
				next_repeat++;
			}
			if (!found_mask)
			{
				printf("%s\n", clean_buf.c_str());
			}
		}
	}
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
    
    string map_file_name = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string repeat_file_name = argv[optind++];
	
	    
	FILE* map_file = fopen(map_file_name.c_str(), "r");
	if (map_file == NULL)
	{
		fprintf(stderr, "Error: cannot open SAM map file %s for reading\n",
				map_file_name.c_str());
		exit(1);
	}
	
	FILE* repeat_file = fopen(repeat_file_name.c_str(), "r");
	if (repeat_file == NULL)
	{
		fprintf(stderr, "Error: cannot open SAM map file %s for reading\n",
				repeat_file_name.c_str());
		exit(1);
	}
	
	driver(map_file, repeat_file);
    
    return 0;
}
