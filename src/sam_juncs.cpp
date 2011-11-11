/*
 *  sam_juncs.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 7/5/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif

#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "junctions.h"


void get_junctions_from_hitstream(HitStream& hitstream,
								  ReadTable& it, 
								  JunctionSet& junctions)
{
	HitsForRead curr_hit_group;
	
	
	hitstream.next_read_hits(curr_hit_group);
	
	
	uint32_t curr_obs_order = it.observation_order(curr_hit_group.insert_id);
	
	
	// While we still have unreported hits...
	while(curr_obs_order != VMAXINT32)
	{		
        for (size_t i = 0; i < curr_hit_group.hits.size(); ++i)
        {
            const BowtieHit& bh = curr_hit_group.hits[i];
            junctions_from_alignment(bh, junctions);
        }
        
        //fprintf(stderr, "#Hits = %d\n", curr_hit_group.hits.size()); 
        
        //curr_hit_group = HitsForRead();
        // Get next hit group
        
        hitstream.next_read_hits(curr_hit_group);
        curr_obs_order = it.observation_order(curr_hit_group.insert_id);
	}
	
	hitstream.reset();
}


void driver(FILE* hit_map)
{	
    ReadTable it;
    
    RefSequenceTable rt(sam_header, true);
	
    SAMHitFactory hit_factory(it,rt);
	
	//HitStream hitstream(hit_map, &hit_factory, false, true, true);
	
    JunctionSet junctions;
    
    while (hit_map && !feof(hit_map))
	{
        char bwt_buf[2048];
		if (!fgets(bwt_buf, 2048, hit_map))
        {
            break;
        }
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		// Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		
		if (hit_factory.get_hit_from_buf(bwt_buf, bh, false))
		{
            junctions_from_alignment(bh, junctions);

		}
	}
    
    for (JunctionSet::iterator itr = junctions.begin();
         itr != junctions.end();
         ++itr)
	{
        const char* ref_name = rt.get_name(itr->first.refid);
    
        fprintf(stdout,
                "%s\t%d\t%d\t%c\n",
                ref_name,
                itr->first.left - 1,
                itr->first.right,
                itr->first.antisense ? '-' : '+');
    }
    
	fprintf(stderr, "Extracted %lu junctions\n", junctions.size()); 
}

void print_usage()
{
    fprintf(stderr, "Usage:   sam_juncs <hits.sam>\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

int main(int argc, char** argv)
{
	fprintf(stderr, "sam_juncs v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
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
    
    string map_filename = argv[optind++];

	FILE* map_file = fopen(map_filename.c_str(), "r");
	if (!map_file)
	{
		fprintf(stderr, "Error: cannot open map file %s for reading\n",
				map_filename.c_str());
		exit(1);
	}

    driver(map_file);
    
    return 0;
}
