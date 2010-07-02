#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  wiggles.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
#include <cassert>
#include <stdio.h>
#include <vector>
#include <string>
#include "wiggles.h"

using namespace std;

void print_wiggle_header(FILE* coverage_out)
{
	fprintf(coverage_out, "track type=bedGraph name=\"TopHat - read coverage\"\n");
}

void print_wiggle_for_ref(FILE* coverage_out,
						  const string& ref_name,
						  const vector<unsigned int>& DoC)
{
	unsigned short last_doc = 0; // Last DoC value we wrote to the file
	size_t last_pos = 0; // Postition where the last written DoC came from
	//fprintf(coverage_out, "variableStep chrom=chr%s\n", name.c_str());
	for (size_t i = 0; i < DoC.size(); ++i)
	{
		if (last_doc != DoC[i])
		{
			size_t j = last_pos;
			while (i - j > 10000000)
			{
				fprintf(coverage_out,"%s\t%d\t%d\t%d\n",ref_name.c_str(),(int)j, (int)j + 10000000, last_doc);
				j += 10000000;
			}
			if (i > 0)
			{
				fprintf(coverage_out,"%s\t%d\t%d\t%d\n",ref_name.c_str(),(int)j, (int)i, last_doc);
			}
			last_pos = i;
			last_doc = DoC[i];
		}
	}
	if (last_doc)
	{
		fprintf(coverage_out,"%s\t%d\t%d\t%d\n",ref_name.c_str(),(int)last_pos, (int)DoC.size() -1, last_doc);
	}
}

void driver(FILE* map_file, FILE* coverage_file)
{	
    ReadTable it;
    RefSequenceTable rt(true);
	
    SAMHitFactory hit_factory(it,rt);
	
	char bwt_buf[2048];
	bwt_buf[0] = 0;
	
	uint32_t last_ref = 0;
	
	vector<unsigned int> DoC;
	
	print_wiggle_header(coverage_file);
	
	while (map_file && !feof(map_file))
	{
		fgets(bwt_buf, 2048, map_file);
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		// Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		
		if (hit_factory.get_hit_from_buf(bwt_buf, bh, false))
		{
			if (bh.ref_id() != last_ref)
			{
				if (last_ref != 0)
				{
					// print wiggle
					
					string ref_name = rt.get_name(last_ref);
					print_wiggle_for_ref(coverage_file,
											  ref_name,
										 DoC);
				}
				DoC.clear();
			}
			
			last_ref = bh.ref_id();
			add_hit_to_coverage(bh, DoC);
		}
	}
	
	if (last_ref != 0)
	{
		string ref_name = rt.get_name(last_ref);
		print_wiggle_for_ref(coverage_file,
							 ref_name,
							 DoC);
		
	}
}

void print_usage()
{
    fprintf(stderr, "Usage:   wiggles <accepted_hits.sam> <coverage.wig>\n");
}


int main(int argc, char** argv)
{
	fprintf(stderr, "wiggles v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "---------------------------------------\n");
	
    int parse_ret = parse_options(argc, argv, print_usage);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string map_filename = argv[optind++];
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string coverage_file_name = argv[optind++];
	
	FILE* map_file = fopen(map_filename.c_str(), "r");
	if (!map_file)
	{
		fprintf(stderr, "Error: cannot open map file %s for reading\n",
				map_filename.c_str());
		exit(1);
	}
	
    // Open the approppriate files
	
    FILE* coverage_file = fopen((coverage_file_name).c_str(), "w");
    if (coverage_file == NULL)
    {
        fprintf(stderr, "Error: cannot open WIG file %s for writing\n",
                coverage_file_name.c_str());
        exit(1);
    }

    driver(map_file, coverage_file);
	
    return 0;
}

