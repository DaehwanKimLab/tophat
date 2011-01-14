/*
 *  fix_map_ordering.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/28/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */


#include <queue>
#include <cstring>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <getopt.h>
#include "common.h"
#include "reads.h"
#include "bwt_map.h"

using namespace seqan;
using namespace std;

struct MapOrdering
{
	MapOrdering(ReadTable& IT) : it(IT) {} 
	bool operator()(pair<uint64_t, char*>& lhs, pair<uint64_t, char*>& rhs)
	{
		uint64_t lhs_id = lhs.first;
		uint64_t rhs_id = rhs.first;
		
		return lhs_id > rhs_id;
		
		//return it.observation_order(lhs_id) > it.observation_order(rhs_id);
	}
	
	ReadTable& it;
};

void driver(FILE* reads_file, FILE* map_file)
{
	ReadTable it;
	
	char bwt_buf[4096];
	
	MapOrdering order(it);
	
	priority_queue< pair<uint64_t,char*>, 
					vector<pair<uint64_t, char*> >, 
					MapOrdering >  map_pq(it);
	
	while (fgets(bwt_buf, sizeof(bwt_buf), map_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		if (*bwt_buf == 0)
			continue;
		
		const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";
		static const int buf_size = 2048;
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
			//fprintf(stderr, "Warning: found malformed record, skipping\n");
			continue;
		}

		//char* p1 = strdup(name);
		char* p2 = strdup(bwt_buf);
		uint64_t id = it.get_id(name);
		
		map_pq.push(make_pair(id,p2));
		
		if (map_pq.size() > 1000000)
		{
			const pair<uint64_t, char*>& t = map_pq.top();
			printf("%s\n", t.second);
			//free(t.first);
			free(t.second);
			map_pq.pop();
		}
	}
	
	while (map_pq.size())
	{
		const pair<uint64_t, char*>& t = map_pq.top();
		printf("%s\n", t.second);
		
		free(t.second);
		map_pq.pop();
	}
}

void print_usage()
{
    fprintf(stderr, "Usage:   fix_map_ordering <reads.fa/.fq> <map.bwtout>\n");
}

//const char *short_options = "fq";
//
//static struct option long_options[] = {
//{"fasta",       required_argument,       0,            'f'},
//{"fastq",       required_argument,       0,            'q'},
//{0, 0, 0, 0} // terminator
//};

//int parse_options(int argc, char** argv)
//{
//    int option_index = 0;
//    int next_option;
//    do {
//        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
//        switch (next_option) { 
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
    
    string reads_file_name = argv[optind++];
    
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    string map_file_name = argv[optind++];
	
	FILE* reads_file = fopen(reads_file_name.c_str(), "r");
	if (reads_file == NULL)
	{
		fprintf(stderr, "Error: cannot open reads file %s for reading\n",
				reads_file_name.c_str());
		exit(1);
	}
    
	FILE* map_file = map_file_name == "-" ? stdin : fopen(map_file_name.c_str(), "r");
	if (map_file == NULL)
	{
		fprintf(stderr, "Error: cannot open Bowtie map file %s for reading\n",
				map_file_name.c_str());
		exit(1);
	}
	
    driver(reads_file, map_file);
    
    return 0;
}
