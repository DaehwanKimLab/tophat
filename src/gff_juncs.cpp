/*
 *  gff_juncs.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/15/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <getopt.h>
#include <string>
#include <cstdio>

#include "common.h"
#include "FSA/gff.h"

using namespace std;
using namespace fsa;

void print_usage()
{
    fprintf(stderr, "Usage:   gff_juncs <genes.gff>\n");
}

//const char *short_options = "v";
//
//static struct option long_options[] = {
//{"verbose",		no_argument,	0,	'v'},
//{0, 0, 0, 0} // terminator
//};
//
//int parse_options(int argc, char** argv)
//{
//	int option_index = 0;
//	int next_option; 
//	do { 
//		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);		
//		switch (next_option) {
//			case 'v':
//				verbose = true;
//				break;
//			case -1: /* Done with options. */
//				break;
//			default: 
//				print_usage();
//				return 1;
//		}
//	} while(next_option != -1);
//	
//	return 0;
//}

uint32_t get_junctions_from_gff(const GFF_database& gff_db)
{
	// Table to hold the exons for each transcript, we'll make introns 
	// from these below
	
	typedef map<string, const GFF*> TransTable;
	typedef map<const GFF*, vector<const GFF*> > TransExonTable; 
	TransTable transcripts;
	TransExonTable transcript_exons;
	
    uint32_t num_juncs_reported = 0;
	
	for(GFF_database::const_iterator gff_itr = gff_db.begin();
		gff_itr != gff_db.end();
		++gff_itr)
	{
		const GFF& gff_rec = *gff_itr;
		if (gff_rec.type == "mRNA")
		{
			GFF::AttributeTable::const_iterator att_itr;
			att_itr = gff_rec.attributes.find("ID");
			if (att_itr == gff_rec.attributes.end() ||
				att_itr->second.size() != 1)
			{
				cerr << "Malformed transcript record " << gff_rec << endl; 
				continue;
			}
			const string& id = att_itr->second.front();
			transcripts.insert(make_pair(id, &gff_rec));
			transcript_exons.insert(make_pair(&gff_rec, vector<const GFF*>()));
		}

	}
	
	for(GFF_database::const_iterator gff_itr = gff_db.begin();
		gff_itr != gff_db.end();
		++gff_itr)
	{
		const GFF& gff_rec = *gff_itr;
		if (gff_rec.type == "exon")
		{
			GFF::AttributeTable::const_iterator att_itr;
			att_itr = gff_rec.attributes.find("Parent");
			if (att_itr == gff_rec.attributes.end())
			{
				cerr << "Malformed exon record " << gff_rec << endl; 
				continue;
			}
			vector<string> parent_transcripts = att_itr->second;
			for (vector<string>::iterator par_itr = parent_transcripts.begin();
				 par_itr != parent_transcripts.end();
				 ++par_itr)
			{
				const string& parent_str = *par_itr;
				TransTable::iterator parent_rec = transcripts.find(parent_str);
				if (parent_rec == transcripts.end())
				{
					cerr << "No transcript with id " << parent_str << endl;
					continue;
				}
				TransExonTable::iterator parent_exons = transcript_exons.find(parent_rec->second);
				if (parent_exons == transcript_exons.end())
				{
					cerr << "No exons with for " << *par_itr << endl;
					continue;
				}
				parent_exons->second.push_back(&gff_rec);
			}
		}
	}
	
	
	for (TransExonTable::iterator trans_itr = transcript_exons.begin();
		 trans_itr != transcript_exons.end();
		 ++trans_itr)
	{
		vector<const GFF*>& exons = trans_itr->second;
		if (exons.size() <= 1)
			continue;
		
		vector<const GFF*>::iterator prev_itr = exons.begin();
		vector<const GFF*>::iterator curr_itr = ++(exons.begin());
		while (curr_itr != exons.end())
		{
			const GFF* five_prime_ex = *prev_itr;
			const GFF* three_prime_ex = *curr_itr;
			if (three_prime_ex->start < five_prime_ex->end)
			{
				fprintf(stderr, "Error: bad transcript annotation:\n");
				cerr << *(trans_itr->first);
				fprintf(stderr, "Offending exons overlapped:\n");
				cerr << *five_prime_ex;
				cerr << *three_prime_ex;
				//exit(2);
                break;
			}
			
			fprintf(stdout, "%s\t%d\t%d\t%c\n",
					five_prime_ex->seqid.c_str(),
					five_prime_ex->end - 1,
					three_prime_ex->start - 1,
					five_prime_ex->strand);
            ++num_juncs_reported;
			++curr_itr;
			++prev_itr;
		}
	}
    return num_juncs_reported;
	
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc, argv, print_usage);
	if (parse_ret)
		return parse_ret;
	
    	
	if(optind >= argc) 
	{
		print_usage();
		return 2;
	}
	
	string gff_filename = argv[optind++];
	
	GFF_database gff_db;
	if (gff_filename == "")
	{
		print_usage();
		exit(2);
	}
	
	fprintf(stderr, "gff_juncs v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "---------------------------\n");
	
	gff_db.from_file(gff_filename);
    gff_db.sort_entries();
	
	uint32_t num_juncs_reported = get_junctions_from_gff(gff_db); 
    fprintf(stderr, "Extracted %u junctions from %s\n", 
            num_juncs_reported, gff_filename.c_str());
    if (!num_juncs_reported)
        return 1;
    return 0;
}
