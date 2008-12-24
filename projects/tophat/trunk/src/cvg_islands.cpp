/*
 *  cvg_islands.cpp
 *  CSAMapper
 *
 *  Created by Cole Trapnell on 7/13/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <deque>
#include <stdio.h>
#include <zlib.h>
#include <cassert>
#include "alphabet.h"
#include "islands.h"

#define LINE_LEN 60

using namespace std;



vector<Island> island_metadata;

void add_island_metadata(const char* ref_src,
						 int start_pos,
						 int length,
						 int depth,
						 int id)
{
	island_metadata.push_back(Island(ref_src, start_pos, length, depth, id)); 
}

void compute_island_d_stats(vector<Island>& islands,
							unsigned long long total_map_depth, 
							unsigned long long total_map_length)
{ 
	long double map_d = (long double) total_map_depth;// / (long double) total_map_length;
	
	for (size_t j = 0; j < islands.size(); ++j)
	{
		Island& i = islands[j];
		long double il_d = (float) i.total_depth / (float) i.length;
		
		i.d_stat = il_d / map_d;
	}
}
						 

void print_islands_gff(FILE* fp_gff, vector<Island>& islands) 
{
	long double min_d_stat = 99999999999.0;
	long double max_d_stat = 0.0;
	
	for (size_t j = 0; j < islands.size(); ++j)
	{
		Island& i = islands[j];
		min_d_stat = min(min_d_stat, i.d_stat);
		max_d_stat = max(max_d_stat, i.d_stat);
	}
	
	fprintf(fp_gff,"track name=TopHat_islands description=\"TopHat RNA-Seq islands\" useScore=1\n");
 
	for (size_t j = 0; j < islands.size(); ++j)
	{
		Island& i = islands[j];
		int d_score = (int)((( i.d_stat ) / max_d_stat) * 1000.0);
		fprintf(fp_gff,"%s\tTopHat\tisland\t%d\t%d\t%d\t.\t.\tIL%08d\n", 
				i.ref_src, i.start_pos, i.start_pos + i.length, d_score, i.id);
	}
}

void print_island(FILE* fp_fasta,
				  FILE* fp_gff,
				  const char* name,
				  const string& seq, 
				  int start_pos,
				  int depth,
				  int length, 
				  unsigned int id)
{
	const char* ref_src = name;
	 
	fprintf(fp_fasta, ">IL%08d %s,%d\n", id, ref_src, start_pos);
	unsigned int i = 0;

	for (; i + 60 < seq.length(); i += 60)
		fprintf(fp_fasta, "%s\n", seq.substr(i, 60).c_str());
	fprintf(fp_fasta, "%s\n", seq.substr(i, (seq.length()  - i)).c_str());
	
	add_island_metadata(ref_src, start_pos + 1, length, depth, id);
	
	//fprintf(fp_gff,"%s\tTopHat\tisland\t%d\t%d\t.\t.\tIL%08d\n", 
	//		ref_src, start_pos + 1, (int)(start_pos + seq.length() + 1), (int)id);
}

/** Reports a list of coverage islands from the assembled map.  Islands may be
 partial exons or non-coding RNAs.  This function takes three files as input: a 
 zipped Maq .cns file produced by `maq assemble`, a GFF file describing the 
 islands coordinates and other information, and a FASTA file containing the 
 island sequences.  The function depends on a number of global parameters, 
 defined above.
 */

// FIXME: This function desperately needs refactoring.  It's a mess.  We should be 
// encapsulating the gzipped consensus file in a class that manages the 
// (reference,donor) base-call pairs in a way that lets us get a bunch of them
// out of the buffer at a time, but that doesn't require fseeking in the 
// zipped file.  I'm delaying the refactoring because TopHat may drop Maq in 
// the near future.
void print_islands(FILE* fp_fasta,
				   FILE* fp_gff,
				   gzFile fp, 
				   float min_depth, 
				   int break_length,
				   int extend_exons, 
				   char qual_thresh,
				   int clear_bases,
				   bool use_ref_seq)
{
	int len;
	unsigned int island_id = 1;
	
	// This counter sums the depth of the map at each base for later 
	// so we can later normalize island coverage levels.
	unsigned long long total_map_depth = 0;
	unsigned long long total_map_length = 0;
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		//uint8_t *qual;
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		//qual = (uint8_t*)malloc(len);
		fprintf(stderr, "processing islands in: %s (%dbp)\n", name, len);
		string seq;
		
		vector<unsigned int> covs;
		
		// This buffer keeps values from the zipped consensus for use with exon
		// extension.  Elements are (low, high) consensus value pairs.  The 
		// buffer keeps the most recent (rightmost) elements at the front.
		deque<pair<uint32_t,uint32_t> > raw_cns_vals;
		static const unsigned int raw_cns_buf_size = 10000; 
		
		deque<pair<uint32_t, uint32_t> > lookforward_buf;
		int max_extend = extend_exons;
		
		//string quals;
		bool skipping = true;
		
				
		// A counter to track consecutive N's in the consensus.  
		// break the consensus string when consec_Ns >= break_length
		int consec_Ns = 0;
		
		int island_begin = 0;
		int processed_bases = 0;
		int unzipped_bases = 0;
		unsigned int bases_to_skip = 0; 
		
		// Process individual columns in the consensus stream
		while (processed_bases < len) {
			
			uint32_t low, high;
			
			if (lookforward_buf.empty())
			{
				gzread(fp, &low, sizeof(uint32_t));
				gzread(fp, &high, sizeof(uint32_t));
				++unzipped_bases;
			}
			else
			{
				//assert (processed_bases + lookforward_buf.size() < len);
				low = lookforward_buf.front().first;
				high = lookforward_buf.front().second;
				lookforward_buf.pop_front();
			}
			
			++processed_bases;
			
			// If the consensus buffer is full, pop a character off the back
			if (raw_cns_vals.size() >= raw_cns_buf_size)
				raw_cns_vals.pop_back();
			
			// Add the new (reference, donor) base-call pair to the consensus
			// buffer
			raw_cns_vals.push_front(make_pair<uint32_t, uint32_t>(low, high));
			
			char c = nst_nt16_rev_table[(high>>24) & 0xf];
			
			// Extract the reference base-call from the packed consensus
			// entry at this contig position
			char r = nst_nt16_rev_table[high>>28];
			assert (r != 'X');
			if (c == 'N')
				++consec_Ns;
			else
				consec_Ns = 0;
			
			if (skipping == true)
			{
				if (consec_Ns == 0)
				{
					char s = c;
					char q = ((high>>16&0xff) + 33 > 126) ? 126 : ((high>>16&0xff) + 33);
					unsigned int cov = low & 0xff;
					if (use_ref_seq || q < qual_thresh || cov < min_depth)
						s = r;
					bases_to_skip = clear_bases;
					if (bases_to_skip > 0)
						--bases_to_skip;
					else
					{
						seq.push_back(s);
						covs.push_back(cov);
						total_map_depth += cov;
					}
					
					skipping = false;
					
					// Stay in 0-based coordinates
					island_begin = processed_bases - 1;
				}
			}
			else
			{
				if (consec_Ns >= break_length || processed_bases == (len - 1))
				{
					skipping = true;
					bases_to_skip = 0;
					
					//unsigned int s_len = seq.length() - (consec_Ns - 1);
					// Trim off the accumulated 'N' bases.
					if (processed_bases != (len - 1))
					{
						
						seq.resize(seq.length() - (consec_Ns - 1) - clear_bases);
						covs.resize(covs.size() - (consec_Ns - 1) - clear_bases);
					}
					
					// Compute the average coverage for this island
					unsigned int total_depth = accumulate(covs.begin(), covs.end(), 0);
					float avg_cov = total_depth / (float)covs.size();
					
					string extension;
					
					string extension_fwd;
					if (avg_cov >= min_depth)
					{
						
						if (extend_exons + clear_bases > 0)
						{
							// raw_cns_vals contains consec_Ns 'N' chars at this
							// point, and seq contains none, since they were 
							// trimmed off, above.
							int k = seq.length() + consec_Ns + clear_bases;
							total_map_length += seq.length() + clear_bases;
							size_t k_start = k;
							unsigned int fwd_pos = island_begin;
													
														
							while (k >= 0 && k < (int)raw_cns_vals.size() &&
								  (int)(k - k_start) < max_extend + clear_bases)
							{
								uint8_t ref_bits_k = (raw_cns_vals[k].second >> 28);
								char ref_base_k = nst_nt16_rev_table[ref_bits_k];
								extension.push_back(ref_base_k);
								k++;
							}
							
							reverse(extension.begin(), 
									extension.end());
							
							// extension already contains the left set of 
							// trimmed bases
							extension_fwd = extension + seq;
							
							fwd_pos -= (extension.length() - clear_bases);
							
							
							extension.clear();
							string trim;
							
							// Add back whatever we trimmed off, taking bases
							// from the reference
							k = consec_Ns;
							while (k >= 0 && k < (int)raw_cns_vals.size() &&
								   k < (int)clear_bases + consec_Ns)
							{
								// Extract the reference base from the consensus
								// history buffer
								uint8_t ref_bits_k = (raw_cns_vals[k].second >> 28);
								char ref_base_k = nst_nt16_rev_table[ref_bits_k];
								assert (ref_base_k != 'X');
								trim.push_back(ref_base_k);
								
								k++;
							}
							
							// All bases from the history buffer come in 
							// reverse order, so reverse them before adding to
							// the island
							if (trim.length())
							{
								reverse(trim.begin(), 
										trim.end());
							}
							
							// If we are extending, we need to first take 
							// bases from the reference that we have already
							// pulled out of the gzipped stream.
							k = consec_Ns - min((int)consec_Ns, max_extend);
							k_start = k;
							while (k >= 0 && k < (int)raw_cns_vals.size() &&
								   (int)(k - k_start) < min((int)consec_Ns, max_extend) )
							{
								// Extract the reference base from the consensus
								// history buffer
								uint8_t ref_bits_k = (raw_cns_vals[k].second >> 28);
								char ref_base_k = nst_nt16_rev_table[ref_bits_k];
								assert (ref_base_k != 'X');
								extension.push_back(ref_base_k);
								
								k++;
							}
							
							// All bases from the history buffer come in 
							// reverse order, so reverse them before adding to
							// the island
							if (extension.length())
							{
								reverse(extension.begin(), 
										extension.end());
							}
							
							// Now we have exhausted bases we've already pulled
							// out of the gzipped stream, so if we need to do
							// more extending, we need to upzip more bases.
							int bound = min(len - processed_bases, (int)(max_extend - (k - k_start)));
							for (k = 0; k < bound; ++k)
							{
								if (k >= (int)lookforward_buf.size())
								{
									// the _extensions_ of two otherwise disjoint
									// islands might overlap, so we can't just
									// toss the bases from the stream; we need
									// to save them until we decide they aren't 
									// part of an upcoming island.
									gzread(fp, &low, sizeof(uint32_t));
									gzread(fp, &high, sizeof(uint32_t));
									++unzipped_bases;
									lookforward_buf.push_back(make_pair<uint32_t, uint32_t>(low, high));
								}
								else 
								{
									low = lookforward_buf[k].first;
									high = lookforward_buf[k].second;
								}
								
								uint8_t ref_bits_k = (high >> 28);
								char ref_base_k = nst_nt16_rev_table[ref_bits_k];
								extension.push_back(ref_base_k);
							}
							
							extension_fwd += trim;
							extension_fwd += extension;
							
							//total_map_length += extension_fwd.length();
							print_island(fp_fasta, 
										 fp_gff,
										 name, 
										 extension_fwd, 
										 fwd_pos, 
										 total_depth, 
										 extension_fwd.length(), 
										 island_id++);
						}
						else
						{
							print_island(fp_fasta,
										 fp_gff,
										 name, 
										 seq, 
										 island_begin + 1, 
										 total_depth, 
										 seq.length(), 
										 island_id++);	
							total_map_length += seq.length();
						}
						
						
					}
					seq.clear();
					covs.clear();
					//lookforward_buf.clear();
				}
				else
				{
					char s = c;
					char q = ((high>>16&0xff) + 33 > 126) ? 126 : ((high>>16&0xff) + 33);
					unsigned int cov = low & 0xff;
					if (use_ref_seq || q < qual_thresh || cov < min_depth)
						s = r;
					if (bases_to_skip > 0)
						--bases_to_skip;
					else
					{
						seq.push_back(s);
						covs.push_back(cov);
						total_map_depth += cov;
					}
				}
			}
		} 
		assert (processed_bases == unzipped_bases);
		//free(name);
	}
	compute_island_d_stats(island_metadata, total_map_depth, total_map_length);
	print_islands_gff(fp_gff, island_metadata);
}

static bool verbose = false;

static void print_usage() 
{
	cerr << "Usage: cvg_islands [options] <in.cns> <out.fa> <out.gff>" << endl;
	cerr << "    -v           verbose output [default: off]" << endl;
	cerr << "    -R           Ignore SNPs in the donor - always take the reference base at each position [default: on]" << endl;
	cerr << "    -d  <int>    minimum avg. depth of coverage for an island to be reported [default: 0.0]" << endl;
	cerr << "    -b  <int>    Maximum gap length before splitting one island into two [default: 6bp]" << endl;
	cerr << "    -e  <int>    extend islands by this length on each side [default: 45bp]" << endl;
	cerr << "    -c  <int>    # of bases to trim from donor sequence before extending with reference sequence.  Ignored if -R specified [default: 2]" << endl;   
	cerr << "    -q  <char>   Only use donor base-calls with this Phred quality or higher. Ignored if -R specified [default: 'I']" << endl;   
	
	cerr << "\nPrints FASTA and GFF files of coverage islands from a Maq consensus" << endl;
	cerr << "Fasta header is >REF_SOURCE_{ISLAND_ID,POS_IN_SOURCE,AVG_COV}" << endl;
}


/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			print_usage();
			exit(1);
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	print_usage();
	exit(1);
	return -1;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloat(float lower, const char *errmsg) {
	float l;
	l = atof(optarg);
	
	if (l < lower) {
		cerr << errmsg << endl;
		print_usage();
		exit(1);
	}
	return l;
	
	cerr << errmsg << endl;
	print_usage();
	exit(1);
	return -1;
}

int main(int argc, char** argv)
{
	string cns_filename;
	const char *short_options = "vd:b:e:q:c:R";
	float min_depth = 0.0;
	unsigned int break_length = 4;
	char qual_thresh = 33;
	unsigned int extend_exons = 30;
	unsigned int clear_bases = 0;
	bool use_ref_seq = false;
	int next_option;
	do { 
		next_option = getopt(argc, argv, short_options);
		switch (next_option) {
	   		case 'v': /* verbose */
				verbose = true;
				break;
			case 'd':
				min_depth = parseFloat(0, "-d arg must be at least 0");
				break;
			case 'b':
				break_length = parseInt(0, "-b arg must be at least 0");
				break;
			case 'e':
				extend_exons = parseInt(0, "-e arg must be at least 0");;
				break;
			case 'c':
				clear_bases = parseInt(0, "-c arg must be at least 1");;
				break;
			case 'R':
				use_ref_seq = true;
				break;
			case 'q':
				if (strlen(optarg) > 1)
				{
					cerr << "-q takes a single quality char as its argument" << endl;
					exit(1);
				}
				qual_thresh = *(optarg);
				break;
			case -1: /* Done with options. */
				break;
			default: 
				print_usage();
				return 1;
		}
	} while(next_option != -1);
	
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	cns_filename = argv[optind++];
	gzFile cns_file = gzopen(cns_filename.c_str(), "r");
	if (!cns_file)
	{
		fprintf(stderr, "Error: couldn't open consensus file %s\n", cns_filename.c_str());
		return 1;
	}
	
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	string fasta_file_name = argv[optind++];
	
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	string gff_file_name = argv[optind++];
	
	
	FILE* fasta_file = fopen(fasta_file_name.c_str(), "w");
	
	
	if (fasta_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",fasta_file_name.c_str());
		exit(1);
	}
	
	FILE* gff_file = fopen(gff_file_name.c_str(), "w");
	if (gff_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",gff_file_name.c_str());
		exit(1);
	}
	
	
	//fprintf(stderr, "Beginning run with params: b = %d\n", break_length);
	print_islands(fasta_file,
				  gff_file,
				  cns_file, 
				  min_depth, 
				  break_length, 
				  extend_exons, 
				  qual_thresh,
				  clear_bases,
				  use_ref_seq);
	return 0;
}

