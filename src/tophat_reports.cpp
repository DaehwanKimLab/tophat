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
#endif

#include <cassert>
#include <cstdio>
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
#include "FSA/gff.h"
#include "FSA/sequence.h"

#ifdef PAIRED_END
#include "inserts.h"
#endif

using namespace std;
using namespace seqan;
using fsa::GFF;
using fsa::GFF_database;
using std::set;



static bool filter_junctions = true;
static float min_isoform_fraction = 0.15;
static bool accept_all = true;
string gff_file = "";
string ref_fasta = "";
string output_dir = "tophat_out";
bool verbose = false;

void print_usage()
{
#ifdef PAIRED_END
    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <map1.bwtout> [splice_map1.sbwtout] [map2.bwtout] [splice_map2.sbwtout]\n");
#else
	fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <map1.bwtout> [splice_map1.sbwtout]\n");
#endif
}

#ifdef PAIRED_END
void insert_best_pairings(SequenceTable& rt,
						  HitTable& hits1,
						  HitTable& hits2,
						  BestInsertAlignmentTable& best_pairings)						
{
	for(SequenceTable::const_iterator ci = rt.begin();
		ci != rt.end();
		++ci) 
	{		
		
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
		vector<size_t> map1_singletons;
		vector<size_t> map2_singletons;
		vector<pair<size_t, size_t> > happy_mates;
		
		string name = ci->first;
		uint32_t ref_id = rt.get_id(name);
		HitList* hits1_in_ref = hits1.get_hits(ref_id);
		HitList* hits2_in_ref = hits2.get_hits(ref_id);
		
		if (!hits1_in_ref || !hits2_in_ref)
			continue;
		
		if (verbose)
			fprintf(stderr, "Looking for best insert mappings in %s\n", name.c_str());
		
		best_insert_mappings(ref_id,
							 name, 
							 *hits1_in_ref, 
							 *hits2_in_ref, 
							 best_pairings);
	}	
}
#endif

void fragment_best_alignments(SequenceTable& rt,
							  HitTable& hits1,
							  BestFragmentAlignmentTable& best_alignments)						
{
	for(SequenceTable::const_iterator ci = rt.begin();
		ci != rt.end();
		++ci) 
	{		
		
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
				
		string name = ci->first;
		uint32_t ref_id = rt.get_id(name);
		HitList* hits_in_ref = hits1.get_hits(ref_id);
		
		if (!hits_in_ref)
			continue;
		
		if (verbose)
			fprintf(stderr, "Looking for best alignments in %s\n", name.c_str());
		
		best_fragment_mappings(ref_id,
							 name, 
							 *hits_in_ref, 
							 best_alignments);
	}	
}

typedef map<string, const GFF*> GeneTable;
typedef map<const GFF*, vector<const GFF*> > GeneExonTable; 

uint32_t total_exonic_depth(const GeneExonTable& genes, 
							const vector<short>& DoC,
							fsa::Sequence* ref_str)
{
	uint32_t total_depth = 0;
	for (GeneExonTable::const_iterator gene_itr = genes.begin(); 
		 gene_itr != genes.end();
		 ++gene_itr)
	{
		const GFF& gene_gff = *(gene_itr->first);
		const vector<const GFF*>& exons = gene_itr->second;
		vector<bool> exonic_coords(gene_gff.end - gene_gff.start + 1);
		for (vector<const GFF*>::const_iterator exon_itr = exons.begin();
			 exon_itr != exons.end();
			 ++exon_itr)
		{
			const GFF* exon_gff = *exon_itr;
			assert (exon_gff->start >= gene_gff.start);
			for (uint32_t i = exon_gff->start;
				 i < exon_gff->end;
				 ++i)
			{
				if (i - 1 >= DoC.size())
					break;
				
				// Did we already count this one?
				if (!exonic_coords[i - gene_gff.start] && (!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
				{
					total_depth += DoC[i - 1];
				}
				exonic_coords[i - gene_gff.start] = true;
			}
		}
	}
	
	return total_depth;
}

uint32_t gene_depth(const GFF& gene_gff, 
						   const vector<const GFF*>& exons,
						   const vector<short>& DoC,
						   fsa::Sequence* ref_str)
{
	uint32_t gene_DoC = 0;
	vector<bool> exonic_coords(gene_gff.end - gene_gff.start + 1);
	for (vector<const GFF*>::const_iterator exon_itr = exons.begin();
		 exon_itr != exons.end();
		 ++exon_itr)
	{
		const GFF* exon_gff = *exon_itr;
		assert (exon_gff->start >= gene_gff.start);
		for (uint32_t i = exon_gff->start;
			 i < exon_gff->end;
			 ++i)
		{
			if (i - 1 >= DoC.size())
				break;
			// Did we already count this one?
			if (!exonic_coords[i - gene_gff.start] && (!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
			{
				gene_DoC += DoC[i - 1];
			}
			exonic_coords[i - gene_gff.start] = true;
		}
	}
	
	return gene_DoC;
}

uint32_t gene_length(const GFF& gene_gff, 
					 const vector<const GFF*>& exons,
					 fsa::Sequence* ref_str)
{
	uint32_t gene_exonic_length = 0;
	vector<bool> exonic_coords(gene_gff.end - gene_gff.start + 1);
	for (vector<const GFF*>::const_iterator exon_itr = exons.begin();
		 exon_itr != exons.end();
		 ++exon_itr)
	{
		const GFF* exon_gff = *exon_itr;
		assert (exon_gff->start >= gene_gff.start);
		for (uint32_t i = exon_gff->start;
			 i < exon_gff->end;
			 ++i)
		{
			// Did we already count this one?
			if (!exonic_coords[i - gene_gff.start] && 
				(!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
			{
				gene_exonic_length++;
			}
			exonic_coords[i - gene_gff.start] = true;
		}
	}
	
	return gene_exonic_length;
}

double gene_mend(const GFF& gene_gff, 
				 const vector<const GFF*>& exons,
				 const vector<short>& DoC,
				 uint32_t total_map_depth,
				 fsa::Sequence* ref_str)
{
	vector<bool> exonic_coords(gene_gff.end - gene_gff.start + 1);
	vector<double> gene_depths;
	for (vector<const GFF*>::const_iterator exon_itr = exons.begin();
		 exon_itr != exons.end();
		 ++exon_itr)
	{
		const GFF* exon_gff = *exon_itr;
		assert (exon_gff->start >= gene_gff.start);
		for (uint32_t i = exon_gff->start;
			 i < exon_gff->end;
			 ++i)
		{
			if (i - 1 >= DoC.size())
				break;
			// Did we already count this one?
			if (!exonic_coords[i - gene_gff.start] && 
				(!ref_str || !ref_str->is_pos_hardmasked(i - 1)))
			{
				gene_depths.push_back(DoC[i - 1]);
			}
			exonic_coords[i - gene_gff.start] = true;
		}
	}
	double gene_mend = 0.0;
	if (gene_depths.size())
	{
		sort(gene_depths.begin(), gene_depths.end());
		double median_DoC = gene_depths[gene_depths.size()/2];
		
		gene_mend = 1000000000.0 * (median_DoC / (double)total_map_depth);
	}
	
	return gene_mend;
}

double gene_rpkm(const GFF& gene_gff, 
				 const vector<const GFF*>& exons,
				 const vector<short>& DoC,
				 uint32_t total_map_depth,
				 fsa::Sequence* ref_str)
{
	uint32_t length = gene_length(gene_gff, exons, ref_str);
	uint32_t depth = gene_depth(gene_gff, exons, DoC, ref_str);
	double gene_avg_DoC = depth / (double) length;
	
	double gene_rpkm = 1000000000 * ( gene_avg_DoC / (double)total_map_depth );
	return gene_rpkm;
}

string gene_short_name(const GFF& gene_gff)
{
	GFF::AttributeTable::const_iterator att_itr;
	string gene_name;
	att_itr = gene_gff.attributes.find("Name");
	if (att_itr == gene_gff.attributes.end())
	{
		att_itr = gene_gff.attributes.find("ID");
		if (att_itr == gene_gff.attributes.end())
		{
			cerr << "Malformed gene record " << gene_gff << endl; 
			return "";
		}
		else
		{
			gene_name = att_itr->second.front();
		}
	}
	else
	{
		gene_name = att_itr->second.front();	
	}
	
	return gene_name;
}

void genes_for_ref(const GFF_database& gff_db, 
				   const string& ref_name,
				   GeneExonTable& gene_exons)
{
	
	GFF_database ref_features = gff_db.chromosome_features(ref_name); 
	
	typedef map<string, string> TransTable;
	
	GeneTable genes;
	TransTable transcripts;
	
	for(GFF_database::const_iterator gff_itr = ref_features.begin();
		gff_itr != ref_features.end();
		++gff_itr)
	{
		const GFF* gff_rec = new GFF(*gff_itr);
		if (gff_rec->type == "gene")
		{
			GFF::AttributeTable::const_iterator att_itr;
			att_itr = gff_rec->attributes.find("ID");
			if (att_itr == gff_rec->attributes.end() ||
				att_itr->second.size() != 1)
			{
				cerr << "Malformed gene record " << *gff_rec << endl; 
				continue;
			}
			const string& id = att_itr->second.front();
			genes.insert(make_pair(id, gff_rec));
			gene_exons.insert(make_pair(gff_rec, vector<const GFF*>()));
		}
		
	}
	
	for(GFF_database::const_iterator gff_itr = ref_features.begin();
		gff_itr != ref_features.end();
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
			string id = att_itr->second.front();
			
			att_itr = gff_rec.attributes.find("Parent");
			if (att_itr == gff_rec.attributes.end() ||
				att_itr->second.size() != 1)
			{
				cerr << "Malformed transcript record " << gff_rec << endl; 
				continue;
			}
			string gene = att_itr->second.front();
			
			transcripts.insert(make_pair(id, gene));
		}
		
	}
	
	for(GFF_database::const_iterator gff_itr = ref_features.begin();
		gff_itr != ref_features.end();
		++gff_itr)
	{
		const GFF* gff_rec = new GFF(*gff_itr);
		if (gff_rec->type == "exon")
		{
			GFF::AttributeTable::const_iterator att_itr;
			att_itr = gff_rec->attributes.find("Parent");
			if (att_itr == gff_rec->attributes.end())
			{
				cerr << "Malformed exon record " << gff_rec << endl; 
				continue;
			}
			vector<string> parent_transcripts = att_itr->second;
			for (vector<string>::iterator par_itr = parent_transcripts.begin();
				 par_itr != parent_transcripts.end();
				 ++par_itr)
			{
				// Get a valid transcript for this exon
				const string& transcript_str = *par_itr;
				TransTable::iterator transcript_itr = transcripts.find(transcript_str);
				if (transcript_itr == transcripts.end())
				{
					cerr << "No transcript with id " << transcript_str << endl;
					continue;
				}
				
				// Now find the gene GFF record for that transcript
				GeneTable::iterator gene_itr = genes.find(transcript_itr->second);
				if (gene_itr == genes.end())
				{
					cerr << "No gene associated with transcript " << transcript_str << endl;
					continue;
				}
				
				
				GeneExonTable::iterator gene_exon_itr = gene_exons.find(gene_itr->second);
				if (gene_exon_itr == gene_exons.end())
				{
					cerr << "No exons for " << gene_itr->second << endl;
					continue;
				}
				
				const GFF* gene_for_exon = gene_exon_itr->first;
				if (gene_for_exon->start <= gff_rec->start && gene_for_exon->end >= gff_rec->end)
				{
					gene_exon_itr->second.push_back(gff_rec);
				}
				else
				{
					cerr << "Exon falls outside of its enclosing gene" << endl;
					cerr << "Offending exon:" << endl << gff_rec;
					cerr << "In transcript: " << transcript_str << endl;
					cerr << "Enclosing gene:" << endl << *gene_for_exon;
				}
				break;
				
			}
		}
	}
}


void mend_for_genes(FILE* quant_expression_out,
					const vector<short>& DoC,
					const GeneExonTable& gene_exons,
					uint32_t map_depth,
					fsa::Sequence* ref_str)
{
	
	for (GeneExonTable::const_iterator gene_itr = gene_exons.begin(); 
		 gene_itr != gene_exons.end();
		 ++gene_itr)
	{
		//const string& gene_name = att_itr->second.front();
		string gene_name = gene_short_name(*(gene_itr->first));
		double mend = gene_mend(*(gene_itr->first), gene_itr->second, DoC, map_depth, ref_str);
		fprintf(quant_expression_out, "%s\t%lf\n", gene_name.c_str(), mend);
	}
}

void rpkm_for_genes(FILE* quant_expression_out,
					const vector<short>& DoC,
					const GeneExonTable& gene_exons,
					uint32_t map_depth,
					fsa::Sequence* ref_str)
{
	
	for (GeneExonTable::const_iterator gene_itr = gene_exons.begin(); 
		 gene_itr != gene_exons.end();
		 ++gene_itr)
	{
		//const string& gene_name = att_itr->second.front();
		double rpkm = gene_rpkm(*(gene_itr->first), gene_itr->second, DoC, map_depth, ref_str);
		string gene_name = gene_short_name(*(gene_itr->first));
		fprintf(quant_expression_out, "%s\t%lf\n", gene_name.c_str(), rpkm);
	}
}

struct GeneExpression
{
	GeneExpression(double R = 0.0, double M = 0.0) : rpkm(R), mend(M) {}
	double rpkm;
	double mend;
};

typedef map<string, GeneExpression> GETable;

void calculate_gene_expression(GETable& gene_expression,
								 const vector<short>& DoC,
								 const GeneExonTable& gene_exons,
								 uint32_t map_depth,
								 fsa::Sequence* ref_str)
{
	for (GeneExonTable::const_iterator gene_itr = gene_exons.begin(); 
		 gene_itr != gene_exons.end();
		 ++gene_itr)
	{
		string gene_name = gene_short_name(*(gene_itr->first));
		if (gene_name == "pxb")
			int a = 5;
		//const string& gene_name = att_itr->second.front();
		double rpkm = gene_rpkm(*(gene_itr->first), gene_itr->second, DoC, map_depth, ref_str);
		double mend = gene_mend(*(gene_itr->first), gene_itr->second, DoC, map_depth, ref_str);
		gene_expression[gene_name] = GeneExpression(rpkm, mend);
	}
}

void print_gene_expression(FILE* expr_out, const GETable& gene_expression)
{
	for (GETable::const_iterator itr = gene_expression.begin();
		 itr != gene_expression.end();
		 ++itr)
	{
		fprintf(expr_out, "%s\t%lf\t%lf\n", itr->first.c_str(), itr->second.rpkm, itr->second.mend);
	}
}

bool is_masked_char(char c)
{
	if (c == 'N' || c == 'n')
		return true;
	return false;
}

void driver(FILE* map1, 
			FILE* splice_map1,
			FILE* map2, 
			FILE* splice_map2,
			FILE* coverage_out,
			FILE* junctions_out,
			GFF_database gff_db,
			FILE* quant_expression_out,
			ifstream& ref_stream)
{
	typedef String<char> Reference;
	map<string, ifstream::pos_type> ref_file_offsets;
	
	gff_db.sort_entries();
	
	uint32_t total_map_depth = 0;
	
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;
		readMeta(ref_stream, name, Fasta());
		ifstream::pos_type offset = ref_stream.tellg();
		ref_file_offsets[name] = offset;
		read(ref_stream, ref_str, Fasta());
	}

	ref_stream.clear();
	ref_stream.seekg(0, ios::beg);
	
	bool paired_end = map2 && splice_map2;
	SequenceTable it(paired_end);
	SequenceTable rt(true);
    HitTable hits1(it,rt);
    HitTable hits2(it,rt);
    
    get_mapped_reads(map1, hits1, false);
	if (splice_map1)
		get_mapped_reads(splice_map1, hits1, true);
	
	JunctionSet junctions;
	
	print_wiggle_header(coverage_out);
	
	if (!paired_end)
	{
		fprintf(stderr, "Finished reading alignments\n");
		
		if (accept_all)
		{
			accept_all_hits(hits1);
		}
		else
		{
			BestFragmentAlignmentTable best_alignments(it.size());
			fragment_best_alignments(rt, hits1, best_alignments);
			
			accept_unique_hits(best_alignments);
		}
		junctions_from_alignments(hits1, junctions);
		
		for (SequenceTable::const_iterator ci = rt.begin();
			 ci != rt.end();
			 ++ci)
		{
			vector<short> DoC;
			fsa::Sequence* ref_seq = NULL;
			
			if (ref_stream.good())
			{	
				map<string,ifstream::pos_type>::iterator pi = ref_file_offsets.find(ci->first);
				if (pi != ref_file_offsets.end())
				{	
					ifstream::pos_type pos = pi->second;
					ref_stream.seekg(pos, ios::beg);
					Reference ref;
					read(ref_stream, ref, Fasta());
					ref_seq = new fsa::Sequence("", string(toCString(ref)));
					ref_seq->init_hardmasking(1, is_masked_char);
					//get_mask(string(toCString(ref_seq)), mask_ref);
				}
			}
			
			const HitList* h1 = hits1.get_hits(ci->second);
			if (h1)
				add_hits_to_coverage(*h1, DoC);
			
			if (filter_junctions)
				accept_valid_junctions(junctions, ci->second, DoC, min_isoform_fraction);
			else
				accept_all_junctions(junctions, ci->second);
			
			print_wiggle_for_ref(coverage_out, ci->first, DoC);
			
			GeneExonTable gene_exons;
			genes_for_ref(gff_db, ci->first, gene_exons);
			total_map_depth += total_exonic_depth(gene_exons,DoC, ref_seq);
			delete ref_seq;
		}
		
		ref_stream.clear();
		ref_stream.seekg(0, ios::beg);
		
		GETable gene_expression;
		
		for (SequenceTable::const_iterator ci = rt.begin();
			 ci != rt.end();
			 ++ci)
		{
			vector<short> DoC;
			fsa::Sequence* ref_seq = NULL;
			
			if (ref_stream.good())
			{	
				map<string,ifstream::pos_type>::iterator pi = ref_file_offsets.find(ci->first);
				if (pi != ref_file_offsets.end())
				{	
					ifstream::pos_type pos = pi->second;
					ref_stream.seekg(pos, ios::beg);
					Reference ref;
					read(ref_stream, ref, Fasta());
					ref_seq = new fsa::Sequence("", string(toCString(ref)));
					ref_seq->init_hardmasking(1, is_masked_char);
					//get_mask(string(toCString(ref_seq)), mask_ref);
				}
			}
			
			const HitList* h1 = hits1.get_hits(ci->second);
			if (h1)
				add_hits_to_coverage(*h1, DoC);
			
			if (filter_junctions)
				accept_valid_junctions(junctions, ci->second, DoC, min_isoform_fraction);
			else
				accept_all_junctions(junctions, ci->second);
			
			GeneExonTable gene_exons;
			genes_for_ref(gff_db, ci->first, gene_exons);
			
			if (quant_expression_out)
			{
				calculate_gene_expression(gene_expression, 
										  DoC, 
										  gene_exons, 
										  total_map_depth, 
										  ref_seq);
			}
			delete ref_seq;
		}
		if (quant_expression_out)
			print_gene_expression(quant_expression_out, gene_expression);
		print_junctions(junctions_out, junctions, rt);
	}
#ifdef PAIRED_END	
	else
	{
		get_mapped_reads(map2, hits2, false);
		get_mapped_reads(splice_map2, hits2, true);
		
		BestInsertAlignmentTable best_pairings(it.size());
		insert_best_pairings(rt, hits1, hits2, best_pairings);
		
		accept_valid_hits(best_pairings);
		
		junctions_from_alignments(hits1, junctions);
		junctions_from_alignments(hits2, junctions);
		
		print_junctions(junctions_out, junctions, rt);
		
		for (SequenceTable::const_iterator ci = rt.begin();
			 ci != rt.end();
			 ++ci)
		{
			vector<short> DoC;
			const HitList* h1 = hits1.get_hits(ci->second);
			if (h1)
				add_hits_to_coverage(*h1, DoC);
			
			const HitList* h2 = hits2.get_hits(ci->second);
			if (h2)
				add_hits_to_coverage(*h2, DoC);
			
			print_wiggle_for_ref(coverage_out, ci->first, DoC);
		}
    }
#endif
	
	uint32_t accepted_junctions = 0;
	for (JunctionSet::iterator itr = junctions.begin(); itr != junctions.end(); ++itr)
	{
		if(itr->second.accepted)
			accepted_junctions++;
	}
	
    fprintf(stderr, "Found %d junctions from happy spliced reads\n", accepted_junctions);
}

const char *short_options = "r:I:d:s:va:AF:G:o:R:";

#define USE_RPKM 260

static struct option long_options[] = {
{"insert-len",      required_argument,       0,            'I'},
{"insert-stddev",      required_argument,       0,            's'},
{"read-len",       required_argument,       0,            'r'},
{"max-dist",       required_argument,       0,            'd'},
{"min-anchor",       required_argument,       0,            'a'},
{"min-intron",       required_argument,       0,            'i'},
{"min-isoform-fraction",       required_argument,       0,            'F'},
{"verbose",		no_argument,	0,							'v'},
{"accept-all-hits",      no_argument,       0,            'A'},
{"gff-annotations",      no_argument,       0,            'G'},
{"ref-fasta",      no_argument,       0,            'R'},
{"output-dir",      no_argument,       0,            'o'},
{0, 0, 0, 0} // terminator
};

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
static float parseFloat(float lower, float upper, const char *errmsg) {
	float l;
	l = atof(optarg);
	
	if (l < lower) {
		cerr << errmsg << endl;
		print_usage();
		exit(1);
	}
	
	if (l > upper)
	{
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

int parse_options(int argc, char** argv)
{
	int option_index = 0;
	int next_option; 
	do { 
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);		
		switch (next_option) {
			case 'd':
	   			max_mate_inner_dist = (uint32_t)parseInt(0, "-d/--max-dist arg must be at least 0");
	   			break;
			case 'I':
	   			insert_len = (uint32_t)parseInt(1, "-I/--insert-len arg must be at least 1");
	   			break;
			case 's':
	   			insert_len_std_dev = (uint32_t)parseInt(1, "-s/--insert-stddev arg must be at least 1");
	   			break;
			case 'v':
				verbose = true;
				break;
			case 'a':
				min_anchor_len = (uint32_t)parseInt(4, "-a/--min-anchor arg must be at least 4");
				break;
			case 'i':
				min_intron_length = (uint32_t)parseInt(1, "-a/--min-intron arg must be at least 1");
				break;
			case 'F':
				min_isoform_fraction = parseFloat(0.0, 1.0, "-a/--min-isoform-fraction arg must be [0.0,1.0]");
				if (min_isoform_fraction == 0.0)
				{
					filter_junctions = false;
				}
				break;
			case 'A':
				accept_all = true;
				break;
			case 'G':
				gff_file = optarg;
				break;
			case 'R':
				ref_fasta = optarg;
				break;
			case 'o':
				output_dir = optarg;
				break;
			case -1: /* Done with options. */
				break;
			default: 
				print_usage();
				return 1;
		}
	} while(next_option != -1);
	
	return 0;
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
	if (parse_ret)
		return parse_ret;
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}
	
	string coverage_file_name = argv[optind++];
	
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
	
	string map1_file_name = argv[optind++];
	
	string splice_map1_file_name;
	string map2_file_name;
	string splice_map2_file_name;
	
	if (optind < argc) 
		splice_map1_file_name = argv[optind++];
	
	if (optind < argc)
		map2_file_name = argv[optind++];

	if (optind < argc) 
		splice_map2_file_name = argv[optind++];
		
	// Open the approppriate files

	FILE* coverage_file = fopen((output_dir + "/" + coverage_file_name).c_str(), "w");
	if (coverage_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for writing\n",
				coverage_file_name.c_str());
		exit(1);
	}
	
	FILE* junctions_file = fopen((output_dir + "/" + junctions_file_name).c_str(), "w");
	if (junctions_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for writing\n",
				junctions_file_name.c_str());
		exit(1);
	}
	
    FILE* map1_file = fopen(map1_file_name.c_str(), "r");
	if (map1_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map1_file_name.c_str());
		exit(1);
	}
	
	FILE* splice_map1_file = NULL;
	
	if (!splice_map1_file_name.empty())
	{
		splice_map1_file = fopen(splice_map1_file_name.c_str(), "r");
		if (splice_map1_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map1_file_name.c_str());
			exit(1);
		}
	}
	
	if ((splice_map2_file_name.empty() && !map2_file_name.empty()) ||
		(!splice_map2_file_name.empty() && map2_file_name.empty()))
	{
		fprintf(stderr, "Error: please specify both a contiguous and a spliced map file for paired-end reports\n");
		exit(1);
	}
	
	FILE* map2_file = NULL;
	FILE* splice_map2_file = NULL;
		
	if (!splice_map2_file_name.empty() && !map2_file_name.empty())
	{
		map2_file = fopen(map2_file_name.c_str(), "r");
		if (map2_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					map2_file_name.c_str());
			exit(1);
		}
		
		splice_map2_file = fopen(splice_map2_file_name.c_str(), "r");
		if (splice_map2_file == NULL)
		{
			fprintf(stderr, "Error: cannot open %s for reading\n",
					splice_map2_file_name.c_str());
			exit(1);
		}
	}
	
	GFF_database gff_db;
	FILE* quant_expression_out = NULL;
	if (gff_file != "")
	{
		gff_db.from_file(gff_file);
		string::size_type slash = gff_file.rfind("/");
		string::size_type dotGFF = gff_file.rfind(".gff");
		
		string expr_out_filename;
		if (slash == string::npos)
			slash = 0;
				
		expr_out_filename = gff_file.substr(slash, dotGFF) + ".expr";
		
		quant_expression_out = fopen((output_dir + "/" + expr_out_filename).c_str(), "w");
	}
	
	ifstream ref_stream(ref_fasta.c_str(), ifstream::in);;
	
    driver(map1_file, 
		   splice_map1_file, 
		   map2_file,
		   splice_map2_file, 
		   coverage_file, 
		   junctions_file, 
		   gff_db,
		   quant_expression_out,
		   ref_stream);
    
    return 0;
}


