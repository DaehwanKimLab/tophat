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
#include "tokenize.h"
#include "genes.h"
#include "FSA/gff.h"
#include "reads.h"
#include "FSA/sequence.h"

//#include "inserts.h"

using namespace std;
using namespace seqan;
using fsa::GFF;
using fsa::GFF_database;
using std::set;



static bool filter_junctions = true;
static float min_isoform_fraction = 0.15f;
static bool accept_all = true;
string gff_file = "";
string ref_fasta = "";
string gene_filter = "";
string output_dir = "tophat_out";
bool verbose = false;


void fragment_best_alignments(RefSequenceTable& rt,
							  ReadTable& it,
                              HitTable& hits1,
                              BestFragmentAlignmentTable& best_alignments)
{
    for(RefSequenceTable::const_iterator ci = rt.begin();
        ci != rt.end();
        ++ci)
    {

        // Tracks the number of singleton ALIGNMENTS, not the number of singleton
        // READS in each Bowtie map.

        string name = ci->second.name;
        uint64_t ref_id = ci->first;
        HitList* hits_in_ref = hits1.get_hits(ref_id);

        if (!hits_in_ref)
            continue;

        if (verbose)
            fprintf(stderr, "Looking for best alignments in %s\n", name.c_str());

        best_fragment_mappings(ref_id,
                               name,
                               *hits_in_ref,
							   it,
                               best_alignments);
    }
}

bool is_masked_char(char c)
{
    if (c == 'N' || c == 'n')
        return true;
    return false;
}

void load_hits(const vector<string>& filenames,
			   RefSequenceTable& rt,
			   ReadTable& it,
               HitTable& hits)
{
    for (size_t i = 0; i < filenames.size(); ++i)
    {
        const string& filename = filenames[i];
		HitFactory* hit_factory = NULL;
		
        if (filename.rfind(".sbwtout") != string::npos)
		{
			hit_factory = new SplicedBowtieHitFactory(it,rt, i == 0 || i == filenames.size() - 1);
		}
		else if (filename.rfind(".sam") != string::npos)
		{
			hit_factory = new SAMHitFactory(it,rt);
		}
		else
		{
			hit_factory = new BowtieHitFactory(it,rt);
		}
		
        FILE* map = fopen(filename.c_str(), "r");

        if (map == NULL)
        {
            fprintf(stderr, "Error: could not open %s\n", filename.c_str());
            exit(1);
        }
        fprintf(stderr, "Loading hits from %s\n", filename.c_str());
        size_t num_hits_before_load = hits.total_hits();
        get_mapped_reads(map, hits, *hit_factory, true);
        fprintf(stderr, "Loaded %d hits from %s\n", 
                (int)hits.total_hits() - (int)num_hits_before_load, 
                filename.c_str());
		delete hit_factory;
    }
}

void load_hits(RefSequenceTable& rt,
			   ReadTable& it,
			   const string& left_read_maplist,
               HitTable& left_hits,
               const string* right_read_maplist,
               HitTable* right_hits)
{
    vector<string> left_filenames;
    tokenize(left_read_maplist, ",", left_filenames);
    load_hits(left_filenames, rt, it, left_hits);
    if (right_read_maplist && right_hits)
    {
        vector<string> right_filenames;
        tokenize(*right_read_maplist, ",", right_filenames);
        load_hits(right_filenames, rt, it, *right_hits);
    }
}

void print_sam_for_hit(FILE* fout,
					   char* alt_name,
					   const char* bwt_buf, 
					   bool spliced)
{
	const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";
	static const int buf_size = 256;
	char orientation;
	char read_name[buf_size];
	int bwtf_ret = 0;
	//uint32_t seqid = 0;
	char text_name[buf_size];
	unsigned int text_offset;
	char sequence[buf_size];
	
	uint32_t sam_flag = 0;
	uint32_t sam_pos = 0;
	uint32_t map_quality = 255;
	char cigar[256];
	string mate_ref_name = "*";
	uint32_t mate_pos = 0;
	uint32_t insert_size = 0;
	char qualities[buf_size];
	unsigned int other_occs;
	char mismatches[buf_size];
	memset(mismatches, 0, sizeof(mismatches));
	// Get a new record from the tab-delimited Bowtie map
	bwtf_ret = sscanf(bwt_buf,
					  bwt_fmt_str,
					  read_name,
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
		return;
	}
	
	// Stripping the slash and number following it gives the insert name
	char* slash = strrchr(read_name, '/');
	if (slash)
	{
		*slash = 0;
	}
	int read_len = (int)strlen(sequence);
	
	// Add this alignment to the table of hits for this half of the
	// Bowtie map
	if (spliced)
	{
		// Parse the text_name field to recover the splice coords
		vector<string> toks;
		
		tokenize_strict(text_name, "|", toks);
		
		int num_extra_toks = (int)toks.size() - 6;
		
		if (num_extra_toks >= 0)
		{
			static const uint8_t left_window_edge_field = 1;
			static const uint8_t splice_field = 2;
			//static const uint8_t right_window_edge_field = 3;
			//static const uint8_t junction_type_field = 4;
			//static const uint8_t strand_field = 5;
			
			string contig = toks[0];
			for (int t = 1; t <= num_extra_toks; ++t)
			{
				contig += "|";
				contig += toks[t];
			}
			
			vector<string> splice_toks;
			tokenize(toks[num_extra_toks + splice_field], "-", splice_toks);

			
			uint32_t left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
			
			sam_pos = left + 1;
			
			uint32_t spliced_read_len = (uint32_t)strlen(sequence);
			uint32_t left_splice_coord = atoi(splice_toks[0].c_str());
			uint32_t right_splice_coord = atoi(splice_toks[1].c_str());
			int8_t left_splice_overhang = left_splice_coord - left + 1;
			int8_t right_splice_overhang = spliced_read_len - left_splice_overhang;
			
			//uint32_t right = atoi(splice_toks[1].c_str()) + right_splice_overhang;
			
			sprintf(cigar, 
					"%dM%dN%dM", 
					left_splice_overhang,
					right_splice_coord - left_splice_coord - 1,
					right_splice_overhang);
			
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
					if ((orientation == '+' && abs(mismatch_pos - left_splice_overhang) < 5) ||
						(orientation == '-' && abs(((int)spliced_read_len - left_splice_overhang + 1) - mismatch_pos)) < 5)
						mismatch_in_anchor = true;
				}
				//mismatch_toks.push_back(pch);
				pch = strtok (NULL, ",");
			}
			strcpy(text_name, contig.c_str());
			//strcpy(sequence, "*");
			//strcpy(qualities,"*");
		}
	}
	else
	{
		sam_pos = text_offset + 1;
		sprintf(cigar, "%dM", read_len);
	}
	
	if (orientation == '-')
		sam_flag |= 0x0010; // BAM_FREVERSE
	
	fprintf(fout,
			"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
			alt_name ? alt_name : read_name,
			sam_flag,
			text_name,
			sam_pos,
			map_quality,
			cigar,
			mate_ref_name.c_str(),
			mate_pos,
			insert_size,
			sequence,
			qualities);
}

void print_sam_for_accepted_hits(FILE* fout,
								 RefSequenceTable& rt,
								 ReadTable& it,
								 const vector<string>& filenames,
								 HitTable& hits,
								 FILE* reads_file)
{
	HitFactory* hit_factory = NULL;
    for (size_t i = 0; i < filenames.size(); ++i)
    {
		rewind(reads_file);
		
        const string& filename = filenames[i];
		bool sam = (filename.rfind(".sam") != string::npos);
		bool spliced_bowtie = (filename.rfind(".sbwtout") != string::npos); 
        if (spliced_bowtie)
		{
			hit_factory = new SplicedBowtieHitFactory(it, rt, i == 0 || i == filenames.size() - 1);
		}
		else if (sam)
		{
			hit_factory = new SAMHitFactory(it, rt);
		}
		else
		{
			hit_factory = new BowtieHitFactory(it, rt);
		}
		
        FILE* map = fopen(filename.c_str(), "r");
		
        if (map == NULL)
        {
            fprintf(stderr, "Error: could not open %s\n", filename.c_str());
            exit(1);
        }

        //fprintf(stderr, "Loading hits from %s\n", filename.c_str());
        //size_t num_hits_before_load = hits.total_hits();
		
		static int buf_size = 2048;
		char bwt_buf[buf_size];
		uint32_t reads_extracted = 0;
		
		char read_name[buf_size];
		char read_seq[buf_size];
		char read_alt_name[buf_size];
		char read_quals[buf_size];
		
		ReadID last_id = 0;
		while (fgets(bwt_buf, 2048, map))
		{
			// Chomp the newline
			char* nl = strrchr(bwt_buf, '\n');
			if (nl) *nl = 0;
			
			if (*bwt_buf == 0)
				continue;
			
			string clean_buf = bwt_buf;
			// Get a new record from the tab-delimited Bowtie map
			BowtieHit bh;
			
			if (hit_factory->get_hit_from_buf(bwt_buf, bh, true))
			{
				if (bh.insert_id() != last_id)
				{
					get_read_from_stream(bh.insert_id(), 
										 it, 
										 reads_file,
										 reads_format,
										 false,
										 read_name, 
										 read_seq,
										 read_alt_name,
										 read_quals);
					last_id = bh.insert_id();
				}
				// Only check uniqueness if these hits are spliced
				// hits.add_hit(bh, spliced);
				const HitList* hit_list = hits.get_hits(bh.ref_id());
				if (!hit_list ||
					!binary_search(hit_list->begin(), hit_list->end(), bh, hit_insert_id_lt))
					continue;
				//fprintf(fout, "%s\n", clean_buf.c_str());
				
				if (sam)
				{
					char* first_tab = strchr(bwt_buf, '\t');
					if (!first_tab)
						continue;
					
					fprintf(fout, "%s\t%s\n", read_alt_name, first_tab);
				}
				else
				{
					print_sam_for_hit(fout, read_alt_name, clean_buf.c_str(), spliced_bowtie);
				}
			}
			reads_extracted++;
		}
    }
}

void print_sam_for_accepted_hits(FILE* fout, 
								 RefSequenceTable& rt,
								 ReadTable& it,
								 const string& left_read_maplist,
								 HitTable& left_hits,
								 FILE* left_reads_file,
								 const string* right_read_maplist,
								 HitTable* right_hits,
								 FILE* right_reads_file)
{
    vector<string> left_filenames;
    tokenize(left_read_maplist, ",", left_filenames);

    print_sam_for_accepted_hits(fout, rt, it, left_filenames, left_hits, left_reads_file);

    if (right_read_maplist && right_hits)
    {
        vector<string> right_filenames;
        tokenize(*right_read_maplist, ",", right_filenames);

        print_sam_for_accepted_hits(fout, rt, it, right_filenames, *right_hits, right_reads_file);
    }
}

void driver(const string& left_maps,
			FILE* left_reads_file,
            const string* right_maps,
			FILE* right_reads_file,
            FILE* coverage_out,
            FILE* junctions_out,
			FILE* accepted_hits_out,
            GFF_database gff_db,
            FILE* quant_expression_out,
            ifstream& ref_stream,
            FILE* gene_filter_file)
{
    typedef String<char> Reference;
    map<string, fsa::Sequence*> masks;
    gff_db.sort_entries();

    uint32_t total_map_depth = 0;


	
//    while(ref_stream.good() &&
//          !ref_stream.eof())
//    {
//        Reference ref_str;
//        string name;
//        readMeta(ref_stream, name, Fasta());
////        ifstream::pos_type offset = ref_stream.tellg();
////        ref_file_offsets[name] = offset;
//        read(ref_stream, ref_str, Fasta());
//        fsa::Sequence* ref_seq = new fsa::Sequence("", string(toCString(ref_str)));
//        ref_seq->init_hardmasking(1, is_masked_char);
//        ref_seq->seq.clear();
//        masks[name] = ref_seq;
//    }
    
    std::set<string>* gene_filter = NULL;
    if (gene_filter_file)
    {
        gene_filter = new std::set<string>();
        char filter_buf[2048];
        while(!feof(gene_filter_file) &&
              fgets(filter_buf, 2048, gene_filter_file))
        {
            string short_name = fsa::Util::trim(filter_buf, " \t\n");
            gene_filter->insert(short_name);
        }
    }

    ref_stream.clear();
    ref_stream.seekg(0, ios::beg);

    // Load the set of left maps, and if provided, the set of right maps

    ReadTable it;
    RefSequenceTable rt(true);
	
	//HitFactory hit_factory(it,rt);
	
    HitTable left_hits;

    HitTable* right_hits = NULL;
    if (right_maps != NULL)
    {
        right_hits = new HitTable();
    }

    load_hits(rt, it, left_maps, left_hits, right_maps, right_hits);

    JunctionSet junctions;

    print_wiggle_header(coverage_out);

    GeneFactory gene_factory;
    map<string, Expression*> gene_expression;

	fprintf(stderr, "Finished reading alignments\n");

	if (accept_all)
	{
		accept_all_hits(left_hits);
	}
	else
	{
		BestFragmentAlignmentTable best_alignments(it.size());
		fragment_best_alignments(rt, it, left_hits, best_alignments);
	
		accept_unique_hits(best_alignments);
	}
	junctions_from_alignments(left_hits, junctions);
	
	
	fprintf (stderr, "Building coverage wiggles...");
	
	for (RefSequenceTable::const_iterator ci = rt.begin();
		 ci != rt.end();
		 ++ci)
	{
		vector<short> DoC;
		fsa::Sequence* ref_seq = NULL;
		map<string, fsa::Sequence*>::iterator seq_itr = masks.find(ci->second.name);
		if (seq_itr != masks.end())
			ref_seq = seq_itr->second;
		
		const HitList* h1 = left_hits.get_hits(ci->first);
		if (h1)
			add_hits_to_coverage(*h1, DoC);
		
		if (filter_junctions)
			accept_valid_junctions(junctions, ci->first, DoC, min_isoform_fraction);
		else
			accept_all_junctions(junctions, ci->first);

		print_wiggle_for_ref(coverage_out, ci->second.name, DoC);

		//GeneExonTable gene_exons;
		GeneTable ref_genes;
		gene_factory.get_genes(ci->second.name, gff_db, ref_genes, gene_filter);
		total_map_depth += total_exonic_depth(ref_genes,DoC, ref_seq);
	}

	fprintf (stderr, "done\n");
	fprintf (stderr, "Calculating RPKMs...");
	for (RefSequenceTable::const_iterator ci = rt.begin();
		 ci != rt.end();
		 ++ci)
	{
		vector<short> DoC;
		fsa::Sequence* ref_seq = NULL;


		map<string, fsa::Sequence*>::iterator seq_itr = masks.find(ci->second.name);
		if (seq_itr != masks.end())
			ref_seq = seq_itr->second;

		const HitList* h1 = left_hits.get_hits(ci->first);
		if (h1)
			add_hits_to_coverage(*h1, DoC);

		if (quant_expression_out)
		{
			GeneTable ref_genes;
			gene_factory.get_genes(ci->second.name, gff_db, ref_genes, gene_filter);
			calculate_gene_expression(ref_genes,
									  DoC,
									  ref_seq,
									  total_map_depth,
									  gene_expression);
		}
		delete ref_seq;
	}
	
	
	if (quant_expression_out)
		print_gene_expression(quant_expression_out, gene_expression);
	
	fprintf (stderr, "done\n");
	
	fprintf (stderr, "Printing junction BED track...");
	print_junctions(junctions_out, junctions, rt);
	fprintf (stderr, "done\n");

    uint32_t accepted_junctions = 0;
    for (JunctionSet::iterator itr = junctions.begin(); itr != junctions.end(); ++itr)
    {
        if(itr->second.accepted)
		{
            accepted_junctions++;
		}
		else
		{
			JunctionStats& s = itr->second;
			for (std::set<BowtieHit*>::iterator hi = s.supporting_hits.begin();
				 hi != s.supporting_hits.end();
				 ++hi)
			{
				(*hi)->accepted(false);
			}
		}
    }

	fprintf (stderr, "Reporting final accepted alignments...");
	print_sam_for_accepted_hits(accepted_hits_out, 
								rt, 
								it,
								left_maps, 
								left_hits,
								left_reads_file,
								right_maps, 
								right_hits,
								right_reads_file);
	
	fprintf (stderr, "done\n");
    fprintf(stderr, "Found %d junctions from happy spliced reads\n", accepted_junctions);
}

const char *short_options = "r:I:d:s:va:AF:G:o:R:f:";

#define USE_RPKM 260

static struct option long_options[] = {
    {"insert-stddev",      required_argument,       0,            's'},
    {"insert-len",       required_argument,       0,            'r'},
    {"max-dist",       required_argument,       0,            'd'},
    {"min-anchor",       required_argument,       0,            'a'},
    {"min-intron",       required_argument,       0,            'i'},
    {"min-isoform-fraction",       required_argument,       0,            'F'},
    {"verbose",     no_argument,    0,                          'v'},
    {"accept-all-hits",      no_argument,       0,            'A'},
    {"gff-annotations",      no_argument,       0,            'G'},
    {"ref-fasta",      no_argument,       0,            'R'},
    {"output-dir",      no_argument,       0,            'o'},
    {"gene-filter",      no_argument,       0,            'f'},
    {0, 0, 0, 0} // terminator
};


void print_usage()
{
    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <left_map1,...,left_mapN> <left_reads.fq>\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}


/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static float parseFloat(float lower, float upper, const char *errmsg) {
    float l;
    l = (float)atof(optarg);

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
            max_mate_inner_dist = (uint32_t)parseInt(0, "-d/--max-dist arg must be at least 0", print_usage);
            break;
        case 'r':
            insert_len = (uint32_t)parseInt(1, "-I/--insert-len arg must be at least 1", print_usage);
            break;
        case 's':
            insert_len_std_dev = (uint32_t)parseInt(1, "-s/--insert-stddev arg must be at least 1", print_usage);
            break;
        case 'v':
            verbose = true;
            break;
        case 'a':
            min_anchor_len = (uint32_t)parseInt(4, "-a/--min-anchor arg must be at least 4", print_usage);
            break;
        case 'i':
            min_intron_length = (uint32_t)parseInt(1, "-a/--min-intron arg must be at least 1", print_usage);
            break;
        case 'F':
            min_isoform_fraction = parseFloat(0.0f, 1.0f, "-a/--min-isoform-fraction arg must be [0.0,1.0]");
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
        case 'f':
            gene_filter = optarg;
            break;
        case -1:     /* Done with options. */
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
	reads_format = FASTQ;
	
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
	
	string accepted_hits_file_name = argv[optind++];
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }

    string left_maps = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string left_reads_filename = argv[optind++];
	
    string* right_maps = NULL;
	string* right_reads_filename = NULL;
	FILE* right_reads_file = NULL;

    // Open the approppriate files

    FILE* coverage_file = fopen((output_dir + "/" + coverage_file_name).c_str(), "w");
    if (coverage_file == NULL)
    {
        fprintf(stderr, "Error: cannot open WIG file %s for writing\n",
                coverage_file_name.c_str());
        exit(1);
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
        if (!quant_expression_out)
        {
            fprintf(stderr, "Error: cannot open %s for writing\n",
                    expr_out_filename.c_str());
            exit(1);
        }
    }
    
    FILE* gene_filter_file = NULL;
    if (gene_filter != "")
    {
        gene_filter_file = fopen(gene_filter.c_str(), "r");
        
        if (!gene_filter_file)
        {
            fprintf(stderr, "Error: cannot open %s for reading\n",
                    gene_filter.c_str());
            exit(1);
        }
    }

    ifstream ref_stream(ref_fasta.c_str(), ifstream::in);;

    driver(left_maps,
		   left_reads_file,
           right_maps,
		   right_reads_file,
           coverage_file,
           junctions_file,
		   accepted_hits_file,
           gff_db,
           quant_expression_out,
           ref_stream,
           gene_filter_file);

    return 0;
}


