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
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif

#include <getopt.h>
#include <string>
#include <cstdio>
#include <set>
#include "GList.hh"
#include "gtf_tracking.h"

#include "common.h"
#include "bwt_map.h"

using namespace std;


void print_usage()
{
    fprintf(stderr, "Usage:   gtf_juncs <transcripts.gtf>\n");
}

uint32_t get_junctions_from_gff(FILE* ref_mRNA_file,
                                RefSequenceTable& rt)
{
    
    uint32_t num_juncs_reported = 0;
    
    GList<GSeqData> ref_rnas;
	
	if (ref_mRNA_file)
	{
		//read_mRNAs(ref_mRNA_file, false, ref_rnas, ref_rnas, NULL, -1, false);
		read_mRNAs(ref_mRNA_file, ref_rnas);
	}
    
    set<pair<string, pair<int, int> > > uniq_juncs;
    
    // Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			char* name = GffObj::names->gseqs.getName(ref_rnas[j]->gseq_id);
			uint32_t ref_id = rt.get_id(name, NULL, 0);
			for (int i = 0; i < ref_rnas[j]->mrnas_f.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_f[i]);
				
				
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
				
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
                        
                        fprintf(stdout, "%s\t%d\t%d\t%c\n",
                                name,
                                ex.end - 1,
                                next_ex.start - 1,
                                '+');
                        uniq_juncs.insert(make_pair(name, make_pair(ex.end - 1, next_ex.start - 1))); 
					}
				}
				
			}
			
			for (int i = 0; i < ref_rnas[j]->mrnas_r.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_r[i]);
				
				
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
					
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
						//ops.push_back(AugmentedCuffOp(CUFF_INTRON, ex.end, next_ex.start - ex.end - 1));
                        fprintf(stdout, "%s\t%d\t%d\t%c\n",
                                name,
                                ex.end - 1,
                                next_ex.start - 1,
                                '-');
                        uniq_juncs.insert(make_pair(name, make_pair(ex.end - 1, next_ex.start - 1)));
					}
				}
			}
		}
    }
    
    return uniq_juncs.size();
	
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
	
	string gtf_filename = argv[optind++];
	
	//GFF_database gff_db;
	if (gtf_filename == "")
	{
		print_usage();
		exit(2);
	}
    
    FILE* ref_gtf = fopen(gtf_filename.c_str(), "r");
    if (!ref_gtf)
    {
        fprintf (stderr, "Error: could not open GTF file %s for reading\n", gtf_filename.c_str());
        exit(1);
    }
	
	fprintf(stderr, "gtf_juncs v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
	fprintf(stderr, "---------------------------\n");
	
	//gff_db.from_file(gff_filename);
//    gff_db.sort_entries();
//	
    
    RefSequenceTable rt(true);
	uint32_t num_juncs_reported = get_junctions_from_gff(ref_gtf, rt);
    
    
    //uint32_t num_juncs_reported = 0;
    
    fprintf(stderr, "Extracted %u junctions from %s\n", 
            num_juncs_reported, gtf_filename.c_str());
    if (!num_juncs_reported)
        return 1;
    return 0;
}
