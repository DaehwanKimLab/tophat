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
#include "gff.h"

#include "common.h"
#include "bwt_map.h"

using namespace std;


void print_usage()
{
    fprintf(stderr, "Usage:   gtf_juncs <transcripts.gtf>\n");
}

void read_transcripts(FILE* f, GffReader& gffr) { 
  //assume gffr was just created but not initialized
  gffr.init(f, true, true); //(gffile, mRNA-only, sortByLoc)
  gffr.showWarnings(verbose);
  //(keepAttr,   mergeCloseExons,  noExonAttr)
  gffr.readAll(false, true, true); 
  //now all parsed GffObjs are in gffr.gflst, grouped by genomic sequence
  }

uint32_t get_junctions_from_gff(FILE* ref_mRNA_file,
                                RefSequenceTable& rt)
{
	GffReader gff_reader(ref_mRNA_file, true); //only recognizable transcript features, sort them by locus
	if (ref_mRNA_file)
	{
		read_transcripts(ref_mRNA_file, gff_reader);
	}
	
	set<pair<string, pair<int, int> > > uniq_juncs;
	
	//if any ref data was loaded
	int last_gseqid=-1;
	const char* gseqname=NULL;
	for (int i=0;i<gff_reader.gflst.Count();i++) {
		//ref data is grouped by genomic sequence
		GffObj& rna = *(gff_reader.gflst[i]);
		uint tlen=rna.len();
		if (rna.hasErrors() || (tlen+500>GFF_MAX_LOCUS)) { 
			//if (verbose) 
			GMessage("Warning: transcript %s discarded (structural errors found, length=%d).\n", rna.getID(), tlen);
			continue;
			}
		if (rna.isDiscarded()) {
			//discarded generic "gene" or "locus" features with no other detailed subfeatures
			continue;
			}
		if (rna.exons.Count()==0) {
				//if (verbose)
				// GMessage("Warning: %s %s found without exon segments (adding default exon).\n",rna.getFeatureName(), rna.getID());
				rna.addExon(rna.start,rna.end);
				}

		if (rna.gseq_id!=last_gseqid) {
		    gseqname=rna.getGSeqName();
		    rt.get_id(gseqname, NULL, 0);
		    last_gseqid=rna.gseq_id;
		    }
		for (int e = 1; e < rna.exons.Count(); ++e) {
		    GffExon& ex = *(rna.exons[e]);
		    GffExon& prex = *(rna.exons[e-1]);
		    fprintf(stdout, "%s\t%d\t%d\t%c\n",
		                 gseqname,
		                 prex.end-1, ex.start-1, rna.strand);
		    uniq_juncs.insert(make_pair(gseqname, make_pair(prex.end - 1, ex.start - 1)));
		    }
		} //for each loaded GFF record
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
