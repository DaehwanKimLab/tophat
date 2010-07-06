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
#include "Glist.hh"
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
                        num_juncs_reported++;
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
                        num_juncs_reported++;
					}
				}
			}
		}
    }
    
    /*
     fprintf(stdout, "%s\t%d\t%d\t%c\n",
     five_prime_ex->seqid.c_str(),
     five_prime_ex->end - 1,
     three_prime_ex->start - 1,
     five_prime_ex->strand);
     */
    
    
//	// Table to hold the exons for each transcript, we'll make introns 
//	// from these below
//	
//	typedef map<string, const GFF*> TransTable;
//	typedef map<const GFF*, vector<const GFF*> > TransExonTable; 
//	TransTable transcripts;
//	TransExonTable transcript_exons;
//	
//    uint32_t num_juncs_reported = 0;
//	
//	for(GFF_database::const_iterator gff_itr = gff_db.begin();
//		gff_itr != gff_db.end();
//		++gff_itr)
//	{
//		const GFF& gff_rec = *gff_itr;
//		if (gff_rec.type == "mRNA")
//		{
//			GFF::AttributeTable::const_iterator att_itr;
//			att_itr = gff_rec.attributes.find("ID");
//			if (att_itr == gff_rec.attributes.end() ||
//				att_itr->second.size() != 1)
//			{
//				cerr << "Malformed transcript record " << gff_rec << endl; 
//				continue;
//			}
//			const string& id = att_itr->second.front();
//			transcripts.insert(make_pair(id, &gff_rec));
//			transcript_exons.insert(make_pair(&gff_rec, vector<const GFF*>()));
//		}
//
//	}
//	
//	for(GFF_database::const_iterator gff_itr = gff_db.begin();
//		gff_itr != gff_db.end();
//		++gff_itr)
//	{
//		const GFF& gff_rec = *gff_itr;
//		if (gff_rec.type == "exon")
//		{
//			GFF::AttributeTable::const_iterator att_itr;
//			att_itr = gff_rec.attributes.find("Parent");
//			if (att_itr == gff_rec.attributes.end())
//			{
//				cerr << "Malformed exon record " << gff_rec << endl; 
//				continue;
//			}
//			vector<string> parent_transcripts = att_itr->second;
//			for (vector<string>::iterator par_itr = parent_transcripts.begin();
//				 par_itr != parent_transcripts.end();
//				 ++par_itr)
//			{
//				const string& parent_str = *par_itr;
//				TransTable::iterator parent_rec = transcripts.find(parent_str);
//				if (parent_rec == transcripts.end())
//				{
//					cerr << "No transcript with id " << parent_str << endl;
//					continue;
//				}
//				TransExonTable::iterator parent_exons = transcript_exons.find(parent_rec->second);
//				if (parent_exons == transcript_exons.end())
//				{
//					cerr << "No exons with for " << *par_itr << endl;
//					continue;
//				}
//				parent_exons->second.push_back(&gff_rec);
//			}
//		}
//	}
//	
//	
//	for (TransExonTable::iterator trans_itr = transcript_exons.begin();
//		 trans_itr != transcript_exons.end();
//		 ++trans_itr)
//	{
//		vector<const GFF*>& exons = trans_itr->second;
//		if (exons.size() <= 1)
//			continue;
//		
//		vector<const GFF*>::iterator prev_itr = exons.begin();
//		vector<const GFF*>::iterator curr_itr = ++(exons.begin());
//		while (curr_itr != exons.end())
//		{
//			const GFF* five_prime_ex = *prev_itr;
//			const GFF* three_prime_ex = *curr_itr;
//			if (three_prime_ex->start < five_prime_ex->end)
//			{
//				fprintf(stderr, "Error: bad transcript annotation:\n");
//				cerr << *(trans_itr->first);
//				fprintf(stderr, "Offending exons overlapped:\n");
//				cerr << *five_prime_ex;
//				cerr << *three_prime_ex;
//				//exit(2);
//                break;
//			}
//			
//			fprintf(stdout, "%s\t%d\t%d\t%c\n",
//					five_prime_ex->seqid.c_str(),
//					five_prime_ex->end - 1,
//					three_prime_ex->start - 1,
//					five_prime_ex->strand);
//            ++num_juncs_reported;
//			++curr_itr;
//			++prev_itr;
//		}
//	}
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
	
	fprintf(stderr, "gff_juncs v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
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
