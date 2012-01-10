/*
 *  insertions.cpp
 *  TopHat
 *
 *  Created by Ryan Kelley on 11/04/2010.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#define SVN_REVISION "XXX"
#endif


#include <cassert>
#include <cstdio>
#include <cstring>
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
#include "insertions.h"
#include "fragments.h"
#include "wiggles.h"
#include "tokenize.h"
#include "reads.h"

#include "inserts.h"

/**
 * Add insertions from an alignment to an InsertionSet.
 * This will look for insertion in the alignment specified by bh. If the 
 * insertion is already in insertions, it will updated the count. Otherwise,
 * it will add the insertion to the set and initialize the count to 1.
 * @param bh The bowtie hit to be used to specify alignment infromation.
 * @param insertions The InsertionSet that will be updated with the insertion information from teh alignment.
 */
void insertions_from_alignment(const BowtieHit& bh, InsertionSet& insertions){

	vector<Insertion> new_insertions;
	insertions_from_spliced_hit(bh, new_insertions);

	for(size_t i = 0; i < new_insertions.size(); ++i){
		Insertion insertion = new_insertions[i];
		InsertionSet::iterator itr = insertions.find(insertion);
		if (itr != insertions.end()){
			itr->second += 1;
		}
		else{
			assert(insertion.refid != VMAXINT32);
			insertions[insertion] = 1;
		}
	}
	return;
}

/**
 * Print insertions in BED format.
 * Note, as per the BED-standard (http://genome.ucsc.edu/FAQ/FAQformat)
 *   -The coordinates should be 0-based
 *   -The chromEnd field should not include the actual feature
 *   -The name will be the inserted sequence
 *   -The score will be the number of supporing counts, which is capped at 1,000
 * By (my) convention, the chromStart will be the last genome postion
 * before hte insertio.
 * 
 * <chrom>\t<left>\t<left>\t<inserted sequence>\t<read count>\n
 * @param insertions_out The output file
 * @param insertions Maps from insertions to number of supporting reads
 * @param ref_sequences The table of reference sequences 
 */
void print_insertions(FILE* insertions_out, const InsertionSet& insertions, RefSequenceTable& ref_sequences){
	fprintf(insertions_out, "track name=insertions description=\"TopHat insertions\"\n");
	for(InsertionSet::const_iterator i = insertions.begin(); i != insertions.end(); ++i){
		int counts = i->second;
		if(counts > 1000){
			counts = 1000;
		}
		fprintf(insertions_out, "%s\t%d\t%d\t%s\t%d\n",
			ref_sequences.get_name(i->first.refid),
			i->first.left,
			i->first.left,
			(i->first.sequence).c_str(),
			counts);
	}
}

/**
 * Extract a list of insertions from a bowtie hit.
 * Given a bowtie hit, extract a vector of insertions.  
 * @param bh The bowtie hit to use for alignment information.
 * @param insertions Used to store the resultant vector of insertions.
 */
void insertions_from_spliced_hit(const BowtieHit& bh, vector<Insertion>& insertions){
	const vector<CigarOp>& cigar = bh.cigar();
	unsigned int positionInGenome = bh.left();
	unsigned int positionInRead = 0;

	bool bSawFusion = false;
	for(size_t c = 0; c < cigar.size(); ++c){
	  Insertion insertion;
	  switch(cigar[c].opcode){
	  case REF_SKIP:
	    positionInGenome += cigar[c].length;
	    break;
	  case rEF_SKIP:
	    positionInGenome -= cigar[c].length;
	    break;
	  case MATCH:
	  case mATCH:
	    if (cigar[c].opcode == MATCH)
	      positionInGenome += cigar[c].length;
	    else
	      positionInGenome -= cigar[c].length;
	    positionInRead += cigar[c].length;
	    break;
	  case DEL:
	    positionInGenome += cigar[c].length;
	    break;
	  case dEL:
	    positionInGenome -= cigar[c].length;
	    break;
	  case INS:
	  case iNS:
	    /*
	     * Note that the reported position in the genome from the SAM
	     * alignment is 1-based, since the insertion object is expecting
	     * a 0-based co-ordinate, we need to subtract 1
	     */
	    if (bSawFusion)
	      insertion.refid = bh.ref_id2();
	    else
	      insertion.refid = bh.ref_id();

	    if (cigar[c].opcode == INS)
	      insertion.left = positionInGenome - 1;
	    else
	      insertion.left = positionInGenome + 1;
	      
	    insertion.sequence = bh.seq().substr(positionInRead, cigar[c].length);
	    
	    insertions.push_back(insertion);
	    positionInRead += cigar[c].length;
	    break;
	  case FUSION_FF:
	  case FUSION_FR:
	  case FUSION_RF:
	    bSawFusion = true;
	    positionInGenome = cigar[c].length;
	    break;
	  default:
	    break;
		}	
	}	
	return;
} 

void merge_with(InsertionSet& insertions, const InsertionSet& other)
{
  for (InsertionSet::const_iterator insertion = other.begin(); insertion != other.end(); ++insertion)
    {
      InsertionSet::iterator itr = insertions.find(insertion->first);
      if (itr != insertions.end())
	{
	  itr->second += insertion->second;
	}
      else
	{
	  insertions[insertion->first] = insertion->second;
	}
    }
}
