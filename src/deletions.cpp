/*
 *  deletions.cpp
 *  TopHat
 *
 *  Created by Ryan Kelley on 10/09/2010.
 *
 */
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include "common.h"
#include "deletions.h"





/*
 * Print deletions in BED format
 * As per the BED-standard (http://genome.ucsc.edu/FAQ/FAQformat)
 *	-The coordinates should be 0-based
 *	-The chromEnd field should not contain the actual feature
 *	-The name will be "-"
 *	-The score will be count of supporting reads (max of 1,000)
 *
 * chromStart refers to the position of the first deleted based 
 * <chrom>\t<left>\t<right>\t-\t<read count>\n
 * @param deletions_out The output file
 * @param deletions Maps from deletions to number of supporting reads
 * @pram ref_sequences The table of reference sequences
 *	
 */
void print_deletions(FILE* deletions_out, const DeletionSet& deletions, RefSequenceTable& ref_sequences){
	fprintf(deletions_out, "track name=deletions description=\"TopHat deletions\"\n");
	for(DeletionSet::const_iterator i = deletions.begin(); i != deletions.end(); ++i){
		fprintf(deletions_out, "%s\t%d\t%d\t%s\t%d\n",
			ref_sequences.get_name(i->first.refid),
			i->first.left + 1,
			i->first.right,
			"-",
			i->second);
	}
}

/**
 * Add deletions from an alignment to an DeletionSet.
 * This will look for deletion in the alignment specified by bh. If the 
 * deletion is already in deletions, it will updated the count. Otherwise,
 * it will add the deletion to the set and initialize the count to 1.
 * @param bh The bowtie hit to be used to specify alignment infromation.
 * @param deletions The DeletionSet that will be updated with the deletion information from teh alignment.
 */
void deletions_from_alignment(const BowtieHit& bh, DeletionSet& deletions){
	vector<Deletion> new_deletions;
	deletions_from_spliced_hit(bh, new_deletions);
	
	for(size_t i = 0; i < new_deletions.size(); ++i){
		Deletion deletion = new_deletions[i];
		DeletionSet::iterator itr = deletions.find(deletion);
		if (itr != deletions.end()){
			itr->second += 1;
		}
		else{
			deletions[deletion] = 1;
		}
	}
	return;
}




/**
 * Extract a list of deletions from a bowtie hit.
 * Given a bowtie hit, extract a vector of deletions.  
 * @param bh The bowtie hit to use for alignment information.
 * @param insertions Used to store the resultant vector of deletions.
 */
void deletions_from_spliced_hit(const BowtieHit& bh, vector<Deletion>& deletions){
  const vector<CigarOp>& cigar = bh.cigar();
  unsigned int positionInGenome = bh.left();
  unsigned int positionInRead = 0;

  bool bSawFusion = false;
  for(size_t c = 0; c < cigar.size(); ++c){
    Junction deletion;
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
    case dEL:
      if (bSawFusion)
	deletion.refid = bh.ref_id2();
      else
	deletion.refid = bh.ref_id();

      if (cigar[c].opcode == DEL)
	{
	  deletion.left = positionInGenome - 1;
	  deletion.right = positionInGenome + cigar[c].length;
	}
      else
	{
	  deletion.left = positionInGenome - cigar[c].length;
	  deletion.right = positionInGenome + 1;
	}
      
      deletions.push_back(deletion);
      positionInGenome += cigar[c].length;
      break;
    case INS:
    case iNS:
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

void merge_with(DeletionSet& deletions, const DeletionSet& other)
{
  for (DeletionSet::const_iterator deletion = other.begin(); deletion != other.end(); ++deletion)
    {
      DeletionSet::iterator itr = deletions.find(deletion->first);
      if (itr != deletions.end())
	{
	  itr->second += deletion->second;
	}
      else
	{
	  deletions[deletion->first] = deletion->second;
	}
    }
}
