/*
 *  fusions.cpp
 *  TopHat
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
#include <seqan/modifier.h>
#include <seqan/align.h>
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "junctions.h"
#include "fusions.h"
#include "fragments.h"
#include "wiggles.h"
#include "tokenize.h"
#include "reads.h"

#include "inserts.h"

// daehwan - replace this with that of SeqAn
int difference(const string& first, const string& second)
{
  int len = seqan::length(first);
  if (len != (int)seqan::length(second))
    return 0;

  int min_value = 10000;
  short value1[1024] = {0,};
  short value2[1024] = {0,};

  short* curr = value1;
  short* prev = value2;
  
  for (int j = 0; j < len; ++j)
    {
      for (int i = 0; i < len; ++i)
	{
	  int value = 10000;
	  int match = first[i] == second[j] ? 0 : 1;

	  // right
	  if (i == 0)
	    value = j * 2 + match;
	  else if (j > 0)
	    value = prev[i] + 2;

	  int temp_value = 10000;

	  // down
	  if (j == 0)
	    temp_value = i * 2 + match;
	  else if (i > 0)
	    temp_value = curr[i-1] + 2;

	  if (temp_value < value)
	    value = temp_value;

	  // match
	  if (i > 0 && j > 0)
	    temp_value = prev[i-1] + match;

	  if (temp_value < value)
	    value = temp_value;

	  curr[i] = value;

	  if ((i == len - 1 || j == len - 1) && value < min_value)
	    min_value = value;
	}

      short* temp = prev;
      prev = curr;
      curr = temp;
   }

  return min_value;
}


/**
 * Add fusions from an alignment to an FusionSet.
 * This will look for fusion in the alignment specified by bh.
 * @param bh The bowtie hit to be used to specify alignment infromation.
 * @param fusions The FusionSet that will be updated with the insertion information from teh alignment.
 */
void fusions_from_alignment(const BowtieHit& bh,
			    FusionSet& fusions,
			    RefSequenceTable& rt,
			    bool update_stat)
{
  RefSequenceTable::Sequence* ref_str1 = rt.get_seq(bh.ref_id());
  RefSequenceTable::Sequence* ref_str2 = rt.get_seq(bh.ref_id2());

  if (!ref_str1 || !ref_str2)
    return;

  vector<Fusion> new_fusions;
  fusions_from_spliced_hit(bh, new_fusions);

  for(size_t i = 0; i < new_fusions.size(); ++i)
    {
      Fusion fusion = new_fusions[i];
      const vector<CigarOp>& cigars = bh.cigar();

      /*
       * Assume read is in the same direction as fusion.
       */
     
      // Find the position of Fusion.
      size_t fusion_pos = 0;
      for (; fusion_pos < cigars.size(); ++fusion_pos)
	{
	  CigarOpCode opcode = cigars[fusion_pos].opcode;
	  if (opcode == FUSION_FF || opcode == FUSION_FR || opcode == FUSION_RF || opcode == FUSION_RR)
	    break;
	}

      if (fusion_pos <= 0 || fusion_pos + 1 >= cigars.size())
	continue;

      // For left bases,
      size_t left_pos = 0;
      for (int j = (int)fusion_pos - 1; j >= 0; --j)
	{
	  const CigarOp& cigar = cigars[j];
	  switch (cigar.opcode)
	    {
	    case MATCH:
	    case mATCH:
	    case REF_SKIP:
	    case rEF_SKIP:
	    case DEL:
	    case dEL:
	      left_pos += cigar.length;
	      break;
	    default:
	      break;
	    }
	}

      // For right bases,
      size_t right_pos = 0;
      for (size_t j = fusion_pos + 1; j < cigars.size(); ++j)
	{
	  const CigarOp& cigar = cigars[j];
	  switch (cigar.opcode)
	    {
	    case MATCH: 
	    case mATCH:
	    case REF_SKIP:
	    case rEF_SKIP:
	    case DEL:
	    case dEL:
	      right_pos += cigar.length;
	      break;
	    default:
	      break;
	    }
	}

      if (left_pos < fusion_anchor_length || right_pos < fusion_anchor_length)
	continue;

      FusionSet::iterator itr = fusions.find(fusion);
      if (itr != fusions.end())
	{
	  itr->second.count += 1;
	}
      else
	{
	  assert(fusion.refid1 != 0xFFFFFFFF);
	  FusionStat fusionStat;
	  fusionStat.count = 1;
	  fusions[fusion] = fusionStat;
	  
	  itr = fusions.find(fusion);

	  if (!update_stat)
	    {
	      /*
	       * make a reversed fusion.
	       * this is necessary to detect reads that contradict the fusion.
	       */

	      FusionStat fusionStat_rev;
	      fusionStat_rev.count = 1;
	      fusionStat_rev.reversed = true;

	      Fusion fusion_rev(fusion.refid2, fusion.refid1, fusion.right, fusion.left, fusion.dir);
	      fusions[fusion_rev] = fusionStat_rev;
	    }
	}
      
      if (update_stat)
	{
	  if (itr->second.chr1_seq.length() <= 0)
	    {
	      size_t len = 100;
	      size_t half_len = len / 2;
	      size_t increase = 20;

	      if (fusion.left >= half_len && fusion.left + half_len <= seqan::length(*ref_str1) &&
		  fusion.right >= half_len && fusion.right + half_len <= seqan::length(*ref_str2))
		{
		  seqan::Dna5String left, right;

		  if (fusion.dir == FUSION_RF || fusion.dir == FUSION_RR)
		    {
		      left = seqan::infix(*ref_str1, fusion.left - half_len, fusion.left + half_len);
		      seqan::reverseComplement(left);
		    }
		  else
		    left = seqan::infix(*ref_str1, fusion.left - half_len + 1, fusion.left + half_len + 1);
		  
		  if (fusion.dir == FUSION_FR || fusion.dir == FUSION_RR)
		    {
		      right = seqan::infix(*ref_str2, fusion.right - half_len + 1, fusion.right + half_len + 1);
		      seqan::reverseComplement(right);
		    }
		  else
		    right = seqan::infix(*ref_str2, fusion.right - half_len, fusion.right + half_len);
		  
		  itr->second.chr1_seq = DnaString_to_string(left);
		  itr->second.chr2_seq = DnaString_to_string(right);

		  for (size_t j = 0; j < 5; ++j)
		    {
		      size_t pos = (5 - j - 1) * increase / 2;
		      const string& left_sub = itr->second.chr1_seq.substr(pos, (j+1) * increase);
		      const string& right_sub = itr->second.chr2_seq.substr(pos, (j+1) * increase);

		      itr->second.diffs.push_back(difference(left_sub, right_sub));
		    }
		}
	    }
		  
	  assert (bh.ref_id() == itr->first.refid1);

	  itr->second.left_ext = max((size_t)itr->second.left_ext, left_pos);
	  itr->second.right_ext = max((size_t)itr->second.right_ext, right_pos);
	  
	  for (size_t k = 0; k < left_pos && k < itr->second.left_bases.size(); ++k)
	    {
	      ++(itr->second.left_bases[k]);
	    }
	  
	  for (size_t k = 0; k < right_pos && k < itr->second.right_bases.size(); ++k)
	    {
	      ++(itr->second.right_bases[k]);
	    }
	}
    }
}

void unsupport_fusions(const BowtieHit& bh, FusionSet& fusions, const FusionSet& fusions_ref)
{
  if (bh.fusion_opcode() != FUSION_NOTHING || bh.is_spliced() || bh.read_len() < 40)
    return;

  FusionSet::const_iterator lb, ub;
  
  uint32_t left = bh.left() + 20;
  uint32_t right = bh.right() - 20;

  lb = fusions_ref.upper_bound(Fusion(0u, 0u, left, 0));
  ub = fusions_ref.lower_bound(Fusion(0xffffffffu, 0xffffffffu, right, 0xffffffffu));
  while (lb != ub && lb != fusions_ref.end())
    {
      if (lb->first.refid1 == bh.ref_id())
	{
	  // daehwan
#if 0
	  // MCF-7	RPS6KB1	17:57970443-58027925:1	TMEM49	17:57784863-57917950:1
	  if ((lb->first.left == 57992061 && lb->first.right == 57917126) ||
	      (lb->first.left == 57917126 && lb->first.right == 57992061))
	    {
	      const char* dir_str = "ff";
	      if (lb->first.dir == FUSION_FR)
		dir_str = "fr";
	      else if (lb->first.dir == FUSION_RF)
		dir_str = "rf";

	      cout << "fusion: " << lb->first.left << "-" << lb->first.right << endl;
	      cout << dir_str << endl;
	      cout << bh.insert_id() << ": " << bh.left() << "-" << bh.right() << "\t";
	      cout << print_cigar(bh.cigar()) << " " << (int)bh.edit_dist() << endl;
	      cout << bh.seq() << endl;
	      cout << bh.qual() << endl;
	      cout << endl;
	    }
#endif
	  
	  FusionSet::iterator itr;
	  if (lb->second.reversed)
	    itr = fusions.find(Fusion(lb->first.refid2, lb->first.refid1, lb->first.right, lb->first.left, lb->first.dir));
	  else
	    itr = fusions.find(lb->first);

	  if (itr == fusions.end())
	    {
	      FusionStat fusionStat;
	      fusionStat.unsupport_count = 1;
	      fusions[lb->first] = fusionStat;
	    }
	  else
	    ++(itr->second.unsupport_count);
	}

      ++lb;
    }
}

/**
 */
void print_fusions(FILE* fusions_out, FusionSet& fusions, RefSequenceTable& ref_sequences)
{
  // fprintf(fusions_out, "track name=fusions description=\"TopHat fusions\"\n");

  vector<Fusion> vFusion;
  for (FusionSet::iterator itr = fusions.begin(); itr != fusions.end(); ++itr)
    {
      vFusion.push_back(itr->first);
    }
  sort(vFusion.begin(), vFusion.end());
  
  for (size_t i = 0; i < vFusion.size(); ++i)
    {
      FusionSet::iterator itr = fusions.find(vFusion[i]);
      int counts = itr->second.count;
      if (counts <= 0 || itr->second.reversed)
	continue;
      
      const char* dir = "";
      if (itr->first.dir == FUSION_FF)
	dir = "ff";
      else if(itr->first.dir == FUSION_FR)
	dir = "fr";
      else if(itr->first.dir == FUSION_RF)
	dir = "rf";
      else
	dir = "rr";

      assert (itr->second.left_bases.size() == itr->second.right_bases.size());

      float symm = 0.0f;
      for (uint32_t i = 0; i < itr->second.left_bases.size(); ++i)
	{
	  float term = ((int)itr->second.left_bases[i] - (int)itr->second.right_bases[i]) / (float)counts;
	  symm += (term * term);
	}

      fprintf(fusions_out, "%s-%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f",
	      ref_sequences.get_name(itr->first.refid1),
	      ref_sequences.get_name(itr->first.refid2),
	      itr->first.left,
	      itr->first.right,
	      dir,
	      counts,
	      itr->second.pair_count,
	      itr->second.pair_count_fusion,
	      itr->second.unsupport_count,
	      itr->second.left_ext,
	      itr->second.right_ext,
	      symm);

      fprintf(fusions_out, "\t@\t");

      for (uint32_t i = 0; i < itr->second.diffs.size(); ++i)
	{
	  fprintf(fusions_out, "%d ", itr->second.diffs[i]);
	}

      fprintf(fusions_out, "\t@\t");

      uint32_t half_length = itr->second.chr1_seq.length() / 2;
      fprintf(fusions_out, "%s %s\t@\t", itr->second.chr1_seq.substr(0, half_length).c_str(), itr->second.chr1_seq.substr(half_length).c_str());
      fprintf(fusions_out, "%s %s\t@\t", itr->second.chr2_seq.substr(0, half_length).c_str(), itr->second.chr2_seq.substr(half_length).c_str());

      for (uint32_t i = 0; i < itr->second.left_bases.size(); ++i)
	{
	  fprintf(fusions_out, "%d ", itr->second.left_bases[i]);
	}

      fprintf(fusions_out, "\t@\t");

      for (uint32_t i = 0; i < itr->second.right_bases.size(); ++i)
	{
	  fprintf(fusions_out, "%d ", itr->second.right_bases[i]);
	}

      fprintf(fusions_out, "\t@\t");

      sort(itr->second.vPairSupport.begin(), itr->second.vPairSupport.end());
      for (uint32_t i = 0; i < min((size_t)200, itr->second.vPairSupport.size()); ++i)
	{
	  fprintf(fusions_out, "%d:%d ", itr->second.vPairSupport[i].ldist, itr->second.vPairSupport[i].rdist);
	}

      fprintf(fusions_out, "\n");
    }
}

/**
 * Extract a list of fusions from a bowtie hit.
 * Given a bowtie hit, extract a vector of insertions.  
 * @param bh The bowtie hit to use for alignment information.
 * @param insertions Used to store the resultant vector of insertions.
 */
void fusions_from_spliced_hit(const BowtieHit& bh, vector<Fusion>& fusions, bool auto_sort)
{
  const vector<CigarOp>& cigar = bh.cigar();
  unsigned int positionInGenome = bh.left();

  for(size_t c = 0; c < cigar.size(); ++c)
    {
      Fusion fusion;
      switch(cigar[c].opcode)
	{
	case REF_SKIP:
	case MATCH:
	case DEL:
	  positionInGenome += cigar[c].length;
	  break;
	case rEF_SKIP:
	case mATCH:
	case dEL:
	  positionInGenome -= cigar[c].length;
	  break;
	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  fusion.dir = cigar[c].opcode;
	  if (fusion.dir == FUSION_RF || fusion.dir == FUSION_RR)
	    positionInGenome = positionInGenome + 1;
	  else
	    positionInGenome = positionInGenome - 1;
	  
	  if (bh.ref_id() < bh.ref_id2() ||
	     (bh.ref_id() == bh.ref_id2() && positionInGenome < cigar[c].length) ||
	     !auto_sort)
	    {
	      fusion.refid1 = bh.ref_id();
	      fusion.refid2 = bh.ref_id2();
	      fusion.left = positionInGenome;
	      fusion.right = cigar[c].length;
	    }
	  else
	    {
	      assert (auto_sort);
	      fusion.refid1 = bh.ref_id2();
	      fusion.refid2 = bh.ref_id();
	      fusion.left = cigar[c].length;
	      fusion.right = positionInGenome;
	    }
	  
	  fusions.push_back(fusion);
	  break;
	default:
	  break;
	}	
    }	
} 

void pair_support(const vector<pair<BowtieHit, BowtieHit> >& best_hits, FusionSet& fusions, FusionSet& fusions_ref)
{
  if (best_hits.size() > fusion_multipairs)
    return;

  for (size_t i = 0; i < best_hits.size(); ++i)
    {
      const BowtieHit& lh = best_hits[i].first;
      const BowtieHit& rh = best_hits[i].second;

      bool left_fusionSpanned = lh.fusion_opcode() != FUSION_NOTHING;
      bool right_fusionSpanned = rh.fusion_opcode() != FUSION_NOTHING;
      if (left_fusionSpanned && right_fusionSpanned)
	continue;
      
      bool fusionSpanned = left_fusionSpanned || right_fusionSpanned;
      bool fusion_leftSide = false;
      
      uint32_t ref_id1 = lh.ref_id2();
      uint32_t ref_id2 = rh.ref_id();
      
      int dir = FUSION_FF;
      if (fusionSpanned)
	{
	  if (left_fusionSpanned)
	    dir = lh.fusion_opcode();
	  else
	    dir = rh.fusion_opcode();
	}
      else
	{
	  if (!lh.antisense_align() && !rh.antisense_align())
	    dir = FUSION_FR;
	  else if (lh.antisense_align() && rh.antisense_align())
	    dir = FUSION_RF;
	  else
	    {
	      if (lh.ref_id() == rh.ref_id())
		{
		  if ((lh.antisense_align() && lh.left() > rh.left()) ||
		      (!lh.antisense_align() && lh.left() < rh.left()))
		    dir = FUSION_FF;
		  else
		    dir = FUSION_RR;
		}
	      else
		{
		  if ((lh.antisense_align() && lh.ref_id() > rh.ref_id()) ||
		      (!lh.antisense_align() && lh.ref_id() < rh.ref_id()))
		    dir = FUSION_FF;
		  else
		    dir = FUSION_RR;
		}
	    }
	}
      
      FusionSet::iterator lb, ub;
      bool unsupport = false;
      
      // int inner_dist = max_report_intron_length * 2;

      int inner_dist = 10000;
      int outer_dist = 10000 * 2;
      int max_dist = 10000 * 2;
      int left1 = 0, left2 = 0, right1 = 0, right2 = 0;
      
      if (fusionSpanned)
	{
	  vector<Fusion> new_fusions;
	  if (left_fusionSpanned)
	    fusions_from_spliced_hit(lh, new_fusions);
	  else
	    fusions_from_spliced_hit(rh, new_fusions);
	  
	  Fusion& fusion = new_fusions[0];
	  dir = fusion.dir;
	  ref_id1 = fusion.refid1;
	  ref_id2 = fusion.refid2;

	  if (left_fusionSpanned && lh.ref_id() != ref_id2 && lh.ref_id2() != ref_id2)
	    unsupport = true;
	  
	  if (right_fusionSpanned && rh.ref_id() != ref_id1 && rh.ref_id2() != ref_id1)
	    unsupport = true;

	  int temp1, temp2;
	  const BowtieHit* fusionHit;
	  const BowtieHit* otherHit;
	  
	  if (left_fusionSpanned)
	    {
	      fusionHit = &lh;
	      otherHit = &rh;
	    }
	  else
	    {
	      fusionHit = &rh;
	      otherHit = &lh;
	    }

	  // check the normal hit (otherHit) is on the left hand side of the fusion.
	  if (ref_id1 == ref_id2)
	    {
	      if (dir == FUSION_FF || dir == FUSION_FR)
		fusion_leftSide = otherHit->left() < fusionHit->left();
	      else if (dir == FUSION_RF)
		fusion_leftSide = otherHit->left() < fusionHit->right();
	      else
		fusion_leftSide = otherHit->right() > fusionHit->left();
	    }
	  else
	    fusion_leftSide = fusionHit->ref_id() == otherHit->ref_id();

	  // *****
	  // daehwan - make sure what '+' and '-' really mean for FF, FR, RF, RR cases
	  // *****
	  if ((dir != FUSION_RF && dir != FUSION_RR) && fusionHit->antisense_align() && (!fusion_leftSide || otherHit->antisense_align()))
	    unsupport = true;

	  if ((dir != FUSION_FR && dir != FUSION_RR) && !fusionHit->antisense_align() && (fusion_leftSide || !otherHit->antisense_align()))
	    unsupport = true;

	  if ((dir == FUSION_FR || dir == FUSION_RR) && !fusionHit->antisense_align() && (fusion_leftSide || otherHit->antisense_align()))
	    unsupport = true;

	  if ((dir == FUSION_RF || dir == FUSION_RR) && fusionHit->antisense_align() && (!fusion_leftSide || !otherHit->antisense_align()))
	    unsupport = true;

	  temp1 = otherHit->left();
	  temp2 = otherHit->right();

	  if ((fusion_leftSide && dir == FUSION_RF) || (!fusion_leftSide && dir != FUSION_FR))
	    {
	      temp2 = temp1 + inner_dist;
	      if (temp1 > outer_dist)
		temp1 = temp1 - outer_dist;
	      else
		temp1 = 0;
	    }
	  else
	    {
	      if (temp2 >= inner_dist)
		temp1 = temp2 - inner_dist;
	      else
		temp1 = 0;
	      
	      temp2 = temp2 + outer_dist;
	    }
	  
	  if (fusion_leftSide)
	    {
	      left1 = temp1;
	      left2 = temp2;		  
	    }
	  else
	    {
	      right1 = temp1;
	      right2 = temp2;
	    }
	  
	  lb = fusions_ref.find(fusion);
	  ub = fusions_ref.end();
	  
	  // daehwan - debug
#if 0
	  if (fusion.left == 6994359 && fusion.right == 17581683)
	    {
	      cout << "daehwan - test - pair_with_fusion: " << lh.insert_id() << endl;
	      cout << "edit dist: " << (int)lh.edit_dist() << ", " << (int)rh.edit_dist() << endl;
	      const char* dir_str = "ff";
	      if (dir == FUSION_FR)
		dir_str = "fr";
	      else if (dir == FUSION_RF)
		dir_str = "rf";
	      else if (dir == FUSION_RR)
		dir_str = "rr";
	      
	      cout << dir_str << " : " << (lh.antisense_align() ? "-" : "+") << " " << (rh.antisense_align() ? "-" : "+") <<  endl;
	      cout << lh.ref_id() << ": " << lh.left() << "-" << lh.right() << endl;
	      cout << lh.ref_id2() << ": " << print_cigar(lh.cigar()) << endl;
	      cout << rh.ref_id() << ": " << rh.left() << "-" << rh.right() << endl;
	      cout << rh.ref_id2() << ": " << print_cigar(rh.cigar()) << endl;
	      cout << "found: " << (lb == ub ? "no" : "yes") << endl;
	      cout << "unsupport: " << (unsupport ? "yes" : "no") << endl;
	      cout << "fusion_left: " << (fusion_leftSide ? "yes" : "no") << endl;
	      cout << endl;
	    }
#endif
	}
      else
	{
	  if (dir == FUSION_FF)
	    {
	      if (lh.antisense_align())
		{
		  if (rh.right() >= inner_dist)
		    right1 = rh.right() - inner_dist;
		  right2 = right1 + outer_dist;
		  
		  left2 = lh.left() + inner_dist;
		  if (left2 > outer_dist)
		    left1 = left2 - outer_dist;
		}
	      else
		{
		  if (lh.right() >= inner_dist)
		    left1 = lh.right() - inner_dist;
		  left2 = left1 + outer_dist;
		  
		  right2 = rh.left() + inner_dist;
		  if (right2 > outer_dist)
		    right1 = right2 - outer_dist;
		}
	    }
	  else if (dir == FUSION_FR)
	    {
	      if (lh.right() >= inner_dist)
		left1 = lh.right() - inner_dist;
	      left2 = left1 + outer_dist;

	      if (rh.right() >= inner_dist)
		right1 = rh.right() - inner_dist;
	      right2 = right1 + outer_dist;
	    }
	  else if (dir == FUSION_RF)
	    {
	      left2 = lh.left() + inner_dist;
	      right2 = rh.left() + inner_dist;
	      
	      if (left2 > outer_dist)
		left1 = left2 - outer_dist;
	      
	      if (right2 > outer_dist)
		right1 = right2 - outer_dist;
	    }
	  else // if (dir == FUSION_RR)
	    {
	      if (lh.antisense_align())
		{
		  left2 = lh.left() + inner_dist;
		  if (left2 > outer_dist)
		    left1 = left2 - outer_dist;

		  if (rh.right() >= inner_dist)
		    right1 = rh.right() - inner_dist;
		  right2 = right1 + outer_dist;
		}
	      else
		{
		  if (lh.right() >= inner_dist)
		    left1 = lh.right() - inner_dist;
		  left2 = left1 + outer_dist;

		  right2 = rh.left() + inner_dist;
		  if (right2 > outer_dist)
		    right1 = right2 - outer_dist;
		}
	    }
	  
	  // daehwan - debug
#if 0
	  if (fusion.left == 6994359 && fusion.right == 17581683)
	    {
	      const char* dir_str = "ff";
	      if (dir == FUSION_FR)
		dir_str = "fr";
	      else if (dir == FUSION_RF)
		dir_str = "rf";
	      else if (dir == FUSION_RR)
		dir_str = "rr";
	      
	      cout << "paired-end from two chromosomes" << endl;
	      cout << "insert id: " << lh.insert_id() << endl;
	      cout << dir_str << " : " << (lh.antisense_align() ? "-" : "+") << " " << (rh.antisense_align() ? "-" : "+") <<  endl;
	      cout << lh.ref_id() << ": " << lh.left() << "-" << lh.right() << endl;
	      cout << lh.ref_id2() << " " << print_cigar(lh.cigar()) << endl;
	      cout << rh.ref_id() << ": " << rh.left() << "-" << rh.right() << endl;
	      cout << rh.ref_id2() << " " << print_cigar(rh.cigar()) << endl;
	      cout << "left: " << left1 << "-" << left2 << endl;
	      cout << "right: " << right1 << "-" << right2 << endl;
	      cout << endl;
	    }
#endif

	  if (ref_id1 > ref_id2 || (ref_id1 == ref_id2 && lh.left() > rh.left()))
	    {
	      uint32_t temp = ref_id1;
	      ref_id1 = ref_id2;
	      ref_id2 = temp;
	      
	      temp = left1;
	      left1 = right1;
	      right1 = temp;
	      
	      temp = left2;
	      left2 = right2;
	      right2 = temp;
	    }
	  
	  lb = fusions_ref.upper_bound(Fusion(ref_id1, ref_id2, left1, right1));
	  ub = fusions_ref.lower_bound(Fusion(ref_id1, ref_id2, left2, right2));

	  // daehwan
#if 0
	  static const uint32_t chr_id1 = RefSequenceTable::hash_string("chr2");
	  static const uint32_t chr_id2 = RefSequenceTable::hash_string("chr3");

	  if ((lh.ref_id() == chr_id1 && rh.ref_id() == chr_id2) ||
	      (lh.ref_id() == chr_id2 && rh.ref_id() == chr_id1))
	    {
	      // KPL-4   BSG     19:571325-583492:1      NFIX    19:13106584-13209610:1
	      // const uint32_t left1 = 571325, right1 = 583492, left2 = 13106584, right2 = 13209610;
	      
	      // SK-BR-3	DHX35	20:37590942-37668366:1	ITCH	20:32951041-33099198:1
	      // const uint32_t left1 = 37590942, right1 = 37668366, left2 = 32951041, right2 = 33099198;
	      
	      // SK-BR-3	NFS1	20:34220262-34287287:-1	PREX1	20:47240790-47444420:-1
	      // const uint32_t left1 = 34220262, right1 = 34287287, left2 = 47240790, right2 = 47444420;

	      // VCaP	TIA1	2:70436576-70475792:-1	DIRC2	3:122513642-122599986:1
	      uint32_t left1 = 70436576, right1 = 70475792, left2 = 122513642, right2 = 122599986;
	      if (lh.ref_id() != chr_id1)
		{
		  uint32_t temp = left1;
		  left1 = left2;
		  left2 = temp;

		  temp = right1;
		  right1 = right2;
		  right2 = temp;
		}
	      
	      if ((lh.left() >= left1 && lh.left() <= right1 && rh.left() >= left2 && rh.left() <= right2) ||
		  (lh.left() >= left2 && lh.left() <= right2 && rh.left() >= left1 && rh.left() <= right1))
		{
		  for (size_t t = 0; t < left.size(); ++t)
		    {
		      const BowtieHit& lh = left[t];
		      const BowtieHit& rh = right[t];
		      
		      const char* dir_str = "ff";
		      if (dir == FUSION_FR)
			dir_str = "fr";
		      else if (dir == FUSION_RF)
			dir_str = "rf";
		      else if (dir == FUSION_RR)
			dir_str = "rr";
		      
		      cout << "paired-end from two chromosomes" << endl;
		      cout << "insert id: " << lh.insert_id() << endl;
		      cout << dir_str << " : " << (lh.antisense_align() ? "-" : "+") << " " << (rh.antisense_align() ? "-" : "+") <<  endl;
		      cout << lh.ref_id() << ": " << lh.left() << "-" << lh.right() << endl;
		      cout << lh.ref_id2() << " " << print_cigar(lh.cigar()) << endl;
		      cout << rh.ref_id() << ": " << rh.left() << "-" << rh.right() << endl;
		      cout << rh.ref_id2() << " " << print_cigar(rh.cigar()) << endl;
		      cout << "left: " << left1 << "-" << left2 << endl;
		      cout << "right: " << right1 << "-" << right2 << endl;
		      cout << endl;
		    }
		}
	    }
#endif
	}
      
      while (lb != ub && lb != fusions_ref.end())
	{
	  if (lb->first.dir == dir &&
	      lb->first.refid1 == ref_id1 &&
	      lb->first.refid2 == ref_id2 &&
	      ((!fusionSpanned && lb->first.right >= right1 && lb->first.right <= right2) || fusionSpanned) &&
	      !lb->second.reversed)
	    {
	      int dist = 0, left_dist = 0, right_dist = 0;

	      // daehwan - check this out
	      if (!fusionSpanned || fusion_leftSide)
		{
		  if (dir == FUSION_RF || dir == FUSION_RR)
		    left_dist = (left2 - inner_dist) - (int)lb->first.left;
		  else
		    left_dist = (int)lb->first.left - (left1 + inner_dist);
		}
	      
	      if (!fusionSpanned || !fusion_leftSide)
		{
		  if (dir == FUSION_FR || dir == FUSION_RR)
		    right_dist = (int)lb->first.right - (right1 + inner_dist);
		  else
		    right_dist = (right2 - inner_dist) - (int)lb->first.right;
		}

	      dist = abs(left_dist) + abs(right_dist);

	      // daehwan - fix this later
	      if (dist < 0 && fusionSpanned)
	      	unsupport = true;
	      bool pass = dist < max_dist;
	      
	      // daehwan
#if 0
	      // if(!fusionSpanned)
	      // if (pass && !unsupport)
	      if (lb->first.left == 6994359 && lb->first.right == 17581683)
		{
		  const char* dir_str = "ff";
		  if (dir == FUSION_FR)
		    dir_str = "fr";
		  else if (dir == FUSION_RF)
		    dir_str = "rf";
		  else if (dir == FUSION_RR)
		    dir_str = "rr";
		  
		  cout << dir_str << endl;
		  cout << "dist: " << dist << endl;
		  cout << lb->first.refid1 << " " << lb->first.refid2 << endl;
		  cout << lb->first.left << "-" << lb->first.right << endl;
		  cout << "unsupport: " << (unsupport ? "yes" : "no") << endl;
		  cout << "pass: " << (pass ? "yes" : "no") << endl;
		  cout << "ids: " << lh.insert_id() << " : " << rh.insert_id() << endl;
		  cout << endl << endl;
		}
#endif
		  
	      FusionSet::iterator itr = fusions.find(lb->first);
	      if (itr != fusions.end())
		{
		  if (unsupport)
		    ++(itr->second.unsupport_count_pair);
		  else if (pass)
		    {
		      if (fusionSpanned)
			++(itr->second.pair_count_fusion);
		      else
			{
			  itr->second.vPairSupport.push_back(FusionPairSupport(left_dist, right_dist));
			  ++(itr->second.pair_count);
			  
			  if (itr->second.vPairSupport.size() >= 300)
			    {
			      sort(itr->second.vPairSupport.begin(), itr->second.vPairSupport.end());
			      itr->second.vPairSupport.erase(itr->second.vPairSupport.begin() + 200, itr->second.vPairSupport.end());
			    }
			}
		    }
		}
	      else
		{
		  FusionStat fusionStat;
		  if (unsupport)
		    fusionStat.unsupport_count_pair = 1;
		  else if (pass)
		    {
		      if (fusionSpanned)
			fusionStat.pair_count_fusion = 1;
		      else
			{
			  fusionStat.vPairSupport.push_back(FusionPairSupport(left_dist, right_dist));
			  fusionStat.pair_count = 1;
			}
		    }
		  
		  fusions[lb->first] = fusionStat;			    
		}
	    }
	  
	  if (fusionSpanned)
	    break;
	  
	  ++lb;
	}
    }
}

void merge_with(FusionSimpleSet& fusions, const FusionSimpleSet& other_fusions)
{
  for (FusionSimpleSet::const_iterator other_itr = other_fusions.begin(); other_itr != other_fusions.end(); ++other_itr)
    {
      FusionSimpleSet::iterator itr = fusions.find(other_itr->first);
      if (itr != fusions.end())
	{
	  FusionSimpleStat& curr  = itr->second;
	  curr.merge_with(other_itr->second);
	}
      else
	{
	  fusions[other_itr->first] = other_itr->second;
	}
    }
}

void merge_with(FusionSet& fusions, const FusionSet& other_fusions)
{
  for (FusionSet::const_iterator other_itr = other_fusions.begin(); other_itr != other_fusions.end(); ++other_itr)
    {
      FusionSet::iterator itr = fusions.find(other_itr->first);
      if (itr != fusions.end())
	{
	  FusionStat& curr  = itr->second;
	  curr.merge_with(other_itr->second);
	}
      else
	{
	  fusions[other_itr->first] = other_itr->second;
	}
    }
}
