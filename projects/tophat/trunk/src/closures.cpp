/*
 *  closures.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 1/15/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "bwt_map.h"
#include "inserts.h"
#include "closures.h"
#include "reads.h"
#include "tokenize.h"
#include "fusions.h"

using namespace seqan;
using namespace std;


bool possible_cotranscript(const BowtieHit& h1, const BowtieHit& h2, bool check_strand = true)
{
  if (h1.insert_id() != h2.insert_id()) 
    return false;
  
  InsertAlignmentGrade grade(h1, h2);
  return (!grade.too_far && !grade.too_close && grade.opposite_strands);
}


void check_mates(const HitList& hits1_in_ref,
		 const HitList& hits2_in_ref,
		 vector<pair<size_t, size_t> >& happy_mates,
		 vector<size_t>& map1_singletons,
		 vector<size_t>& map2_singletons)
{
  std::set<size_t> marked;
  // TODO: if this shows up on the profile, replace it with a linear
  // time algorithm.  This one is 2*n*lg(n).
  HitList::const_iterator last_good = hits2_in_ref.begin();
  
  for (size_t i = 0; i < hits1_in_ref.size(); ++i)
    {
      pair<HitList::const_iterator, HitList::const_iterator> range_pair;
      range_pair = equal_range(last_good, hits2_in_ref.end(),
			       hits1_in_ref[i], hit_insert_id_lt);
      bool found_hit = false;
      if (range_pair.first != range_pair.second)
	last_good = range_pair.first;
      for (HitList::const_iterator f = range_pair.first;
	   f != range_pair.second;
	   ++f)
	{
	  if (possible_cotranscript(hits1_in_ref[i], *f))
	    {
	      happy_mates.push_back(make_pair(i,f - hits2_in_ref.begin()));
	      marked.insert(f - hits2_in_ref.begin());
	      found_hit = true;
	    }
	}
      if (!found_hit)
	map1_singletons.push_back(i);
    }
  
  for (size_t i = 0; i < hits2_in_ref.size(); ++i)
    {
      if (marked.find(i) == marked.end())
	{
	  map2_singletons.push_back(i);
	}
    }	
}

void find_fusion_closure(HitsForRead& left_hits,
			 HitsForRead& right_hits,
			 std::set<Fusion>& fusions)
{
  if (left_hits.hits.size() > 3 || right_hits.hits.size() > 3)
    return;

  for (size_t left_index = 0; left_index < left_hits.hits.size(); ++left_index)
    {
      for (size_t right_index = 0; right_index < right_hits.hits.size(); ++right_index)
	{
	  BowtieHit* leftHit = &left_hits.hits[left_index];
	  BowtieHit* rightHit = &right_hits.hits[right_index];

	  uint32_t dir = FUSION_FF;
	  if (leftHit->antisense_align() != rightHit->antisense_align())
	    dir = FUSION_FF;
	  else if (!leftHit->antisense_align())
	    dir = FUSION_FR;
	  else
	    dir = FUSION_RF;

	  if (dir == FUSION_FF && leftHit->antisense_align())
	    {
	      BowtieHit * tmp = leftHit;
	      leftHit = rightHit;
	      rightHit = tmp;
	    }
	  
	  uint32_t left, right;
	  if (dir == FUSION_FF || dir == FUSION_FR)
	    left = leftHit->right() - 4;
	  else
	    left = leftHit->left() + 4;
	  
	  if (dir == FUSION_FF || dir == FUSION_RF)
	    right = rightHit->left() + 4;
	  else
	    right = rightHit->right() - 4;

	  if (leftHit->ref_id() == rightHit->ref_id() && dir == FUSION_FF)
	    {
	      int dist = 0;
	      dist = rightHit->left() - leftHit->right();
	      
	      if (dist >= 0 && dist <= (int)fusion_min_dist)
		continue;
	    }
	  
	  uint32_t ref_id1 = leftHit->ref_id();
	  uint32_t ref_id2 = rightHit->ref_id();
	  
	  if (dir == FUSION_FR || dir == FUSION_RF)
	    {
	      if ((ref_id2 < ref_id1) || (ref_id1 == ref_id2 && left > right))
		{
		  uint32_t temp = ref_id1;
		  ref_id1 = ref_id2;
		  ref_id2 = temp;
		  
		  temp = left;
		  left = right;
		  right = temp;
		}
	    }

	  fusions.insert(Fusion(ref_id1, ref_id2, left, right, dir));
	}
    }
}

bool prefer_shorter_pairs = true;

typedef pair<InsertAlignmentGrade, vector<InsertAlignment> > BestPairingHits; 

template<typename Visitor>
void visit_best_pairing(HitsForRead& left_hit_group,
			HitsForRead& right_hit_group,
			Visitor& visitor)
{
  BestPairingHits insert_best;
	
  for (size_t i = 0; i < left_hit_group.hits.size(); ++i)
    {
      BowtieHit& h1 = left_hit_group.hits[i];
      
      for (size_t j = 0; j < right_hit_group.hits.size(); j++)
	{
	  BowtieHit& h2 = right_hit_group.hits[j];
	  if (h1.ref_id() != h2.ref_id())
	    continue;
	  
	  uint32_t refid = h1.ref_id();
	  InsertAlignmentGrade s(h1, h2);
	  
	  //pair<InsertAlignmentGrade, vector<InsertAlignment> >& insert_best
	  //					= best_status_for_inserts[curr_left_obs_order];
	  InsertAlignmentGrade& current = insert_best.first;
	  // Is the new status better than the current best one?
	  if (current < s)
	    {
	      insert_best.second.clear();
	      current = s;
	      insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
	    }
	  else if (! (s < current))
	    {
	      if (prefer_shorter_pairs && current.num_mapped == 2)
		{
		  pair<int, int> dc = pair_distances(*(insert_best.second[0].left_alignment), *(insert_best.second[0].right_alignment));
		  pair<int, int> ds = pair_distances(h1,h2);
		  if (ds.second < dc.second)
		    {
		      //chucked_for_shorter_pair += insert_best.second.size();
		      insert_best.second.clear();
		      current = s;
		      insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
		      
		    }
		}
	      else
		{
		  insert_best.second.push_back(InsertAlignment(refid, &h1, &h2));
		}
	    }
	}
    }
  
  visitor.visit(insert_best);
}

class CovMapIntronFinder
{	
  std::set<SpliceMotif> _tmp_hits;
  vector<SpliceMotif> _hits;
  RefSequenceTable::Sequence* _ref_seq;
public:
  CovMapIntronFinder() : _ref_seq(NULL) {}
  CovMapIntronFinder(RefSequenceTable::Sequence* ref_seq) : 
    _ref_seq(ref_seq), 
    _lms(vector<bool>(length(*ref_seq), false)),
    _rms(vector<bool>(length(*ref_seq), false)){}
  
  void add_motifs_in_window(int left, 
			    int right,
			    const string& left_motif, 
			    const string& right_motif)
  {
    size_t l_start = max(0,left); 
    for (int i = l_start; 
	 i < min(right, (int)length(*_ref_seq) - 2); 
	 ++i)
      {
	seqan::Infix<RefSequenceTable::Sequence>::Type curr
	  = seqan::infix(*_ref_seq,i, i + 2);
	if (curr == left_motif)
	  _lms[i] = true;
	else if (curr == right_motif)
	  _rms[i] = true;
      }	
  }
  
  void finalize()
  {
    int pos = 0;
    for (vector<bool>::iterator itr = _lms.begin();
	 itr != _lms.end();
	 ++itr, ++pos)
      {
	if (*itr)
	  _hits.push_back(SpliceMotif(pos, true));
      }
    
    pos = 0;
    for (vector<bool>::iterator itr = _rms.begin();
	 itr != _rms.end();
	 ++itr, ++pos)
      {
	if (*itr)
	  _hits.push_back(SpliceMotif(pos, false));
      }
    
    sort(_hits.begin(), _hits.end());
  }
	
  size_t seq_len() const { return length(*_ref_seq); }
  const vector<SpliceMotif>& hits() const { return _hits; }
  vector<bool> _lms;
  vector<bool> _rms;
};

typedef CovMapIntronFinder CIF;

struct RefCIF
{
	RefCIF() : ref_seq(NULL) {}
	RefCIF(const CIF& f, const CIF& r, RefSequenceTable::Sequence* s) :
		fwd_cif(f), rev_cif(r), ref_seq(s) {}
	CIF fwd_cif;
	CIF rev_cif;
	RefSequenceTable::Sequence* ref_seq;
};

class CoverageMapVisitor
{
public:
  map<uint32_t, RefCIF> finders;
  
  CoverageMapVisitor(istream& ref_stream, 
		     RefSequenceTable& rt)
  {
    
    while(ref_stream.good() && 
	  !ref_stream.eof()) 
      {
	RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence;
	string name;
	readMeta(ref_stream, name, Fasta());
	string::size_type space_pos = name.find_first_of(" \t\r");
	if (space_pos != string::npos)
	  {
	    name.resize(space_pos);
	  }
	read(ref_stream, *ref_str, Fasta());
	
	uint32_t ref_id = rt.get_id(name, NULL, 0);
	finders[ref_id] = RefCIF(CIF(ref_str), CIF(ref_str), ref_str);
      }
  }
  
  void visit(BestPairingHits& pairings)
  {
    if (!pairings.first.num_mapped == 2)
      return;

    static string fwd_lm("GT");
    static string fwd_rm("AG");
    static string rev_lm("CT"); 
    static string rev_rm("AC");
    
    for (size_t i = 0; 
	 i < pairings.second.size(); 
	 ++i)
      {
	
	InsertAlignment& al = pairings.second[i];
	
	BowtieHit& bh_left = *(al.left_alignment);
	BowtieHit& bh_right = *(al.right_alignment);
	
	assert (bh_left.ref_id() == bh_right.ref_id());
	
	map<uint32_t, RefCIF >::iterator if_itr = finders.find(bh_left.ref_id());
	assert (if_itr != finders.end());
	
	RefCIF& ref_finders = if_itr->second;
	
	pair<int, int> ds = pair_distances(bh_left, bh_right);
	
	int minor_hit_start, major_hit_start;
	int minor_hit_end, major_hit_end;
	if (bh_left.left() < bh_right.left())
	  {
	    minor_hit_start = (int)bh_left.left();
	    minor_hit_end = (int)bh_left.right();
	    major_hit_start = (int)bh_right.left();
	    major_hit_end = (int)bh_right.right();
	  }
	else
	  {
	    minor_hit_start = (int)bh_right.left();
	    minor_hit_end = (int)bh_right.right();
	    major_hit_start = (int)bh_left.left();
	    major_hit_end = (int)bh_left.right();
	  }
	
	int inner_dist = major_hit_start - minor_hit_end;

	bool skip_fwd = false;
	bool skip_rev = false;

	if (library_type == FR_FIRSTSTRAND)
	  {
	    if (bh_left.antisense_align()) skip_rev = true;
	    else skip_fwd = true;
	  }

	if (library_type == FR_SECONDSTRAND)
	  {
	    if (bh_left.antisense_align()) skip_fwd = true;
	    else skip_rev = true;
	  }

	if (inner_dist > inner_dist_mean + 2 * inner_dist_std_dev)
	  {
	    // Forward strand
	    if (!skip_fwd)
	      {
		ref_finders.fwd_cif.add_motifs_in_window((int)bh_left.left() - 10, 
							 bh_left.right() + 10, 
							 fwd_lm, 
							 fwd_rm);
	    
		ref_finders.fwd_cif.add_motifs_in_window((int)bh_right.left() - 10, 
							 bh_right.right() + 10, 
							 fwd_lm, 
							 fwd_rm);
	      }
	    
	    // Reverse strand
	    if (!skip_rev)
	      {
		ref_finders.rev_cif.add_motifs_in_window((int)bh_left.left() - 10, 
							 bh_left.right() + 10, 
							 rev_lm, 
							 rev_rm);
		
		ref_finders.rev_cif.add_motifs_in_window((int)bh_right.left() - 10, 
							 bh_right.right() + 10, 
							 rev_lm, 
							 rev_rm);
	      }
	  }
      }
  }
  
  void finalize()
  {
    for (map<uint32_t, RefCIF >::iterator itr = finders.begin();
	 itr != finders.end(); 
	 ++itr)
      {
	itr->second.fwd_cif.finalize();
	itr->second.rev_cif.finalize();
	delete itr->second.ref_seq;
	itr->second.ref_seq = NULL;
      }	
  }
};


class JunctionMapVisitor
{
public:
  typedef JunctionFinder<RefSequenceTable::Sequence, CIF> JF;
  
  struct JunctionTable
  {
    JunctionTable() : jf(NULL), possible_splices(NULL) {}
    JunctionTable(ClosureJunctionSet* ps, JF* _jf, bool as) 
      :  jf(_jf), possible_splices(ps), antisense(as) {}
    
    JF* jf;
    ClosureJunctionSet* possible_splices;
    bool antisense;
  };
  
  
  map<uint32_t, pair<JunctionTable, JunctionTable> > _finders;
  
  JunctionMapVisitor(ClosureJunctionSet& fwd_splices, 
		     ClosureJunctionSet& rev_splices, 
		     map<uint32_t, RefCIF >& finders) 
  {
    static const uint32_t bowtie_padding = 5;
    for (map<uint32_t, RefCIF >::iterator itr = finders.begin();
	 itr != finders.end();
	 ++itr)
      {
	JF* fwd_jf = new JF(itr->second.fwd_cif,
			    inner_dist_mean,
			    inner_dist_std_dev,
			    min_closure_intron_length,
			    min_closure_exon_length,
			    bowtie_padding,
			    island_extension);
	JF* rev_jf = new JF(itr->second.rev_cif,
			    inner_dist_mean,
			    inner_dist_std_dev,
			    min_closure_intron_length,
			    min_closure_exon_length,
			    bowtie_padding,
			    island_extension);
	
	_finders[itr->first] = make_pair(JunctionTable(&fwd_splices, fwd_jf,false),
					 JunctionTable(&rev_splices, rev_jf, true));
      }
  }
  
  
  void visit(BestPairingHits& pairings)
  {
    if (!pairings.first.num_mapped == 2)
      return;
    for (size_t i = 0; 
	 i < pairings.second.size(); 
	 ++i)
      {
	
	InsertAlignment& al = pairings.second[i];
	
	BowtieHit& bh_left = *(al.left_alignment);
	BowtieHit& bh_right = *(al.right_alignment);
	
	assert (bh_left.ref_id() == bh_right.ref_id());
	
	map<uint32_t, pair<JunctionTable, JunctionTable> >::iterator if_itr = _finders.find(bh_left.ref_id());
	assert (if_itr != _finders.end());
	
	JF& fwd_jf = *(if_itr->second.first.jf);
	fwd_jf.possible_junctions(*(if_itr->second.first.possible_splices), bh_left, bh_right);
	
	JF& rev_jf = *(if_itr->second.second.jf);
	rev_jf.possible_junctions(*(if_itr->second.second.possible_splices), bh_left, bh_right);
      }
  }
};

void closure_driver(vector<FZPipe>& map1, 
		    vector<FZPipe>& map2, 
		    ifstream& ref_stream, 
		    FILE* juncs_file,
		    FILE* fusions_out)
{
  typedef RefSequenceTable::Sequence Reference;
  
  ReadTable it;
  RefSequenceTable rt(true);

  BowtieHitFactory hit_factory(it, rt);

  std::set<Fusion> fusions;
  
  fprintf (stderr, "Finding near-covered motifs...");
  CoverageMapVisitor cov_map_visitor(ref_stream, rt);
  uint32_t coverage_attempts = 0;
  
  assert(map1.size() == map2.size());
  for (size_t num = 0; num < map1.size(); ++num)
    {
      HitStream left_hs(map1[num].file, &hit_factory, false, true, false);
      HitStream right_hs(map2[num].file, &hit_factory, false, true, false);
      
      HitsForRead curr_left_hit_group;
      HitsForRead curr_right_hit_group;
      
      left_hs.next_read_hits(curr_left_hit_group);
      right_hs.next_read_hits(curr_right_hit_group);
      
      uint32_t curr_right_obs_order = it.observation_order(curr_left_hit_group.insert_id);
      uint32_t curr_left_obs_order = it.observation_order(curr_right_hit_group.insert_id);
      
      while(curr_left_obs_order != VMAXINT32 &&
	    curr_right_obs_order != VMAXINT32)
	{
	  while (curr_left_obs_order < curr_right_obs_order&&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {
	      // Get hit group
	      
	      left_hs.next_read_hits(curr_left_hit_group);
	      curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	    }
	  
	  while (curr_left_obs_order > curr_right_obs_order &&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {
	      // Get hit group
	      
	      right_hs.next_read_hits(curr_right_hit_group);
	      curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	    }
	  
	  while (curr_left_obs_order == curr_right_obs_order &&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {
	      if (num == 0)
		find_fusion_closure(curr_left_hit_group, curr_right_hit_group, fusions);
	      
	      if (coverage_attempts++ % 10000 == 0)
		fprintf (stderr, "Adding covered motifs from pair %d\n", coverage_attempts);

	      visit_best_pairing(curr_left_hit_group, curr_right_hit_group, cov_map_visitor);
	      
	      left_hs.next_read_hits(curr_left_hit_group);
	      curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		    
	      right_hs.next_read_hits(curr_right_hit_group);
	      curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	    }
	}
    }
  
  cov_map_visitor.finalize();
  fprintf (stderr, "done\n");
  
  ClosureJunctionSet fwd_splices;
  ClosureJunctionSet rev_splices;
  
  JunctionMapVisitor junc_map_visitor(fwd_splices, rev_splices, cov_map_visitor.finders);
  fprintf (stderr, "Searching for closures...");
  uint32_t closure_attempts = 0;
  
  for (size_t num = 0; num < map1.size(); ++num)
    {
      map1[num].rewind();
      map2[num].rewind();
      
      HitStream left_hs = HitStream(map1[num].file, &hit_factory, false, true, false);
      HitStream right_hs = HitStream(map2[num].file, &hit_factory, false, true, false);
      
      HitsForRead curr_left_hit_group;
      HitsForRead curr_right_hit_group;
      
      left_hs.next_read_hits(curr_left_hit_group);
      right_hs.next_read_hits(curr_right_hit_group);
      
      uint32_t curr_right_obs_order = it.observation_order(curr_left_hit_group.insert_id);
      uint32_t curr_left_obs_order = it.observation_order(curr_right_hit_group.insert_id);
      
      while(curr_left_obs_order != VMAXINT32 &&
	    curr_right_obs_order != VMAXINT32)
	{
	  while (curr_left_obs_order < curr_right_obs_order &&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {
	      // Get hit group
	      
	      left_hs.next_read_hits(curr_left_hit_group);
	      curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	    }
	  
	  while (curr_left_obs_order > curr_right_obs_order &&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {
	      // Get hit group
	      
	      right_hs.next_read_hits(curr_right_hit_group);
	      curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	    }
	  
	  while (curr_left_obs_order == curr_right_obs_order &&
		 curr_left_obs_order != VMAXINT32 && curr_right_obs_order != VMAXINT32)
	    {	
	      if (closure_attempts++ % 10000 == 0)
		fprintf (stderr, "Trying to close pair %d\n", closure_attempts);

	      visit_best_pairing(curr_left_hit_group, curr_right_hit_group, junc_map_visitor);
	      left_hs.next_read_hits(curr_left_hit_group);
	      curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	      
	      right_hs.next_read_hits(curr_right_hit_group);
	      curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	    }
	}
    }

  for (size_t num = 0; num < map1.size(); ++num)
    {
      map1[num].close();
      map2[num].close();
    }
  
  fprintf(stderr, "%lu Forward strand splices\n", fwd_splices.size());
  fprintf(stderr, "%lu Reverse strand splices\n", rev_splices.size());
  
  fprintf (stderr, "done\n");
  uint32_t num_potential_splices = 0;
  fprintf (stderr, "Reporting possible junctions...");
  map<uint32_t, pair<JunctionMapVisitor::JunctionTable, JunctionMapVisitor::JunctionTable> >::iterator f_itr;
  f_itr = junc_map_visitor._finders.begin();
  
  ClosureJunctionSet::iterator j_itr;
  j_itr = fwd_splices.begin();
  while (j_itr != fwd_splices.end())
    {
      fprintf (juncs_file,"%s\t%u\t%u\t%c\n",
	       rt.get_name(j_itr->refid),
	       j_itr->left,j_itr->right,'+');
      ++num_potential_splices;
      ++j_itr;
    }
  
  j_itr = rev_splices.begin();
  while (j_itr != rev_splices.end())
    {
      fprintf (juncs_file,"%s\t%u\t%u\t%c\n",
	       rt.get_name(j_itr->refid),
	       j_itr->left,j_itr->right,'-');
      ++num_potential_splices;
      ++j_itr;
    }
  
  //accept_all_best_hits(best_status_for_inserts);
  fprintf(stderr, "done\n");
  fprintf(stderr, "Searched for closures between %d pairs\n", searched);
  fprintf(stderr, "Successfully closed %d pairs\n", closed);
  
  fprintf(stderr, "Found %d total possible splices\n", num_potential_splices);

  // daehwan
#if 0
  fprintf (stderr, "Reporting potential fusions...\n");
  if(fusions_out){
    for(std::set<Fusion>::iterator itr = fusions.begin(); itr != fusions.end(); ++itr){
      const char* ref_name1 = rt.get_name(itr->refid1);
      const char* ref_name2 = rt.get_name(itr->refid2);
      
      const char* dir = "";
      if (itr->dir == FUSION_FR)
	dir = "fr";
      else if(itr->dir == FUSION_RF)
	dir = "rf";
      else
	dir = "ff";
      
      fprintf(fusions_out,
	      "%s\t%d\t%s\t%d\t%s\n",
	      ref_name1,
	      itr->left,
	      ref_name2,
	      itr->right,
	      dir);
    }
    fclose(fusions_out);
  }else{
    fprintf(stderr, "Failed to open fusions file for writing\n");
  }
#endif
}

void print_usage()
{
  fprintf(stderr, "Usage:   closure_juncs <closure.juncs> <ref.fa> <left_map.bwtout>  <right_map.bwtout>\n");
}

int main(int argc, char** argv)
{
  fprintf(stderr, "closure_juncs v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
  fprintf(stderr, "---------------------------\n");
  
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  
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

  string fusions_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string ref_fasta = argv[optind++];

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string left_file_list = argv[optind++];
  vector<string> left_file_names;
  vector<FZPipe> left_files;
  tokenize(left_file_list, ",", left_file_names);

  string unzcmd = getUnpackCmd(left_file_names[0], false);
  for (size_t i = 0; i < left_file_names.size(); ++i)
    {
      FZPipe seg_file(left_file_names[i], unzcmd);
      if (seg_file.file == NULL)
        {
	  fprintf(stderr, "Error: cannot open file %s for reading\n",
		  left_file_names[i].c_str());
	  exit(1);
        }
      left_files.push_back(seg_file);
    }

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string right_file_list = argv[optind++];
  vector<string> right_file_names;
  vector<FZPipe> right_files;
  tokenize(right_file_list, ",", right_file_names);
  for (size_t i = 0; i < right_file_names.size(); ++i)
    {
      FZPipe seg_file(right_file_names[i], unzcmd);
      if (seg_file.file == NULL)
	{
	  fprintf(stderr, "Error: cannot open %s for reading\n",
		  right_file_names[i].c_str());
	  exit(1);
	}
      right_files.push_back(seg_file);
    }

  ifstream ref_stream(ref_fasta.c_str(), ifstream::in);

  FILE* splice_db = fopen(junctions_file_name.c_str(), "w");
  if (splice_db == NULL)
    {
      fprintf(stderr, "Error: cannot open junctions file %s for writing\n",
	      junctions_file_name.c_str());
      exit(1);
    }

  FILE* fusion_db = fopen(fusions_file_name.c_str(), "w");
  if (splice_db == NULL)
    {
      fprintf(stderr, "Error: cannot open fusions file %s for writing\n",
	      fusions_file_name.c_str());
      exit(1);
    }

  closure_driver(left_files,
		 right_files,
		 ref_stream,
		 splice_db,
		 fusion_db);
  
  return 0;
}
