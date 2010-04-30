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
#include "bwt_map.h"
#include "inserts.h"
#include "closures.h"
#include "reads.h"
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;


bool possible_cotranscript(const BowtieHit& h1, const BowtieHit& h2, bool check_strand = true)
{
	if (h1.insert_id() != h2.insert_id()) 
		return false;
	int min_mate_inner_dist = inner_dist_mean - inner_dist_std_dev;
	if (max_mate_inner_dist == -1)
	{
		max_mate_inner_dist = inner_dist_mean + inner_dist_std_dev;
	}
	
	InsertAlignmentGrade grade(h1,h2, min_mate_inner_dist, max_mate_inner_dist);
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
			
			int min_mate_inner_dist = inner_dist_mean - inner_dist_std_dev;
			if (max_mate_inner_dist == -1)
			{
				max_mate_inner_dist = inner_dist_mean + inner_dist_std_dev;
			}
			InsertAlignmentGrade s(h1, h2, min_mate_inner_dist, max_mate_inner_dist);
			
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
		
		//fprintf (stderr, "Loading reference stream\n");
		//ofstream splice_db = cout;
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
			
			uint32_t ref_id = rt.get_id(name, NULL);
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
			
			//int exp_inner_dist = inner_dist_mean;
			
			if (bh_right.left() - bh_left.right() > inner_dist_mean + 2 * inner_dist_std_dev)
			{
			
				// Forward strand
				ref_finders.fwd_cif.add_motifs_in_window((int)bh_left.left() - 10, 
														 bh_left.right() + 10, 
														 fwd_lm, 
														 fwd_rm);

				ref_finders.fwd_cif.add_motifs_in_window((int)bh_right.left() - 10, 
														 bh_right.right() + 10, 
														 fwd_lm, 
														 fwd_rm);
				
				// Reverse strand
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
			:  jf(_jf), possible_splices(ps), antisense(as) 
		{
			
		}
		
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

void closure_driver(FILE* map1, 
					FILE* map2, 
					ifstream& ref_stream, 
					FILE* juncs_file)
{
	typedef RefSequenceTable::Sequence Reference;
	
	ReadTable it;
	RefSequenceTable rt(true);

    BowtieHitFactory hit_factory(it,rt);
	
	HitStream left_hs(map1, &hit_factory, false, true, false);
	HitStream right_hs(map2, &hit_factory, false, true, false);
	
	HitsForRead curr_left_hit_group;
	HitsForRead curr_right_hit_group;
	
	left_hs.next_read_hits(curr_left_hit_group);
	right_hs.next_read_hits(curr_right_hit_group);
	
	uint32_t curr_right_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	uint32_t curr_left_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	
	
	fprintf (stderr, "Finding near-covered motifs...");
	CoverageMapVisitor cov_map_visitor(ref_stream, rt);
	uint32_t coverage_attempts = 0;
	while(curr_left_obs_order != 0xFFFFFFFF && 
		  curr_right_obs_order != 0xFFFFFFFF)
	{
		//uint64_t next_left_id = left_hs.next_group_id();
		//uint64_t next_right_id = right_hs.next_group_id();
		
		//uint32_t next_left_order = unmapped_reads.observation_order(next_left_id);
		//uint32_t next_left_order = unmapped_reads.observation_order(next_right_id);
		
		while (curr_left_obs_order < curr_right_obs_order&&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{
			// Get hit group
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		}
		
		while (curr_left_obs_order > curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{
			// Get hit group
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
			
		}
		
		while (curr_left_obs_order == curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{			
			if (coverage_attempts++ % 1000 == 0)
				fprintf (stderr, "Adding covered motifs from pair %d\n", coverage_attempts); 
			visit_best_pairing(curr_left_hit_group, curr_right_hit_group, cov_map_visitor);
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
	}
	
	cov_map_visitor.finalize();
	fprintf (stderr, "done\n");
	
	rewind(map1);
	rewind(map2);
	
	
	left_hs = HitStream(map1, &hit_factory, false, true, false);
	right_hs = HitStream(map2, &hit_factory, false, true, false);
	
	left_hs.next_read_hits(curr_left_hit_group);
	right_hs.next_read_hits(curr_right_hit_group);
	
	curr_right_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	curr_left_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	
	ClosureJunctionSet fwd_splices;
	ClosureJunctionSet rev_splices;
	
	JunctionMapVisitor junc_map_visitor(fwd_splices, rev_splices, cov_map_visitor.finders);
	fprintf (stderr, "Searching for closures...");
	uint32_t closure_attempts = 0;
	while(curr_left_obs_order != 0xFFFFFFFF && 
		  curr_right_obs_order != 0xFFFFFFFF)
	{
		//uint64_t next_left_id = left_hs.next_group_id();
		//uint64_t next_right_id = right_hs.next_group_id();
		
		//uint32_t next_left_order = unmapped_reads.observation_order(next_left_id);
		//uint32_t next_left_order = unmapped_reads.observation_order(next_right_id);
		
		while (curr_left_obs_order < curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{
			// Get hit group
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		}
		
		while (curr_left_obs_order > curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{
			// Get hit group
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
			
		}
		
		while (curr_left_obs_order == curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{	
			if (closure_attempts++ % 1000 == 0)
				fprintf (stderr, "Trying to close pair %d\n", closure_attempts); 
			visit_best_pairing(curr_left_hit_group, curr_right_hit_group, junc_map_visitor);
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
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
}

void print_usage()
{
    fprintf(stderr, "Usage:   closure_juncs <closure.juncs> <ref.fa> <left_map.bwtout>  <right_map.bwtout>\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

//const char *short_options = "i:I:e:r:s:fq";
//
//static struct option long_options[] = {
//{"min-anchor",       required_argument,       0,            'a'},
//{"insert-len",       required_argument,       0,            'r'},
//{"insert-len-std-dev",       required_argument,       0,    's'},
//{"min-anchor",       required_argument,       0,            'a'},
//{"min-intron",       required_argument,       0,            'i'},
//{"island-extension",       required_argument,       0,            'e'},
//{0, 0, 0, 0} // terminator
//};

//int parse_options(int argc, char** argv)
//{
//    int option_index = 0;
//    int next_option;
//    do {
//        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
//        switch (next_option) {
//			case 'r':
//	   			insert_len = (uint32_t)parseInt(1, "-I/--insert-len arg must be at least 1", print_usage);
//	   			break;
//			case 's':
//	   			insert_len_std_dev = (uint32_t)parseInt(1, "-s/--insert-stddev arg must be at least 1", print_usage);
//	   			break;
//			case 'i':
//				min_intron_length = (uint32_t)parseInt(1, "-a/--min-intron arg must be at least 1", print_usage);
//				break;
//			case 'I':
//				max_intron_length = parseInt(1,"-I arg must be at least 1", print_usage);
//				break;
//			case 'e':
//				island_extension = parseInt(0, "-e arg must be at least 1", print_usage);
//				break;
//			case 'f':
//				reads_format = FASTA;
//				break;
//			case 'q':
//				reads_format = FASTQ;
//				break; 
//            case -1:     /* Done with options. */
//                break;
//            default:
//                print_usage();
//                return 1;
//        }
//    } while(next_option != -1);
//    
//	max_mate_inner_dist = max_intron_length + insert_len;
//	
//    return 0;
//}

int main(int argc, char** argv)
{
	fprintf(stderr, "closure_juncs v%s\n", PACKAGE_VERSION); 
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
	
    string ref_fasta = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string left_map_name = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string right_map_name = argv[optind++];
	
    FILE* left_map = fopen(left_map_name.c_str(), "r");
    if (left_map == NULL)
    {
        fprintf(stderr, "Error: cannot open %s for reading\n",
                left_map_name.c_str());
        exit(1);
    }
	
	FILE* right_map = fopen(right_map_name.c_str(), "r");
    if (right_map == NULL)
    {
        fprintf(stderr, "Error: cannot open %s for reading\n",
                right_map_name.c_str());
        exit(1);
    }
	
    ifstream ref_stream(ref_fasta.c_str(), ifstream::in);
	
	FILE* splice_db = fopen(junctions_file_name.c_str(), "w");
    if (splice_db == NULL)
    {
        fprintf(stderr, "Error: cannot open junctions file %s for writing\n",
                junctions_file_name.c_str());
        exit(1);
    }
	
    closure_driver(left_map,
				   right_map,
				   ref_stream,
				   splice_db);
	
    return 0;
}
