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

using namespace std;
using namespace seqan;
using std::set;


// Length of the outer dimension of a single insert from the paired-end library
static int min_anchor_len = 5;
static int min_intron_len = 40;
static int min_exon_len = 50;

static bool verbose = false;

void print_usage()
{
    fprintf(stderr, "Usage:   juncs_db <map1.bwtout> <map2.bwtout> <ref_ebwt_basename>\n");
}

typedef vector<string> Mapped;

struct r_to_l_lexcompare
{
 bool operator()(const string& s1, const string& s2) const
 {
     return lexicographical_compare(s1.rbegin(), s1.rend(),
                                    s2.rbegin(), s2.rend());
 }
};


bool possible_cotranscript(const BowtieHit& h1, const BowtieHit& h2, bool check_strand = true)
{
	if (h1.insert_id != h2.insert_id) 
		return false;
	int min_mate_inner_dist = insert_len - h1.read_len() - 
		h2.read_len() - insert_len_std_dev;
	if (max_mate_inner_dist == -1)
	{
		max_mate_inner_dist = insert_len - h1.read_len() - 
		h2.read_len() + insert_len_std_dev;
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
	
	// Sanity checking to verify hits are sorted by insert id.
#if !NDEBUG		
	for (size_t i = 1; i < hits1_in_ref.size(); ++i)
	{
		const BowtieHit& h1 = hits1_in_ref[i];
		const BowtieHit& h2 = hits1_in_ref[i-1];
		assert(h1.insert_id >= h2.insert_id);
	}
	
	for (size_t i = 1; i < hits2_in_ref.size(); ++i)
	{
		const BowtieHit& h1 = hits2_in_ref[i];
		const BowtieHit& h2 = hits2_in_ref[i-1];
		assert(h1.insert_id >= h2.insert_id);
	}
#endif
	
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

bool splice_junc_lt(const pair<size_t, size_t>& lhs, 
					const pair<size_t, size_t>& rhs)
{
	if (lhs.first < rhs.first)
		return true;
	else
		return lhs.second < rhs.second;
}




/* The following computation identifies possible splice sites.
 
 |||||||||||||||-----------------------------------------||||||||||||||||
    read_len                  inner_dist                      read_len    
 Where r is the length of a read and x is the *internal* distance between 
 mates in the *genomic* coordinate space.  We first check to see whether
 x is too long to be an insert from a single contiguous substring of the 
 genome.  That is, if x is longer than you'd expect from a unspliced 
 molecule, than the insert likely spans a splice junction.
 
 The expected value of x when the insert doesn't span a splice is the
 library insert size mean (I) minus twice the read length, plus or minus 
 the library insert size standard deviation.  If x is longer than that,
 the insert probably spans a splice.
 
 For an insert that spans a splice junction, that junction must fall between
 the ends of the mates.  Let the right side of the alignment of the left
 read be called minor_hit_end, and the left side of the alignment of the right 
 read be major_hit_start.  Then an insert's inner_dist = major_hit_start - 
 minor_hit_end.  Let the expected inner distance for a insert that doesn't 
 cross a junction be Insert_mean - 2* read_len. Then a splice-crossing insert 
 must have inner_dist >= expected_inner_dist.  
 
 Let the actual splice position = (splice_left, splice_right). Then 
 (splice_left - minor_hit_end) + (major_hit_start - splice_right) = 
 expected_inner_dist +/- std_dev.
 */


// Hops that are (right_motif,left_motif) represent exonic links, and
// increase the transcriptomic distance.

// Hops that are (left_motif, right_motif) represent intronic links,
// and don't increase the transcriptomic distance, but we may emit them as
// potential splice sites, if they are on a closure between the ends of an
// insert.

template<typename TStr>
void print_splices(const std::set<pair<size_t, size_t> >& potential_splices,
				   int read_len,
				   const string& tag,
				   TStr& ref_str,
				   const string& ref_name,
				   ostream& splice_db)
{
	int half_splice_len = read_len - min_anchor_len;
	for (std::set<pair<size_t, size_t> >::iterator j = potential_splices.begin();
		 j != potential_splices.end();
		 ++j)
	{
		size_t left_start, right_start;
		uint32_t left_end, right_end;
		
		left_start = (int)j->first - half_splice_len + 1 >= 0 ? (int)j->first - half_splice_len + 1 : 0;
		left_end = left_start + half_splice_len;
		
		right_start = j->second;
		right_end = right_start + half_splice_len < length(ref_str) ? right_start + half_splice_len : length(ref_str) - right_start;
		
		
		
		Infix<String<Dna5, Alloc<> > >::Type left_splice = infix(ref_str,
																 left_start, 
																 left_end);
		Infix<String<Dna5, Alloc<> > >::Type right_splice = infix(ref_str, 
																  right_start, 
																  right_end);
		
		splice_db << ">" << ref_name << "|" << left_start << "|" << j->first <<
			"-" << j->second << "|" << right_end << "|" << tag << endl;

		splice_db << left_splice << right_splice << endl;
		
	}	
}


void driver(FILE* map1, FILE* map2, ifstream& ref_stream)
{
	typedef String< Dna5, Alloc<> > Reference;

	SequenceTable it;
	SequenceTable rt;
    HitTable hits1(it,rt);
    HitTable hits2(it,rt);
    
    get_mapped_reads(map1, hits1, false);
    get_mapped_reads(map2, hits2, false);
    
	uint32_t num_happy = 0;

	static const uint32_t bowtie_padding = 5;
	
	uint32_t num_potential_splices = 0;

	ofstream splice_db("possible_junc_db.fa");
	//ofstream splice_db = cout;
	while(ref_stream.good() && 
		  !ref_stream.eof()) 
	{
		Reference ref_str;
		string name;
		readMeta(ref_stream, name, Fasta());
		read(ref_stream, ref_str, Fasta());
		
		// Tracks the number of singleton ALIGNMENTS, not the number of singleton
		// READS in each Bowtie map.
		vector<size_t> map1_singletons;
		vector<size_t> map2_singletons;
		vector<pair<size_t, size_t> > happy_mates;
		
		
		uint32_t ref_id = rt.get_id(name);
		const HitList* p_hits1_in_ref = hits1.get_hits(ref_id);
		const HitList* p_hits2_in_ref = hits2.get_hits(ref_id);
		
		if (!p_hits1_in_ref || !p_hits2_in_ref)
			continue;
		const HitList& hits1_in_ref = *p_hits1_in_ref;
		const HitList& hits2_in_ref = *p_hits2_in_ref;
		
		if (verbose)
			fprintf(stderr, "Looking for happy mates in %s\n", name.c_str());
		
		check_mates(hits1_in_ref,
					hits2_in_ref,
					happy_mates,
					map1_singletons,
					map2_singletons);
		//print_mate_distances(made_dist, happy_mates, hits1_in_ref, hits2_in_ref);
		
		
		std::set<pair<size_t, size_t> > fwd_splices;
		std::set<pair<size_t, size_t> > rev_splices;
		
		int min_read_len = 9999999;
		for (HitList::const_iterator ci = hits1_in_ref.begin(); 
			 ci != hits1_in_ref.end();
			 ++ci)
		{
			if (min_read_len > ci->read_len())
				min_read_len = ci->read_len();
		}
		
		for (HitList::const_iterator ci = hits2_in_ref.begin(); 
			 ci != hits2_in_ref.end();
			 ++ci)
		{
			if (min_read_len > ci->read_len())
				min_read_len = ci->read_len();
		}
		
		if (verbose)
		{
			fprintf(stderr, "\tFound %d happy mates\n", (int)happy_mates.size());
			fprintf(stderr, "Looking for possible splices in %s, with a minimum read length = %d\n", name.c_str(), min_read_len);
		}

		vector<const HitList*> all_hits;
		all_hits.push_back(p_hits1_in_ref);
		all_hits.push_back(p_hits2_in_ref);
		
		uint32_t max_gap = 1.5 * min_read_len;

		
		typedef MappedIntronFinder<Reference> IF;
		typedef JunctionFinder<Reference, MappedIntronFinder<Reference> > JF;
		IF fwd_intron_finder(ref_str, all_hits, "GT", "AG", max_gap);  // Very aggressive allowance for finding donors and acceptors
		
		JF fwd_finder(ref_str, fwd_intron_finder, insert_len, insert_len_std_dev, min_intron_len, min_exon_len, bowtie_padding, max_gap);
		
		if (verbose)
		{
			fprintf(stderr, "\tSearching for forward closures...");
		}
		fwd_finder.possible_junctions(fwd_splices,
									  happy_mates,
									  hits1_in_ref,
									  hits2_in_ref);
		if (verbose)
		{
			fprintf(stderr, "done\n");
		}
		num_potential_splices += fwd_splices.size();
				
		print_splices(fwd_splices,min_read_len, "GTAG|fwd", ref_str, name, splice_db);
		
		IF rev_intron_finder(ref_str, all_hits, "CT", "AC", max_gap);  // Very aggressive allowance for finding donors and acceptors
		
		JF rev_finder(ref_str, rev_intron_finder, insert_len, insert_len_std_dev, min_intron_len, min_exon_len, bowtie_padding, max_gap);
		
		if (verbose)
		{
			fprintf(stderr, "\tSearching for reverse closures...");
		}
		rev_finder.possible_junctions(rev_splices,
									  happy_mates,
									  hits1_in_ref,
									  hits2_in_ref);
		if (verbose)
		{
			fprintf(stderr, "done\n");
		}
		
		print_splices(rev_splices, min_read_len, "GTAG|rev", ref_str, name, splice_db);
		
		num_potential_splices += rev_splices.size();
		
		if (verbose)
		{
			fprintf(stderr, "\tFound %d possible splices\n", (int)rev_splices.size() + (int)fwd_splices.size());
		}
		num_happy += happy_mates.size();
	}
	
    fprintf(stderr, "Found %d happy mates\n", num_happy);
	fprintf(stderr, "Found %d total possible splices\n", num_potential_splices);
}


const char *short_options = "r:I:d:s:va:e:i:";
static struct option long_options[] = {
{"insert-len",      required_argument,       0,            'I'},
{"insert-stddev",      required_argument,       0,            's'},
{"read-len",       required_argument,       0,            'r'},
{"max-dist",       required_argument,       0,            'd'},
{"min-anchor",       required_argument,       0,            'a'},
{"min-intron",       required_argument,       0,            'i'},
{"min-exon",       required_argument,       0,            'e'},
{"verbose",		no_argument,	0,							'v'},
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
				min_intron_len = (uint32_t)parseInt(20, "-i/--min-intron arg must be at least 40");
				break;
			case 'e':
				min_exon_len = (uint32_t)parseInt(20, "-e/--min-exon arg must be at least 20");
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
	
	string map1_file_name = argv[optind++];
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}

	string map2_file_name = argv[optind++];
	
	if(optind >= argc) 
	{
		print_usage();
		return 1;
	}
	
	string ref_base_name = argv[optind++];
	
    FILE* map1_file = fopen(map1_file_name.c_str(), "r");
	if (map1_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map1_file_name.c_str());
		exit(1);
	}

	FILE* map2_file = fopen(map2_file_name.c_str(), "r");
	if (map2_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				map2_file_name.c_str());
		exit(1);
	}
	
	// TODO: probably should check environment for this file, and generate it
	// with bowtie-inspect if the EBWT's there and the FASTA isn't.
	
	// It's really annoying that I'm mixing stream types, but I'd rather use
	// the (well-tested) FASTA reading code from SeqAn.
	string ref_file_name = ref_base_name;
	ifstream ref_stream(ref_file_name.c_str());
	
	if (!ref_stream.good())
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				ref_file_name.c_str());
		exit(1);
	}
    
    driver(map1_file, map2_file, ref_stream);
    
    return 0;
}
