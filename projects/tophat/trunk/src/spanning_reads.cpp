/*
 *  spanning_reads.c
 *  CSAMapper
 *
 *  Created by Cole Trapnell on 8/4/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>

#include <algorithm>
#include <numeric>
#include <cassert>
#include <bitset>
#include "common.h"
#include "reads.h"
#include "alphabet.h"
#include "timer.h"
#include "islands.h"

using namespace std;

typedef string RefContig; 
typedef unsigned int RefID;

static unsigned int seq_key_len = 4;
static unsigned int max_span_mismatches = 0;
static unsigned int seed_size = 28;

static unsigned int max_memory_megs = 1024;
static int single_island_juncs_above = 9999; // score should max out at 1000, so the default excludes all single island juncs

static bool verbose = false;
//static const unsigned int seq_table_size = 16777216;
//static const unsigned int seq_table_size = 1048576;
//static const unsigned int seq_table_size = max_intron_length;



static unsigned int reads_processed = 0;

struct Exon
{
	string short_name;
	string long_name;
	RefID reference_ctg_id;
	unsigned int pos_in_ref;
	int score;
	string seq;
};

typedef map<RefID, vector<Exon*> > EXONS_FOR_REFS;
EXONS_FOR_REFS exons_for_reference;


map<string, Exon*> names_to_exons;
map<string, RefID> names_to_refctg_ids;

map<RefID, RefContig> ids_to_refctgs;
RefID next_ref_id = 0u;
RefID id_for_refctg(const string& ref_name)
{
	
	map<string, RefID>::iterator i = names_to_refctg_ids.find(ref_name);
	if (i == names_to_refctg_ids.end())
	{	
		names_to_refctg_ids[ref_name] = next_ref_id;
		ids_to_refctgs[next_ref_id] = ref_name;
		return next_ref_id++;
	}
	return i->second;
}


static void print_usage() 
{
	fprintf(stderr, "Usage: spanning_reads [options] <islands.fa> <islands.gff> <unmapped_reads.fa>\n");
	fprintf(stderr, "    -v       verbose output\n");
	fprintf(stderr, "    -s <int> Seed length                                [default: 28]\n");
	fprintf(stderr, "    -a <int> anchor length                              [default: 5]\n");
	fprintf(stderr, "    -m <int> mismatches allowed in extension            [default: 0]\n");
	fprintf(stderr, "    -I <int> maximum intron length                      [default: 20000]\n");
	fprintf(stderr, "    -i <int> minimum intron length                      [default: 70]\n");
	fprintf(stderr, "    -S       allow junctions within individual islands \n");
	fprintf(stderr, "    -M <int> maximum memory footprint (MB)              [default: 2048]\n");
	
}


void read_exons(FILE* exon_fasta, FILE* exon_gff)
{
	int total_exon_length = 0;
	int num_exons = 0;
	while(!feof(exon_fasta))
	{
		Exon* e = new Exon;
		if(!next_fasta_record(exon_fasta, e->long_name, e->seq))
		{	
			delete e;
			break;	
		}
		
		total_exon_length += e->seq.length();
		num_exons++;
		
		RefID ref = 0u;
		unsigned int pos_in_ref = 0;
		
		char* long_name = strdup(e->long_name.c_str());
		char* ln = long_name;
		
		char* tok = strsep(&long_name, " ");
		if (!*tok)
			continue;
		string short_name = tok;
		
		for(int i = 0; (tok = strsep(&long_name, ","));)
		{
			if (*tok)
			{
				switch(i)
				{
					
					case 0: ref = id_for_refctg(tok); break;
					case 1: pos_in_ref = atoi(tok); break;
					default: break;
				}
				++i;
			}
		}
		
		e->reference_ctg_id = ref;
		e->pos_in_ref = pos_in_ref;
		e->short_name = short_name;
		e->score = 0;
		
		names_to_exons[e->short_name] = e;
		
		EXONS_FOR_REFS::iterator i = exons_for_reference.find(ref);
		if (i == exons_for_reference.end())
			exons_for_reference[ref] = vector<Exon*>();
		exons_for_reference[ref].push_back(e);
		
		free(ln);
	}
	
	char buf[2048];
	
	// Throw out the track header
	fgets(buf, sizeof(buf), exon_gff);
	int single_junc_islands = 0;
	int single_junc_island_length = 0;
	while(!feof(exon_gff))
	{
		char name[2048];
		char source[64];
		char type[64];
		int start;
		int stop;
		int score;
		char c[64];
		char d[64];
		char ID[256];
	
		fscanf(exon_gff, "%s%s%s%d%d%d%s%s%s",
			   name, source, type, &start, &stop, &score, c, d, ID);
		map<string, Exon*>::iterator itr = names_to_exons.find(ID);
		if (itr != names_to_exons.end())
		{
			itr->second->score = score;
			if (score > single_island_juncs_above)
			{
				single_junc_island_length += itr->second->seq.length();
				single_junc_islands++;
			}
		}
		//Island i(name, start, stop - start, score, 0);
	}
	
	if (verbose)
	{
		fprintf(stderr, "Loaded %d bp in %d islands\n", total_exon_length, num_exons);
		fprintf(stderr, "Will examine %d islands for single island junctions, %d bp total\n",
				single_junc_islands, single_junc_island_length);
	}
}



unsigned int num_pair_hits = 0;
unsigned int num_valid_pair_hits = 0;
unsigned int num_donors = 0;
unsigned int num_acceptors = 0;
unsigned int num_perfect_splits = 0;
unsigned int num_single_hash_hits = 0;
unsigned int num_double_hash_hits = 0;
unsigned int num_spanning_reads = 0;
unsigned int num_pairs_without_donor_acceptor = 0;

struct SortExonsByRefPos
{
	bool operator()(const Exon* lhs, const Exon* rhs) const
	{
		return lhs->pos_in_ref < rhs->pos_in_ref;
	}
};

struct Splice
{
	Splice(Exon* L, Exon* R, size_t PL, size_t PR)
		: l(L), r(R), pos_in_l(PL), pos_in_r(PR) {}
	Exon* l;
	Exon* r;
	size_t pos_in_l;
	size_t pos_in_r;
};

struct PackedSplice
{
	PackedSplice() : left(0u), seed(0u), right(0u) {}
	PackedSplice(uint32_t l, uint64_t s, uint32_t r) : left(l), seed(s), right(r) {}
	uint32_t left;
	uint64_t seed;
	uint32_t right;
};


// The second element of these pairs is the left (or right) side of a splice
// seed from a possible junction.  The first element is the sequence flanking
// that seed
typedef pair<uint32_t, uint32_t> PackedSpliceHalf;

PackedSpliceHalf pack_left_splice_half(Exon* left, 
									   uint32_t pos_in_l, 
									   unsigned int key_len)
{
	const char* l = left->seq.c_str(); 
	l += pos_in_l;
	
	const char* left_end = l;
	l -= 16;
	
	assert (l + seq_key_len < left->seq.c_str() + left->seq.length());
	
	PackedSpliceHalf packed_half = make_pair(0u,0u);
	
	// Pack up to 32 bits (16 bases) of sequence into left
	if (l < left->seq.c_str())
		l = left->seq.c_str();
	while (l < left_end)
	{
		packed_half.first <<= 2;
		packed_half.first |= (0x3 & charToDna5[(size_t)*l]);
		++l;
	}
	
	// Pack up the seed bits
	for (unsigned int i = 0; i < key_len; ++i)
	{
		packed_half.second <<= 2;
		packed_half.second |= (0x3 & charToDna5[(size_t)*(l + i)]);
	}
	return packed_half;
}


PackedSpliceHalf pack_right_splice_half(Exon* right, 
										uint32_t pos_in_r, 
										unsigned int key_len)
{
	const char* r = right->seq.c_str(); 
	r += pos_in_r - seq_key_len;
		
	PackedSpliceHalf packed_half = make_pair(0u,0u);
	
	// Pack the seed bits
	for (unsigned int i = 0; i < key_len; ++i)
	{
		packed_half.second <<= 2;
		packed_half.second |= (0x3 & charToDna5[(size_t)*(r + i)]);
	}
	
	r += seq_key_len;
	// Now pack 32 bits (16 bases) of sequence into left
	const char* right_end = r + 16;
	
	if ((size_t)(right_end - right->seq.c_str()) > right->seq.length())
		right_end = right->seq.c_str() + right->seq.length();
	while (r < right_end)
	{
		packed_half.first <<= 2;
		packed_half.first |= (0x3 & charToDna5[(size_t)*r]);
		++r;
	}	
	
	return packed_half;
}

PackedSplice combine_splice_halves(const PackedSpliceHalf& left_half, 
								   const PackedSpliceHalf& right_half)
{
	uint64_t seed = left_half.second << (seq_key_len << 1) | right_half.second;
	return PackedSplice(left_half.first,seed, right_half.first);
}

PackedSplice pack_splice(const Splice& s, unsigned int key_len)
{
	const char* l = s.l->seq.c_str(); // l points to beginning of left exon sequence
	l += s.pos_in_l; 

	assert (l + seq_key_len < s.l->seq.c_str() + s.l->seq.length());
	
	const char* r = s.r->seq.c_str(); // r points to beginning of right exon sequence
	r += s.pos_in_r - seq_key_len; 
	//r += 2; // r points to right side of junction;
	
	uint64_t seed = 0;
	uint64_t left = 0;
	uint64_t right = 0;
	
	// Pack up the seed bits
	for (unsigned int i = 0; i < key_len; ++i)
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*(l + i)]);
	}
	
	for (unsigned int i = 0; i < key_len; ++i)
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*(r + i)]);
	}
	
	// Now pack 32 bits (16 bases) of sequence into left
	const char* left_end = l;
	l -= 16;
	if (l < s.l->seq.c_str())
		l = s.l->seq.c_str();
	while (l < left_end)
	{
		left <<= 2;
		left |= (0x3 & charToDna5[(size_t)*l]);
		++l;
	}
	
	r += seq_key_len;
	// Now pack 32 bits (16 bases) of sequence into left
	const char* right_end = r + 16;
	
	if ((size_t)(right_end - s.r->seq.c_str()) > s.r->seq.length())
		right_end = s.r->seq.c_str() + s.r->seq.length();
	while (r < right_end)
	{
		right <<= 2;
		right |= (0x3 & charToDna5[(size_t)*r]);
		++r;
	}
	
	return PackedSplice(left,seed,right);
}

/** Returns the number of characters in strings w1 and w2 that match,
 *  starting from right to left
 */
int get_matching_chars(uint32_t w1, uint32_t w2)
{
	//find the least significant mismatching bit between w1 and w2
	int mismatch_bit = ffs(w1 ^ w2);
	
	// If there is no mismatching bit, the words are equal
	if (!mismatch_bit)
		return -1;
	
	// Given the mismatching bit, determine where the mismatching base is
	mismatch_bit -= 1;
	mismatch_bit -= ((mismatch_bit) & 1);
	mismatch_bit >>= 1;
	
	// Return the number of matching characters.
	return mismatch_bit;
}


/* Represents a hit between a splice seed and a read. */
// TODO: consider packing pos and meta into a single 32-bit int.
struct ReadHit
{
	ReadHit(uint32_t l, uint32_t r, uint32_t p, uint32_t m, bool rc) 
		: left(l), right(r), pos(p), meta(m), reverse_complement(rc) {}
	uint32_t left; // 2-bits per base rep of the left remainder
	uint32_t right; //2-bits per base rep of the right remainder
	uint32_t pos; // position of the seed within the read
	uint32_t meta : 31;
	bool reverse_complement : 1;
} __attribute__((packed));

/* Holds information about a read that we don't want to throw away,
   and only really need for printing */
struct ReadMetadata
{
	ReadMetadata(const string& n, const string& s) : name(n)/*, seq(s)*/{}
	string name;
	//string seq;
};

// read_metadata stores information we don't need for alignment, but that 
// we might want for output, like the name of a read
typedef vector<ReadMetadata> MetadataTable;
MetadataTable read_metadata;

// A MerTable maps k-mers to hits in indexed reads.  See the comment for 
// mer_table
typedef vector<ReadHit> ReadHitList;
typedef vector<ReadHitList> MerTable;

// This table is indexed by k-mer stores positions in indexed reads where those
// k-mers may be found.
MerTable mer_table;

// This just tracks the total number of entries (hits) in the set of indexed
// reads, and is used only in verbose output
unsigned int read_hits = 0;

/** Add a single read to the global mer_table.  The caller can pass either the
 *  forward or the reverse strand of the read - and is responsible for taking
 *  the correct seed region of that read.  Returns the increase in index size
 *  as a result of adding this read
 */

size_t index_read(const string& seq, unsigned int read_num, bool reverse_complement)
{
//	assert (seq.length() <= 32);
//	if (2 * seq_key_len > 32)
//		return 0;
	
	// h is will hold the 2-bit-per-base representation of the k-mer seeds for
	// this read.
	
	uint64_t seed = 0;
	bitset<256> left = 0;
	bitset<256> right = 0;
	const char* p = seq.c_str();
	
	unsigned int seq_len = seq.length();
	const char* seq_end = p + seq_len;
	
	// Build the first seed
	while (p < seq.c_str() + (2 * seq_key_len))
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*p]);
		++p;
	}
	
	// Build the rest of them with a sliding window, adding successive bases
	// to the "right-remainder" word.  
	while (p < seq_end)
	{
		
		right <<= 2;
		right |= (0x3 & charToDna5[(size_t)*p]);		
		++p;
	}	
	
	size_t cap_increase = 0;
	
	// At this point, seed contains the 5'-most 2*seq_key_len bases of the 
	// read, and right contains everthing else on the 3' end.
	
	// This loop will construct successive seed, along with 32-bit words 
	// containing the left and right remainders for each seed
	size_t i = 0;
	size_t new_hits = 0;
	do 
	{
		// Let's not make an out-of-bounds write, if this fails the global 
		// mer_table is too small
		assert (seed < mer_table.size());
		
		read_hits++;
		
		// How many base pairs exist in the right remainder beyond what we have
		// space for ?
		int extra_right_bp = (seq.length() - (i + 2*seq_key_len)) - 16;
		
		uint32_t hit_right;
		if (extra_right_bp > 0)
		{
			//bitset<32> tmp_hit_right = (right >> (extra_right_bp << 1));
			hit_right = (right >> (extra_right_bp << 1)).to_ulong();
		}
		else
		{
			hit_right = right.to_ulong();
		}
		
		uint32_t hit_left = ((left << (256 - 32)) >> (256 - 32)).to_ulong();
		
		size_t prev_cap = mer_table[seed].capacity();
		
		mer_table[seed].push_back(ReadHit(hit_left, hit_right,i, read_num, reverse_complement));
		cap_increase += (mer_table[seed].capacity() - prev_cap) * sizeof (ReadHit);
		new_hits++;
		
		// Take the leftmost base of the seed and stick it into bp
		uint64_t bp = seed & (0x3uLL << ((seq_key_len << 2) - 2));
		
		// Move that base down to the least significant bits of bp
		bp >>= ((seq_key_len << 2) - 2);
		
		// And tack it onto the left remainder of the read
		left <<= 2;
		left |= bp;
		
		// Now take the leftmost base of the right remainder and stick it into 
		// the rightmost position of the seed
		uint32_t right_len = seq_len - (i + seq_key_len * 2);
		//bp = right & (0x3uLL << ((right_len - 1) << 1));
		
		seed <<= 2;
		//cout << right << endl;
		bitset<256> tmp_right = (right >> ((right_len - 1) << 1));
		//cout <<tmp_right << endl;
		seed |= tmp_right.to_ulong();
		seed &=  ~(0xFFFFFFFFFFFFFFFFuLL << (seq_key_len << 2));
			
		//Now remove that leftmost base of the right remainder
		//right &=  ~(0xFFFFFFFFFFFFFFFFuLL << (right_len - 1 << 1));
		if (right_len)
		{
			right.set((right_len - 1 << 1), 0);
			right.set((right_len - 1  << 1) + 1, 0);
		}
		++i;
		
	}while(i <= (size_t)(seq_end - seq.c_str()) - (2 * seq_key_len));
	return cap_increase;
}

/** Adds up to max_reads, both forward and reverse, to the global mer_table
 */
void index_reads(FILE* reads_file, long long max_memory_megs)
{
	string defline, seq;
	
	mer_table.clear();
	
	// Reserve an entry for each k-mer we might see
	size_t mer_table_size = 1 << ((seq_key_len << 1)<<1);
	
	mer_table.resize(mer_table_size);
	
	// We want to reserve some space in the hit vector for each k-mer in the 
	// lookup table, to avoid excessive resizing.
	
	// Estimate the hits per k-mer as the length of the seed region of each read
	//  minus the length of the k-mer plus one.  
	size_t mers_per_read = (seed_size - (2*seq_key_len) + 1);

	long long max_memory = (max_memory_megs * 1024 * 1024);
	long long reads_per_batch = max_memory / (sizeof(ReadHit) * 2 * (mers_per_read));
	
	if (verbose)
	{
		fprintf(stderr, "Index will not consume more than %lld bytes\n", (max_memory) );
		fprintf(stderr, "Using %ld mers per read, %ld bytes per mer, %ld bytes per read\n", (int)2*mers_per_read, sizeof(ReadHit), (sizeof(ReadHit) * 2 * (mers_per_read)));
		fprintf(stderr, "TopHat will load at most %lld reads at a time\n", reads_per_batch);
	}

	
	// Multiply by the total number of reads (including their RCs)
	// and then divide by the number of distinct mers we could see
	// to get the memory needed per mer in the table
	size_t est_hits_per_mer = reads_per_batch * mers_per_read / mer_table_size;
	
	long long  memory_used = 0;
	if (verbose)
	{
		fprintf(stderr, "Estimating %d hits per %d-mer, for a total of %ld reserved entries\n",
				(int)est_hits_per_mer, (2*seq_key_len),est_hits_per_mer * mer_table_size);
	}
	
	for (size_t i = 0; i < mer_table.size(); ++i)
	{
		mer_table[i].reserve(est_hits_per_mer);
		memory_used += mer_table[i].capacity() * sizeof(ReadHit);
	}
	
	if (verbose)
	{
		fprintf(stderr, "Initialized empty table, starting capacity of %lld bytes\n", memory_used);
	}
	
	
	// read_metadata stores the name of each read, and other information we 
	// don't need for alignment, but can't throw out.
	read_metadata.clear();
	
	int read_num = 0;
	
	// Get up to max_reads from the file, adding them to the lookup table
	// as we go
	while(!feof(reads_file) && memory_used <= max_memory)
	{
		seq.clear();
		defline.clear();
		
		// Get a new FASTA record
		// TODO: support FASTQ here?
		if (!next_fasta_record(reads_file, defline, seq))
			return;
		
		string seed;
		if (seq.length() > seed_size) 
			seed = seq.substr(0, seed_size);
		else
			seed = seq;
		read_metadata.push_back(ReadMetadata(defline,seed));
		
		// Add the forward version of the read to the index
		memory_used += index_read(seed, read_num, false);
		
		reverse_complement(seq);
				
		// When we take the RC, we need to take a suffix, rather than a prefix,
		// as the seed. The reads always are given 5' to 3' in the file, so the
		// 5' becomes a suffix under reverse complementing.
		if (seq.length() > seed_size)
			seed = seq.substr(seq.length() - seed_size);
		else
			seed = seq;
		
		// Uncomment the two lines below to get a separate entry
		// in the read_metadata table for a read and its reverse complement
		//read_num++;
		//read_metadata.push_back(ReadMetadata(defline,seed));
		
		// Add the reverse version of the read to the index
		memory_used += index_read(seed, read_num, true);
		
		reads_processed++;
		read_num++;
		
		if (verbose && read_num % 100000 == 0)
		{
			fprintf(stderr, "Processed %d reads, index occupies %lld bytes\n", read_num, memory_used);
		}
		
	} // end while loop over the reads
}

/** 
 *  Computes the Hamming distance between two 16bp words, up to a specified
 *  maximum number of mismatches.
 */
uint32_t mismatching_bases(uint32_t w1_word, 
					  uint32_t w2_word, 
					  int len,
					  uint32_t max_mis)
{
	uint32_t diffs = 0;
	
	int shift = 0;
	int L = len;
	uint32_t misses = 0;
	
	// While we haven't yet exceeded the maximum allowable mismatches,
	// and there are still unaligned bases, keep shift-anding
	while (shift < len && misses <= max_mis)
	{
		int match_chars = 0;
		
		// Get the number of characters matching on the right sides of 
		// both words
		match_chars = get_matching_chars(w1_word, w2_word);
		
		// If they are equal for this shift, we are done,
		// the loop will stop at the next iteration
		if (match_chars == -1)
		{
			match_chars = len;
			shift = len;
		}
		else
		{
			// If there is a mismatch in the remaining words
			// decide how much to shift by and try again
			match_chars = min(len, match_chars);
			int shift_chars = (match_chars + 1);
			
			L -= shift_chars;
			
			int shift_bits = shift_chars << 1;
			
			// Shift right past the matching part and the first mismatching base
			w1_word >>= (shift_bits);
			w2_word >>= (shift_bits);
			
			shift += shift_chars;
			diffs++;
			misses++;
		}
	}
	return diffs;
}


void lookup_splice_in_read_index(RefID ref_ctg_id,
								 const PackedSplice& p, 
								 Exon* left,
								 Exon* right,
								 size_t pos_in_l,
								 size_t pos_in_r,
								 bool antisense)
{	
	//PackedSplice p = pack_splice(s, seq_key_len);
	assert(p.seed < mer_table.size());
	ReadHitList& hl = mer_table[p.seed];
	
	bool filter_debug = false;//left->pos_in_ref == 19961027 && pos_in_l == 201 && right->pos_in_ref == 19962790 && pos_in_r == 268;
	
	string ref;
	
	for (size_t hit = 0; hit < hl.size(); ++hit)
	{
		ReadHit& rh = hl[hit];
		uint32_t pos = rh.pos;
		uint32_t left_mismatches = 0;
		
		if (pos)
		{
			uint32_t read_left = rh.left;
			uint32_t splice_left = p.left;
			if (pos < 16)
			{
				uint32_t left_mask = ~(0xFFFFFFFFu << (pos << 1));
				splice_left &= left_mask;
			}

			
			left_mismatches = mismatching_bases(read_left, 
												splice_left, 
												pos,
												max_span_mismatches);
			
			
			if (filter_debug)
			{
				//cerr << u32ToDna(read_left, pos) << " " << u32ToDna(splice_left, pos)  << " " << left_mismatches << endl;
				ref += u32ToDna(splice_left, min(16, (int)pos)) + " ";
				ref += u32ToDna(p.seed, 2*seq_key_len) + " ";
			}
			
			if (left_mismatches > max_span_mismatches)
			{
				continue;
			}

		}
		 
		if (pos < seed_size - 2 * seq_key_len)
		{

			assert ((int)seed_size - 2 * seq_key_len - pos > 0);
			uint32_t right_read_bases = min(16u,seed_size - 2 * seq_key_len - pos);
			uint32_t right_bits = right_read_bases << 1;
			//uint32_t right_bits = (32ul - (2*seq_key_len) - ((pos + 2 * seq_key_len) - 1)) << 1;
			uint32_t read_right = rh.right;
			uint32_t right_mask = (0xFFFFFFFFu << (32ul - right_bits));
			uint32_t splice_right = p.right & right_mask;
			splice_right >>= ((16 - right_read_bases) << 1);
			
			
			
			uint32_t right_mismatches = mismatching_bases(read_right, 
															splice_right, 
															right_read_bases,
															max_span_mismatches - left_mismatches);
			
			if (filter_debug)
			{
				//cerr << u32ToDna(read_right, right_read_bases) << " " << u32ToDna(splice_right, right_read_bases) << " " << right_mismatches << endl;
				ref += u32ToDna(splice_right, right_read_bases);

			}
			
			if (right_mismatches > max_span_mismatches - left_mismatches)
			{
				continue;
			}
		}
		
		ReadMetadata& rm = read_metadata[rh.meta];
		
		//fprintf(stderr, "%d, %s, %d\n", pos, rm.name.c_str(), rh.meta);
//		fprintf(stdout, "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n", 
//				rm.name.c_str(),
//				ids_to_refctgs[ref_ctg_id].c_str(),
//				antisense ? "-" : "+",
//				left->short_name.c_str(),
//				(int)left->pos_in_ref,
//				(int)pos_in_l + seq_key_len - 1,
//				right->short_name.c_str(),
//				(int)right->pos_in_ref,
//				(int)pos_in_r - seq_key_len);
		
		string left_str = u32ToDna(rh.left, min((int)pos, 16));
		string seed = u32ToDna(p.seed, 2 * seq_key_len);
		string right_str = u32ToDna(rh.right, min(16, (int)seed_size - (int)(2 * seq_key_len) - (int)pos));
		string seq =  left_str + seed + right_str;
		
		if (filter_debug)
		{
			//cerr << u32ToDna(read_right, right_read_bases) << " " << u32ToDna(splice_right, right_read_bases) << " " << right_mismatches << endl;
			//ref += u32ToDna(splice_right, right_read_bases);
			cerr << ref << endl;
			cerr << left_str + " " + seed + " " + right_str << endl;
		}
		
		int left_splice_coord = (int)left->pos_in_ref + (int)pos_in_l + seq_key_len - 1;
		int right_splice_coord = (int)right->pos_in_ref + (int)pos_in_r - seq_key_len;
		int window_length = seq.length() - seq_key_len;
		int left_window_edge = left_splice_coord - window_length;
		int right_window_edge = right_splice_coord + window_length;
		int align_pos = window_length - (left_str.length() + seq_key_len);
		assert (align_pos >= 0 && align_pos <= window_length);
		// TODO: report mismatches in the spliced bowtie output
		fprintf(stdout, "%s\t%s\t%s|%d|%d-%d|%d|GTAG|%s\t%d\t%s\t%s\t1\n", 
				rm.name.c_str(), // read name
				rh.reverse_complement ? "-" : "+", // sense of the alignment
				ids_to_refctgs[ref_ctg_id].c_str(), //reference contig
				left_window_edge,
				left_splice_coord,
				right_splice_coord,
				right_window_edge,
				antisense ? "rev" : "fwd",
				align_pos,
				seq.c_str(),
				string(seq.length(), 'I').c_str()
				);
	}
}



unsigned long long total_splices = 0;

void map_possible_exon_juncs(FILE* reads_file)
{
	total_splices = 0;
	unsigned int forward_splices = 0;
	unsigned int reverse_splices = 0;
	
	Timer splice_map_timer(cerr, "", false);
	typedef pair<size_t, PackedSpliceHalf> SpliceHalf;
	for (EXONS_FOR_REFS::iterator j = exons_for_reference.begin();
		 j != exons_for_reference.end();
		 ++j)
	{
 
		unsigned int contig_fwd_splices = 0;
		unsigned int contig_rev_splices = 0;
		
		vector<SpliceHalf> forward_donors;
		vector<SpliceHalf> forward_acceptors;
		vector<SpliceHalf> reverse_donors;
		vector<SpliceHalf> reverse_acceptors;
		
		// Parallel arrays for convenience
		vector<Exon*> forward_donor_exons;
		vector<Exon*> forward_acceptor_exons;
		vector<Exon*> reverse_donor_exons;
		vector<Exon*> reverse_acceptor_exons;
		
		for (vector<Exon*>::iterator e = j->second.begin();
			e != j->second.end();
			++e)
		{
			string& seq = (*e)->seq;
			Exon* exon = *e;
			unsigned int pos = (*e)->pos_in_ref;
			for (size_t z = seq_key_len + 1; z < seq.length() - seq_key_len - 2; ++z)
			{
				char l = seq[z - 1];
				char r = seq[z];
				if (l == 'G' && r == 'T')
				{
					size_t donor_pos = pos + z - 1;
					size_t s = donor_pos - exon->pos_in_ref - seq_key_len;
					PackedSpliceHalf p = pack_left_splice_half(exon, s, seq_key_len);
					forward_donors.push_back(make_pair(donor_pos,p));
					forward_donor_exons.push_back(exon);
				}
				if (l == 'A' && r == 'G')
				{
					size_t acceptor_pos = pos + z - 1;
					size_t s = acceptor_pos - exon->pos_in_ref + 2 + seq_key_len;
					PackedSpliceHalf p = pack_right_splice_half(exon, s, seq_key_len);
					forward_acceptors.push_back(make_pair(acceptor_pos,p));
					forward_acceptor_exons.push_back(exon);
				}
				if (l == 'C' && r == 'T')
				{
					size_t acceptor_pos = pos + z - 1;
					size_t s = acceptor_pos - exon->pos_in_ref - seq_key_len;
					PackedSpliceHalf p = pack_left_splice_half(exon, s, seq_key_len);
					reverse_acceptors.push_back(make_pair(pos + z - 1,p));
					
					reverse_acceptor_exons.push_back(exon);
				}
				if (l == 'A' && r == 'C')
				{
					size_t donor_pos = pos + z - 1;
					size_t s = donor_pos - exon->pos_in_ref + 2 + seq_key_len;
					PackedSpliceHalf p = pack_right_splice_half(exon, s, seq_key_len);
					reverse_donors.push_back(make_pair(donor_pos,p));
					reverse_donor_exons.push_back(exon);
				}
			}
			//exon->seq.resize(0);
			exon->long_name.resize(0);
		}
		
		for (size_t d = 0; d < forward_donors.size(); ++d)
		{
			bool broke_out = false;
			
			// start pos is a lower bound on downstream acceptor positions
			// to consider
			size_t start_pos = forward_donors[d].first + min_intron_length;

			SpliceHalf dummy = make_pair(start_pos,PackedSpliceHalf());
			vector<SpliceHalf>::iterator lb = upper_bound(forward_acceptors.begin(), 
														  forward_acceptors.end(),
														  dummy);
			
			if (lb == forward_acceptors.end())
				break;
			for (size_t a = lb - forward_acceptors.begin(); 
				 a < forward_acceptors.size(); 
				 ++a)
			{
				if (forward_acceptors[a].first - forward_donors[d].first > max_intron_length)
				{
					broke_out = true;
					break;
				}
				
				// Don't allow splices within exons
				if (forward_acceptor_exons[a] == forward_donor_exons[d] &&
					forward_acceptor_exons[a]->score < single_island_juncs_above)
					continue;
				
				forward_splices++;
				contig_fwd_splices++;
				total_splices++;
				
				Exon* l = forward_donor_exons[d];
				Exon* r = forward_acceptor_exons[a];
				size_t pos_in_l = forward_donors[d].first  - l->pos_in_ref - seq_key_len;
				size_t pos_in_r = forward_acceptors[a].first - r->pos_in_ref + 2 + seq_key_len;
				PackedSplice p = combine_splice_halves(forward_donors[d].second,
													   forward_acceptors[a].second);
				lookup_splice_in_read_index(l->reference_ctg_id,
											p, 
											l, 
											r, 
											pos_in_l, 
											pos_in_r, 
											false);
				

			}
			if (broke_out)
				continue;
		}
		 
		if (verbose)
		{
			fprintf(stderr, 
					"Sense strand for %s, has %d total %d-mer splices\n", 
					ids_to_refctgs[j->first].c_str(),
					contig_fwd_splices,
					2*seq_key_len);
		}
		
		for (size_t a = 0; a < reverse_acceptors.size(); ++a)
		{
			bool broke_out = false;
			// start pos is a lower bound on downstream donor positions
			// to consider
			size_t start_pos = reverse_acceptors[a].first + min_intron_length;
			SpliceHalf dummy = make_pair(start_pos,PackedSpliceHalf());
			
			
			vector<SpliceHalf>::iterator lb = upper_bound(reverse_donors.begin(), 
														  reverse_donors.end(),
														  dummy);
			if (lb == reverse_donors.end())
				break;
			for (size_t d = lb - reverse_donors.begin(); 
				 d < reverse_donors.size(); 
				 ++d)
			{
				if (reverse_donors[d].first - reverse_acceptors[a].first > max_intron_length)
				{
					broke_out = true;
					break;
				}
				
				// Don't allow splices within exons
				if (reverse_donor_exons[d] == reverse_acceptor_exons[a] &&
					reverse_donor_exons[d]->score < single_island_juncs_above)
					continue;
					
				
				reverse_splices++;
				contig_rev_splices++;
				total_splices++;
				
				Exon* l = reverse_acceptor_exons[a];
				Exon* r = reverse_donor_exons[d];
				size_t pos_in_l = reverse_acceptors[a].first - l->pos_in_ref - seq_key_len;
				size_t pos_in_r = reverse_donors[d].first - r->pos_in_ref + 2 + seq_key_len;
				
				PackedSplice p = combine_splice_halves(reverse_acceptors[a].second,
													   reverse_donors[d].second);
				lookup_splice_in_read_index(l->reference_ctg_id,
											p, 
											l, 
											r, 
											pos_in_l, 
											pos_in_r,
											true);
				 
			}
			if (broke_out)
				continue;
		}
		
		if (verbose)
		{
			fprintf(stderr, 
					"Antisense strand for %s, has %d total %d-mer splices\n", 
					ids_to_refctgs[j->first].c_str(),
					contig_rev_splices,
					2*seq_key_len);
		}
		
		if (verbose) 
		{
			//fprintf(stderr, "%d possible splices\n", total_splices);
			//fprintf(stderr, "%ul\n", (int)dummy);
		}
		
		if (verbose)
		{
			fprintf(stderr, "finished %lld splices in %d seconds ~(%f splices/sec)\n", 
					total_splices, (int)splice_map_timer.elapsed(), total_splices/(float)(splice_map_timer.elapsed()));
			
		}
	}
}

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
			fprintf(stderr,"%s\n",errmsg);
			print_usage();
			exit(1);
		}
		return (int32_t)l;
	}
	fprintf(stderr,"%s\n",errmsg);
	print_usage();
	exit(1);
	return -1;
}

void destroy_exons()
{
	for (EXONS_FOR_REFS::iterator j = exons_for_reference.begin();
		 j != exons_for_reference.end();
		 ++j)
	{
		for (vector<Exon*>::iterator e = j->second.begin();
			 e != j->second.end();
			 ++e)
		{
			Exon* exon = *e;
			delete exon;
		}
	}		
}

int main(int argc, char** argv)
{
	const char *short_options = "va:m:s:i:I:S:M:";
	int next_option; 
	do { 
		next_option = getopt(argc, argv, short_options);
		switch (next_option) {
	   		case 'v': /* verbose */
				verbose = true;
				break;
			case 'a':
				seq_key_len = parseInt(1, "-a arg must be at least 1");
				break;
			case 'm':
				max_span_mismatches = parseInt(0,"-m arg must be at least 0");
				break;
			case 's':
				seed_size = parseInt(1,"-s arg must be at least 1");
				break;
			case 'I':
				max_intron_length = parseInt(1,"-I arg must be at least 1");
				break;
			case 'i':
				min_intron_length = parseInt(1,"-i arg must be at least 1");
				break;
			case 'M':
				max_memory_megs = parseInt(50,"-M arg must be at least 50");
				break;
			case 'S':
				single_island_juncs_above = parseInt(0,"-S arg must be at least 0");
				break;
			case -1: /* Done with options. */
				break;
			default: 
				print_usage();
				return 1;
		}
	} while(next_option != -1);

	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	string exon_fasta_file_name = argv[optind++];
	
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	string exon_gff_file_name = argv[optind++];
	
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	
	string read_file_name = argv[optind++];
	
		
	FILE* read_file = NULL;
	if (read_file_name == "-")
		read_file = stdin;
	else
		read_file = fopen(read_file_name.c_str(), "r");
	
	if (read_file == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",read_file_name.c_str());
		exit(1);
	}
	  
	// Open the FASTA file containing the exon records scraped by cvg_islands
	FILE* exon_fasta = fopen(exon_fasta_file_name.c_str(), "r");
	if (exon_fasta == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				exon_fasta_file_name.c_str());
		exit(1);
	}
	 
	// Open the GFF file containing the exon metadata scraped by cvg_islands
	FILE* exon_gff = fopen(exon_gff_file_name.c_str(), "r");
	if (exon_gff == NULL)
	{
		fprintf(stderr, "Error: cannot open %s for reading\n",
				exon_gff_file_name.c_str());
		exit(1);
	}
	
	if (verbose)
		fprintf(stderr, "Files OK, starting TopHat mapping run...\n");
	
	//static const int MAX_MEM = 2000000000; //2 GB
	//static const int MAX_MEM = 300000000; //300 MB
	
	if (verbose)
		fprintf(stderr, "Loading exons\n");
	read_exons(exon_fasta, exon_gff);
	
		
	while (!feof(read_file))
	{ 
		
		index_reads(read_file, max_memory_megs);
		// Get the scraped exons
		map_possible_exon_juncs(read_file);
	}
	
	if (verbose)
	{
		fprintf(stderr,"TopHat examined %lld possible splices\n", total_splices);
	}
	
	destroy_exons();
	
	//match_reads_against_contigs(read_file);
	if (verbose)
	{
		fprintf(stderr, "TopHat examined %d IUM reads\n", reads_processed);
		fprintf(stderr, "\t %d total splice seeds\n", read_hits);
	}
	//while (1) sleep(60);
	return 0;
}

