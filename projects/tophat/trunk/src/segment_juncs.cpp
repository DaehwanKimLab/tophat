#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  segment_juncs.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/5/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <cassert>
#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
//#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <bitset>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "segments.h"
#include "reads.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"

using namespace seqan;
using namespace std;
using namespace __gnu_cxx;

// daehwan //geo
//#define B_DEBUG 1

void print_usage()
{
    fprintf(stderr, "Usage:   segment_juncs <ref.fa> <segment.juncs> <segment.insertions> <segment.deletions> <left_reads.fq> <left_reads.bwtout> <left_seg1.bwtout,...,segN.bwtout> [right_reads.fq right_reads.bwtout right_seg1.bwtout,...,right_segN.bwtout]\n");
}

// This is the maximum number of bowtie mismatches allower per segment hit
static const int num_bowtie_mismatches = 2;

static const int max_cov_juncs = 5000000;
// static const int max_cov_juncs = std::numeric_limits<int>::max();
static const int max_seg_juncs = 10000000;
int max_microexon_stretch = 2000;
int butterfly_overhang = 6;
int min_cov_length = 20;

void get_seqs(istream& ref_stream,
			  RefSequenceTable& rt,
			  bool keep_seqs = true,
			  bool strip_slash = false)
{    
    while(ref_stream.good() &&
          !ref_stream.eof())
    {
		RefSequenceTable::Sequence* ref_str = new RefSequenceTable::Sequence();
        string name;
        readMeta(ref_stream, name, Fasta());
		string::size_type space_pos = name.find_first_of(" \t\r");
		if (space_pos != string::npos)
		{
			name.resize(space_pos);
		}
		//fprintf(stderr, "\tLoading %s...", name.c_str());
		seqan::read(ref_stream, *ref_str, Fasta());
		//fprintf(stderr, "done\n");
        rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
		if (!keep_seqs)
			delete ref_str;
    }	
}


RefSeg seg_from_bowtie_hit(const BowtieHit& T)
{
  RefSeg r_seg(T.ref_id(), POINT_DIR_DONTCARE, T.antisense_align(), READ_DONTCARE, 0, 0);
	
  if (T.antisense_align())
    {
      r_seg.left = max(0, T.right() - 2);
      r_seg.right = T.right() + (T.right() - T.left() + 1); // num allowed bowtie mismatches 
      r_seg.points_where = POINT_DIR_RIGHT;
    }
  else
    {
      r_seg.left = max(0, T.left() - (T.right() - T.left() + 1));
      r_seg.right = T.left() + 2; // num allowed bowtie mismatches
      r_seg.points_where = POINT_DIR_LEFT;
    }
  
  return r_seg;
}

pair<RefSeg, RefSeg> segs_from_bowtie_hits(const BowtieHit& T,
					   const BowtieHit& H)
{
  pair<RefSeg, RefSeg> seg_pair;
  if (H.antisense_align() == false && 
      abs((H.right() + 1) - T.left()) < (int)max_segment_intron_length)
    {
      RefSeg left_seg(H.ref_id(), POINT_DIR_RIGHT, H.antisense_align(), READ_DONTCARE, 0, 0);
      left_seg.left = max(0, H.right() - 2);
      left_seg.right = H.right() + (H.right() - H.left() + 1); // num allowed bowtie mismatches 
      
      RefSeg right_seg(T.ref_id(), POINT_DIR_LEFT, T.antisense_align(), READ_DONTCARE, 0, 0);
      right_seg.left = max(0, T.left() - (T.right() - T.left() + 1));
      right_seg.right = T.left() + 2; // num allowed bowtie mismatches
      
      seg_pair = make_pair(left_seg, right_seg);
    }
  else if (H.antisense_align() == true &&
	   abs((T.right() + 1) - H.left()) < (int)max_segment_intron_length)
    {
      RefSeg left_seg(T.ref_id(), POINT_DIR_RIGHT, T.antisense_align(), READ_DONTCARE, 0, 0);
      left_seg.left = max(0, T.right() - 2);
      left_seg.right = T.right() + (T.right() - T.left() + 1); // num allowed bowtie mismatches 
      
      RefSeg right_seg(H.ref_id(), POINT_DIR_LEFT, H.antisense_align(), READ_DONTCARE, 0, 0);
      right_seg.left = max(0, H.left() - (H.right() - H.left() + 1));
      right_seg.right = H.left() + 2; // num allowed bowtie mismatches
      
      seg_pair = make_pair(left_seg, right_seg);
    }
  return seg_pair;
}

//static const size_t half_splice_mer_len = 6;
//static const size_t splice_mer_len = 2 * half_splice_mer_len;

struct MerExtension 
{
	static const int MAX_EXTENSION_BP = 14;
	uint32_t left_dna_str : 28;     // up to 14bp encoded in 2-bits-per-base
	uint8_t  left_ext_len : 4;  // how many bases in this extension
	
	uint32_t right_dna_str : 28;     // up to 14bp encoded in 2-bits-per-base
	uint8_t  right_ext_len : 4;  // how many bases in this extension
	
	MerExtension() : left_dna_str(0), left_ext_len(0), right_dna_str(0), right_ext_len(0) {}
	
	bool operator<(const MerExtension& rhs) const
	{
		if (left_dna_str != rhs.left_dna_str)
			return left_dna_str < rhs.left_dna_str;
		if (left_ext_len != rhs.left_ext_len)
			return left_ext_len < rhs.left_ext_len;
		if (right_dna_str != rhs.right_dna_str)
			return right_dna_str < rhs.right_dna_str;
		if (right_ext_len != rhs.right_ext_len)
			return right_ext_len < rhs.right_ext_len;
		return false;
	}
	
	bool operator==(const MerExtension& rhs) const
	{
		bool eq = left_dna_str == rhs.left_dna_str &&
				  right_dna_str == rhs.right_dna_str &&
				  left_ext_len == rhs.left_ext_len &&
				  right_ext_len == rhs.right_ext_len;

		return eq;
	}
};


/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
uint8_t charToDna5[] = {
/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*  48 */ 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
/*    A     C           G                    N */
/*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*             T */
/*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
/*    a     c           g                    n */
/* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/*             t */
/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};


typedef vector<vector<MerExtension> > MerExtensionTable;
typedef vector<uint32_t> MerExtensionCounts;

MerExtensionTable extensions;
MerExtensionCounts extension_counts;

uint64_t dna5str_to_idx(const string& str)
{
	uint64_t idx = 0;
	
	for (size_t i = 0; i < str.length(); ++i)
	{
		idx <<=2;
		char c = toupper(str[i]);
		idx |= (0x3 & charToDna5[(size_t)c]);
	}
	return idx;
}

uint64_t colorstr_to_idx(const string& str)
{
	uint64_t idx = 0;
	
	for (size_t i = 0; i < str.length(); ++i)
	{
		idx <<=2;
		char c = str[i];
		idx |= (0x3 & charToDna5[(size_t)c]);
	}
	return idx;
}

void store_read_extensions(MerExtensionTable& ext_table,
			   int seq_key_len,
			   int min_ext_len,
			   const string& seq,
			   bool use_precount_table)
{
	// h is will hold the 2-bit-per-base representation of the k-mer seeds for
	// this read.
	
	uint64_t seed = 0;
	bitset<256> left = 0;
	bitset<256> right = 0;
	const char* p = seq.c_str();
	
	unsigned int seq_len = (int)seq.length();
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
	
	// This loop will construct successive seed, along with 32-bit words 
	// containing the left and right remainders for each seed
	uint32_t i = 0;
	size_t new_hits = 0;
	do 
	{
		// How many base pairs exist in the right remainder beyond what we have
		// space for ?
		int extra_right_bp = ((int)seq.length() - 
			(i + 2 * seq_key_len)) - MerExtension::MAX_EXTENSION_BP;
		
		uint32_t hit_right = 0;
		if (extra_right_bp > 0)
		{
			//bitset<32> tmp_hit_right = (right >> (extra_right_bp << 1));
			hit_right = (uint32_t)(right >> (extra_right_bp << 1)).to_ulong();
		}
		else
		{
			hit_right = (uint32_t)right.to_ulong();
		}
		
		uint32_t hit_left = (uint32_t)((left << (256 - 32)) >> (256 - 32)).to_ulong();

		//size_t prev_cap = (*mer_table)[seed].capacity();
		//(*mer_table)[seed].push_back(ReadHit(hit_left, hit_right,i, read_num, reverse_complement));
		//cap_increase += ((*mer_table)[seed].capacity() - prev_cap) * sizeof (ReadHit);
		
		MerExtension ext;

		ext.right_dna_str = hit_right;
		ext.right_ext_len = min(seq_len - (2 * seq_key_len) - i,
							  (unsigned int)MerExtension::MAX_EXTENSION_BP);

		ext.left_dna_str = hit_left;
		ext.left_ext_len = min(i, (unsigned int)MerExtension::MAX_EXTENSION_BP);
		if (use_precount_table)
		{
			int curr_seed = --extension_counts[seed]; 
			if (curr_seed < 0 || curr_seed > (int)ext_table[seed].size())
			{
				fprintf(stderr, "Error: curr_seed is %d, max is %lu\n", curr_seed, (long unsigned int)ext_table[seed].size());
			}
			
			ext_table[seed][curr_seed] = ext;
		}
		else
		{
			ext_table[seed].push_back(ext);
		}
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
			right.set(((right_len - 1) << 1), 0);
			right.set(((right_len - 1)  << 1) + 1, 0);
		}
		++i;
		
	}while(i <= (size_t)(seq_end - seq.c_str()) - (2 * seq_key_len));
	return;
}

void count_read_extensions(MerExtensionCounts& ext_counts,
						   int seq_key_len,
						   int min_ext_len,
						   const string& seq)
{
	// h is will hold the 2-bit-per-base representation of the k-mer seeds for
	// this read.
	
	uint64_t seed = 0;
	bitset<256> left = 0;
	bitset<256> right = 0;
	const char* p = seq.c_str();
	
	unsigned int seq_len = (int)seq.length();
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
	
	// This loop will construct successive seed, along with 32-bit words 
	// containing the left and right remainders for each seed
	uint32_t i = 0;
	size_t new_hits = 0;
	do 
	{
		ext_counts[seed]++;

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
			right.set(((right_len - 1) << 1), 0);
			right.set(((right_len - 1)  << 1) + 1, 0);
		}
		++i;
		
	}while(i <= (size_t)(seq_end - seq.c_str()) - (2 * seq_key_len));
	return;
}

//void count_read_mers(FILE* reads_file, size_t half_splice_mer_len)
void count_read_mers(FZPipe& reads_file, size_t half_splice_mer_len)
{
	Read read;
	size_t splice_mer_len = 2 * half_splice_mer_len;
	size_t mer_table_size = 1 << ((splice_mer_len)<<1);
	extension_counts.resize(mer_table_size);
    FLineReader fr(reads_file);
	//while(!feof(reads_file))
    while (!fr.isEof())
	{
		read.clear();
		
		// Get the next read from the file
		if (!next_fasta_record(fr, read.name, read.seq, reads_format))
		  break;
		if (reads_format == FASTQ)
		{
		  if (!next_fastq_record(fr, read.seq, read.alt_name, read.qual, reads_format))
		    break;
		}

		if (color)
		  // erase the primer and the adjacent color
		  read.seq.erase(0, 2);

		if (read.seq.size() > 32)
			read.seq.resize(32);
		
		count_read_extensions(extension_counts,
							  half_splice_mer_len,
							  half_splice_mer_len,
							  read.seq);
	}	
	
	//rewind(reads_file);
    reads_file.rewind();
}

void compact_extension_table()
{
	for (size_t i = 0; i < extensions.size(); ++i)
	{
		vector<MerExtension>& exts = extensions[i];
		sort(exts.begin(), exts.end());
		vector<MerExtension>::iterator new_end  = unique(exts.begin(), exts.end());
		exts.erase(new_end, exts.end());
		vector<MerExtension>(exts).swap(exts);
	}
}	

void prune_extension_table(uint8_t max_extension_bp)
{
	uint32_t mask =  ~(0xFFFFFFFFuLL << (max_extension_bp << 1));
	
	for (size_t i = 0; i < extensions.size(); ++i)
	{
		vector<MerExtension>& exts = extensions[i];
		for (size_t j = 0; j < exts.size(); ++j)
		{
			MerExtension& ex = exts[j];
			if (ex.left_ext_len > max_extension_bp)
			{
				ex.left_ext_len = max_extension_bp;
				ex.left_dna_str &= mask;
			}
			
			if (ex.right_ext_len > max_extension_bp)
			{
			  ex.right_dna_str >>= ((ex.right_ext_len - max_extension_bp) << 1);
			  ex.right_ext_len = max_extension_bp; 
			}
		}
	}
}

//void store_read_mers(FILE* reads_file, size_t half_splice_mer_len)
void store_read_mers(FZPipe& reads_file, size_t half_splice_mer_len)
{
	Read read;
	size_t splice_mer_len = 2 * half_splice_mer_len;
	
	size_t mer_table_size = 1 << ((splice_mer_len)<<1);
	extensions.resize(mer_table_size);
	
	size_t num_indexed_reads = 0;
	FLineReader fr(reads_file);
	//while(!feof(reads_file))
	while(!fr.isEof())
	{
		read.clear();
		
		// Get the next read from the file
		if (!next_fasta_record(fr, read.name, read.seq, reads_format))
		  break;
		if (reads_format == FASTQ)
		{
		  if (!next_fastq_record(fr, read.seq, read.alt_name, read.qual, reads_format))
		    break;
		}

		if (color)
		  // erase the primer and the adjacent color
		  read.seq.erase(0, 2);
		
		if (read.seq.size() > 32)
			read.seq.resize(32);
		
		store_read_extensions(extensions,
				      half_splice_mer_len,
				      half_splice_mer_len,
				      read.seq, 
				      true);
		
		// Do NOT index the reverse of the reads
		
		num_indexed_reads++;
		if (num_indexed_reads % 1000000 == 0)
		{
			//fprintf(stderr, "Indexed %lu reads, compacting extension table\n", num_indexed_reads);
			//compact_extension_table();
		}
	}	
	
	//fprintf(stderr, "Indexed %lu reads, compacting extension table\n", num_indexed_reads)
	
	uint64_t num_extensions = 0;
	for (size_t i = 0; i < extensions.size(); ++i)
	{
		num_extensions += extensions[i].size();
	}
	//fprintf (stderr, "Total extensions: %lu\n", (long unsigned int)num_extensions);
  reads_file.rewind();
}

//void index_read_mers(vector<FILE*> reads_files,
void index_read_mers(vector<FZPipe>& reads_files,
					 size_t half_splice_mer_len)
{
	extensions.clear();
	for (size_t i = 0; i < reads_files.size(); ++i)
	{
		count_read_mers(reads_files[i], half_splice_mer_len);
	}
	
	extensions.resize(extension_counts.size());
	
	for (size_t i = 0; i < extension_counts.size(); ++i)
	{
		extensions[i].resize(extension_counts[i]);
	}
	
	for (size_t i = 0; i < reads_files.size(); ++i)
	{
		store_read_mers(reads_files[i], half_splice_mer_len);
	}
	
	compact_extension_table();
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
		if (match_chars == -1 || match_chars >= len)
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

uint64_t rc_dna_str(uint64_t dna_str)
{
	dna_str = ~dna_str;
	uint64_t rc = 0;
	for (int i = 0; i < 32; i++)
	{
		rc <<= 2;
		rc |= (dna_str & 0x3);
		dna_str >>= 2;
	}
	return rc;
}

uint64_t rc_color_str(uint64_t color_str)
{
  uint64_t rc = 0;
  for (int i = 0; i < 32; ++i)
    {
      rc <<= 2;
      rc |= (color_str & 0x3);
      color_str >>= 2;
    }
  return rc;
}

struct DnaSpliceStrings
{
  DnaSpliceStrings(uint64_t f, uint64_t r) : fwd_string(f), rev_string(r), first_in_string('N'), last_in_string('N') {}
  uint64_t fwd_string;
  uint64_t rev_string;

  // for color-space purposes
  char first_in_string;
  char last_in_string;
  
  bool operator<(const DnaSpliceStrings& rhs) const
  {
    if (fwd_string != rhs.fwd_string)
      return fwd_string < rhs.fwd_string;
    if (rev_string != rhs.rev_string)
      return rev_string < rhs.rev_string;
    return false;
  }
  
  bool operator==(const DnaSpliceStrings& rhs) const
  {
    return fwd_string == rhs.fwd_string && rev_string == rhs.rev_string;
  }
};

struct IntronMotifs
{
  IntronMotifs(uint32_t rid) : ref_id(rid) {}
  uint32_t ref_id;
	
  vector<pair<size_t, DnaSpliceStrings> > fwd_donors;
  vector<pair<size_t, DnaSpliceStrings> > fwd_acceptors;
  vector<pair<size_t, DnaSpliceStrings> > rev_donors;
  vector<pair<size_t, DnaSpliceStrings> > rev_acceptors;	
  
  void unique(vector<pair<size_t, DnaSpliceStrings> >& f)
  {
    sort(f.begin(), f.end());
    vector<pair<size_t, DnaSpliceStrings> >::iterator i = std::unique(f.begin(), f.end());
    f.erase(i, f.end());
  }
  
  void unique()
  {
    unique(fwd_donors);
    unique(fwd_acceptors);
    unique(rev_donors);
    unique(rev_acceptors);
  }
  
  void attach_mers(RefSequenceTable::Sequence& ref_str)
  {
    attach_upstream_mers(ref_str, fwd_donors);
    attach_upstream_mers(ref_str, rev_acceptors);
    
    attach_downstream_mers(ref_str, rev_donors);
    attach_downstream_mers(ref_str, fwd_acceptors);		
  }
  
  void attach_upstream_mers(RefSequenceTable::Sequence& ref_str,
			    vector<pair<size_t, DnaSpliceStrings> >& dinucs)
  {
    for (size_t i = 0; i < dinucs.size(); ++i)
      {
	size_t pos = dinucs[i].first;
	int half_splice_mer_len = 32;
	
	if (color)
	  {
	    if (pos <= (size_t)half_splice_mer_len+1 || pos >= length(ref_str))
	      continue; 
	    
	    Dna5String seg_str = seqan::infixWithLength(ref_str,
							pos - half_splice_mer_len - 1,
							half_splice_mer_len + 1);
	    stringstream ss(stringstream::in | stringstream::out);
	    string s;
	    ss << seg_str;
	    ss >> s;
	    
	    string col_seg_str = convert_bp_to_color(s, true);
	    uint64_t idx = colorstr_to_idx(col_seg_str);
	    
	    dinucs[i].second.fwd_string = idx;
	    dinucs[i].second.rev_string = rc_color_str(idx);
	    dinucs[i].second.first_in_string = s[1];
	    dinucs[i].second.last_in_string = s[half_splice_mer_len];
	  }
	else
	  {
	    if (pos <= (size_t)half_splice_mer_len || pos >= length(ref_str))
	      continue; 
	    
	    Dna5String seg_str = seqan::infixWithLength(ref_str,
							pos - half_splice_mer_len,
							half_splice_mer_len);
	    
	    stringstream ss(stringstream::in | stringstream::out);
	    string s;
	    ss << seg_str;
	    ss >> s;
	    uint64_t idx = dna5str_to_idx(s);
	    dinucs[i].second.fwd_string = idx;
	    dinucs[i].second.rev_string = rc_dna_str(idx);
	  }
      }
  }
  
  
  void attach_downstream_mers(RefSequenceTable::Sequence& ref_str,
			      vector<pair<size_t, DnaSpliceStrings> >& dinucs)
  {
    for (size_t i = 0; i < dinucs.size(); ++i)
      {
	size_t pos = dinucs[i].first;
	
	int half_splice_mer_len = 32;
	
	if (pos + 2 + half_splice_mer_len >= length(ref_str))
	  continue; 
	
	if (color)
	  {
	    Dna5String seg_str = seqan::infixWithLength(ref_str,
							pos + 2 - 1,
							half_splice_mer_len + 1);
	    stringstream ss(stringstream::in | stringstream::out);
	    string s;
	    ss << seg_str;
	    ss >> s;
	    
	    string col_seg_str = convert_bp_to_color(s, true);
	    uint64_t idx = colorstr_to_idx(col_seg_str);
	    
	    dinucs[i].second.fwd_string = idx;
	    dinucs[i].second.rev_string = rc_color_str(idx);
	    dinucs[i].second.first_in_string = s[1];
	    dinucs[i].second.last_in_string = s[half_splice_mer_len];
	  }
	else
	  {
	    Dna5String seg_str = seqan::infixWithLength(ref_str,
							pos + 2,
							half_splice_mer_len);
	    
	    stringstream ss(stringstream::in | stringstream::out);
	    string s;
	    ss << seg_str;
	    ss >> s;
	    uint64_t idx = dna5str_to_idx(s);
	    
	    dinucs[i].second.fwd_string = idx;
	    dinucs[i].second.rev_string = rc_dna_str(idx);
	  }
      }
  }
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

static inline std::string u32ToDna(uint32_t a, int len) {
	char buf[17];
	assert(len <= 16);
	for(int i = 0; i < len; i++) {
		buf[len-i-1] = "ACGT"[a & 3];
		a >>= 2;
	}
	buf[len] = '\0';
	return std::string(buf);
}


PackedSpliceHalf pack_left_splice_half(const string& seq, 
									   uint32_t pos_in_l, 
									   unsigned int seq_key_len)
{
	const char* l = seq.c_str(); 
	l += pos_in_l;
	
	const char* left_end = l;
	l -= 16;
	
	assert (l + seq_key_len < seq.c_str() + seq.length());
	
	PackedSpliceHalf packed_half = make_pair(0u,0u);
	
	// Pack up to 32 bits (16 bases) of sequence into left
	if (l < seq.c_str())
		l = seq.c_str();
	while (l < left_end)
	{
		packed_half.first <<= 2;
		packed_half.first |= (0x3 & charToDna5[(size_t)*l]);
		++l;
	}
	
	// Pack up the seed bits
	for (unsigned int i = 0; i < seq_key_len; ++i)
	{
		packed_half.second <<= 2;
		packed_half.second |= (0x3 & charToDna5[(size_t)*(l + i)]);
	}
	return packed_half;
}


PackedSpliceHalf pack_right_splice_half(const string& seq, 
										uint32_t pos, 
										unsigned int seq_key_len)
{
	const char* r = seq.c_str(); 
	r += pos - seq_key_len;
	
	PackedSpliceHalf packed_half = make_pair(0u,0u);
	
	// Pack the seed bits
	for (unsigned int i = 0; i < seq_key_len; ++i)
	{
		packed_half.second <<= 2;
		packed_half.second |= (0x3 & charToDna5[(size_t)*(r + i)]);
	}
	
	r += seq_key_len;
	// Now pack 32 bits (16 bases) of sequence into left
	const char* right_end = r + 16;
	
	if ((size_t)(right_end - seq.c_str()) > seq.length())
		right_end = seq.c_str() + seq.length();
	while (r < right_end)
	{
		packed_half.first <<= 2;
		packed_half.first |= (0x3 & charToDna5[(size_t)*r]);
		++r;
	}	
	
	return packed_half;
}

PackedSplice combine_splice_halves(const PackedSpliceHalf& left_half, 
								   const PackedSpliceHalf& right_half,
								   int seq_key_len)
{
	uint64_t seed = left_half.second << (seq_key_len << 1) | right_half.second;
	return PackedSplice(left_half.first,seed, right_half.first);
}

PackedSplice pack_splice(const string& seq,
						 int l_pos_in_seq,
						 int r_pos_in_seq,
						 unsigned int seq_key_len)
{
	const char* l = seq.c_str(); // l points to beginning of left exon sequence
	l += l_pos_in_seq; 
	
	assert (l + seq_key_len < seq.c_str() + seq.length());
	
	const char* r = seq.c_str(); // r points to beginning of right exon sequence
	r += r_pos_in_seq - seq_key_len; 
	//r += 2; // r points to right side of junction;
	
	uint64_t seed = 0;
	uint64_t left = 0;
	uint64_t right = 0;
	
	// Pack up the seed bits
	for (unsigned int i = 0; i < seq_key_len; ++i)
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*(l + i)]);
	}
	
	for (unsigned int i = 0; i < seq_key_len; ++i)
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*(r + i)]);
	}
	
	// Now pack 32 bits (16 bases) of sequence into left
	const char* left_end = l;
	l -= 16;
	if (l < seq.c_str())
		l = seq.c_str();
	while (l < left_end)
	{
		left <<= 2;
		left |= (0x3 & charToDna5[(size_t)*l]);
		++l;
	}
	
	r += seq_key_len;
	// Now pack 32 bits (16 bases) of sequence into left
	const char* right_end = r + 16;
	
	if ((size_t)(right_end - seq.c_str()) > seq.length())
		right_end = seq.c_str() + seq.length();
	while (r < right_end)
	{
		right <<= 2;
		right |= (0x3 & charToDna5[(size_t)*r]);
		++r;
	}
	
	return PackedSplice((uint32_t)left,seed,(uint32_t)right);
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

// A MerTable maps k-mers to hits in indexed reads.  See the comment for 
// mer_table
typedef vector<ReadHit> ReadHitList;
typedef vector<ReadHitList> MerTable;

size_t index_read(MerTable* mer_table, 
				  int seq_key_len,
				  const string& seq, 
				  unsigned int read_num,
				  bool reverse_complement,
				  vector<uint64_t>& seeds)
{
	// h is will hold the 2-bit-per-base representation of the k-mer seeds for
	// this read.
	
	uint64_t seed = 0;
	bitset<256> left = 0;
	bitset<256> right = 0;
	const char* p = seq.c_str();
	
	unsigned int seq_len = (int)seq.length();
	const char* seq_end = p + seq_len;
	
	// Build the first seed
	while (p < seq.c_str() + (2 * seq_key_len))
	{
		seed <<= 2;
		seed |= (0x3 & charToDna5[(size_t)*p]);
		++p;
	}
	
	seeds.push_back(seed);
	
	// Build the rest of them with a sliding window, adding successive bases
	// to the "right-remainder" word.  
	while (p < seq_end)
	{
		
		right <<= 2;
		right |= (0x3 & charToDna5[(size_t)*p]);		
		++p;
	}	
	
	size_t cap_increase = 0;
	
	// At this point, seed contains the 5'-most 2*min_anchor_len bases of the 
	// read, and right contains everthing else on the 3' end.
	
	// This loop will construct successive seed, along with 32-bit words 
	// containing the left and right remainders for each seed
	uint32_t i = 0;
	size_t new_hits = 0;
	do 
	{
		// Let's not make an out-of-bounds write, if this fails the global 
		// mer_table is too small
		assert (!mer_table || seed < mer_table->size());
		
		// How many base pairs exist in the right remainder beyond what we have
		// space for ?
		int extra_right_bp = ((int)seq.length() - (i + 2 * seq_key_len)) - 16;
		
		uint32_t hit_right;
		if (extra_right_bp > 0)
		{
			//bitset<32> tmp_hit_right = (right >> (extra_right_bp << 1));
			hit_right = (uint32_t)(right >> (extra_right_bp << 1)).to_ulong();
		}
		else
		{
			hit_right = (uint32_t)right.to_ulong();
		}
		
		uint32_t hit_left = (uint32_t)((left << (256 - 32)) >> (256 - 32)).to_ulong();
		
		if (mer_table)
		{
			size_t prev_cap = (*mer_table)[seed].capacity();
			(*mer_table)[seed].push_back(ReadHit(hit_left, hit_right,i, read_num, reverse_complement));
			cap_increase += ((*mer_table)[seed].capacity() - prev_cap) * sizeof (ReadHit);
		}
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
		
		seeds.push_back(seed);
		
		//Now remove that leftmost base of the right remainder
		//right &=  ~(0xFFFFFFFFFFFFFFFFuLL << (right_len - 1 << 1));
		if (right_len)
		{
			right.set(((right_len - 1) << 1), 0);
			right.set(((right_len - 1)  << 1) + 1, 0);
		}
		++i;
		
	}while(i <= (size_t)(seq_end - seq.c_str()) - (2 * seq_key_len));
	return cap_increase;
}



struct SeedExtension
{
	SeedExtension(int lp,
				  int rp,
				  int rd_p,
				  int le,
				  int re,
				  int mm) : 
	l_pos_in_ref(lp),
	r_pos_in_ref(rp),
	read_pos(rd_p),
	left_extent(le),
	right_extent(re),
	mismatches(mm){}
	
	int l_pos_in_ref;
	int r_pos_in_ref;
	int read_pos;
	int left_extent;
	int right_extent;
	int mismatches;
};

pair<string::size_type, int> left_extend(const string& ref,
										 const string& read,
										 int ref_pos,
										 int read_pos,
										 int num_mismatches)
{
	string::size_type ext = 0;
	int mm_encountered = 0;
	while(ref_pos >= 0 &&
		  read_pos >= 0)
	{
		//char ref_char = ref[ref_pos];
		//char read_char = read[read_pos];
		if (ref[ref_pos] != read[read_pos])
		{
			if (mm_encountered + 1 > num_mismatches)
				return make_pair(ext, mm_encountered);
			mm_encountered++;
		}
		ext++;
		
		--ref_pos;
		--read_pos;
	}
	
	return make_pair(ext, mm_encountered);
}

pair<string::size_type, int>  right_extend(const string& ref,
										   const string& read,
										   int ref_pos,
										   int read_pos,
										   int num_mismatches)
{
	string::size_type ext = 0;
	int mm_encountered = 0;
	while(ref_pos < (int)ref.size() &&
		  read_pos < (int)read.size())
	{
		if (ref[ref_pos] != read[read_pos])
		{
			if (mm_encountered + 1 > num_mismatches)
				return make_pair(ext, mm_encountered);
			mm_encountered++;
		}
		
		ext++;
		
		++ref_pos;
		++read_pos;
	}
	
	return make_pair(ext, mm_encountered);
}


void extend_from_seeds(vector<SeedExtension>& extensions,
					   const PackedSplice& p,
					   const MerTable& mer_table,
					   const string& ref, 
					   const string& read,
					   size_t l_pos_in_ref,
					   size_t r_pos_in_ref,
					   int seq_key_len)
{	
	assert(p.seed < mer_table.size());
	const ReadHitList& hl = mer_table[p.seed];
	
	for (size_t hit = 0; hit < hl.size(); ++hit)
	{
		const ReadHit& rh = hl[hit];
		uint32_t pos = rh.pos;
		
		pair<string::size_type, int> left_extension;
		pair<string::size_type, int> right_extension;
		left_extension = left_extend(ref, read, l_pos_in_ref - seq_key_len + 1, pos, 2);
		right_extension = right_extend(ref, read, r_pos_in_ref + seq_key_len, pos + 2 * seq_key_len, 2);
		extensions.push_back(SeedExtension(l_pos_in_ref,
										   r_pos_in_ref,
										   pos + seq_key_len, 
										   left_extension.first,
										   right_extension.first,
										   left_extension.second + right_extension.second));
		
	}
}

typedef pair<size_t, PackedSpliceHalf> SpliceHalf;

void get_seed_extensions(const string& ref, 
						 const string& read,
						 int seq_key_len,
						 MerTable& mer_table,
						 vector<SpliceHalf>& donors,
						 vector<SpliceHalf>& acceptors,
						 vector<SeedExtension>& extensions)
{
	for (size_t d = 0; d < donors.size(); ++d)
	{
		bool broke_out = false;
		
		// start pos is a lower bound on downstream acceptor positions
		// to consider
		size_t start_pos = donors[d].first + min_report_intron_length;
		
		SpliceHalf dummy = make_pair(start_pos,PackedSpliceHalf());
		vector<SpliceHalf>::iterator lb = upper_bound(acceptors.begin(), 
													  acceptors.end(),
													  dummy);
		
		if (lb == acceptors.end())
			break;
		for (size_t a = lb - acceptors.begin(); 
			 a < acceptors.size(); 
			 ++a)
		{
			if (acceptors[a].first - donors[d].first > (size_t)max_microexon_stretch)
			{
				broke_out = true;
				break;
			}
			
			size_t l_pos_in_ref = donors[d].first - 1;
			size_t r_pos_in_ref = acceptors[a].first + 2;
			PackedSplice p = combine_splice_halves(donors[d].second,
												   acceptors[a].second,
												   seq_key_len);
			extend_from_seeds(extensions,
							  p,
							  mer_table,
							  ref,
							  read,
							  l_pos_in_ref,
							  r_pos_in_ref,
							  seq_key_len);
		}
		if (broke_out)
			continue;
	}
	
}

void hits_from_seed_extension(uint32_t ref_id,
			      int ref_offset,
			      uint64_t insert_id,
			      bool antisense,
			      vector<SeedExtension>& extensions,
			      vector<BowtieHit>& hits_out,
			      int left_read_edge,
			      int right_read_edge,
			      int seq_key_len)
{
	for (size_t i = 0; i < extensions.size(); ++i)
	{
		SeedExtension& s = extensions[i];
		if (s.read_pos >= right_read_edge || 
			s.read_pos < left_read_edge)
			continue;
		if (s.read_pos - seq_key_len - s.left_extent <= left_read_edge &&
			s.read_pos + seq_key_len + s.right_extent >= right_read_edge  && s.mismatches <= 2 )
		{
			vector<CigarOp> cigar;
			int off_adjust;
			if (antisense)
			{
				CigarOp m1 = CigarOp(MATCH, s.read_pos - left_read_edge);
				CigarOp skip = CigarOp(REF_SKIP, s.r_pos_in_ref - s.l_pos_in_ref);
				CigarOp m2 = CigarOp(MATCH, right_read_edge - s.read_pos);
				cigar.push_back(m1);
				cigar.push_back(skip);
				cigar.push_back(m2); 
				off_adjust = m1.length;
			}
			else
			{
				CigarOp m1 = CigarOp(MATCH, s.read_pos - left_read_edge + 1);
				CigarOp skip = CigarOp(REF_SKIP, s.r_pos_in_ref - s.l_pos_in_ref);
				CigarOp m2 = CigarOp(MATCH, right_read_edge - s.read_pos - 1);
				cigar.push_back(m1);
				cigar.push_back(skip);
				cigar.push_back(m2); 
				off_adjust = m1.length;
			}
			// daehwan - check this
			bool end = false;
			BowtieHit bh(ref_id, 
				     insert_id,
				     ref_offset + s.l_pos_in_ref - off_adjust + 1,
				     cigar,
				     antisense,
				     false,
				     s.mismatches,
				     0,
				     end);
			hits_out.push_back(bh);
		}
	}
}

void align(uint32_t ref_id,
		   uint64_t insert_id,
		   bool antisense,
		   const string& ref, 
		   const string& read, 
		   int ref_offset,
		   int left_read_edge,
		   int right_read_edge,
		   MerTable& mer_table,
		   int seq_key_len,
		   vector<BowtieHit>& hits_out)
{
	// Reserve an entry for each k-mer we might see
	size_t mer_table_size = 1 << ((seq_key_len << 1)<<1);
	
	mer_table.resize(mer_table_size);
	
	vector<uint64_t> seeds;
	index_read(&mer_table, seq_key_len, read, 0, antisense, seeds);
	
	vector<SpliceHalf> forward_donors;
	vector<SpliceHalf> forward_acceptors;
	vector<SpliceHalf> reverse_donors;
	vector<SpliceHalf> reverse_acceptors;
	
	const string& seq = ref;
	
	unsigned int pos = 0;
	
	
	for (size_t z = seq_key_len + 1; z < seq.length() - seq_key_len - 2; ++z)
	{
		char l = seq[z - 1];
		char r = seq[z];
		if (l == 'G' && r == 'T')
		{
			size_t donor_pos = pos + z - 1;
			size_t s = donor_pos - seq_key_len;
			PackedSpliceHalf p = pack_left_splice_half(seq, s, seq_key_len);
			forward_donors.push_back(make_pair(donor_pos,p));
		}
		if (l == 'A' && r == 'G')
		{
			size_t acceptor_pos = pos + z - 1;
			size_t s = acceptor_pos + 2 + seq_key_len;
			PackedSpliceHalf p = pack_right_splice_half(seq, s, seq_key_len);
			forward_acceptors.push_back(make_pair(acceptor_pos,p));
		}
		if (l == 'C' && r == 'T')
		{
			size_t acceptor_pos = pos + z - 1;
			size_t s = acceptor_pos - seq_key_len;
			PackedSpliceHalf p = pack_left_splice_half(seq, s, seq_key_len);
			reverse_acceptors.push_back(make_pair(pos + z - 1,p));
		}
		if (l == 'A' && r == 'C')
		{
			size_t donor_pos = pos + z - 1;
			size_t s = donor_pos + 2 + seq_key_len;
			PackedSpliceHalf p = pack_right_splice_half(seq, s, seq_key_len);
			reverse_donors.push_back(make_pair(donor_pos,p));
		}
	}
	
	vector<SeedExtension> fwd_extensions;
	get_seed_extensions(seq, 
						read, 
						seq_key_len, 
						mer_table, 
						forward_donors, 
						forward_acceptors, 
						fwd_extensions);
	
	hits_from_seed_extension(ref_id, 
							 ref_offset, 
							 insert_id, 
							 antisense,
							 fwd_extensions, 
							 hits_out, 
							 left_read_edge, 
							 right_read_edge, 
							 seq_key_len);
	
	//fprintf(stderr, "Found %d seed hits\n", fwd_extensions.size());
	
	vector<SeedExtension> rev_extensions;
	get_seed_extensions(seq, 
						read,
						seq_key_len, 
						mer_table, 
						reverse_donors, 
						reverse_acceptors, 
						rev_extensions);
	
	hits_from_seed_extension(ref_id, 
							 ref_offset, 
							 insert_id, 
							 antisense,
							 rev_extensions, 
							 hits_out, 
							 left_read_edge, 
							 right_read_edge, 
							 seq_key_len);
	
	for (size_t i = 0; i < seeds.size(); ++i)
		mer_table[seeds[i]].clear();
}

int extension_mismatches = 0;


bool left_extendable_junction(uint64_t upstream_dna_str,
							  size_t key,
							  size_t splice_mer_len,
							  size_t min_ext_len)
{
	vector<MerExtension>& exts = extensions[key];
	for (size_t i = 0; i < exts.size(); ++i)
	{
		const MerExtension& ext = exts[i];
		if (ext.left_ext_len < min_ext_len)
			continue;
		uint64_t upstream = upstream_dna_str & ~(0xFFFFFFFFFFFFFFFFull << (ext.left_ext_len << 1));
		int mism = mismatching_bases(ext.left_dna_str, upstream, ext.left_ext_len, extension_mismatches);
		if (mism <= extension_mismatches)
			return true;
	}
	return false;
}

bool right_extendable_junction(uint64_t downstream_dna_str,
							   size_t key,
							   size_t splice_mer_len,
							   size_t min_ext_len)
{
	vector<MerExtension>& exts = extensions[key];
	for (size_t i = 0; i < exts.size(); ++i)
	{
		const MerExtension& ext = exts[i];
		if (ext.right_ext_len < min_ext_len)
			continue;
		uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull >> (ext.right_ext_len << 1));
		uint64_t downstream = downstream_dna_str & mask;
		downstream >>= ((32 - ext.right_ext_len) << 1);
		int mism = mismatching_bases(ext.right_dna_str, downstream, ext.right_ext_len, extension_mismatches);

		if (mism <= extension_mismatches)
			return true;
	}
	return false;
}

uint32_t junction_key(uint64_t upstream_dna_str,
					  uint64_t downstream_dna_str,
					  size_t splice_mer_len)
{
	uint64_t upstream_mask = ~(0xFFFFFFFFFFFFFFFFull << (splice_mer_len));
	uint64_t upstream_key_half = upstream_dna_str & upstream_mask;
	uint64_t downstream_mask =  ~(0xFFFFFFFFFFFFFFFFull >> (splice_mer_len));
	uint64_t downstream_key_half = (downstream_dna_str & downstream_mask) >> (64 - splice_mer_len);
	uint32_t key = ((uint32_t)upstream_key_half << splice_mer_len) | (uint32_t)downstream_key_half;
	return key;
}

bool extendable_junction(uint64_t upstream_dna_str,
			 uint64_t downstream_dna_str,
			 size_t splice_mer_len,
			 size_t min_ext_len,
			 bool reverse,
			 char last_in_upstream = 'N',
			 char first_in_downstream = 'N')
{
  if (color)
    {
      string two_bp; two_bp.push_back(last_in_upstream); two_bp.push_back(first_in_downstream);
      string color = convert_bp_to_color(two_bp, true);
      char num = (color[0] - '0')  & 0x3;

      if (reverse)
	{
	  upstream_dna_str = (upstream_dna_str >> 2) << 2;
	  upstream_dna_str |= (uint64_t)num;
	}
      else
	{
	  downstream_dna_str = (downstream_dna_str << 2) >> 2;
	  downstream_dna_str |= ((uint64_t)num << 62);
	}
    }
  
  uint32_t key = junction_key(upstream_dna_str, 
			      downstream_dna_str,
			      splice_mer_len);

  upstream_dna_str >>= splice_mer_len;
  downstream_dna_str <<= splice_mer_len;
  
  bool extendable = (left_extendable_junction(upstream_dna_str,
					      key, splice_mer_len, min_ext_len) || 
		     right_extendable_junction(downstream_dna_str,
					       key, splice_mer_len, min_ext_len));
  return extendable;
}

typedef std::set<Junction, skip_count_lt> PotentialJuncs;

struct RecordExtendableJuncs
{
  void record(uint32_t ref_id,
	      const vector<pair<size_t, DnaSpliceStrings> >& left_sites,
	      const vector<pair<size_t, DnaSpliceStrings> >& right_sites,
	      bool antisense,
	      PotentialJuncs& juncs,
	      int min_intron,
	      int max_intron,
	      size_t max_juncs,
	      size_t half_splice_mer_len)
  {
    
    size_t splice_mer_len = 2 * half_splice_mer_len;
    
    size_t curr_R = 0;
    for (size_t L = 0; L < left_sites.size(); ++L)
      {
	while (curr_R < right_sites.size() && 
	       right_sites[curr_R].first < left_sites[L].first + min_intron)
	  {
	    curr_R++;
	  }
	
	size_t left_pos = left_sites[L].first;
	size_t max_right_pos = left_pos + max_intron;
	for (size_t R = curr_R; 
	     R < right_sites.size() && right_sites[R].first < max_right_pos; ++R)
	  {
	    uint64_t upstream_dna_str = left_sites[L].second.fwd_string;
	    char last_in_upstream = left_sites[L].second.last_in_string;
	    uint64_t downstream_dna_str = right_sites[R].second.fwd_string;
	    char first_in_downstream = right_sites[R].second.first_in_string;
	    uint64_t rc_upstream_dna_str = left_sites[L].second.rev_string;
	    uint64_t rc_downstream_dna_str = right_sites[R].second.rev_string;

	    if (extendable_junction(upstream_dna_str,
				    downstream_dna_str, splice_mer_len, 7, false,
				    last_in_upstream, first_in_downstream) ||
		extendable_junction(rc_downstream_dna_str,
				    rc_upstream_dna_str, splice_mer_len, 7, true,
				    last_in_upstream, first_in_downstream))
	      {
		juncs.insert(Junction(ref_id,
				      left_sites[L].first - 1,
				      right_sites[R].first + 2,
				      antisense,
				      R - curr_R));
	      }
	    if (juncs.size() > max_juncs)
		juncs.erase(*(juncs.rbegin()));
	      }			
	  }
      }	
};

struct RecordAllJuncs
{
  void record(uint32_t ref_id,
	      const vector<pair<size_t, DnaSpliceStrings> >& left_sites,
	      const vector<pair<size_t, DnaSpliceStrings> >& right_sites,
	      bool antisense,
	      PotentialJuncs& juncs,
	      int min_intron,
	      int max_intron,
	      size_t max_juncs,
	      size_t half_splice_mer_len)
  {
    size_t curr_R = 0;
    for (size_t L = 0; L < left_sites.size(); ++L)
      {
	while (curr_R < right_sites.size() && 
	       right_sites[curr_R].first < left_sites[L].first + min_intron)
	  {
	    curr_R++;
	  }
	
	size_t left_pos = left_sites[L].first;
	size_t max_right_pos = left_pos + max_intron;
	for (size_t R = curr_R; 
	     R < right_sites.size() && right_sites[R].first < max_right_pos; ++R)
	  {
	    Junction j(ref_id,
		       left_sites[L].first - 1,
		       right_sites[R].first + 2,
		       antisense,
		       R - curr_R);
	    
	    juncs.insert(j);
	    
	    if (juncs.size() > max_juncs)
	      juncs.erase(*(juncs.rbegin()));
	  }
      }
  }
};

struct RecordSegmentJuncs
{
  void record(uint32_t ref_id,
	      const vector<pair<size_t, DnaSpliceStrings> >& left_sites,
	      const vector<pair<size_t, DnaSpliceStrings> >& right_sites,
	      bool antisense,
	      PotentialJuncs& juncs,
	      int min_intron,
	      int max_intron,
	      size_t max_juncs,
	      size_t half_splice_mer_len)
  {
    if (left_sites.size() != right_sites.size())
	return;

    for (size_t i = 0; i < left_sites.size(); ++i)
      {
	Junction j(ref_id,
		   left_sites[i].first - 1,
		   right_sites[i].first + 2,
		   antisense);
	
	juncs.insert(j);
	
	if (juncs.size() > max_juncs)
	  juncs.erase(*(juncs.rbegin()));
      }
  }
};

struct ButterflyKey 
{
	uint32_t pos;
	uint32_t key;
	
	ButterflyKey(uint32_t p, uint32_t k) : pos(p), key(k) {}
	
	bool operator<(const ButterflyKey& rhs) const
	{
		if (key != rhs.key)
			return key < rhs.key;
		if (pos != rhs.pos)
			return pos < rhs.pos;
		return false;
	}
	
	bool operator==(const ButterflyKey& rhs) const
	{
		return pos == rhs.pos && key == rhs.key;
	}
};

uint32_t get_left_butterfly_key(uint64_t upstream_key, 
								const MerExtension& ext,
								size_t half_splice_mer_len)
{
	uint64_t key = ext.right_dna_str >> ((ext.right_ext_len - half_splice_mer_len) << 1);
	uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull << (half_splice_mer_len << 1));
	uint64_t top_half = upstream_key & mask;
	key |= (top_half << (half_splice_mer_len << 1));
	return (uint32_t)key;
}

uint32_t get_right_butterfly_key(uint64_t downstream_key, 
								const MerExtension& ext,
								size_t half_splice_mer_len)
{
	uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull << (half_splice_mer_len << 1));
	uint64_t key = (ext.left_dna_str & mask) << (half_splice_mer_len << 1);
	uint64_t bottom_half = (downstream_key >> (half_splice_mer_len << 1));
	key |= bottom_half;
	return (uint32_t)key;
}

struct RecordButterflyJuncs
{
  void record(uint32_t ref_id,
	      const vector<pair<size_t, DnaSpliceStrings> >& all_left_sites,
	      const vector<pair<size_t, DnaSpliceStrings> >& all_right_sites,
	      bool antisense,
	      PotentialJuncs& juncs,
	      int min_intron,
	      int max_intron,
	      size_t max_juncs,
	      size_t half_splice_mer_len)
  {
    size_t key_length = 2 * half_splice_mer_len;
    size_t extension_length = butterfly_overhang;
    uint64_t bottom_bit_mask = ~(0xFFFFFFFFFFFFFFFFull << (key_length<<1));
    uint64_t top_bit_mask =  ~(0xFFFFFFFFFFFFFFFFull >> (key_length<<1));
		
    if (all_left_sites.empty() || all_right_sites.empty())
      return;
    
    size_t last_site = max(all_left_sites.back().first, 
			   all_right_sites.back().first);
    
    size_t curr_left_site = 0;
    size_t curr_right_site = 0;
    for (size_t window_left_edge = 0; 
	 window_left_edge < last_site; 
	 window_left_edge += max_intron)
      {
	//fprintf(stderr, "\twindow %lu - %lu\n", window_left_edge, window_left_edge + 2 * max_intron); 
	vector<pair<size_t, DnaSpliceStrings> > left_sites;
	vector<pair<size_t, DnaSpliceStrings> > right_sites;
	
	while(curr_left_site < all_left_sites.size() &&
	      all_left_sites[curr_left_site].first < window_left_edge)
	  {
	    curr_left_site++;
	  }
	
	while(curr_right_site < all_right_sites.size() &&
	      all_right_sites[curr_right_site].first < window_left_edge)
	  {
	    curr_right_site++;
	  }
	
	for (size_t ls = curr_left_site; ls < all_left_sites.size(); ++ls)
	  {
	    if (all_left_sites[ls].first < window_left_edge + 2 * max_intron)
	      {
		left_sites.push_back(all_left_sites[ls]);
	      }
	  }
	
	for (size_t rs = curr_right_site; rs < all_right_sites.size(); ++rs)
	  {
	    if (all_right_sites[rs].first < window_left_edge + 2 * max_intron)
	      {
		right_sites.push_back(all_right_sites[rs]);
	      }
	  }
	
	vector<ButterflyKey> left_keys;
	for (size_t L = 0; L < left_sites.size(); ++L)
	  {
	    uint64_t fwd_upstream_dna_str = left_sites[L].second.fwd_string;
	    uint64_t fwd_upstream_key = fwd_upstream_dna_str & bottom_bit_mask;
	    
	    assert (fwd_upstream_key < extensions.size());
	    
	    vector<MerExtension>& fwd_exts = extensions[fwd_upstream_key];
	    for (size_t i = 0; i < fwd_exts.size(); ++i)
	      {
		const MerExtension& ext = fwd_exts[i];
		if (ext.right_ext_len < extension_length)
		  continue;
		
		/*
		  < f_u_key ><ext.right>
		  NNNNNNNNNN  GT
		*/
		
		// take the top bits of the right extension
		uint64_t key = ext.right_dna_str >> ((ext.right_ext_len - extension_length) << 1);
		
		// and the bottom bits of the site key
		uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull << (extension_length << 1));
		uint64_t top_half = fwd_upstream_key & mask;
		
		// and cat them together
		key |= (top_half << (extension_length << 1));
		left_keys.push_back(ButterflyKey((uint32_t)left_sites[L].first, key));
	      }
	    
	    uint64_t rev_upstream_dna_str = left_sites[L].second.rev_string;
	    uint64_t rev_upstream_key = (rev_upstream_dna_str & top_bit_mask) >> (64 - (key_length<<1));
	    
	    assert (rev_upstream_key < extensions.size());
	    
	    vector<MerExtension>& rev_exts = extensions[rev_upstream_key];
	    for (size_t i = 0; i < rev_exts.size(); ++i)
	      {
		const MerExtension& ext = rev_exts[i];
		if (ext.left_ext_len < extension_length)
		  continue;
		
		/*
		  < r_u_key ><ext.left>
		  NNNNNNNNNN  GT
		*/
				
		// reverse complement the left extension, and we will need 
		// what were the bottom bits.  these become the top bits in the 
		// rc.
		uint64_t ext_str = color ? rc_color_str(ext.left_dna_str) : rc_dna_str(ext.left_dna_str);
		ext_str >>= 64 - (ext.left_ext_len << 1);
		
		// now take the top bits of the rc, make them the bottom of 
		// the key
		uint64_t key = ext_str >> ((ext.left_ext_len - extension_length) << 1);
		
		// now add in the seed key bottom bits, making them the top of 
		// the key
		uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull << (extension_length << 1));
		uint64_t top_half = fwd_upstream_key & mask;
		key |= (top_half << (extension_length << 1));
		left_keys.push_back(ButterflyKey((uint32_t)left_sites[L].first, key));
	      }
	  }
	sort (left_keys.begin(), left_keys.end());
	vector<ButterflyKey>::iterator new_end = unique(left_keys.begin(), left_keys.end());
	left_keys.erase(new_end, left_keys.end());
	
	vector<ButterflyKey> right_keys;
	for (size_t R = 0; R < right_sites.size(); ++R)
	  {
	    uint64_t fwd_downstream_dna_str = right_sites[R].second.fwd_string;
	    uint64_t fwd_downstream_key = (fwd_downstream_dna_str & top_bit_mask) >> (64 - (key_length<<1));
	    
	    assert (fwd_downstream_key < extensions.size());
	    
	    vector<uint64_t> fwd_downstream_keys;
	    if (color)
	      {
		for(size_t color_value = 0; color_value < 4; ++color_value)
		  {
		    uint64_t tmp_key = (fwd_downstream_key << 2) >> 2 | (color_value << ((key_length - 1) << 1));
		    fwd_downstream_keys.push_back(tmp_key);
		  }
	      }
	    else
	      {
		fwd_downstream_keys.push_back(fwd_downstream_key);
	      }
	    
	    for(size_t key = 0; key < fwd_downstream_keys.size(); ++key)
	      {
		uint64_t tmp_fwd_downstream_key = fwd_downstream_keys[key];
		vector<MerExtension>& fwd_exts = extensions[tmp_fwd_downstream_key];
		for (size_t i = 0; i < fwd_exts.size(); ++i)
		  {
		    const MerExtension& ext = fwd_exts[i];
		    if (ext.left_ext_len < extension_length)
		      continue;
		    
		    /*
		      <ext.left>< f_d_key >
		      AG NNNNNNNNNN
		    */
		    
		    // take the bottom bits of the left extension, making them the
		    // top of the key.
		    uint64_t mask = ~(0xFFFFFFFFFFFFFFFFull << (extension_length << 1));
		    uint64_t key = (ext.left_dna_str & mask) << (extension_length << 1);
		    
		    // add in the top bits of the seed key, making them the bottom bits
		    // of the key.
		    uint64_t bottom_half = (tmp_fwd_downstream_key >> ((key_length - extension_length) << 1));
		    key |= bottom_half;
		    right_keys.push_back(ButterflyKey((uint32_t)right_sites[R].first, key));
		  }
	      }

	    uint64_t rev_downstream_dna_str = right_sites[R].second.rev_string;
	    uint64_t rev_downstream_key = rev_downstream_dna_str & bottom_bit_mask;
	    
	    assert (rev_downstream_key < extensions.size());
	    
	    vector<uint64_t> rev_downstream_keys;
	    if (color)
	      {
		for(size_t color_value = 0; color_value < 4; ++color_value)
		  {
		    uint64_t tmp_key = (rev_downstream_key >> 2) << 2 | color_value;
		    rev_downstream_keys.push_back(tmp_key);
		  }
	      }
	    else
	      {
		rev_downstream_keys.push_back(rev_downstream_key);
	      }
				
	    for(size_t key = 0; key < rev_downstream_keys.size(); ++key)
	      {
		uint64_t tmp_rev_downstream_key = rev_downstream_keys[key];
		uint64_t tmp_fwd_downstream_key = fwd_downstream_key;
		if (color)
		  {
		    tmp_fwd_downstream_key = rc_color_str(tmp_rev_downstream_key) >> (64 - (key_length << 1));
		  }
		
		vector<MerExtension>& rev_exts = extensions[tmp_rev_downstream_key];
		for (size_t i = 0; i < rev_exts.size(); ++i)
		  {
		    const MerExtension& ext = rev_exts[i];
		    if (ext.right_ext_len < extension_length)
		      continue;
		    
		    /*
		      <ext.right>< r_d_key >
		      AG NNNNNNNNNN
		    */
		    
		    // reverse complement the right_extension.  we want the 
		    // top bits of the extension, but these become the bottom bits
		    // under the rc.
		    uint64_t ext_str = color ? rc_color_str(ext.right_dna_str) : rc_dna_str(ext.right_dna_str);
		    ext_str >>= 64 - (ext.right_ext_len << 1);
		    
		    // take the bottom bits of the rc and make it the top of the key
		    uint64_t key = ext_str << (extension_length << 1);
		    
		    // take the top bits of the seed key and make them the bottom
		    // of the key.
		    uint64_t bottom_half = (tmp_fwd_downstream_key >> ((key_length - extension_length) << 1));
		    key |= bottom_half;
		    
		    right_keys.push_back(ButterflyKey((uint32_t)right_sites[R].first, key));
		  }
	      }
	  }
	
	sort (right_keys.begin(), right_keys.end());
	new_end = unique(right_keys.begin(), right_keys.end());
	right_keys.erase(new_end, right_keys.end());
			
	size_t lk = 0;
	size_t rk = 0;
	
	while (lk < left_keys.size() && rk < right_keys.size())
	  {
	    while (lk < left_keys.size() &&
		   left_keys[lk].key < right_keys[rk].key) { ++lk; }
	    
	    if (lk == left_keys.size())
	      break;
	    
	    while (rk < right_keys.size() &&
		   right_keys[rk].key < left_keys[lk].key) { ++rk; }
				
	    if (rk == right_keys.size())
	      break;
	    
	    if (lk < left_keys.size() && rk < right_keys.size() &&
		right_keys[rk].key == left_keys[lk].key)
	      {
					
		size_t k = right_keys[rk].key;
		size_t lk_end = lk;
		size_t rk_end = rk;
		while (rk_end < right_keys.size() && right_keys[rk_end].key == k) {++rk_end;}
		while (lk_end < left_keys.size() && left_keys[lk_end].key == k) {++lk_end;}
		
		size_t tmp_lk = lk;
					
		while (tmp_lk < lk_end)
		  {
		    size_t tmp_rk = rk;
		    while (tmp_rk < rk_end)
		      {
			int donor = (int)left_keys[tmp_lk].pos - 1;
			int acceptor = (int)right_keys[tmp_rk].pos + 2;
							
			if (acceptor - donor > min_intron && acceptor - donor < max_intron)
			  {
			    Junction j(ref_id,
				       donor,
				       acceptor,
				       antisense,
				       acceptor - donor); // just prefer shorter introns
			    juncs.insert(j);
			    if (juncs.size() > max_juncs)
			      {
				juncs.erase(*(juncs.rbegin()));
			      }
			  }
			++tmp_rk;
		      }
		    
		    ++tmp_lk;
		  }
		
		lk = lk_end;
		rk = rk_end;
	      }
	  }
      }
  }
};

template <class JunctionRecorder>
void juncs_from_ref_segs(RefSequenceTable& rt,
                         vector<RefSeg>& expected_don_acc_windows,
                         PotentialJuncs& juncs,
                         const DnaString& donor_dinuc,
                         const DnaString& acceptor_dinuc,
                         int max_intron,
                         int min_intron,
                         size_t max_juncs,
                         bool talkative,
                         size_t half_splice_mer_len)
{	
	
    typedef map<uint32_t, IntronMotifs> MotifMap;
    
    MotifMap ims;
	
    seqan::DnaStringReverseComplement rev_donor_dinuc(donor_dinuc);
    seqan::DnaStringReverseComplement rev_acceptor_dinuc(acceptor_dinuc);
    
    if (talkative)
        fprintf(stderr, "Collecting potential splice sites in islands\n");

    // 
    bool all_both = true;
    for (size_t r = 0; r < expected_don_acc_windows.size(); ++r)
    {
        const RefSeg& seg = expected_don_acc_windows[r];

	if (seg.points_where != POINT_DIR_BOTH)
	  all_both = false;
        
        RefSequenceTable::Sequence* ref_str = rt.get_seq(seg.ref_id);
        
        if (!ref_str)
            continue;

        bool skip_fwd = false;
        bool skip_rev = false;
        
        if (library_type == FR_FIRSTSTRAND)
        {
            if (seg.read == READ_LEFT)
            {
                if (seg.antisense) skip_rev = true;
                else skip_fwd = true;
            }
            else if(seg.read == READ_RIGHT)
            {
                if (seg.antisense) skip_fwd = true;
                else skip_rev = true;
            }
        }
        
        if (library_type == FR_SECONDSTRAND)
        {
            if (seg.read == READ_LEFT)
            {
                if (seg.antisense) skip_fwd = true;
                else skip_rev = true;
            }
            else if(seg.read == READ_RIGHT)
            {
                if (seg.antisense) skip_rev = true;
                else skip_fwd = true;
            }
        }
        
        pair<map<uint32_t, IntronMotifs>::iterator, bool> ret = 
        ims.insert(make_pair(seg.ref_id, IntronMotifs(seg.ref_id)));
		
        IntronMotifs& motifs = ret.first->second;

	int left_color_offset = 0, right_color_offset = 0;
	if (color)
	  {
	    if (seg.antisense)
	      right_color_offset = 1;
	    else
	      left_color_offset = -1;
	  }
	
        if (seg.left + left_color_offset < 0 || seg.right + right_color_offset >= (int)length(*ref_str) - 1)
            continue;
        
        DnaString org_seg_str = seqan::infix(*ref_str, seg.left + left_color_offset, seg.right + right_color_offset);
	String<char> seg_str;
	assign(seg_str, org_seg_str);

#ifdef B_DEBUG2
   cout << "coord: " << seg.left << " " << seg.right << endl;
	    //<< "seg_str: " << seg_str << endl;
#endif
	if (color)
	  {
	    bool remove_primer = true;
	    seg_str = convert_bp_to_color(org_seg_str, remove_primer);
	  }

	size_t to = 0;
	size_t seg_len = length(seg_str);
	size_t read_len = seg.support_read.size();
	if (read_len <= 0)
	  to = seg_len - 2;
	else
	  to = read_len - 2;

	const size_t max_segment_len = 128;
	uint8_t left_mismatches[max_segment_len] = {0,};
	uint8_t right_mismatches[max_segment_len] = {0,};

	if (max_segment_len < read_len)
	  {
	    fprintf(stderr, "Error: read len(%d) is greater than %d\n", (int)read_len, (int)max_segment_len);
	    exit(-1);
	  }
	
	if (read_len == seg_len || seg.points_where == POINT_DIR_BOTH)
	  {
	    if (seg.points_where == POINT_DIR_RIGHT || seg.points_where == POINT_DIR_BOTH)
	      {
		size_t num_mismatches = 0;
		for (size_t i = 0; i < read_len - 1; ++i)
		  {
		    if (seg_str[i] != seg.support_read[i])
		      ++num_mismatches;

		    left_mismatches[i] = num_mismatches;
		    if (num_mismatches > 2)
		      {
			to = i;
			break;
		      }
		  }
	      }

	    if (seg.points_where == POINT_DIR_LEFT || seg.points_where == POINT_DIR_BOTH)
	      {
		size_t num_mismatches = 0;
		for (int i = read_len - 1; i >= 0; --i)
		  {
		    if (seg_str[i + (seg_len - read_len)] != seg.support_read[i])
		      ++num_mismatches;

		    right_mismatches[i] = num_mismatches;
		    if (num_mismatches > 2)
			break;
		  }
	      }

	    // daehwan
#ifdef B_DEBUG2
		cout << "antisense: " << (seg.antisense ? "-" : "+") << endl
		     << seqan::infix(seg_str, 0, segment_length) << " "
		     << seqan::infix(seg_str, length(seg_str) - segment_length, length(seg_str)) << endl
		     << seg.support_read << endl
		     << 0 << " - " << to << endl;

		for (unsigned int i = 0; i < read_len; ++i)
		  cout << (int)left_mismatches[i];
		cout << "\t";

		for (unsigned int i = 0; i < read_len; ++i)
		  cout << (int)right_mismatches[i];
		cout << endl;		
#endif
	  }

	if (seg.points_where == POINT_DIR_BOTH)
	  {
	    for (size_t i = 0; i <= to; ++i)
	      {
                // Look at a slice of the reference without creating a copy.
                DnaString curr = seqan::infix(org_seg_str, i - left_color_offset, i + 2 - left_color_offset);

		if ((!skip_fwd && curr == donor_dinuc) || (!skip_rev && curr == rev_acceptor_dinuc))
		  {
		    DnaString partner;
		    if (curr == donor_dinuc)
		      partner = acceptor_dinuc;
		    else
		      partner = rev_donor_dinuc;

		    uint8_t left_mismatch = 0;
		    if (i > 0)
		      left_mismatch = left_mismatches[i-1];
		    
		    // daehwan
#ifdef B_DEBUG2
			cout << "i: " << i << endl
			     << "mismatches: " << (int)left_mismatch
			     << " - " << (int)right_mismatches[i] << endl;
#endif
		    if (left_mismatch + right_mismatches[i] <= 2)
		      {
			size_t pos = length(seg_str) - (read_len - i) - 2;
			if (partner == seqan::infix(org_seg_str, pos - left_color_offset, pos + 2 - left_color_offset))
			  {
			    if (curr == donor_dinuc)
			      {
				motifs.fwd_donors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
				motifs.fwd_acceptors.push_back(make_pair(seg.left + pos, DnaSpliceStrings(0,0)));
			      }
			    else
			      {
				motifs.rev_acceptors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
				motifs.rev_donors.push_back(make_pair(seg.left + pos, DnaSpliceStrings(0,0)));
			      }

			    // daehwan
        #ifdef B_DEBUG2
				cout << curr << ":" << partner << " added" << endl;
        #endif
			  }
		      }
		  }
	      }
	  }

	else if (seg.points_where == POINT_DIR_LEFT)
	  {
            // A ref segment that "points left" is one that was flanked 
            // on the right by a partial bowtie hit, indicating that we
            // should be looking for an intron to the left of the hit
	    // In this seg, that means either an "AG" or an "AC"
            for (size_t i = 0; i <= to; ++i)
	      {
                // Look at a slice of the reference without creating a copy.
                DnaString curr = seqan::infix(org_seg_str, i - left_color_offset, i + 2 - left_color_offset);
                if (curr == acceptor_dinuc && !skip_fwd)
		  motifs.fwd_acceptors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
                else if (curr == rev_donor_dinuc && !skip_rev)
		  motifs.rev_donors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
	      }
	  }
        else
	  {
            // A right pointing ref seg wants either a "GT" or a "CT"
            for (size_t i = 0; i <= to; ++i)
	      {
                // Look at a slice of the reference without creating a copy.
                DnaString curr = seqan::infix(org_seg_str, i - left_color_offset, i + 2 - left_color_offset);
                if (curr == donor_dinuc && !skip_fwd)
		  motifs.fwd_donors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
                else if (curr == rev_acceptor_dinuc && !skip_rev)
		  motifs.rev_acceptors.push_back(make_pair(seg.left + i, DnaSpliceStrings(0,0)));
	      }
	  }
    }
    
    if (talkative)
    {
        fprintf(stderr, "reporting synthetic splice junctions...\n");
    }
    for (MotifMap::iterator motif_itr = ims.begin(); motif_itr != ims.end(); ++motif_itr)
    {
        uint32_t ref_id = motif_itr->first;
        
        RefSequenceTable::Sequence* ref_str = rt.get_seq(ref_id);
        
        if (!ref_str)
            err_die("Error: couldn't get ref string for %u\n", ref_id);

        if (talkative)
            fprintf(stderr, "Examining donor-acceptor pairings in %s\n", rt.get_name(ref_id));
        IntronMotifs& motifs = motif_itr->second;

	if (!all_both)
	  motifs.unique();
	
        //motifs.attach_mer_counts(*ref_str);
        motifs.attach_mers(*ref_str);
        
        vector<pair<size_t, DnaSpliceStrings> >& fwd_donors = motifs.fwd_donors;
        vector<pair<size_t, DnaSpliceStrings> >& fwd_acceptors = motifs.fwd_acceptors;
        vector<pair<size_t, DnaSpliceStrings> >& rev_acceptors = motifs.rev_acceptors;
        vector<pair<size_t, DnaSpliceStrings> >& rev_donors = motifs.rev_donors;

        //const char* ref_name = rt.get_name(motif_itr->second.ref_id);
        
        JunctionRecorder recorder;
        recorder.record(ref_id,
                        fwd_donors, 
                        fwd_acceptors, 
                        false, 
                        juncs,
                        min_intron,
                        max_intron, 
                        max_juncs,
                        half_splice_mer_len);
        
        recorder.record(ref_id, 
                        rev_acceptors, 
                        rev_donors, 
                        true,
                        juncs,
                        min_intron,
                        max_intron, 
                        max_juncs,
                        half_splice_mer_len);
    }
    //fprintf(stderr, "Found %d total splices\n", num_juncs);
}


/**
 * Performs a simple global alignment.
 * This function will perform a restricted global alignment. The restriction is that only one insertion/deletion
 * is allowed in the final alignment.
 * @param shortSequence The short sequence to be aligned.
 * @param leftReference The left end of the reference to be aligned, must be exactly as long as the short sequence
 * @param leftReference The right end of the reference to be aligned, must be exactly as long as the short sequenc
 * @param insertPosition This will contain the 0-based index of the first position in the shorter sequence after the insertion/deletion. A value of -1 indicates that the alignment could not be performed.
 * @param mismatchCount This will contain the number of mismatches in the optimal restricted global alignment. The number and length of insertions/deletions is fixed. A value of -1 indicates that the alignment could not be performed.
 */
void simpleSplitAlignment(seqan::String<char>& shorterSequence,
			  seqan::String<char>& leftReference,
			  seqan::String<char>& rightReference,
			  int& insertPosition,
			  int& mismatchCount)
{
			/*
			 * In this restricted alignment, we already know the length and number (1) of insertions/deletions.
			 * We simply need to know where to put it. Do a linear scan through sequence counting the number of induced
			 * errors before and after putting the insertion at each sequence.
			 */

			/*
			 * Note that we could have a case, where both the alignment and the read have the unknonw
			 * nucleotide ('N') and we don't want to reward cases where these characters match
			 */
			vector<unsigned short> beforeErrors(seqan::length(shorterSequence));
			for(int idx = seqan::length(shorterSequence) - 1; idx >= 0; idx -= 1){
				unsigned short prevCount = 0;
				/*
				 * We guarentee idx >= 0, so cast to hide the compiler
				 * warning here
				 */
				if(((size_t)idx) < seqan::length(shorterSequence) - 1){
					prevCount = beforeErrors.at(idx + 1);
				}
				unsigned short currentMismatch = 0;
				if(rightReference[idx] == 'N' || shorterSequence[idx] == 'N' || rightReference[idx] != shorterSequence[idx]){
					currentMismatch = 1;
				}
				beforeErrors.at(idx) = prevCount + currentMismatch;
			}


			vector<unsigned short> afterErrors(seqan::length(shorterSequence));
			for(size_t idx = 0; idx < seqan::length(shorterSequence) ; idx += 1){
				unsigned short prevCount = 0;
				if(idx > 0){
					prevCount = afterErrors.at(idx - 1);
				}
				unsigned short currentMismatch = 0;
				if(leftReference[idx] == 'N' || shorterSequence[idx] == 'N' || leftReference[idx] != shorterSequence[idx]){
					currentMismatch = 1;
				}
				afterErrors.at(idx) = prevCount + currentMismatch;
			}


			mismatchCount = seqan::length(shorterSequence) + 1;
			insertPosition = -1;

			/*
			 * Technically, we could allow the insert position to be at the end or beginning of the sequence,
			 * but we are disallowing it here
			 */
			for(size_t currentInsertPosition = 1; currentInsertPosition < seqan::length(shorterSequence); currentInsertPosition += 1){
				size_t errorCount = beforeErrors.at(currentInsertPosition) + afterErrors.at(currentInsertPosition - 1);
				if(((int)errorCount) < mismatchCount){
					mismatchCount = (int)errorCount;
					insertPosition = currentInsertPosition;
				}
			}
			return;
}


/**
 * Try to detect a small insertion.
 * This code will try to identify a small insertion based on the ungapped alignment of two neighboring
 * segments. The general idea is to try to realign the local region, and see if we can reduce the
 * number of errors. Note that the function makes use of the global parameter "max_insertion_length" to limit the maximum
 * size of a detected insertion.
 * @param rt Sequence table used to lookup sequence information
 * @param leftHit The alignment of the left segment. Note that the leftHit must have a left position less than that of the right hit.
 * @param rightHit The alignment of the right segment. Note that the rightHit must have a left position greater than that of the left hit.
 * @param insertions If an insertion is sucessfully detected, it will be added to this set
 */
void detect_small_insertion(RefSequenceTable& rt,
		seqan::String<char>& read_sequence,
		BowtieHit& leftHit,
		BowtieHit& rightHit,
		std::set<Insertion>& insertions)
{

	RefSequenceTable::Sequence* ref_str = rt.get_seq(leftHit.ref_id());
	if(!ref_str){
		fprintf(stderr, "Error accessing sequence record\n");
	}else{
		size_t read_length = seqan::length(read_sequence);
		int begin_offset = 0;
		int end_offset = 0;

		if(color){
		  if(leftHit.antisense_align())
		    end_offset = 1;
		  else
		    begin_offset = -1;
		}
		
		if(leftHit.left() + begin_offset < 0)
		  return;

		/*
		 * If there is in fact a deletion, we are expecting the genomic sequence to be shorter than
		 * the actual read sequence
		 */
		int discrepancy = read_length - (rightHit.right() - leftHit.left());
		DnaString genomic_sequence_temp = seqan::infix(*ref_str, leftHit.left() + begin_offset, rightHit.right() + end_offset);
		String<char> genomic_sequence;
		assign(genomic_sequence, genomic_sequence_temp);

		if(color)
		  genomic_sequence = convert_bp_to_color(genomic_sequence, true);

		String<char> left_read_sequence = seqan::infix(read_sequence, 0, 0 + seqan::length(genomic_sequence));
		String<char> right_read_sequence = seqan::infix(read_sequence, read_length - seqan::length(genomic_sequence), read_length);

		int bestInsertPosition = -1;
		int minErrors = -1;
		simpleSplitAlignment(genomic_sequence, left_read_sequence, right_read_sequence, bestInsertPosition, minErrors);

		/*
		 * Need to decide if the insertion is suitably improves the alignment
		 */
		/*
		 * If these two segments anchors constitue the entire read, then we require
		 * that this alignment actually improve the number of errors observed in the alignment
		 * Otherwise, it is OK as long as the number of errors doesn't increase.
		 */
		int adjustment = 0;
		if(leftHit.read_len() + rightHit.read_len() >= (int)read_length){
			adjustment = -1;
		}
		if(minErrors <= (leftHit.edit_dist()+rightHit.edit_dist()+adjustment)){
			String<char> insertedSequence = seqan::infix(left_read_sequence, bestInsertPosition, bestInsertPosition + discrepancy);
			if(color)
			  insertedSequence = convert_color_to_bp(genomic_sequence_temp[bestInsertPosition - begin_offset + end_offset - 1], insertedSequence);
			
			insertions.insert(Insertion(leftHit.ref_id(),
					leftHit.left() + bestInsertPosition - 1 + end_offset,
					seqan::toCString(insertedSequence)));
		}
	}
	return;
}

/**
 * Try to detect a small deletion.
 * This code will try to identify a small deletion based on the ungapped alignment of two neighboring
 * segments. The general idea is to try to realign the local region, and see if we can reduce the
 * number of errors. Note that the function makes use of the global parameter "max_deletion_length" to limit the maximum
 * size of a detected deletion.
 * @param rt Sequence table used to lookup sequence information
 * @param leftHit The alignment of the left segment. Note that the leftHit must have a left position less than that of the right hit.
 * @param rightHit The alignment of the right segment. Note that the rightHit must have a left position greater than that of the left hit.
 * @param deletion_juncs If a deletion is sucessfully detected, it will be added to this set
 */
void detect_small_deletion(RefSequenceTable& rt,
		seqan::String<char>& read_sequence,
		BowtieHit& leftHit,
		BowtieHit& rightHit,
		std::set<Deletion>& deletions)
{

	RefSequenceTable::Sequence* ref_str = rt.get_seq(leftHit.ref_id());
	if(!ref_str){
		fprintf(stderr, "Error accessing sequence record\n");
	}else{
		int begin_offset = 0;
		int end_offset = 0;

		if(color){
		  if(leftHit.antisense_align())
		    end_offset = 1;
		  else
		    begin_offset = -1;
		}

		if(leftHit.left() + begin_offset < 0)
		  return;

		size_t read_length = seqan::length(read_sequence);
		if(rightHit.right() + read_length + begin_offset < 0 )
		  return;

		int discrepancy = (rightHit.right() - leftHit.left()) - read_length;
		Dna5String leftGenomicSequence_temp = seqan::infix(*ref_str, leftHit.left() + begin_offset, leftHit.left() + read_length + end_offset);
		Dna5String rightGenomicSequence_temp = seqan::infix(*ref_str, rightHit.right() - read_length + begin_offset, rightHit.right() + end_offset);

		if (length(leftGenomicSequence_temp) < read_length || length(rightGenomicSequence_temp) < read_length)
		  return;

		String<char> leftGenomicSequence;
		assign(leftGenomicSequence, leftGenomicSequence_temp);

		String<char> rightGenomicSequence;
		assign(rightGenomicSequence, rightGenomicSequence_temp);

		if(color){
		  leftGenomicSequence = convert_bp_to_color(leftGenomicSequence, true);
		  rightGenomicSequence = convert_bp_to_color(rightGenomicSequence, true);
		}

		int bestInsertPosition = -1;
		int minErrors = -1;

		simpleSplitAlignment(read_sequence, leftGenomicSequence, rightGenomicSequence, bestInsertPosition, minErrors);

		/*
		 * Need to decide if the deletion is suitably improves the alignment
		 */
		int adjustment = 0;

		/*
		 * If these two segments anchors constitue the entire read, then we require
		 * that this alignment actually improve the number of errors observed in the alignment
		 * Otherwise, it is OK as long as the number of errors doesn't increase.
		 */
		if(leftHit.read_len() + rightHit.read_len() >= (int)read_length){
			adjustment = -1;
		}
		if(minErrors <= (leftHit.edit_dist()+rightHit.edit_dist()+adjustment)){
			deletions.insert(Deletion(leftHit.ref_id(),
					leftHit.left() + bestInsertPosition - 1 + end_offset,
					leftHit.left() + bestInsertPosition + discrepancy + end_offset,
					false));
		}
	}
	return;
}


void find_insertions_and_deletions(RefSequenceTable& rt,
		ReadStream& reads_file,
		vector<HitsForRead>& hits_for_read,
		std::set<Deletion>& deletions,
		std::set<Insertion>& insertions){

	if(hits_for_read.empty()){
			return;
	}

	size_t last_segment = hits_for_read.size()-1;
	size_t first_segment = 0;
	if(last_segment == first_segment){
		return;
	}

	/*
	 * We can check up front whether the first or last element is empty
	 * and avoid doing any more work. Note that the following code requires
	 * that there be at least one elment in each
	 */
	if(hits_for_read[first_segment].hits.empty() || hits_for_read[last_segment].hits.empty()){
		return;
	}

	/*
	 * Need to identify the appropriate insert id for this group of reads
	 */
  Read read;
	/*bool got_read = get_read_from_stream(hits_for_read.back().insert_id,
			reads_file,
			FASTQ,
			false,
      read);*/
  bool got_read  = reads_file.getRead(hits_for_read.back().insert_id, read);
	if(!got_read){
	  err_die("Error: could not get read# %d from stream!",
	      (int)hits_for_read.back().insert_id);
		//return;
	  }

	/*
	 * Work through all combinations of mappings for the first and last segment to see if any are indicative
	 * of a small insertions or deletion
	 */
	HitsForRead& left_segment_hits = hits_for_read[first_segment];
	HitsForRead& right_segment_hits = hits_for_read[last_segment];

	/*
	 * If either of the segment match lists is empty, we could try
	 * to be smarter and work our way in until we find good a segment
	 * match; however, won't do that for noe.
	 */
	if(left_segment_hits.hits.empty() || right_segment_hits.hits.empty()){
		return;
	}

	seqan::String<char> fullRead, rcRead;
	if(color){
	  fullRead = read.seq.c_str() + 1;
	  rcRead = fullRead;
	  seqan::reverseInPlace(rcRead);
	}else{
	  fullRead = read.seq;
	  rcRead = read.seq;
	  seqan::convertInPlace(rcRead, seqan::FunctorComplement<Dna>());
	  seqan::reverseInPlace(rcRead);
	}

	size_t read_length = seqan::length(fullRead);
	for(size_t left_segment_index = 0; left_segment_index < left_segment_hits.hits.size(); left_segment_index++){
		for(size_t right_segment_index = 0; right_segment_index < right_segment_hits.hits.size(); right_segment_index++){
			BowtieHit* leftHit = &left_segment_hits.hits[left_segment_index];
			BowtieHit* rightHit = &right_segment_hits.hits[right_segment_index];
			/*
			 * Now we have found a pair of segment hits to investigate. Need to ensure
			 * that
			 * 1. the alignment orientation is consistent
			 * 2. the distance separation is in the appropriate range
			 * 3. Both hits are aligned to the same contig
			 */
			if(leftHit->ref_id() != rightHit->ref_id()){
				continue;
			}

			if(leftHit->antisense_align() != rightHit->antisense_align()){
				continue;
			}

			seqan::String<char>* modifiedRead = &fullRead;
			/*
			 * If we are dealing with an antisense alignment, then the left
			 * read will actually be on the right, fix this now, to simplify
			 * the rest of the logic, in addition, we will need to use the reverse
			 * complement of the read sequence
			 */
			if(leftHit->antisense_align()){
				BowtieHit * tmp = leftHit;
				leftHit = rightHit;
				rightHit = tmp;
				modifiedRead = &rcRead;
			}

			size_t apparent_length = rightHit->right() - leftHit->left();
			int length_discrepancy = apparent_length - read_length;
			if(length_discrepancy > 0 && length_discrepancy <= (int)max_deletion_length){
				/*
				 * Search for a deletion
				 */
				detect_small_deletion(rt, *modifiedRead, *leftHit, *rightHit, deletions);
			}
			if(length_discrepancy < 0 && length_discrepancy >= -(int)max_insertion_length){
				/*
				 * Search for an insertion
				 */
				detect_small_insertion(rt, *modifiedRead, *leftHit, *rightHit, insertions);
			}
		}
	}
}

/*
 */
int map_read_to_contig(const String<char>& contig, const String<char>& read)
{
  int contig_len = length(contig);
  int read_len = length(read);

  int pos = -1;
  int mismatch = 3;

  for (int i = 0; i < contig_len - read_len; ++i)
    {
      int temp_mismatch = 0;
      for (int j = 0; j < read_len; ++j)
	{
	  if (contig[i+j] != read[j])
	    ++temp_mismatch;

	  if (temp_mismatch >= mismatch)
	    break;
	}

      if (temp_mismatch < mismatch)
	{
	  pos = i;
	  mismatch = temp_mismatch;
	}
    }

  return pos;
}

void find_gaps(RefSequenceTable& rt,
	       ReadStream& reads_file,
	       vector<HitsForRead>& hits_for_read,
	       HitStream& partner_hit_stream,
	       HitStream& seg_partner_hit_stream,
	       std::set<Junction, skip_count_lt>& seg_juncs,
	       eREAD read_side)
{
  if (hits_for_read.empty())
    return;

  size_t last_segment = hits_for_read.size() - 1;
  while (last_segment > 0)
    {
      if (!hits_for_read[last_segment].hits.empty())
	break;
      
      --last_segment;
    }

  hits_for_read.resize(last_segment + 1);

  size_t first_segment = 0;
  if (last_segment == first_segment &&
      (hits_for_read[first_segment].hits.empty() || hits_for_read[first_segment].hits[0].end()))
    return;

  uint32_t insert_id = hits_for_read[last_segment].insert_id;
  
  HitsForRead partner_hit_group;
  uint32_t next_order = partner_hit_stream.next_group_id();

  bool has_partner = false;
  while (insert_id >= next_order && next_order != 0)
    {
      partner_hit_stream.next_read_hits(partner_hit_group);
      next_order = partner_hit_stream.next_group_id();
    }
  
  has_partner = insert_id == partner_hit_group.insert_id;

  if (!has_partner)
    {
      next_order = seg_partner_hit_stream.next_group_id();
      while (insert_id >= next_order && next_order != 0)
	{
	  seg_partner_hit_stream.next_read_hits(partner_hit_group);
	  next_order = seg_partner_hit_stream.next_group_id();
	}

      has_partner = insert_id == partner_hit_group.insert_id;
    }
  
  Read read;
  /*
  bool got_read = get_read_from_stream(hits_for_read[last_segment].insert_id, 
				       reads_file,
				       FASTQ,
				       false,
				       read); */
  bool got_read = reads_file.getRead(hits_for_read[last_segment].insert_id, read);
  if (!got_read) {
    err_die("Error: could not get read# %d from stream!",
        (int)hits_for_read[last_segment].insert_id);
    return;
    }

  HitsForRead partner_hits_for_read;
  if (first_segment != last_segment)
    partner_hits_for_read = hits_for_read[last_segment];
  
  HitsForRead& left_segment_hits = hits_for_read[first_segment];
  HitsForRead& right_segment_hits = partner_hits_for_read;

  bool check_partner = true;
  if (first_segment != last_segment)
    {
      for (size_t i = 0; i < left_segment_hits.hits.size(); ++i)
      {
      BowtieHit& leftHit = left_segment_hits.hits[i];
      for (size_t j = 0; j < right_segment_hits.hits.size(); ++j)
        {
          BowtieHit& rightHit = right_segment_hits.hits[j];

          if (leftHit.ref_id() == rightHit.ref_id() && leftHit.antisense_align() == rightHit.antisense_align())
          {
            int dist = 0;
            if (leftHit.antisense_align())
              dist = leftHit.left() - rightHit.right();
            else
              dist = rightHit.left() - leftHit.right();

            if (dist >= min_segment_intron_length && dist < (int)max_segment_intron_length)
              {
                check_partner = false;
                break;
              }
          }
        }

      if (!check_partner)
        break;
      }
    }

  if (check_partner && has_partner)
    {
      // empty hits from 1 to num_segments - 1
      for (size_t i = first_segment + 1; i < hits_for_read.size(); ++i)
        {
          hits_for_read[i].hits.clear();
        }
      
      seqan::String<char> fullRead, rcRead;
      fullRead = read.seq;
      rcRead = read.seq;
      seqan::convertInPlace(rcRead, seqan::FunctorComplement<Dna>());
      seqan::reverseInPlace(rcRead);
      size_t read_length =  read.seq.length();

      for (size_t l = 0; l < left_segment_hits.hits.size(); ++l)
      {
        BowtieHit& leftHit = left_segment_hits.hits[l];
        for (size_t r = 0; r < partner_hit_group.hits.size(); ++r)
          {
            BowtieHit& rightHit = partner_hit_group.hits[r];
            if (leftHit.ref_id() != rightHit.ref_id() || leftHit.antisense_align() == rightHit.antisense_align())
        continue;

            int dist = 0;
            if (leftHit.antisense_align())
        dist = leftHit.left() - rightHit.right();
            else
        dist = rightHit.left() - leftHit.right();

            if (dist < min_segment_intron_length && dist >= (int)max_segment_intron_length)
        continue;

            RefSequenceTable::Sequence* ref_str = rt.get_seq(rightHit.ref_id());
            const size_t part_seq_len = inner_dist_std_dev > inner_dist_mean ? inner_dist_std_dev - inner_dist_mean : 0;
            const size_t flanking_seq_len = inner_dist_mean + inner_dist_std_dev;

            Dna5String right_flanking_seq;
            size_t left = 0;
            if (rightHit.antisense_align())
        {
          if (flanking_seq_len <= rightHit.left())
            {
              left = rightHit.left() - flanking_seq_len;
              right_flanking_seq = seqan::infix(*ref_str, left, left + flanking_seq_len + part_seq_len);
            }
          else
            break;
        }
            else
        {
          if (part_seq_len <= rightHit.right())
            {
              left = rightHit.right() - part_seq_len;
              right_flanking_seq = seqan::infix(*ref_str, left, left + flanking_seq_len + part_seq_len);
            }
          else
            break;
        }

            const size_t check_read_len = min(15, segment_length - segment_mismatches - 3);
            seqan::String<char> fwd_read = infix(fullRead, read_length - check_read_len, read_length);
            seqan::String<char> rev_read = infix(rcRead, 0, check_read_len);

            int fwd_pos = map_read_to_contig(right_flanking_seq, fwd_read);
            if (fwd_pos >= 0)
        {
          BowtieHit hit(rightHit.ref_id(), rightHit.insert_id(),
              left + fwd_pos, check_read_len, false, 0, true);

          hits_for_read[last_segment].hits.push_back(hit);
        }

            int rev_pos = map_read_to_contig(right_flanking_seq, rev_read);

            if (rev_pos >= 0)
        {
          BowtieHit hit(rightHit.ref_id(), rightHit.insert_id(),
              left + rev_pos, check_read_len, true, 0, true);

          hits_for_read[last_segment].hits.push_back(hit);
        }

            // daehwan - for debug purposes
    #if B_DEBUG
            cerr << "daehwan!!!" << endl
           << "insert id: " << rightHit.insert_id() << endl
           << "first segment: " << first_segment << ", last_segment: " << last_segment << endl
           << right_flanking_seq << " : " << seqan::length(right_flanking_seq) << endl
           << fwd_read << " : " << fwd_pos << endl
           << rev_read << " : " << rev_pos << endl
           << "left: " << leftHit.left() << "-" << leftHit.right() << (leftHit.antisense_align() ? " -" : " +") << endl
           << "right: " << rightHit.left() << "-" << rightHit.right() << (rightHit.antisense_align() ? " -" : " +") << endl;
            if (fwd_pos >= 0 || rev_pos >= 0)
        {
          const BowtieHit& hit = hits_for_read[last_segment].hits.back();
          cerr << "back: " << hit.left() << "-" << hit.right() << (hit.antisense_align() ? " -" : " +") << endl;
        }
    #endif
          }
      }
    }

  vector<RefSeg> expected_don_acc_windows;
  string seq(read.seq);

  for (size_t s = 0; s < hits_for_read.size(); ++s)
    {
      HitsForRead& curr = hits_for_read[s];
      for (size_t h = 0; h < curr.hits.size(); ++h)
	{
	  bool found_right_seg_partner = s == hits_for_read.size() - 1;
	  BowtieHit& bh = curr.hits[h];

	  // "drs" is distant seg right partner
	  // "rrs" is right of right seg partner
	  vector<BowtieHit*> drs_bhs;
	  vector<BowtieHit*> rrs_bhs;

	  if (s < hits_for_read.size() - 1)
	    {
	      // Look for a right partner for the current hit
	      HitsForRead& right = hits_for_read[s + 1];
	      
	      for (size_t r = 0; r < right.hits.size(); ++r)
		{
		  BowtieHit& rh = right.hits[r];
		  if (bh.antisense_align() != rh.antisense_align() || bh.ref_id() != rh.ref_id())
		    continue;
		  
		  if ((bh.antisense_align() && rh.right() == bh.left()) ||
		      (!bh.antisense_align() && bh.right() == rh.left() ))
		    {
		      found_right_seg_partner = true;
		      break;
		    }

		  int dist = 0;
		  if (bh.antisense_align())
		    dist = bh.left() - rh.right();
		  else
		    dist = rh.left() - bh.right();
		  
		  if (dist >= min_segment_intron_length && dist < (int)max_segment_intron_length)
		    drs_bhs.push_back(&rh);
		}
	    }

	  if (!found_right_seg_partner && s < hits_for_read.size() - 2)
	    {
	      // Look for a right of right partner for the current hit
	      HitsForRead& right_right = hits_for_read[s + 2];
	      
	      for (size_t r = 0; r < right_right.hits.size(); ++r)
		{
		  BowtieHit& rrh = right_right.hits[r];
		  if (bh.antisense_align() != rrh.antisense_align() || bh.ref_id() != rrh.ref_id())
		    continue;

		  int dist = 0;
		  if (bh.antisense_align())
		    dist = bh.left() - rrh.right();
		  else
		    dist = rrh.left() - bh.right();
		  
		  if (dist >= min_segment_intron_length + segment_length && dist < (int)max_segment_intron_length + segment_length)
		    rrs_bhs.push_back(&rrh);
		}
	    }
	  
	  if (!found_right_seg_partner && (drs_bhs.size() > 0 || rrs_bhs.size() > 0))
	    {
	      const int look_bp = 8;
	      const size_t color_offset = color ? 1 : 0;

	      vector<BowtieHit*> d_bhs = rrs_bhs.size() > 0 ? rrs_bhs : drs_bhs;
	      for (size_t r = 0; r < d_bhs.size(); ++r)
		{
		  string support_read;
		  if (rrs_bhs.size() <= 0)
		    support_read = seq.substr(color_offset + (s+1) * segment_length - look_bp, look_bp * 2);
		  else
		    support_read = seq.substr(color_offset + (s+1) * segment_length - look_bp, segment_length + look_bp * 2);
		  
		  BowtieHit& d_bh = *(d_bhs[r]);
		  if (!bh.antisense_align())
		    {
		      RefSeg right_seg(bh.ref_id(), POINT_DIR_BOTH, bh.antisense_align(), read_side, 0, 0, support_read);
		      right_seg.left = max(0, bh.right() - look_bp);
		      right_seg.right = d_bh.left() + look_bp;
		      expected_don_acc_windows.push_back(right_seg);
		    }
		  else
		    {
		      if (color)
			reverse(support_read.begin(), support_read.end());
		      else
			reverse_complement(support_read);
		      
		      RefSeg left_seg(bh.ref_id(), POINT_DIR_BOTH, bh.antisense_align(), read_side, 0, 0, support_read);
		      left_seg.left = d_bh.right() - look_bp;
		      left_seg.right = bh.left() + look_bp;
		      expected_don_acc_windows.push_back(left_seg);	
		    }

		  // daehwan
#ifdef B_DEBUG2
		  cout << "insert id: " << bh.insert_id() << endl
		       << (bh.antisense_align() ? "-" : "+") << endl
		       << seq << endl
		       << "(" << s << ") - " << support_read << endl;
#endif
		}
	    }
	}
    } //for each hits_for_read

  juncs_from_ref_segs<RecordSegmentJuncs>(rt, 
				      expected_don_acc_windows, 
				      seg_juncs, 
				      "GT", 
				      "AG",  
				      max_segment_intron_length, 
				      min_segment_intron_length, 
				      max_seg_juncs,
				      false,
				      0);

  juncs_from_ref_segs<RecordSegmentJuncs>(rt, 
				      expected_don_acc_windows, 
				      seg_juncs, 
				      "GC", 
				      "AG", 
				      max_segment_intron_length, 
				      min_segment_intron_length, 
				      max_seg_juncs,
				      false,
				      0);
  
  juncs_from_ref_segs<RecordSegmentJuncs>(rt, 
				      expected_don_acc_windows, 
				      seg_juncs, 
				      "AT", 
				      "AC",  
				      max_segment_intron_length, 
				      min_segment_intron_length, 
				      max_seg_juncs,
				      false,
				      0);
}


MerTable mer_table;

int seed_alignments = 0;
int microaligned_segs = 0;

map<RefSeg, vector<string>* > microexon_windows;

bool overlap_in_genome(int ll, int lr, int rl, int rr)
{
	if (ll >= rl && ll < rr)
		return true;
	if (lr > rl && lr < rr)
		return true;
	if (rl >= ll && rl < lr)
		return true;
	if (rr > ll && rr < lr)
		return true;
	return false;
}

void add_to_microexon_windows(uint32_t ref_id, 
			      int left_boundary, 
			      int right_boundary, 
			      const string& dna_str,
			      eREAD read)
{
  RefSeg left_dummy(ref_id, POINT_DIR_DONTCARE, false, read, left_boundary, right_boundary);
  RefSeg right_dummy(ref_id, POINT_DIR_DONTCARE, false, read, right_boundary, right_boundary + 1);
  
  map<RefSeg, vector<string>* >::iterator lb = microexon_windows.lower_bound(left_dummy);
  map<RefSeg, vector<string>* >::iterator ub = microexon_windows.lower_bound(right_dummy);
  vector<string>* new_vec = NULL;
  if (lb == microexon_windows.end())
    {
      microexon_windows.insert(make_pair(left_dummy, new vector<string>(1, dna_str)));
      return;
    }
  
  map<RefSeg, vector<string>* >::iterator first_to_be_erased = microexon_windows.end();
  map<RefSeg, vector<string>* >::iterator last_to_be_erased = ub;
  
  while (lb != ub)
    {
      // everyone in this range that overlaps with the new interval needs to 
      // be merged together.
      if (overlap_in_genome(lb->first.left, lb->first.right, left_boundary, right_boundary))
	{
	  if (!new_vec)
	    new_vec = new vector<string>();
	  if (first_to_be_erased == microexon_windows.end())
	    first_to_be_erased = lb;
	  left_dummy.left = min(lb->first.left, left_boundary);
	  left_dummy.right = max(lb->first.right, right_boundary);
	  
	  new_vec->insert(new_vec->end(), lb->second->begin(), lb->second->end()); 
	  delete lb->second;
	}
      else if (first_to_be_erased != microexon_windows.end())
	{
	  last_to_be_erased = lb;
	}
      
      ++lb;
    }
  
  if (first_to_be_erased != microexon_windows.end())
    {
      microexon_windows.erase(first_to_be_erased, last_to_be_erased);
    }
  
  if (!new_vec)
    {
      // never found an overlapping window, so just add this one and bail
      microexon_windows.insert(make_pair(left_dummy, new vector<string>(1, dna_str)));
      return;
    }
  else
    {
      new_vec->push_back(dna_str);
      microexon_windows.insert(make_pair(left_dummy, new_vec));
      return;
    }
  
}

void align_microexon_segs(RefSequenceTable& rt,
			  std::set<Junction, skip_count_lt>& juncs,
			  int max_juncs,
			  int half_splice_mer_len)
{
	int num_segments = 0;
	for (map<RefSeg, vector<string>* >::iterator itr = microexon_windows.begin(); 
		 itr != microexon_windows.end(); ++itr)
	{
		vector<string>& unaligned_segments = *itr->second;
		num_segments += unaligned_segments.size();
	}
	
	fprintf(stderr, "Aligning %d microexon segments in %lu windows\n",
	       num_segments, (long unsigned int)microexon_windows.size());
	
	extensions.clear();
	
	size_t splice_mer_len = 2 * half_splice_mer_len;
	size_t mer_table_size = 1 << ((splice_mer_len)<<1);
	
	extensions.resize(mer_table_size);

	int window_num = 0;
	for (map<RefSeg, vector<string>* >::iterator itr = microexon_windows.begin(); 
		 itr != microexon_windows.end(); ++itr)
	{
		window_num++;
		if ((window_num % 100) == 0)
			fprintf(stderr, "\twindow %d\n",window_num);
		
		
		stringstream ss(stringstream::in | stringstream::out);

		for (size_t j = 0; j < extensions.size(); ++j)
		{
			extensions[j].clear();
		}
		
		vector<string>& unaligned_segments = *itr->second;
		
		for (size_t j = 0; j < unaligned_segments.size(); ++j)
		{
			stringstream ss(stringstream::in | stringstream::out);
			string s;
			//cerr << w.unaligned_segments[j];
			ss << unaligned_segments[j];
			ss >> s;

			store_read_extensions(extensions,
					      half_splice_mer_len,
					      half_splice_mer_len,
					      s,
					      false);
		}
		
		vector<RefSeg> segs;
		segs.push_back(itr->first);
		RefSeg r = itr->first;
		r.points_where = POINT_DIR_LEFT;
		segs.push_back(r);
		
		juncs_from_ref_segs<RecordExtendableJuncs>(rt, 
							   segs, 
							   juncs, 
							   "GT", 
							   "AG", 
							   max_microexon_stretch, 
							   min_coverage_intron_length, 
							   max_juncs,
							   false,
							   half_splice_mer_len);
		num_segments += unaligned_segments.size();
		delete itr->second;
	}
	fprintf(stderr, "Checked %d segments against %lu windows for microexon junctions\n",
	               num_segments, (long unsigned int)microexon_windows.size());
	fprintf(stderr, "Found %ld potential microexon junctions\n", (long int)juncs.size());
}

/*
 * Easy guys ... this function puts its pants on just like the rest of you -- one leg
 * at a time. Except, once its pants are on, it makes gold records. Alright, here we
 * go.
 */
 
void look_for_hit_group(RefSequenceTable& rt,
    ReadStream& readstream,
    ReadStream& readstream_for_segment_search,
    ReadStream& readstream_for_indel_discovery,
			ReadTable& unmapped_reads,
			vector<HitStream>& seg_files,
			int curr_file,
			uint64_t insert_id,
			vector<HitsForRead>& hits_for_read, //<-- collecting segment hits for this read
			HitStream& partner_hit_stream,
			HitStream& seg_partner_hit_stream,
			std::set<Junction, skip_count_lt>& juncs,
			std::set<Deletion>& deletions,
			std::set<Insertion>& insertions,
			eREAD read_side)

{
  HitStream& hit_stream = seg_files[curr_file];
  HitsForRead hit_group;
  uint32_t order = unmapped_reads.observation_order(insert_id);
  int seq_key_len = min((int)min_anchor_len, 6);
  while(true) {
	uint64_t next_group_id = hit_stream.next_group_id();
	uint32_t next_order = unmapped_reads.observation_order(next_group_id);
	  // If we would have seen the hits by now, stop looking in this stream,
	  // but forward the search to the next (lower) segment if possible.
	if (order < next_order || next_group_id == 0) {
		 if (curr_file > 0) {
		   //look for next (lower) segment mappings for this read
			 look_for_hit_group(rt,
			  readstream,
			  readstream_for_segment_search,
			  readstream_for_indel_discovery,
			  unmapped_reads,
			  seg_files,
			  curr_file - 1,
			  insert_id,
			  hits_for_read,
					    partner_hit_stream,
					    seg_partner_hit_stream,
			  juncs,
			  deletions,
			  insertions,
			  read_side);
			  }
		 else
		   if (insert_id && !no_microexon_search) {
			  //microexon search
		     Read read;
			  // The hits are missing for the leftmost segment, which means
			  // we should try looking for junctions via seed and extend
			  // using it (helps find junctions to microexons).
			  /* bool got_read = get_read_from_stream(insert_id,
					readstream.file,
					FASTQ,
					false,
					read); */
        bool got_read = readstream.getRead(insert_id, read);
			  if (!got_read) {
				  //fprintf(stderr, "Warning: could not get read with insert_id=%d\n", (int)insert_id);
          //break; //exit loop
			    err_die("Error: could not get read with insert_id=%d from file %s\n",
			        (int)insert_id, readstream.filename());
				  }
			  string fwd_read = read.seq;
			  if (color) // remove the primer and the adjacent color
				 fwd_read.erase(0, 2);
			  // make sure there are hits for all the other segs, all the
			  // way to the root (right-most) one.
			  int empty_seg = 0;
			  for (size_t h = 1; h < hits_for_read.size(); ++h) {
					   if (hits_for_read[h].hits.empty())
						empty_seg = h;
						}
			  // if not, no microexon alignment for this segment
			  if (empty_seg != 0)
				  break;
			  fwd_read = fwd_read.substr(0, segment_length);
			  string rev_read = fwd_read;
			  //check the reverse
			  if (color) reverse(rev_read.begin(), rev_read.end());
					else reverse_complement(rev_read);
			  for (size_t h = 0; h < hits_for_read[empty_seg + 1].hits.size(); ++h) {
				  const BowtieHit& bh = hits_for_read[empty_seg + 1].hits[h];
				  RefSequenceTable::Sequence* ref_str = rt.get_seq(bh.ref_id());
				  if (ref_str == NULL)
					 continue;
				  int ref_len = length(*ref_str);
				  int left_boundary;
				  int right_boundary;
				  bool antisense = bh.antisense_align();
				  vector<BowtieHit> empty_seg_hits;
				  seed_alignments++;
				  if (antisense) {
					  left_boundary = max(0, bh.right() - (int)min_anchor_len);
					  right_boundary = min(ref_len - 2,
						 left_boundary +  max_microexon_stretch);
					  if (right_boundary - left_boundary < 2 * seq_key_len)
						  continue;
					  microaligned_segs++;
					  add_to_microexon_windows(bh.ref_id(), left_boundary, right_boundary, rev_read, read_side); 
					  }
				  else {
					  right_boundary = min(ref_len - 2,
						 bh.left() + (int)min_anchor_len);
					  left_boundary = max(0, right_boundary -  max_microexon_stretch);
					  if (right_boundary - left_boundary < 2 * seq_key_len)
						  continue;
					  microaligned_segs++;
					  add_to_microexon_windows(bh.ref_id(), left_boundary, right_boundary, fwd_read, read_side);  
					  }
				  } //for h
		   } // !no_microexon_search
		 break;
		 }
	else if (hit_stream.next_read_hits(hit_group)) {
		 // if we found hits for the target group in the left stream,
		 // add them to the accumulating vector and continue the search
		 if (hit_group.insert_id == insert_id) {
			hits_for_read[curr_file] = hit_group;
			if (curr_file > 0)
			  // we need to look left (recursively) for the group we 
			  // just read for this stream.
			  look_for_hit_group(rt,
					 readstream,
					 readstream_for_segment_search,
					 readstream_for_indel_discovery,
					 unmapped_reads,
					 seg_files,
					 curr_file - 1,
					 insert_id,
					 hits_for_read,
					     partner_hit_stream,
					     seg_partner_hit_stream,
					 juncs,
					 deletions,
					 insertions,
					 read_side);
			break;
			} //same insert_id (group)
		 else if (curr_file > 0) {
			// different group, we need to start a whole new 
			// search for it, with a whole new vector of hits.
			vector<HitsForRead> hits_for_new_read(seg_files.size());
			hits_for_new_read[curr_file] = hit_group;
			look_for_hit_group(rt,
			   readstream,
			   readstream_for_segment_search,
			   readstream_for_indel_discovery,
			   unmapped_reads,
			   seg_files,
			   curr_file - 1,
			   hit_group.insert_id,
			   hits_for_new_read,
					   partner_hit_stream,
					   seg_partner_hit_stream,
			   juncs,
			   deletions,
			   insertions,
			   read_side);
			find_insertions_and_deletions(rt,
						      readstream_for_indel_discovery,
						      hits_for_new_read,
						      deletions,
						      insertions);
			find_gaps(rt,
				  readstream_for_segment_search,
				  hits_for_new_read,
				  partner_hit_stream,
				  seg_partner_hit_stream,
				  juncs,
				  read_side);
			} //different group
		 }//got next group
	} //while loop
}



bool process_next_hit_group(RefSequenceTable& rt,
       ReadStream&  readstream,
       ReadStream&  readstream_for_segment_search,
       ReadStream&  readstream_for_indel_discovery,
			    ReadTable& unmapped_reads,
			    vector<HitStream>& seg_files, 
			    size_t last_file_idx,
			    HitStream& partner_hit_stream,
			    HitStream& seg_partner_hit_stream,
			    std::set<Junction, skip_count_lt>& juncs,
			    std::set<Deletion>& deletions,
			    std::set<Insertion>& insertions,
			    eREAD read)
{
  HitStream& last_segmap_hitstream = seg_files[last_file_idx];
  HitsForRead hit_group;
  bool result = last_segmap_hitstream.next_read_hits(hit_group);
  vector<HitsForRead> hits_for_read(seg_files.size());
  hits_for_read.back() = hit_group;

  look_for_hit_group(rt,
		     readstream,
		     readstream_for_segment_search,
		     readstream_for_indel_discovery,
		     unmapped_reads, 
		     seg_files, 
		     (int)last_file_idx - 1,
		     hit_group.insert_id,
		     hits_for_read,
		     partner_hit_stream,
		     seg_partner_hit_stream,
		     juncs,
		     deletions,
		     insertions,
		     read);
  
  if (result) {
    find_insertions_and_deletions(rt,
				  readstream_for_indel_discovery,
				  hits_for_read,
				  deletions,
				  insertions);
    
    find_gaps(rt,
	      readstream_for_segment_search,
	      hits_for_read,
	      partner_hit_stream,
	      seg_partner_hit_stream,
	      juncs,
	      read);
    }
  return result;
}

static const int UNCOVERED = 0;
static const int LOOK_LEFT = 1;
static const int LOOK_RIGHT = 2;

uint8_t get_cov(const vector<uint8_t>& cov, uint32_t c)
{
	uint32_t b = c >> 2;
	uint32_t r = c & 0x3;
	uint8_t  s = (r << 1);
	uint8_t  v = cov[b];
	v &= (0x3 << s);
	v >>= s;
	return v;
}

void build_coverage_map(ReadTable& it, RefSequenceTable& rt,
    vector<FZPipe*>& seg_files, map<uint32_t, vector<bool> >& coverage_map) {
 if (!coverage_map.empty()) return;
 BowtieHitFactory hit_factory(it,rt);

 for (size_t f = 0; f < seg_files.size(); ++f)
   {
     //fprintf(stderr, "Adding hits from segment file %d to coverage map\n", (int)f);
     seg_files[f]->rewind();
     FILE* fp = seg_files[f]->file;
     HitStream hs(fp, &hit_factory, false, false, false);
     HitsForRead hit_group;
     while (hs.next_read_hits(hit_group))
     {
       for (size_t h = 0; h < hit_group.hits.size(); ++h)
         {
           BowtieHit& bh = hit_group.hits[h];

           pair<map<uint32_t, vector<bool> >::iterator, bool> ret =
               coverage_map.insert(make_pair(bh.ref_id(), vector<bool>()));
           vector<bool>& ref_cov = ret.first->second;

           size_t right_extent = bh.right();

           if (right_extent >= ref_cov.size())
           {
             ref_cov.resize(right_extent + 1, 0);
           }

           for (uint32_t c = (uint32_t)bh.left(); c < (uint32_t)bh.right(); ++c)
           {
             ref_cov[c] = true;
             //if (ref_cov[c]<255) ref_cov[c]++;
           }
         }
     } //while next_read_hits
   }
}

void pair_covered_sites(ReadTable& it,
			RefSequenceTable& rt,
			vector<FZPipe*>& seg_files,
			std::set<Junction, skip_count_lt>& cov_juncs,
			map<uint32_t, vector<bool> >& coverage_map,
			size_t half_splice_mer_len)
{
  vector<RefSeg> expected_look_left_windows;
  vector<RefSeg> expected_look_right_windows;
  build_coverage_map(it,rt, seg_files, coverage_map);
  
  static const int extend = 45;
  int num_islands = 0;
  
  vector<RefSeg> expected_don_acc_windows;
  
  fprintf(stderr, "Recording coverage islands\n");
  size_t cov_bases = 0;
  for (map<uint32_t, vector<bool> >::iterator itr = coverage_map.begin();
       itr != coverage_map.end();
       ++itr)
    {
      vector<bool>& cov = itr->second;
      
      size_t island_left_edge = 0;
      for (size_t c = 1; c < cov.size(); ++c)
	{
	  if (cov[c])
	    {
	      cov_bases++;
	      if (!cov[c - 1])
		{
		  num_islands += 1;
		  int edge = (int)c - extend;
		  edge = max(edge, 0);
		  island_left_edge = edge;
		}
	    }
	  else
	    {
	      if (cov[c - 1])
		{
		  expected_don_acc_windows.push_back(RefSeg(itr->first,
							    POINT_DIR_LEFT,
							    false, /* not important */
							    READ_DONTCARE,
							    island_left_edge,
							    c + extend));
		  expected_don_acc_windows.push_back(RefSeg(itr->first,
							    POINT_DIR_RIGHT,
							    false, /* not important */
							    READ_DONTCARE,
							    island_left_edge,
							    c + extend));	
		}
	    }
	}
    }
  fprintf(stderr, "Found %d islands covering %ld bases\n", num_islands, (long int)cov_bases);
  
  juncs_from_ref_segs<RecordButterflyJuncs>(rt, 
					    expected_don_acc_windows, 
					    cov_juncs, 
					    "GT", 
					    "AG", 
					    max_coverage_intron_length, 
					    min_coverage_intron_length, 
					    max_cov_juncs,
					    true,
					    half_splice_mer_len);
  fprintf(stderr, "Found %ld potential intra-island junctions\n", (long int)cov_juncs.size());
}

// daehwan
struct ReadInfo
{
  ReadID read_id;
  uint32_t left;
  uint32_t right;

  bool operator<(const ReadInfo& rhs) const 
  {
    if (left != rhs.left)
      return left < rhs.left;
    if (right != rhs.right)
      return right < rhs.right;
    
    return false;
  }
};

void capture_island_ends(ReadTable& it,
			 RefSequenceTable& rt,
			 //vector<FILE*>& seg_files,
			 vector<FZPipe*>& seg_files,
			 std::set<Junction, skip_count_lt>& cov_juncs,
			  map<uint32_t, vector<bool> >& coverage_map,
			 size_t half_splice_mer_len)
{
  //static int island_repeat_tolerance = 10;
  vector<RefSeg> expected_look_left_windows;
  vector<RefSeg> expected_look_right_windows;

  // daehwan
//#define DEBUG_CHECK_EXONS 1
//#define DEBUG_RANGE_ONLY 1

#ifndef DEBUG_CHECK_EXONS
  build_coverage_map(it, rt, seg_files, coverage_map);
#else
  //build coverage map here, so we can debug it
  #ifdef DEBUG_RANGE_ONLY
   static const uint32_t chr14_id = RefSequenceTable::hash_string("chr14");
  #endif
  vector<ReadInfo> hits;
  BowtieHitFactory hit_factory(it,rt);
  
  for (size_t f = 0; f < seg_files.size(); ++f)
	{
	  fprintf(stderr, "Adding hits from segment file %d to coverage map\n", (int)f);
	  seg_files[f]->rewind();
	  FILE* fp = seg_files[f]->file;
	  //rewind(fp);
	  HitStream hs(fp, &hit_factory, false, false, false);
	  
	  HitsForRead hit_group;
	  while (hs.next_read_hits(hit_group))
	  {
	  for (size_t h = 0; h < hit_group.hits.size(); ++h)
		{
		  BowtieHit& bh = hit_group.hits[h];
		  // daehwan
		  //if (check_exons) <-- DEBUG_CHECK_EXONS
		  #ifdef DEBUG_RANGE_ONLY
		     if (bh.ref_id() != chr14_id)
		       continue;

		     // if (bh.left() < 66567028 && bh.right() > 66604392)
		     if (bh.left() < 66400000 || bh.right() > 66700000)
		       continue;

		     ReadInfo read_info;
		     read_info.read_id = bh.insert_id();
		     read_info.left = bh.left();
		     read_info.right = bh.right();
		     hits.push_back(read_info);
		   #endif

		  pair<map<uint32_t, vector<bool> >::iterator, bool> ret = 
		     coverage_map.insert(make_pair(bh.ref_id(), vector<bool>()));
		  vector<bool>& ref_cov = ret.first->second;
		  
		  size_t right_extent = bh.right();
		  
		  if (right_extent >= ref_cov.size())
		   {
		   ref_cov.resize(right_extent + 1, 0);
		   }
		  for (uint32_t c = (uint32_t)bh.left(); c < (uint32_t)bh.right(); ++c)
		   {
		   ref_cov[c] = true;
		   }
		}
	  }
	}

  sort(hits.begin(), hits.end());
#endif 
//  static const int min_cov_length = segment_length + 2;
  long covered_bases = 0;
  int long_enough_bases = 0;
  int left_looking = 0;
  int right_looking = 0;
  static const int extend = 45;
  static const int repeat_tol = 5;
  
  int num_islands = 0;
  
  for (map<uint32_t, vector<bool> >::iterator itr = coverage_map.begin();
       itr != coverage_map.end();
       ++itr)
    {
      #ifdef B_DEBUG
      fprintf (stderr, "Finding pairings in ref seq %s\n", rt.get_name(itr->first));
      #endif
      vector<bool>& cov = itr->second;
      vector<bool> long_enough(cov.size());
      
      size_t last_uncovered = 0;
      
      static const uint8_t min_cov = 1;
      
      for (size_t c = 1; c < cov.size(); ++c)
      {
      uint8_t c_cov = cov[c]; //get_cov(cov, c);
      if (c_cov < min_cov || c == (cov.size()) - 1)
        {
        int putative_exon_length = (int)c - (int)last_uncovered;
        uint32_t last_pos_cov = cov[c - 1]; //get_cov(cov,c - 1);
        if (last_pos_cov >= min_cov && putative_exon_length >= min_cov_length)
          {
          #ifdef B_DEBUG
          fprintf(stderr, "cov. island: %d-%d\n", (int)(last_uncovered + 1), (int)c);
          fprintf(stderr, "\t(putative exon length = %d, min_cov_length=%d)\n",putative_exon_length, min_cov_length);
          #endif
          covered_bases += (c + 1 - last_uncovered);
          for (int l = (int)c; l > (int)last_uncovered; --l)
            {
              long_enough[l] = true;
            }
          }
        last_uncovered = c;
        }
      }
      vector<bool>& ref_cov = long_enough;
      vector<uint8_t> cov_state(ref_cov.size(), UNCOVERED);

      // daehwan - print islands (exons)
      //if (check_exons)
#ifdef DEBUG_CHECK_EXONS
	//{
	  uint32_t left = 0, right = 0;
	  for (size_t c = 1; c < ref_cov.size(); ++c)
		{
		if (ref_cov[c] && !ref_cov[c-1])
		  {
		   left = c;
		   cout << "Exon: " << left << "-";
		  }
		  else if (!ref_cov[c] && ref_cov[c-1])
		  {
		    right = c - 1;
		    cout << right << endl;
		    for (size_t k = 0; k < hits.size(); ++k)
		    {
		      const ReadInfo& hit = hits[k];
		      if (hit.right < left)
		         continue;
		
		      if (hit.left > right)
		         break;
		      
		      cout << "\t" << hit.read_id << " " << hit.left << "-" << hit.right << endl;
		    }
		  }
		}
	//}
#endif
	for (size_t c = 1; c < ref_cov.size(); ++c)
	{
	  if (ref_cov[c])
		{
		  long_enough_bases++;
		  if (!ref_cov[c - 1])
		  {
		  num_islands += 1;
		  
		  for (int r = (int)c - extend; 
		       r >= 0 && r < (int)c + repeat_tol && r < (int)cov_state.size(); 
		       ++r)
		    {
		      cov_state[r] |= LOOK_LEFT;
		    }
		  }
	    }
	  else
	    {
	      if (ref_cov[c - 1])
		{
		  for (int l = (int)c - repeat_tol; 
		       l >= 0 && l < (int)c + extend && l < (int)cov_state.size(); 
		       ++l)
		    {
		      cov_state[l] |= LOOK_RIGHT;
		    }	
		}
	    }
	}
      
      RefSeg* curr_look_left = NULL;
      RefSeg* curr_look_right = NULL;
      
	for (size_t c = 1; c < cov_state.size(); ++c)
	{
	  if (cov_state[c] & LOOK_LEFT)
	    {
	      left_looking++;
	      if (!(cov_state[c-1] & LOOK_LEFT))
		{
		  expected_look_left_windows.push_back(RefSeg(itr->first,
							      POINT_DIR_LEFT,
							      false, /* not important */
							      READ_DONTCARE,
							      c,
							      c + 1));
		  curr_look_left = &(expected_look_left_windows.back());
		}
	      else if (curr_look_left)
		{
		  curr_look_left->right++;
		}
	    }
	  else
	    {
	      if ((cov_state[c-1] & LOOK_LEFT))
		{
		  curr_look_left = NULL;
		}
	    }
	  
	  if (cov_state[c] & LOOK_RIGHT)
	    {
	      right_looking++;
	      if (!(cov_state[c-1] & LOOK_RIGHT))
		  {
		    expected_look_right_windows.push_back(RefSeg(itr->first,
							         POINT_DIR_RIGHT,
							         false, /* not important */
							         READ_DONTCARE,
							         c,
							         c + 1));

		    curr_look_right = &(expected_look_right_windows.back());
		  }
	      else if (curr_look_right)
		   {
		     curr_look_right->right++;
		   }
	    }
	  else
	    {
	      if ((cov_state[c-1] & LOOK_RIGHT))
		   {
		     curr_look_right = NULL;
		   }
	    }
	 }
	}
  
  fprintf(stderr, " Map covers %ld bases\n", covered_bases);
  fprintf(stderr, " Map covers %d bases in sufficiently long segments\n", long_enough_bases);
  fprintf(stderr, " Map contains %d good islands\n", num_islands + 1);
  fprintf(stderr, " %d are left looking bases\n", left_looking);
  fprintf(stderr, " %d are right looking bases\n", right_looking);

  vector<RefSeg> expected_don_acc_windows;
  expected_don_acc_windows.insert(expected_don_acc_windows.end(),
				  expected_look_right_windows.begin(),
				  expected_look_right_windows.end());
  
  expected_don_acc_windows.insert(expected_don_acc_windows.end(),
				  expected_look_left_windows.begin(),
				  expected_look_left_windows.end());
  
  if (!butterfly_search) coverage_map.clear(); //free some memory
  juncs_from_ref_segs<RecordExtendableJuncs>(rt, 
					     expected_don_acc_windows, 
					     cov_juncs, 
					     "GT", 
					     "AG",  
					     max_coverage_intron_length, 
					     min_coverage_intron_length, 
					     max_cov_juncs,
					     true,
					     half_splice_mer_len);
  //fprintf(stderr, "Found %ld potential island-end pairing junctions\n", (long int)cov_juncs.size());
}


void process_segment_hits(RefSequenceTable& rt,
              ReadStream& readstream,
              ReadStream& readstream_for_segment_search,
              ReadStream& readstream_for_indel_discovery,
			  vector<FZPipe>& seg_files,
			  FZPipe& partner_reads_map_file,
			  FZPipe& seg_partner_reads_map_file,
			  BowtieHitFactory& hit_factory,
			  ReadTable& it,
			  std::set<Junction, skip_count_lt>& juncs,
			  std::set<Deletion>& deletions,
			  std::set<Insertion>& insertions,
  			  eREAD read = READ_DONTCARE)
{
  HitStream partner_hit_stream(partner_reads_map_file, &hit_factory, false, false, false);
  HitStream seg_partner_hit_stream(seg_partner_reads_map_file, &hit_factory, false, false, false);
  
  vector<HitStream> hit_streams;
  for (size_t i = 0; i < seg_files.size(); ++i)
    {
      HitStream hs(seg_files[i].file, &hit_factory, false, false, false);
      
      // if the next group id in this stream is zero, it's an empty stream,
      // and we can simply skip it.
      //if (hs.next_group_id() != 0)
      hit_streams.push_back(hs);
    }
  
  int num_group = 0;
  while (process_next_hit_group(rt,
				readstream,
				readstream_for_segment_search,
				readstream_for_indel_discovery,
				it, 
				hit_streams, 
				hit_streams.size() - 1,
				partner_hit_stream,
				seg_partner_hit_stream,
				juncs, deletions, insertions, read) == true)
    {
      num_group++;
      if (num_group % 500000 == 0)
	     fprintf(stderr, "\tProcessed %d root segment groups\n", num_group);
    }
  fprintf(stderr, "Microaligned %d segments\n", microaligned_segs); 
}

void print_juncs(RefSequenceTable& rt, std::set<Junction, skip_count_lt>& juncs, const char* str)
{
  fprintf (stderr, "-- %s --\n", str);
  for(std::set<Junction>::iterator itr = juncs.begin();
      itr != juncs.end();
      ++itr)
    {
      const char* ref_name = rt.get_name(itr->refid);
      
      fprintf(stderr, "%s\t%d\t%d\t%c\n",
	      ref_name,
	      itr->left,
	      itr->right,
	      itr->antisense ? '-' : '+');
    }
  fprintf (stderr, "-- done --\n");
}

void driver(istream& ref_stream,
	    FILE* juncs_out,
	    FILE* insertions_out,
	    FILE* deletions_out,
	    ReadStream& left_readstream,
	    ReadStream& left_readstream_for_segment_search,
	    ReadStream& left_readstream_for_indel_discovery,
	    FZPipe& left_reads_map_file,
	    FZPipe& left_seg_file_for_segment_search,
	    vector<FZPipe>& left_seg_files,
	    ReadStream& right_readstream,
	    ReadStream& right_readstream_for_segment_search,
	    ReadStream& right_readstream_for_indel_discovery,
	    FZPipe& right_reads_map_file,
	    FZPipe& right_seg_file_for_segment_search,
	    vector<FZPipe>& right_seg_files
	    )
{	
  if (left_seg_files.size() == 0)
    {
      fprintf(stderr, "No hits to process, exiting\n");
      exit(0);
    }
  
  std::set<Junction, skip_count_lt> seg_juncs;
  std::set<Junction, skip_count_lt> cov_juncs;
  std::set<Junction, skip_count_lt> butterfly_juncs;
  
  std::set<Junction> juncs;
  
  std::set<Deletion> deletions;
  std::set<Insertion> insertions;
  
  RefSequenceTable rt(true, true);
  
  fprintf (stderr, "Loading reference sequences...\n");
  get_seqs(ref_stream, rt, true, false);
	
  ReadTable it;
  BowtieHitFactory hit_factory(it,rt);
  fprintf(stderr, ">> Performing segment-search:\n");
  
  if (left_seg_files.size() > 1)
    {
      fprintf( stderr, "Loading left segment hits...\n");
      process_segment_hits(rt,
			   left_readstream,
			   left_readstream_for_segment_search,
			   left_readstream_for_indel_discovery,
			   left_seg_files,
			   right_reads_map_file,
			   right_seg_file_for_segment_search,
			   hit_factory,
			   it,
			   seg_juncs,
			   deletions,
			   insertions,
			   READ_LEFT);
      fprintf( stderr, "done.\n");
    }
  
  if (right_seg_files.size() > 1)
    {
      fprintf( stderr, "Loading right segment hits...\n");
      process_segment_hits(rt,
			   right_readstream,
			   right_readstream_for_segment_search,
			   right_readstream_for_indel_discovery,
			   right_seg_files,
			   left_reads_map_file,
			   left_seg_file_for_segment_search,
			   hit_factory,
			   it,
			   seg_juncs,
			   deletions,
			   insertions,
			   READ_RIGHT);
      fprintf( stderr, "done.\n");
    }

  fprintf(stderr, "\tfound %ld potential split-segment junctions\n", (long int)seg_juncs.size());
  fprintf(stderr, "\tfound %ld potential small deletions\n", (long int)deletions.size());
  fprintf(stderr, "\tfound %ld potential small insertions\n", (long int)insertions.size());
  
  //vector<FILE*> all_seg_files;
	vector<FZPipe*> all_seg_files;
  for (vector<FZPipe>::size_type i = 0; i != left_seg_files.size();i++)
    {
       all_seg_files.push_back( &(left_seg_files[i]) );
	   }
  
  for (vector<FZPipe>::size_type i = 0; i != right_seg_files.size();i++)
    {
       all_seg_files.push_back( &(right_seg_files[i]) );
       }
	//copy(left_seg_files.begin(), left_seg_files.end(), back_inserter(all_seg_files));
	//copy(right_seg_files.begin(), right_seg_files.end(), back_inserter(all_seg_files));

#if 0
  // daehwan - check this out as Cole insists on using segments gives better results.
  vector<FZPipe*> all_map_files;
  if (left_seg_files.size() > 1)
    {
      all_map_files.push_back(&left_reads_map_file);
    }

  if (right_seg_files.size() > 1)
    {
      all_map_files.push_back(&right_reads_map_file);
    }
  
  copy(all_seg_files.begin(), all_seg_files.end(), back_inserter(all_map_files));
#endif
  map<uint32_t, vector<bool> > coverage_map;
  if (!no_coverage_search || butterfly_search)
    {
      if (ium_reads != "")
	  {
	    vector<string> ium_read_files;
	    tokenize(ium_reads,",", ium_read_files);
	    vector<FZPipe> iums;
          string unzcmd=getUnpackCmd(ium_read_files[0],false);
	    for (size_t ium = 0; ium < ium_read_files.size(); ++ium)
	      {
              //fprintf (stderr, "Indexing extensions in %s\n", ium_read_files[ium].c_str());
				  FZPipe ium_file(ium_read_files[ium],unzcmd);
				  if (ium_file.file==NULL)
		  {
		    fprintf (stderr, "Can't open file %s for reading, skipping...\n",ium_read_files[ium].c_str());
		    continue;
		  }
	        iums.push_back(ium_file);
	      }

	    index_read_mers(iums, 5);
	  }
      else
	  { //no unmapped reads
	    no_coverage_search = true;
	    butterfly_search = false;
	  }

      if (!no_coverage_search)
	    { //coverage search
	    // looking for junctions by island end pairings
	    fprintf(stderr, ">> Performing coverage-search:\n");
	    capture_island_ends(it,
			      rt,
			      all_seg_files,
			      cov_juncs, 
			      coverage_map,
			      5);
		fprintf(stderr, "\tfound %d potential junctions\n",(int)cov_juncs.size());
	    }
    } //coverage search or butterfly search
  
  if (butterfly_search)
    {
      //looking for junctions between and within islands
      fprintf(stderr, ">> Performing butterfly-search: \n");
      prune_extension_table(butterfly_overhang);
      compact_extension_table();
      pair_covered_sites(it,
			 rt,
			 all_seg_files, 
			 butterfly_juncs,
			 coverage_map,
			 5);
      
      fprintf(stderr, "\tfound %d potential junctions\n",(int)butterfly_juncs.size());
    }
  
  coverage_map.clear();
  
  std::set<Junction, skip_count_lt> microexon_juncs;
  if (!no_microexon_search)
    {
      fprintf(stderr, ">> Performing microexon-search: \n");
      std::set<Junction, skip_count_lt> microexon_juncs;
      align_microexon_segs(rt,
			   microexon_juncs,
			   max_cov_juncs,
			   5);
      fprintf(stderr, "\tfound %d potential junctions\n",(int)microexon_juncs.size());
      juncs.insert(microexon_juncs.begin(), microexon_juncs.end());
    }

  juncs.insert(cov_juncs.begin(), cov_juncs.end());
  juncs.insert(seg_juncs.begin(), seg_juncs.end());
  juncs.insert(butterfly_juncs.begin(), butterfly_juncs.end());
  
  //fprintf(stderr, "Reporting potential splice junctions...");
  for(std::set<Junction>::iterator itr = juncs.begin();
      itr != juncs.end();
      ++itr)
    {
      const char* ref_name = rt.get_name(itr->refid);
      
      fprintf(juncs_out,
	      "%s\t%d\t%d\t%c\n",
	      ref_name,
	      itr->left,
	      itr->right,
	      itr->antisense ? '-' : '+');
    }
	//close all reading pipes, just to exit cleanly

    for (vector<FZPipe*>::size_type i = 0; i != all_seg_files.size();i++) {
       all_seg_files[i]->close();
       }
    left_readstream.close();
    left_readstream_for_indel_discovery.close();
    right_readstream.close();
    right_readstream_for_indel_discovery.close();
  //fprintf (stderr, "done\n");
  fprintf (stderr, "Reported %d total possible splices\n", (int)juncs.size());
  
  fprintf (stderr, "Reporting potential deletions...\n");
	if(deletions_out){
      for(std::set<Deletion>::iterator itr = deletions.begin(); itr != deletions.end(); ++itr){
	const char* ref_name = rt.get_name(itr->refid);
	/*
	 * We fix up the left co-ordinate to reference the first deleted base
	 */
	fprintf(deletions_out,
		"%s\t%d\t%d\n",
		ref_name,
		itr->left + 1,
		itr->right);
      }
      fclose(deletions_out);
	}else{
      fprintf(stderr, "Failed to open deletions file for writing, no deletions reported\n");
    }
  
  fprintf (stderr, "Reporting potential insertions...\n");
	if(insertions_out){
      for(std::set<Insertion>::iterator itr = insertions.begin(); itr != insertions.end(); ++itr){
	const char* ref_name = rt.get_name(itr->refid);
	fprintf(insertions_out,
		"%s\t%d\t%d\t%s\n",
		ref_name,
		itr->left,
		itr->left,
		itr->sequence.c_str());
      }
      fclose(insertions_out);
	}else{
      fprintf(stderr, "Failed to open insertions file for writing, no insertions reported\n");
    }
}

int main(int argc, char** argv)
{
  fprintf(stderr, "segment_juncs v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION);
  fprintf(stderr, "---------------------------\n");
  
  int parse_ret = parse_options(argc, argv, print_usage);
  if (parse_ret)
    return parse_ret;
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string ref_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string juncs_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string insertions_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string deletions_file_name = argv[optind++];
  if(optind >= argc)
    {
      print_usage();
      return 1;
    } 
 
  string left_reads_file_name = argv[optind++];

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_reads_map_file_name = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_segment_file_list = argv[optind++];
    
  string right_segment_file_list; 
  string right_reads_file_name;
  string right_reads_map_file_name;
  if (optind < argc)
    {
      right_reads_file_name = argv[optind++];

      if(optind >= argc)
	{
	  print_usage();
	  return 1;
	}
      
      right_reads_map_file_name = argv[optind++];
      
      if(optind >= argc)
	{
	  print_usage();
	  return 1;
	}
      right_segment_file_list = argv[optind++];
    }
  
  // Open the approppriate files
  
  ifstream ref_stream(ref_file_name.c_str(), ifstream::in);
  if (!ref_stream.good())
	  {
	    fprintf(stderr, "Error: cannot open %s for reading\n",
		    ref_file_name.c_str());
	    exit(1);
	  }
  
  
  FILE* juncs_file = fopen(juncs_file_name.c_str(), "w");
  if (!juncs_file)
    {
      fprintf(stderr, "Error: cannot open %s for writing\n",
	      juncs_file_name.c_str());
      exit(1);
    }

  FILE* insertions_file = fopen(insertions_file_name.c_str(), "w");
  if (!insertions_file)
    {
      fprintf(stderr, "Error: cannot open %s for writing\n",
              insertions_file_name.c_str());
      exit(1);
    }


  FILE* deletions_file = fopen(deletions_file_name.c_str(), "w");
  if (!deletions_file)
    {
      fprintf(stderr, "Error: cannot open %s for writing\n",
              deletions_file_name.c_str());
      exit(1);
    }
	
  //string unzcmd=getUnpackCmd(left_reads_file_name, false);
  ReadStream left_reads_file(left_reads_file_name);
  ReadStream left_reads_file_for_segment_search(left_reads_file_name);
  ReadStream left_reads_file_for_indel_discovery(left_reads_file_name);
  if (left_reads_file.file()==NULL || left_reads_file_for_segment_search.file()==NULL ||
        left_reads_file_for_indel_discovery.file()==NULL)
    {
      fprintf(stderr, "Error: cannot open %s for reading\n",
	      left_reads_file_name.c_str());
      exit(1);
    }

  //FILE* left_reads_map_file = fopen(left_reads_map_file_name.c_str(), "r");
  FZPipe left_reads_map_file(left_reads_map_file_name, true);
  if (left_reads_map_file.file == NULL)
    {
      fprintf(stderr, "Error: cannot open %s for reading\n",
	      left_reads_map_file_name.c_str());
      exit(1);
    }

  vector<string> left_segment_file_names;
  vector<FZPipe> left_segment_files;
  tokenize(left_segment_file_list, ",",left_segment_file_names);
  FZPipe left_segment_file_for_segment_search;
  string unzcmd;
  if (left_segment_file_names.size()>0)
     unzcmd=getUnpackCmd(left_segment_file_names[0], false);
  for (size_t i = 0; i < left_segment_file_names.size(); ++i)
    {
      FZPipe seg_file(left_segment_file_names[i], unzcmd);
      if (seg_file.file == NULL)
        {
        fprintf(stderr, "Error: cannot open segment map file %s for reading\n",
          left_segment_file_names[i].c_str());
        exit(1);
        }
      left_segment_files.push_back(seg_file);

      if (i == left_segment_file_names.size() - 1)
      {
        left_segment_file_for_segment_search = FZPipe(left_segment_file_names[i], unzcmd);
        if (left_segment_file_for_segment_search.file == NULL)
          {
            fprintf(stderr, "Error: cannot open %s for reading\n",
              left_segment_file_names[i].c_str());
            exit(1);
          }
      }
    }
  FZPipe right_reads_map_file;
  vector<FZPipe> right_segment_files;
  ReadStream right_reads_file;
  ReadStream right_reads_file_for_segment_search;
  ReadStream right_reads_file_for_indel_discovery;
  FZPipe right_segment_file_for_segment_search;

  if (right_segment_file_list != "")
    {
      //unzcmd=getUnpackCmd(right_reads_file_name, false);
      right_reads_file.init(right_reads_file_name);
      right_reads_file_for_segment_search.init(right_reads_file_name);
      right_reads_file_for_indel_discovery.init(right_reads_file_name);
      if (right_reads_file.file()==NULL || right_reads_file_for_indel_discovery.file()==NULL)
      {
        fprintf(stderr, "Error: cannot open %s for reading\n",
          right_reads_file_name.c_str());
        exit(1);
      }
      
      right_reads_map_file = FZPipe(right_reads_map_file_name, true);
      if (right_reads_map_file.file == NULL)
      {
        fprintf(stderr, "Error: cannot open %s for reading\n",
          right_reads_map_file_name.c_str());
        exit(1);
      }
            
      vector<string> right_segment_file_names;
      tokenize(right_segment_file_list, ",",right_segment_file_names);
      if (right_segment_file_names.size()>0)
         unzcmd=getUnpackCmd(right_segment_file_names[0], false);

      for (size_t i = 0; i < right_segment_file_names.size(); ++i)
		  {
      FZPipe seg_file(right_segment_file_names[i],unzcmd);
      if (seg_file.file == NULL)
        {
          fprintf(stderr, "Error: cannot open segment map %s for reading\n",
            right_segment_file_names[i].c_str());
          exit(1);
        }
      right_segment_files.push_back(seg_file);

	  if (i == right_segment_file_names.size() - 1)
	    {
	      right_segment_file_for_segment_search = FZPipe(right_segment_file_names[i], unzcmd);
	      if (right_segment_file_for_segment_search.file == NULL)
		{
		  fprintf(stderr, "Error: cannot open %s for reading\n",
			  right_segment_file_names[i].c_str());
		  exit(1);
	}
		  }
	}
    }
  // min_cov_length=20;
  if (min_cov_length>segment_length-2) min_cov_length=segment_length-2;
  
  driver(ref_stream, 
	 juncs_file,
	 insertions_file,
	 deletions_file,
	 left_reads_file,
	 left_reads_file_for_segment_search,
	 left_reads_file_for_indel_discovery, 
	 left_reads_map_file,
	 left_segment_file_for_segment_search,
	 left_segment_files, 
	 right_reads_file,
	 right_reads_file_for_segment_search,
	 right_reads_file_for_indel_discovery,
	 right_reads_map_file,
	 right_segment_file_for_segment_search,
	 right_segment_files);
  
  return 0;
}
