/*
 *  bwt_map.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "reads.h"
using namespace std;

void HitTable::add_hit(const BowtieHit& bh, bool check_uniqueness)
{
	uint32_t reference_id = bh.ref_id();
	
	pair<RefHits::iterator, bool> ret = 
	_hits_for_ref.insert(make_pair(reference_id, HitList()));
	HitList& hl = ret.first->second;

	if (check_uniqueness)
	{
		// Check uniqueness, in case we are adding spliced hits from 
		// several spliced alignment sources (e.g. de novo hashing + Bowtie 
		// against a user-supplied index).  We don't want to count the same
		// alignment twice if it happened to be found by more than one method
		HitList::const_iterator lb = lower_bound(hl.begin(), hl.end(), bh, hit_insert_id_lt);   
		HitList::const_iterator ub = upper_bound(hl.begin(), hl.end(), bh, hit_insert_id_lt); 
		for (; lb != ub && lb != hl.end(); ++lb)
		{
			if (*lb == bh)
			{
				//fprintf(stderr, "Chucking duplicate read %d by identity\n", bh.insert_id());
				return;
			}
			
			if (lb->insert_id() == bh.insert_id() &&
				lb->ref_id() == bh.ref_id() &&
				lb->antisense_align() == bh.antisense_align())
			{
				// If we get here, we may be looking at the same alignment
				// However, spanning_reads may report a shorter, trimmed alignment
				// so not all fields will be equal.  If they just disagree on the 
				// ends, and don't indicate a different junction coord, the 
				// alignments are the same.

				if ((lb->left() <= bh.left() && lb->right() >= bh.right()) ||
					(bh.left() <= lb->left() && bh.right() >= lb->right()))
				{
					vector<pair<int,int> > lb_gaps, bh_gaps;
					lb->gaps(lb_gaps);
					bh.gaps(bh_gaps);
					if (lb_gaps == bh_gaps)
					{
						// One alignment is contained in the other, they agree on 
						// where the gaps, if any, are, and they share an id
						// => this is a redundant aligment, so toss it
						//fprintf(stderr, "Chucking duplicate read %d by gap agreement\n", bh.insert_id());
						return;
					}
				}
			}
		}
	}
    _total_hits++;
	hl.push_back(bh);
}

bool hit_insert_id_lt(const BowtieHit& h1, const BowtieHit& h2)
{
	return h1.insert_id() < h2.insert_id();
}

BowtieHit HitFactory::create_hit(const string& insert_name, 
								 const string& ref_name,
								 int left,
								 const vector<CigarOp>& cigar,
								 bool antisense_aln,
								 bool antisense_splice,
								 unsigned char edit_dist,
								 unsigned char splice_mms)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL);
	
	return BowtieHit(reference_id,
					 insert_id, 
					 left, 
					 cigar, 
					 antisense_aln,
					 antisense_splice,
					 edit_dist,
					 splice_mms);	
}

BowtieHit HitFactory::create_hit(const string& insert_name, 
								 const string& ref_name,
								 uint32_t left,
								 uint32_t read_len,
								 bool antisense_aln,
								 unsigned char edit_dist)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL);
	
	return BowtieHit(reference_id,
					 insert_id, 
					 left,
					 read_len,
					 antisense_aln,
					 edit_dist);	
}

bool BowtieHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
										BowtieHit& bh,
										bool strip_slash,
										char* name_out,
										char* name_tags)
{
	if (!orig_bwt_buf || !*orig_bwt_buf)
		return false;
	
	static const int buf_size = 2048;
	
	char bwt_buf[buf_size];
	strcpy(bwt_buf, orig_bwt_buf);
	//const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";

	char orientation;
	
	//int bwtf_ret = 0;
	//uint32_t seqid = 0;
	
	char text_name[buf_size];
	unsigned int text_offset;
	
	
	//unsigned int other_occs;
	char mismatches[buf_size];
	//memset(mismatches, 0, sizeof(mismatches));
	
	const char* buf = bwt_buf;
	char* name = strsep((char**)&buf,"\t");
	char* orientation_str = strsep((char**)&buf,"\t");
	char* text_name_str = strsep((char**)&buf,"\t");
	
	strcpy(text_name, text_name_str);
	
	char* text_offset_str = strsep((char**)&buf,"\t");
	char* seq_str = strsep((char**)&buf,"\t");
	/*const char* qual_str = */strsep((char**)&buf,"\t");
	/*const char* other_occs_str =*/ strsep((char**)&buf, "\t");
	mismatches[0] = 0;
	char* mismatches_str = strsep((char**)&buf, "\t");
	if (mismatches_str)
		strcpy(mismatches, mismatches_str);
	
	orientation = orientation_str[0];
	text_offset = atoi(text_offset_str);

	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
			strcpy(name_tags, pipe);
		*pipe = 0;
	}
	// Stripping the slash and number following it gives the insert name
	char* slash = strrchr(name, '/');
	if (strip_slash && slash)
		*slash = 0;
	
	size_t read_len = strlen(seq_str);
	
	// Add this alignment to the table of hits for this half of the
	// Bowtie map

	//vector<string> mismatch_toks;
	char* pch = strtok (mismatches,",");
	unsigned char num_mismatches = 0;
	while (pch != NULL)
	{
		char* colon = strchr(pch, ':');
		if (colon) 
		{
			num_mismatches++;
		}
		//mismatch_toks.push_back(pch);
		pch = strtok (NULL, ",");
	}
	
	bh = create_hit(name,
					text_name,
					text_offset, 
					(int)read_len, 
					orientation == '-',
					num_mismatches);
	
	return true;
}

int anchor_mismatch = 0;

bool SplicedBowtieHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
											   BowtieHit& bh,
											   bool strip_slash,
											   char* name_out,
											   char* name_tags)
{
	if (!orig_bwt_buf || !*orig_bwt_buf)
		return false;
	
	static const int buf_size = 2048;
	
	char bwt_buf[buf_size];
	strcpy(bwt_buf, orig_bwt_buf);
	//const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";

	char orientation;
	
	//int bwtf_ret = 0;
	//uint32_t seqid = 0;
	
	char text_name[buf_size];
	unsigned int text_offset;
	
	
	//unsigned int other_occs;
	char mismatches[buf_size];
	//memset(mismatches, 0, sizeof(mismatches));

	char* buf = bwt_buf;
	char* name = strsep((char**)&buf,"\t");
	char* orientation_str = strsep((char**)&buf,"\t");
	char* text_name_str = strsep((char**)&buf,"\t");
	strcpy(text_name, text_name_str);
	
	char* text_offset_str = strsep((char**)&buf,"\t");
	char* seq_str = strsep((char**)&buf,"\t");
	/*const char* qual_str = */strsep((char**)&buf,"\t");
	/*const char* other_occs_str =*/ strsep((char**)&buf, "\t");
	mismatches[0] = 0;
	char* mismatches_str = strsep((char**)&buf, "\t");
	if (mismatches_str)
		strcpy(mismatches, mismatches_str);
	
	orientation = orientation_str[0];
	text_offset = atoi(text_offset_str);
	
	int seg_offset = 0;
	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
			strcpy(name_tags, pipe);
		char* tag_buf = pipe + 1;
		char* colon = strchr(tag_buf, ':');
		if (colon)
		{
			*colon = 0;
			seg_offset = atoi(tag_buf);
		}
		*pipe = 0;
	}
	// Stripping the slash and number following it gives the insert name
	char* slash = strrchr(name, '/');
	if (strip_slash && slash)
		*slash = 0;
	
	//int read_len = strlen(sequence);
	
	// Add this alignment to the table of hits for this half of the
	// Bowtie map

	// Parse the text_name field to recover the splice coords
	vector<string> toks;
	
	tokenize_strict(text_name, "|", toks);
	
	int num_extra_toks = (int)toks.size() - 6;
	
	if (num_extra_toks >= 0)
	{
		static const uint8_t left_window_edge_field = 1;
		static const uint8_t splice_field = 2;
		//static const uint8_t right_window_edge_field = 3;
		//static const uint8_t junction_type_field = 4;
		static const uint8_t strand_field = 5;
		
		string contig = toks[0];
		for (int t = 1; t <= num_extra_toks; ++t)
		{
			contig += "|";
			contig += toks[t];
		}
		
		vector<string> splice_toks;
		tokenize(toks[num_extra_toks + splice_field], "-", splice_toks);
		if (splice_toks.size() != 2)
		{			
			fprintf(stderr, "Warning: found malformed splice record, skipping\n");
			return false;			       
		}
		
		uint32_t left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
		int spliced_read_len = strlen(seq_str);
		int8_t left_splice_pos = atoi(splice_toks[0].c_str()) - left + 1;
		int8_t right_splice_pos = spliced_read_len - left_splice_pos;
		
		if (orientation == '+')
		{
			if (left_splice_pos + seg_offset < _anchor_length)
				return false;
		}
		else
		{
			if (right_splice_pos + seg_offset < _anchor_length)
				return false;
		}
		//uint32_t right = atoi(splice_toks[1].c_str()) + right_splice_pos;
		//atoi(toks[right_window_edge_field].c_str());
		int gap_len = atoi(splice_toks[1].c_str()) - atoi(splice_toks[0].c_str()) - 1;
		
		string junction_strand = toks[num_extra_toks + strand_field];
		if (!(junction_strand == "rev" || junction_strand == "fwd")||
			!(orientation == '-' || orientation == '+'))
		{
			fprintf(stderr, "Warning: found malformed splice record, skipping\n");
			return false;
		}
		
		//vector<string> mismatch_toks;
		char* pch = strtok (mismatches,",");
		int mismatches_in_anchor = 0;
		unsigned char num_mismatches = 0;
		while (pch != NULL)
		{
			char* colon = strchr(pch, ':');
			if (colon) 
			{
				*colon = 0;
				num_mismatches++;
				int mismatch_pos = atoi(pch);
				if ((orientation == '+' && abs(mismatch_pos - left_splice_pos) < (int)min_anchor_len) ||
					(orientation == '-' && abs(((int)spliced_read_len - left_splice_pos + 1) - mismatch_pos)) < (int)min_anchor_len)
					mismatches_in_anchor++;
			}
			//mismatch_toks.push_back(pch);
			pch = strtok (NULL, ",");
		}
		
		// FIXME: we probably should exclude these hits somewhere, but this
		// isn't the right place
		vector<CigarOp> cigar;
		cigar.push_back(CigarOp(MATCH, left_splice_pos));
		cigar.push_back(CigarOp(REF_SKIP, gap_len));
		cigar.push_back(CigarOp(MATCH, right_splice_pos));
		bh = create_hit(name,
						contig,
						left, 
						cigar,
						orientation == '-', 
						junction_strand == "rev",
						num_mismatches,
						mismatches_in_anchor);
		return true;
	}
	else
	{
					fprintf(stderr, "Warning: found malformed splice record, skipping\n");
		//			continue;
		return false;
	}

	return false;
}

bool SAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
									 BowtieHit& bh,
									 bool strip_slash,
									 char* name_out,
									 char* name_tags)
{
	if (!orig_bwt_buf || !*orig_bwt_buf)
		return false;
	char bwt_buf[2048];
	strcpy(bwt_buf, orig_bwt_buf);
	// Are we still in the header region?
	if (bwt_buf[0] == '@')
		return false;
	
	char* buf = bwt_buf;
	char* name = strsep((char**)&buf,"\t");
	char* sam_flag_str = strsep((char**)&buf,"\t");
	char* text_name = strsep((char**)&buf,"\t");
	char* text_offset_str = strsep((char**)&buf,"\t");
	const char* map_qual_str = strsep((char**)&buf,"\t");
	char* cigar_str = strsep((char**)&buf,"\t");
	const char* mate_ref_str =  strsep((char**)&buf,"\t");
	const char* mate_pos_str =  strsep((char**)&buf,"\t");
	const char* inferred_insert_sz_str =  strsep((char**)&buf,"\t");
	const char* seq_str =  strsep((char**)&buf,"\t");
	const char* qual_str =  strsep((char**)&buf,"\t");
	
	if (!name ||
		!sam_flag_str ||
		!text_name ||
		!text_offset_str ||
		!map_qual_str ||
		!cigar_str ||
		!mate_ref_str ||
		!mate_pos_str ||
		!inferred_insert_sz_str ||
		!seq_str ||
		!qual_str)
	{
		// truncated or malformed SAM record
		return false;
	}
	
	
	int sam_flag = atoi(sam_flag_str);
	int text_offset = atoi(text_offset_str);
	
	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
			strcpy(name_tags, pipe);
		*pipe = 0;
	}
	// Stripping the slash and number following it gives the insert name
	char* slash = strrchr(name, '/');
	if (strip_slash && slash)
		*slash = 0;
	
	char* p_cig = cigar_str;
	//int len = strlen(sequence);
	vector<CigarOp> cigar;
	bool spliced_alignment = false;
	// Mostly pilfered direct from the SAM tools:
	while (*p_cig) 
	{
		char* t;
		int length = (int)strtol(p_cig, &t, 10);
		if (length <= 0)
		{
			//fprintf (stderr, "CIGAR op has zero length\n");
			return false;
		}
		char op_char = toupper(*t);
		CigarOpCode opcode;
		if (op_char == 'M') opcode = MATCH;
		else if (op_char == 'I') opcode = INS;
		else if (op_char == 'D') opcode = DEL;
		else if (op_char == 'N')
		{
			if (length > max_report_intron_length)
				return false;
			opcode = REF_SKIP;
			spliced_alignment = true;
		}
		else if (op_char == 'S') opcode = SOFT_CLIP;
		else if (op_char == 'H') opcode = HARD_CLIP;
		else if (op_char == 'P') opcode = PAD;
		else
		{
			fprintf (stderr, "invalid CIGAR operation\n");
			return false;
		}
		p_cig = t + 1;
		//i += length;
		cigar.push_back(CigarOp(opcode, length));
	}
	if (*p_cig)
	{
		fprintf (stderr, "unmatched CIGAR operation\n");
		return false;
	}
	
	//vector<string> attributes;
	//tokenize(tag_buf, " \t",attributes);
	
	bool antisense_splice = false;
	unsigned char num_mismatches = 0;
	unsigned char num_splice_anchor_mismatches = 0;
	const char* tag_buf = buf;
	
//	while((tag_buf = strsep((char**)&buf,"\t")))
//	{
//		vector<string> tuple_fields;
//		tokenize(tag_buf,":", tuple_fields);
//		if (tuple_fields.size() == 3)
//		{
//			if (tuple_fields[0] == "XS")
//			{
//				if (tuple_fields[2] == "-")
//					antisense_splice = true;
//			}
//			else if (tuple_fields[0] == "NM")
//			{
//				num_mismatches = atoi(tuple_fields[2].c_str());
//			}
//			else if (tuple_fields[0] == "NS")
//			{
//				num_splice_anchor_mismatches = atoi(tuple_fields[2].c_str());
//			}
//			else
//			{
//				fprintf(stderr, "%s attribute not supported\n", tuple_fields[0].c_str());
//				return false;
//			}
//		}
//	}
	
	while((tag_buf = strsep((char**)&buf,"\t")))
	{
		vector<string> tuple_fields;
		tokenize(tag_buf,":", tuple_fields);
		if (tuple_fields.size() == 3)
		{
			if (tuple_fields[0] == "XS")
			{
				if (tuple_fields[2] == "-")
					antisense_splice = true;
			}
			else if (tuple_fields[0] == "NM")
			{
				num_mismatches = atoi(tuple_fields[2].c_str());
			}
			else if (tuple_fields[0] == "NS")
			{
				//ignored for now
			}
			else
			{
				//fprintf(stderr, "%s attribute not supported\n", tuple_fields[0].c_str());
				//return false;
			}
		}
	}
	
	
//	vector<string> toks;
//	tokenize(tag_buf, ":", toks);
	if (spliced_alignment)
	{
		bh = create_hit(name,
						text_name,
						text_offset - 1, 
						cigar,
						sam_flag & 0x0010,
						antisense_splice,
						num_mismatches,
						num_splice_anchor_mismatches
						);
		return true;

	}
	else
	{		
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);

		bh = create_hit(name,
						text_name,
						text_offset - 1, // SAM files are 1-indexed 
						cigar[0].length, 
						sam_flag & 0x0010,
						num_mismatches);
		return true;
	}
		
	return false;
}

void get_mapped_reads(FILE* bwtf, 
					  HitTable& hits, 
					  HitFactory& hit_factory,
					  bool strip_slash,
					  bool verbose)
{

	
    char bwt_buf[2048];
	uint32_t reads_extracted = 0;
	
	while (fgets(bwt_buf, 2048, bwtf))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		if (*bwt_buf == 0)
			continue;
		// Get a new record from the tab-delimited Bowtie map
		BowtieHit bh;
		if (hit_factory.get_hit_from_buf(bwt_buf, bh, strip_slash))
		{
			// Only check uniqueness if these hits are spliced
			hits.add_hit(bh, true);
		}
		reads_extracted++;
	}
	
	// This will sort the map by insert id.
	hits.finalize();
	fprintf(stderr, "Extracted %d alignments from Bowtie map\n", reads_extracted);
}



AlignStatus status(const BowtieHit* align)
{
	if (!align)
		return UNALIGNED;
	if (align->contiguous())
		return CONTIGUOUS;
	return SPLICED;
}

void add_hits_to_coverage(const HitList& hits, vector<unsigned short>& DoC)
{
	int max_hit_pos = -1;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		max_hit_pos = max((int)hits[i].right(),max_hit_pos);
	}
	
	if ((int)DoC.size() < max_hit_pos)
		DoC.resize(max_hit_pos);
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		const BowtieHit& bh = hits[i];
		
		// split up the coverage contibution for this reads
		size_t j = bh.left();
		const vector<CigarOp>& cigar = bh.cigar();

		for (size_t c = 0 ; c < cigar.size(); ++c)
		{
			switch(cigar[c].opcode)
			{
				case MATCH:
					for (size_t m = 0; m < cigar[c].length; ++m)
					{
						if (DoC[j + m] < 0xFFFF)
							DoC[j + m]++;
					}
				//fall through this case to REF_SKIP is intentional
				case REF_SKIP:
					j += cigar[c].length;
					break;
				default:
					break;
			}
			 
		}
	}
}

void add_hit_to_coverage(const BowtieHit& bh, vector<unsigned int>& DoC)
{
	if ((int)DoC.size() < bh.right())
		DoC.resize(bh.right());
			
	// split up the coverage contibution for this reads
	size_t j = bh.left();
	const vector<CigarOp>& cigar = bh.cigar();
	
	for (size_t c = 0 ; c < cigar.size(); ++c)
	{
		switch(cigar[c].opcode)
		{
			case MATCH:
				for (size_t m = 0; m < cigar[c].length; ++m)
				{
					if (DoC[j + m] < 0xFFFFFFFF)
						DoC[j + m]++;
				}
				//fall through this case to REF_SKIP is intentional
			case REF_SKIP:
				j += cigar[c].length;
				break;
			default:
				break;
		}
		
	}
}

void print_hit(FILE* fout, 
			   const char* read_name,
			   const BowtieHit& bh,
			   const char* ref_name,
			   const char* sequence,
			   const char* qualities)
{
	//static const int buf_size = 256;
	//char text_name[buf_size];
	//char sequence[buf_size];
	
	string seq;
	string quals;
	if (sequence)
	{
		seq = sequence;
		quals = qualities;
		seq.resize(bh.read_len());
		quals.resize(bh.read_len());
		//if (bh.antisense_align())
		//	reverse_complement(seq);
	}
	else
	{
		seq = "*";
	}
	
	
	if (qualities)
	{
		quals = qualities;
		quals.resize(bh.read_len());
	}
	else
	{
		quals = "*";	
	}
	
	uint32_t sam_flag = 0;
	if (bh.antisense_align())
	{
		sam_flag |= 0x0010; // BAM_FREVERSE
		if (sequence)
		{
			reverse_complement(seq);
			reverse(quals.begin(), quals.end());
		}
	}
	
	uint32_t sam_pos = bh.left() + 1;
	uint32_t map_quality = 255;
	char cigar[256];
	cigar[0] = 0;
	string mate_ref_name = "*";
	uint32_t mate_pos = 0;
	uint32_t insert_size = 0;
	//string qualities = "*";
	
	const vector<CigarOp>& bh_cigar = bh.cigar();
	for (size_t c = 0; c < bh_cigar.size(); ++c)
	{
		char ibuf[64];
		sprintf(ibuf, "%d", bh_cigar[c].length);
		switch(bh_cigar[c].opcode)
		{
			case MATCH:
				strcat(cigar, ibuf);
				strcat(cigar, "M");
				break;
			case REF_SKIP:
				strcat(cigar, ibuf);
				strcat(cigar, "N");
				break;
			default:
				break;
		}
	}
	
	//string q = string(bh.read_len(), '!');
	//string s = string(bh.read_len(), 'N');
	
	fprintf(fout,
			"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			read_name,
			sam_flag,
			ref_name,
			sam_pos,
			map_quality,
			cigar,
			mate_ref_name.c_str(),
			mate_pos,
			insert_size,
			seq.c_str(),
			quals.c_str());
	
	fprintf(fout,
			"\tNM:i:%d",
			bh.edit_dist());
	
	if (!bh.contiguous())
	{
		fprintf(fout,
				"\tXS:A:%c",
				bh.antisense_splice() ? '-' : '+');
	}
	
	fprintf(fout, "\n");
}
