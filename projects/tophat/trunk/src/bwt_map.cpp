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
#include <cmath>

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

void LineHitFactory::openStream(HitStream& hs)
{
 if (hs._hit_file==NULL && !hs._hit_file_name.empty()) {
    //open the file for HitStream here
   hs._hit_file=fopen(hs._hit_file_name.c_str(),"r");
    if (hs._hit_file==NULL)
      err_die("Error opening HitStream file %s\n",hs._hit_file_name.c_str());
    return;
    }
 if (hs._fzpipe!=NULL) {
   hs._hit_file=hs._fzpipe->file;
   }
}

void LineHitFactory::rewind(HitStream& hs)
{
  if (hs._fzpipe!=NULL) {
        hs._fzpipe->rewind();
        hs._hit_file=hs._fzpipe->file;
        }
      else if (hs._hit_file) ::rewind((FILE*)(hs._hit_file));
}

bool LineHitFactory::next_record(HitStream& hs, const char*& buf, size_t& buf_size) {
     FILE* f=(FILE *)(hs._hit_file);
     bool new_rec = (fgets(_hit_buf,  _hit_buf_max_sz - 1, f)!=NULL);
     if (!new_rec || feof(f)) {
             hs._eof=true;
             return false;
             }
     ++_line_num;
     char* nl = strrchr(_hit_buf, '\n');
     if (nl) *nl = 0;
     buf = _hit_buf;
     buf_size = _hit_buf_max_sz - 1;
     return true;
     }

void LineHitFactory::closeStream(HitStream& hs) {
  if (hs._fzpipe!=NULL) {
     hs._fzpipe->close();
     return;
     }
  if (hs._hit_file!=NULL) {
    fclose((FILE*)(hs._hit_file));
    hs._hit_file=NULL;
    }
}
void BAMHitFactory::openStream(HitStream& hs) {
  if (hs._hit_file==NULL) {
     if (hs._hit_file_name.empty())
         //err_die("Error: invalid HitStream set for BAMHitFactory(file name missing)\n");
         return; //invalid stream, could be just a place holder
     //open the file here if not already open
     string fext=getFext(hs._hit_file_name);
     if (fext=="sam")
          hs._hit_file = samopen(hs._hit_file_name.c_str(), "r", 0);
        else
          hs._hit_file = samopen(hs._hit_file_name.c_str(), "rb", 0);

     samfile_t* sam_file=(samfile_t*)(hs._hit_file);

     if (sam_file == NULL)
         err_die("Error opening SAM file %s\n", hs._hit_file_name.c_str());
     if (sam_file->header == NULL)
         err_die("Error: no SAM header found for file %s\n", hs._hit_file_name.c_str());
     memset(&_next_hit, 0, sizeof(_next_hit));
     //_beginning = bgzf_tell(sam_file->x.bam);
     _sam_header=sam_file->header;
     if (inspect_header(hs) == false)
         err_die("Error: invalid SAM header for file %s\n",
                 hs._hit_file_name.c_str());
     }
}

void BAMHitFactory::closeStream(HitStream& hs) {
  if (hs._hit_file) {
        samclose((samfile_t*)(hs._hit_file));
        }
  hs._hit_file=NULL;
  _sam_header=NULL;
}

void BAMHitFactory::rewind(HitStream& hs)
{
  /*
    if (_hit_file && ((samfile_t*)_hit_file)->x.bam)
    {
        bgzf_seek(((samfile_t*)_hit_file)->x.bam, _beginning, SEEK_SET);
         _eof = false;
    }
    */
  this->closeStream(hs);
  this->openStream(hs);
}

string BAMHitFactory::hitfile_rec(HitStream& hs, const char* hit_buf) {
  const bam1_t* bamrec=(const bam1_t*)hit_buf;
  char* tamline=bam_format1(((samfile_t*)(hs._hit_file))->header, bamrec);
  string sam_line(tamline);
  free(tamline);
  return sam_line;
  }

bool BAMHitFactory::next_record(HitStream& hs, const char*& buf, size_t& buf_size) {
  if (_next_hit.data) {
      free(_next_hit.data);
      _next_hit.data = NULL;
      }
  _sam_header=((samfile_t*)(hs._hit_file))->header; //needed by get_hit_from_buf later on
  if (hs.eof() || !hs.ready()) return false;

  //mark_curr_pos();

  memset(&_next_hit, 0, sizeof(_next_hit));

  int bytes_read = samread((samfile_t*)(hs._hit_file), &_next_hit);
  if (bytes_read <= 0) {
      hs._eof = true;
      return false;
      }
  buf = (const char*)&_next_hit;
  buf_size = bytes_read;
  return true;
  }


BowtieHit HitFactory::create_hit(const string& insert_name, 
				 const string& ref_name,
				 int left,
				 const vector<CigarOp>& cigar,
				 bool antisense_aln,
				 bool antisense_splice,
				 unsigned char edit_dist,
				 unsigned char splice_mms,
				 bool end)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL, 0);
	
	return BowtieHit(reference_id,
			 insert_id, 
			 left, 
			 cigar, 
			 antisense_aln,
			 antisense_splice,
			 edit_dist,
			 splice_mms,
			 end);
}

BowtieHit HitFactory::create_hit(const string& insert_name, 
				 const string& ref_name,
				 uint32_t left,
				 uint32_t read_len,
				 bool antisense_aln,
				 unsigned char edit_dist,
				 bool end)
{
	uint64_t insert_id = _insert_table.get_id(insert_name);
	uint32_t reference_id = _ref_table.get_id(ref_name, NULL, 0);
	
	return BowtieHit(reference_id,
			 insert_id, 
			 left,
			 read_len,
			 antisense_aln,
			 edit_dist,
			 end);	
}

bool BowtieHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
					BowtieHit& bh,
					bool strip_slash,
					char* name_out,
					char* name_tags,
					char* seq,
					char* qual)
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
	char* name = get_token((char**)&buf,"\t");
	char* orientation_str = get_token((char**)&buf,"\t");
	char* text_name_str = get_token((char**)&buf,"\t");
	
	strcpy(text_name, text_name_str);
	
	char* text_offset_str = get_token((char**)&buf,"\t");
	char* seq_str = get_token((char**)&buf,"\t");
	if (seq)
	  strcpy(seq, seq_str);
	
	const char* qual_str = get_token((char**)&buf,"\t");
	if (qual)
	  strcpy(qual, qual_str);
	
	/*const char* other_occs_str =*/ get_token((char**)&buf, "\t");
	mismatches[0] = 0;
	char* mismatches_str = get_token((char**)&buf, "\t");
	if (mismatches_str)
		strcpy(mismatches, mismatches_str);
	
	orientation = orientation_str[0];
	text_offset = atoi(text_offset_str);

	bool end = true;
	unsigned int seg_offset = 0;
	unsigned int seg_num = 0;
	unsigned int num_segs = 0;
	
	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
		  strcpy(name_tags, pipe);

		char* tag_buf = pipe + 1;
		if (strchr(tag_buf, ':'))
		  {
		    sscanf(tag_buf, "%u:%u:%u", &seg_offset, &seg_num, &num_segs);
		    if (seg_num + 1 == num_segs)
		      end = true;
		    else
		      end = false;
		  }
		
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
			num_mismatches,
			end);
	
	return true;
}

int anchor_mismatch = 0;

bool SplicedBowtieHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
					       BowtieHit& bh,
					       bool strip_slash,
					       char* name_out,
					       char* name_tags,
					       char* seq,
					       char* qual)
{
	if (!orig_bwt_buf || !*orig_bwt_buf)
		return false;
	
	static const int buf_size = 2048;
	
	char bwt_buf[buf_size];
	strcpy(bwt_buf, orig_bwt_buf);
	//const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";

	char orientation;
	char text_name[buf_size];
	unsigned int text_offset;
	char mismatches[buf_size];
	//memset(mismatches, 0, sizeof(mismatches));

	char* buf = bwt_buf;
	char* name = get_token((char**)&buf,"\t");
	char* orientation_str = get_token((char**)&buf,"\t");
	char* text_name_str = get_token((char**)&buf,"\t");
	strcpy(text_name, text_name_str);
	
	char* text_offset_str = get_token((char**)&buf,"\t");
	char* seq_str = get_token((char**)&buf,"\t");
	if (seq)
	  strcpy(seq, seq_str);
	
	const char* qual_str = get_token((char**)&buf,"\t");
	if (qual)
	  strcpy(qual, qual_str);
	
	/*const char* other_occs_str =*/ get_token((char**)&buf, "\t");
	mismatches[0] = 0;
	char* mismatches_str = get_token((char**)&buf, "\t");
	if (mismatches_str)
		strcpy(mismatches, mismatches_str);
	
	orientation = orientation_str[0];
	text_offset = atoi(text_offset_str);
	
	bool end = true;
	unsigned int seg_offset = 0;
	unsigned int seg_num = 0;
	unsigned int num_segs = 0;

	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
		  strcpy(name_tags, pipe);

		char* tag_buf = pipe + 1;
		if (strchr(tag_buf, ':'))
		  {
		    sscanf(tag_buf, "%u:%u:%u", &seg_offset, &seg_num, &num_segs);
		    if (seg_num + 1 == num_segs)
		      end = true;
		    else
		      end = false;
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
		static const uint8_t junction_type_field = 4;
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
			fprintf(stderr, "Warning: found malformed splice record, skipping:\n");
			//fprintf(stderr, "%s (token: %s)\n", text_name, 
			//        toks[num_extra_toks + splice_field].c_str());
			return false;			
		}

		//
		// check for an insertion hit
		//
		if(toks[num_extra_toks + junction_type_field] == "ins")
		  {
			int8_t spliced_read_len = strlen(seq_str);
			/*
			 * The 0-based position of the left edge of the alignment. Note that this
  			 * value may need to be futher corrected to account for the presence of
			 * of the insertion.  
			*/
			uint32_t left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
			uint32_t right = left + spliced_read_len - 1;

			/*
			 * The 0-based position of the last genomic sequence before the insertion
			 */
			uint32_t left_splice_pos = atoi(splice_toks[0].c_str());
		
			string insertedSequence = splice_toks[1];
			/*
			 * The 0-based position of the first genomic sequence after teh insertion
			 */
			uint32_t right_splice_pos = left_splice_pos + 1; 
			if(left > left_splice_pos){
				/*
				 * The genomic position of the left edge of the alignment needs to be corrected
				 * If the alignment does not extend into the insertion, simply subtract the length
				 * of the inserted sequence, otherwise, just set it equal to the right edge
				 */
				left = left - insertedSequence.length();
				if(left < right_splice_pos){
					left = right_splice_pos;
				}
			}
			if(right > left_splice_pos){
				right = right - insertedSequence.length();
				if(right < left_splice_pos){
					right = left_splice_pos;
				}
			}
			/*
			 * Now, right and left should be properly transformed into genomic coordinates
			 * We should be able to deduce how much the alignment matches the insertion
			 * simply based on the length of the read
			 */
			int left_match_length = 0;
			if(left <= left_splice_pos){
				left_match_length = left_splice_pos - left + 1;
			}
			int right_match_length = 0;
			if(right >= right_splice_pos){
				right_match_length = right - right_splice_pos + 1;
			}
			int insertion_match_length = spliced_read_len - left_match_length - right_match_length;

			if(left_match_length <= 0 || right_match_length <= 0 || insertion_match_length <= 0)
			  return false;

			string junction_strand = toks[num_extra_toks + strand_field];
			if(junction_strand != "rev" && junction_strand != "fwd"){
				fprintf(stderr, "Malformed insertion record\n");
				return false;
			}
			
			char* pch = strtok( mismatches, ",");
			unsigned char num_mismatches = 0;

			/*
			 * remember that text_offset holds the left end of the 
			 * alignment, relative to the start of the contig
			 */

			/*
			 * The 0-based relative position of the left-most character
			 * before the insertion in the contig
			 */
			int relative_splice_pos = left_splice_pos - atoi(toks[num_extra_toks + left_window_edge_field].c_str()); 
			while (pch != NULL)
			{
				char* colon = strchr(pch, ':');
				if (colon) 
				{
					*colon = 0;
					int mismatch_pos = atoi(pch);

					/*
					 * for reversely mapped reads,
					 * find the correct mismatched position.
					 */
					if(orientation == '-'){
					  mismatch_pos = spliced_read_len - mismatch_pos - 1;
					}

					/*
					 * Only count mismatches outside of the insertion region
					 * If there is a mismatch within the insertion,
					 * disallow this hit 
					 */
					if(mismatch_pos + text_offset <= relative_splice_pos || mismatch_pos + text_offset > relative_splice_pos + insertedSequence.length()){
					  num_mismatches++;
					}else{
					  return false; 
					}
				}
				pch = strtok (NULL, ",");
			}
			
			
			vector<CigarOp> cigar;
			cigar.push_back(CigarOp(MATCH, left_match_length));
			cigar.push_back(CigarOp(INS, insertion_match_length)); 
			cigar.push_back(CigarOp(MATCH, right_match_length)); 

			/*
			 * For now, disallow hits that don't span
			 * the insertion. If we allow these types of hits,
			 * then long_spanning.cpp needs to be updated
			 * in order to intelligently merge these kinds
			 * of reads back together
			 * 
			 * Following code has been changed to allow segment that end
			 * in an insertion
			 */
			bh = create_hit(name,
					contig,
					left, 
					cigar,
					orientation == '-', 
					junction_strand == "rev",
					num_mismatches,
					0,
					end);
			return true;
		}	

		else
		  {
		    uint32_t left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
		    int spliced_read_len = strlen(seq_str);
		    int8_t left_splice_pos = atoi(splice_toks[0].c_str()) - left + 1;
		    if(left_splice_pos > spliced_read_len) left_splice_pos = spliced_read_len;		  
		    int8_t right_splice_pos = spliced_read_len - left_splice_pos;
		    
		    if (right_splice_pos <= 0 || left_splice_pos <= 0)
		      return false;
		    
		    if (orientation == '+')
		      {
			if (left_splice_pos + seg_offset < _anchor_length){
			  return false;
			}
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
			//fprintf(stderr, "junction_strand=%s, orientation='%c'\n",
			//           junction_strand.c_str(), orientation);
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
		    if(toks[num_extra_toks + junction_type_field] == "del"){
		      cigar.push_back(CigarOp(DEL, gap_len));
		    }else{
		      cigar.push_back(CigarOp(REF_SKIP, gap_len));
		    }
		    cigar.push_back(CigarOp(MATCH, right_splice_pos));
		    
		    bh = create_hit(name,
				    contig,
				    left, 
				    cigar,
				    orientation == '-', 
				    junction_strand == "rev",
				    num_mismatches,
				    mismatches_in_anchor,
				    end);
		    return true;
		  }
	}
	else
	{
	  fprintf(stderr, "Warning: found malformed splice record, skipping\n");
	  //fprintf(stderr, "%s\n", orig_bwt_buf);
	  //			continue;
		return false;
	}

	return false;
}

bool SAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
				     BowtieHit& bh,
				     bool strip_slash,
				     char* name_out,
				     char* name_tags,
				     char* seq,
				     char* qual)
{
	if (!orig_bwt_buf || !*orig_bwt_buf)
		return false;
	char bwt_buf[2048];
	strcpy(bwt_buf, orig_bwt_buf);
	// Are we still in the header region?
	if (bwt_buf[0] == '@')
		return false;
	
	char* buf = bwt_buf;
	char* name = get_token((char**)&buf,"\t");
	char* sam_flag_str = get_token((char**)&buf,"\t");
	char* text_name = get_token((char**)&buf,"\t");
	char* text_offset_str = get_token((char**)&buf,"\t");
	const char* map_qual_str = get_token((char**)&buf,"\t");
	char* cigar_str = get_token((char**)&buf,"\t");
	const char* mate_ref_str =  get_token((char**)&buf,"\t");
	const char* mate_pos_str =  get_token((char**)&buf,"\t");
	const char* inferred_insert_sz_str =  get_token((char**)&buf,"\t");
	
	const char* seq_str =  get_token((char**)&buf,"\t");
	if (seq)
	  strcpy(seq, seq_str);
	
	const char* qual_str =  get_token((char**)&buf,"\t");
	if (qual)
	  strcpy(qual, qual_str);
	
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

	bool end = true;
	unsigned int seg_offset = 0;
	unsigned int seg_num = 0;
	unsigned int num_segs = 0;

	// Copy the tag out of the name field before we might wipe it out
	char* pipe = strrchr(name, '|');
	if (pipe)
	{
		if (name_tags)
			strcpy(name_tags, pipe);

		char* tag_buf = pipe + 1;
		if (strchr(tag_buf, ':'))
		  {
		    sscanf(tag_buf, "%u:%u:%u", &seg_offset, &seg_num, &num_segs);
		    if (seg_num + 1 == num_segs)
		      end = true;
		    else
		      end = false;
		  }

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
	
//	while((tag_buf = get_token((char**)&buf,"\t")))
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
	
	while((tag_buf = get_token((char**)&buf,"\t")))
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
	

	/*
	 * By convention,the NM field of the SAM record
	 * counts an insertion or deletion. I dont' think
	 * we want the mismatch count in the BowtieHit
	 * record to reflect this. Therefore, subtract out
	 * the mismatches due to in/dels
	 */
	for(vector<CigarOp>::const_iterator itr = cigar.begin(); itr != cigar.end(); ++itr){
		if(itr->opcode == INS){
			num_mismatches -= itr->length;
		}
		if(itr->opcode == DEL){
			num_mismatches -= itr->length;
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
				num_splice_anchor_mismatches,
				end);
		return true;

	}
	else
	{		
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);

		bh = create_hit(name,
				text_name,
				text_offset - 1, // SAM files are 1-indexed 
				cigar,
				sam_flag & 0x0010,
				false,
				num_mismatches,
				0,
				end);
		return true;
	}
		
	return false;
}



bool BAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
				 BowtieHit& bh, bool strip_slash,
				char* name_out, char* name_tags,
				char* seq, char* qual) {
    if (_sam_header==NULL)
      err_die("Error: no SAM header when BAMHitFactory::get_hit_from_buf()!");
	const bam1_t* hit_buf = (const bam1_t*)orig_bwt_buf;


	uint32_t sam_flag = hit_buf->core.flag;
	
	int text_offset = hit_buf->core.pos;
	int text_mate_pos = hit_buf->core.mpos;
	int target_id = hit_buf->core.tid;
	int mate_target_id = hit_buf->core.mtid;
	
	vector<CigarOp> cigar;
	bool spliced_alignment = false;
	int num_hits = 1;
	
	double mapQ = hit_buf->core.qual;
	long double error_prob;
	if (mapQ > 0)
	{
		long double p = (-1.0 * mapQ) / 10.0;
		error_prob = pow(10.0L, p);
	}
	else
	{
		error_prob = 1.0;
	}

	//header->target_name[c->tid]
	
    bool end = true;
    unsigned int seg_offset = 0;
    unsigned int seg_num = 0;
    unsigned int num_segs = 0;
    // Copy the tag out of the name field before we might wipe it out
    char* qname = bam1_qname(hit_buf);
    char* pipe = strrchr(qname, '|');
    if (pipe)
    {
        if (name_tags)
            strcpy(name_tags, pipe);

        char* tag_buf = pipe + 1;
        if (strchr(tag_buf, ':'))
          {
            sscanf(tag_buf, "%u:%u:%u", &seg_offset, &seg_num, &num_segs);
            if (seg_num + 1 == num_segs)
              end = true;
            else
              end = false;
          }

        *pipe = 0;
    }


	if (target_id < 0)	{
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
	    bh = create_hit(qname,
	                "*", //ref_name
	                0, //left coord
	                0, //read_len
	                false, //antisense_aln
	                0, //edit_dist
	                end);
		return true;
	}
	
	//string text_name = ((samfile_t*)(hs._hit_file))->header->target_name[target_id];
	string text_name = _sam_header->target_name[target_id];
	for (int i = 0; i < hit_buf->core.n_cigar; ++i) 
	{
		//char* t;

		int length = bam1_cigar(hit_buf)[i] >> BAM_CIGAR_SHIFT;
		if (length <= 0)
		{
			fprintf (stderr, "BAM error: CIGAR op has zero length\n");
			return false;
		}
		
		CigarOpCode opcode;
		switch(bam1_cigar(hit_buf)[i] & BAM_CIGAR_MASK)
		{
			case BAM_CMATCH: opcode  = MATCH; break; 
			case BAM_CINS: opcode  = INS; break;
			case BAM_CDEL: opcode  = DEL; break; 
			case BAM_CSOFT_CLIP: opcode  = SOFT_CLIP; break;
			case BAM_CHARD_CLIP: opcode  = HARD_CLIP; break;
			case BAM_CPAD: opcode  = PAD; break; 
			case BAM_CREF_SKIP:
				opcode = REF_SKIP;
				spliced_alignment = true;
				if (length > (int)max_report_intron_length)
				{
					//fprintf(stderr, "Encounter REF_SKIP > max_gene_length, skipping\n");
					return false;
				}
				break;
			default:
				fprintf (stderr, "BAM read: invalid CIGAR operation\n");
				return false;
		}
		if (opcode != HARD_CLIP)
			cigar.push_back(CigarOp(opcode, length));
	}
	
	string mrnm;
	if (mate_target_id >= 0) {
		if (mate_target_id == target_id) {
			//mrnm = ((samfile_t*)(hs._hit_file))->header->target_name[mate_target_id];
		    mrnm = _sam_header->target_name[mate_target_id];
		    }
		else {
			//fprintf(stderr, "Trans-spliced mates are not currently supported, skipping\n");
			return false;
		    }
	   }
	else {
	   text_mate_pos = 0;
	   }
	//CuffStrand source_strand = CUFF_STRAND_UNKNOWN;
	bool antisense_splice=false;
	unsigned char num_mismatches = 0;
    unsigned char num_splice_anchor_mismatches = 0;

	uint8_t* ptr = bam_aux_get(hit_buf, "XS");
	if (ptr) {
		char src_strand_char = bam_aux2A(ptr);
		if (src_strand_char == '-')
			antisense_splice = true;
		// else if (src_strand_char == '+')
		//	source_strand = CUFF_FWD;
	    }

	ptr = bam_aux_get(hit_buf, "NM");
	if (ptr) {
		num_mismatches = bam_aux2i(ptr);
		}

	ptr = bam_aux_get(hit_buf, "NH");
	if (ptr) {
		num_hits = bam_aux2i(ptr);
		}
	
    //bool antisense_aln = bam1_strand(hit_buf);
    
    //if (_rg_props.strandedness() == STRANDED_PROTOCOL && source_strand == CUFF_STRAND_UNKNOWN)
	//	source_strand = use_stranded_protocol(sam_flag, antisense_aln, _rg_props.mate_strand_mapping());
	if (spliced_alignment)	{
      //if (source_strand == CUFF_STRAND_UNKNOWN) {
      //  fprintf(stderr, "BAM record error: found spliced alignment without XS attribute\n");
      //  }
      bh = create_hit(qname,
                      text_name,
                      text_offset,  // BAM files are 0-indexed
                      cigar,
                      sam_flag & 0x0010,
                      antisense_splice,
                      num_mismatches,
                      num_splice_anchor_mismatches,
                      end);

		}
	else {
      //assert(_rg_props.strandedness() == STRANDED_PROTOCOL || source_strand == CUFF_STRAND_UNKNOWN);
      //assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
      bh = create_hit(qname,
                        text_name,
                        text_offset,  // BAM files are 0-indexed
                        cigar,
                        sam_flag & 0x0010,
                        false,
                        num_mismatches,
                        0,
                        end);
     }
    if (seq!=NULL) {
       char *bseq = (char*)bam1_seq(hit_buf);
       for(int i=0;i<(hit_buf->core.l_qseq);i++) {
         char v = bam1_seqi(bseq,i);
         seq[i]=bam_nt16_rev_table[v];
         }
       seq[hit_buf->core.l_qseq]=0;
       }
    if (qual!=NULL) {
       char *bq  = (char*)bam1_qual(hit_buf);
       for(int i=0;i<(hit_buf->core.l_qseq);i++) {
          qual[i]=bq[i]+33;
          }
       qual[hit_buf->core.l_qseq]=0;
       }
    return true;
}

bool BAMHitFactory::inspect_header(HitStream& hs)
{
    bam_header_t* header = ((samfile_t*)(hs._hit_file))->header;
    
    if (header == NULL) {
       fprintf(stderr, "Warning: No BAM header\n");
       return false;
       }
    if (header->l_text == 0) {
      fprintf(stderr, "Warning: BAM header has 0 length or is corrupted.  Try using 'samtools reheader'.\n");
      return false;
      }
    return true;
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


/*
AlignStatus status(const BowtieHit* align)
{
	if (!align)
		return UNALIGNED;
	if (align->contiguous())
		return CONTIGUOUS;
	return SPLICED;
}
*/

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
					if (DoC[j + m] < VMAXINT32)
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
	       const char* qualities,
	       bool from_bowtie)
{
	string seq;
	string quals;
	if (sequence)
	{
		seq = sequence;
		quals = qualities;
		seq.resize(bh.read_len());
		quals.resize(bh.read_len());
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
		if (sequence && !from_bowtie)  // if it is from bowtie hit, it's already reversed.
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

	/*
	 * In addition to calculating the cigar string,
	 * we need to figure out how many in/dels are in the 
	 * sequence, so that we can give the correct
	 * value for the NM tag
	 */
	int indel_distance = 0;
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
			case INS:
				strcat(cigar, ibuf);
				strcat(cigar, "I");
				indel_distance += bh_cigar[c].length;
				break;
			case DEL:
				strcat(cigar, ibuf);
				strcat(cigar, "D");
				indel_distance += bh_cigar[c].length;
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
	
    if (!sam_readgroup_id.empty())
    	fprintf(fout, "\tRG:Z:%s", sam_readgroup_id.c_str());
    
	fprintf(fout, "\tNM:i:%d", bh.edit_dist() + indel_distance);

	bool containsSplice = false;
	for(vector<CigarOp>::const_iterator itr = bh.cigar().begin(); itr != bh.cigar().end(); ++itr){
		if(itr->opcode == REF_SKIP){
			containsSplice = true;
			break;
		}
	}

	if (containsSplice)
	  fprintf(fout, "\tXS:A:%c", bh.antisense_splice() ? '-' : '+');

	fprintf(fout, "\n");
}

void print_bamhit(GBamWriter& wbam,
           const char* read_name,
           const BowtieHit& bh,
           const char* ref_name,
           const char* sequence,
           const char* qualities,
           bool from_bowtie)
{
    string seq;
    string quals;
    if (sequence) {
        seq = sequence;
        quals = qualities;
        seq.resize(bh.read_len());
        quals.resize(bh.read_len());
        }
      else {
        seq = "*";
        }
    if (qualities) {
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
        if (sequence && !from_bowtie)  // if it is from bowtie hit, it's already reversed.
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
    /*
     * In addition to calculating the cigar string,
     * we need to figure out how many in/dels are in the
     * sequence, so that we can give the correct
     * value for the NM tag
     */
    int indel_distance = 0;
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
            case INS:
                strcat(cigar, ibuf);
                strcat(cigar, "I");
                indel_distance += bh_cigar[c].length;
                break;
            case DEL:
                strcat(cigar, ibuf);
                strcat(cigar, "D");
                indel_distance += bh_cigar[c].length;
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
    /*
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

    fprintf(fout, "\tNM:i:%d", bh.edit_dist() + indel_distance);
    */
    vector<string> auxdata;
    
    if (!sam_readgroup_id.empty())
    {
        string nm("RG:Z:");
        nm += sam_readgroup_id;
        auxdata.push_back(nm);
    }
    
    string nm("NM:i:");
    str_appendInt(nm, bh.edit_dist() + indel_distance);
    auxdata.push_back(nm);
    bool containsSplice = false;
    for(vector<CigarOp>::const_iterator itr = bh.cigar().begin(); itr != bh.cigar().end(); ++itr)
       if(itr->opcode == REF_SKIP){
            containsSplice = true;
            break;
            }

    if (containsSplice) {
       //fprintf(fout, "\tXS:A:%c", bh.antisense_splice() ? '-' : '+');
       nm="XS:A:";
       nm+=(char)(bh.antisense_splice() ? '-' : '+');
       auxdata.push_back(nm);
       }

    GBamRecord *brec = wbam.new_record(read_name, sam_flag, ref_name, sam_pos, map_quality,
                        cigar, mate_ref_name.c_str(), mate_pos,
                        insert_size, seq.c_str(), quals.c_str(), &auxdata);
                        
    
    
    wbam.write(brec);
    delete brec;
}

/**
 * Print a vector of cigar operations to a file.
 * @param bh_cigar A vector of CigarOps
 * @return a string representation of the cigar string
*/
std::string print_cigar(vector<CigarOp>& bh_cigar){
	char cigar[256];
	cigar[0] = 0;	
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
			case INS:
				strcat(cigar, ibuf);
				strcat(cigar, "I");
				break;
			case DEL:
				strcat(cigar, ibuf);
				strcat(cigar, "D");
				break;
			case REF_SKIP:
				strcat(cigar, ibuf);
				strcat(cigar, "N");
				break;
			default:
				break;
		}
	}
	string result(cigar);
	return result;
}

bool BowtieHit::check_editdist_consistency(const RefSequenceTable& rt)
{
  RefSequenceTable::Sequence* ref_str = rt.get_seq(_ref_id);
  if (!ref_str)
    return false;
  
  const seqan::Dna5String ref_seq = seqan::infix(*ref_str, _left, right());

  size_t pos_seq = 0;
  size_t pos_ref = 0;
  size_t mismatch = 0;
  for (size_t i = 0; i < _cigar.size(); ++i)
    {
      CigarOp cigar = _cigar[i];
      switch(cigar.opcode)
	{
	case MATCH:
	  {
	    for (size_t j = 0; j < cigar.length; ++j)
	      {
		if (_seq[pos_seq++] != ref_seq[pos_ref++])
		  ++mismatch;
	      }
	  }
	  break;
	case INS:
	  {
	    pos_seq += cigar.length;
	  }
	  break;
	  
	case DEL:
	case REF_SKIP:
	  {
	    pos_ref += cigar.length;
	  }
	  break;

	default:
	  break;
	}
    }
  
  return mismatch == _edit_dist;
}
