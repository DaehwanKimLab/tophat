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

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>

#include "common.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "reads.h"
#include "align_status.h"

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

void LineHitFactory::seek(HitStream& hs, int64_t offset)
{
  // daehwan - implement this later
  if (hs._fzpipe != NULL) {
    hs._fzpipe->seek(offset);
    hs._hit_file=hs._fzpipe->file;
  }
  // else if (hs._hit_file) ::seek((FILE*)(hs._hit_file));
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

void BAMHitFactory::seek(HitStream& hs, int64_t offset)
{
  if (hs._hit_file) {
    bgzf_seek(((samfile_t*)hs._hit_file)->x.bam, offset, SEEK_SET);
  }
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
				 const string& ref_name2,
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
  uint32_t reference_id2 = reference_id;
  
  if (ref_name2.length() > 0)
    reference_id2 = _ref_table.get_id(ref_name2, NULL, 0);
  
  return BowtieHit(reference_id,
		   reference_id2,
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
		   reference_id,
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


void parseSegReadName(char* name, char*& name_tags, bool strip_slash,
		bool &end, unsigned int &seg_offset, unsigned int& seg_num,
		                                   unsigned int & num_segs) {
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
}


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
	parseSegReadName(name, name_tags, strip_slash, end, seg_offset, seg_num, num_segs);

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
					"",
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
		    const string& junction_type = toks[num_extra_toks + junction_type_field];
		    string junction_strand = toks[num_extra_toks + strand_field];

		    int spliced_read_len = strlen(seq_str);
		    uint32_t left = atoi(toks[num_extra_toks + left_window_edge_field].c_str());
		    int8_t left_splice_pos = atoi(splice_toks[0].c_str());
		    if (junction_type != "fus" || (junction_strand != "rf" && junction_strand != "rr"))
		      {
			left += text_offset;
			left_splice_pos = left_splice_pos - left + 1;
		      }
		    else
		      {
			left -= text_offset;
			left_splice_pos = left - left_splice_pos + 1;
		      }

		    if(left_splice_pos > spliced_read_len) left_splice_pos = spliced_read_len;		  
		    int8_t right_splice_pos = spliced_read_len - left_splice_pos;

		    int gap_len = 0;
		    if (junction_type == "fus")
		      gap_len = atoi(splice_toks[1].c_str());
		    else
		      gap_len = atoi(splice_toks[1].c_str()) - atoi(splice_toks[0].c_str()) - 1;
	    
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

		    if (!(junction_strand == "ff" || junction_strand == "fr" || junction_strand == "rf" || junction_strand == "rr" || junction_strand == "rev" || junction_strand == "fwd")||
			!(orientation == '-' || orientation == '+'))
		      {
			fprintf(stderr, "Warning: found malformed splice record, skipping\n");
			//fprintf(stderr, "junction_strand=%s, orientation='%c'\n",
			//           junction_strand.c_str(), orientation);
			return false;
		      }

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
			pch = strtok (NULL, ",");
		      }
		    
		    // FIXME: we probably should exclude these hits somewhere, but this
		    // isn't the right place
		    vector<CigarOp> cigar;
		    if (junction_type != "fus" || (junction_strand != "rf" && junction_strand != "rr"))
		      cigar.push_back(CigarOp(MATCH, left_splice_pos));
		    else
		      cigar.push_back(CigarOp(mATCH, left_splice_pos));
		    
		    if(junction_type == "del")
		      cigar.push_back(CigarOp(DEL, gap_len));
		    else if(junction_type == "fus")
		      {
			if (junction_strand == "ff")
			  cigar.push_back(CigarOp(FUSION_FF, gap_len));
			else if (junction_strand == "fr")
			  cigar.push_back(CigarOp(FUSION_FR, gap_len));
			else if (junction_strand == "rf")
			  cigar.push_back(CigarOp(FUSION_RF, gap_len));
			else
			  cigar.push_back(CigarOp(FUSION_RR, gap_len));
		      }
		    else
		      cigar.push_back(CigarOp(REF_SKIP, gap_len));

		    if (junction_type != "fus" || (junction_strand != "fr" && junction_strand != "rr"))
		      cigar.push_back(CigarOp(MATCH, right_splice_pos));
		    else
		      cigar.push_back(CigarOp(mATCH, right_splice_pos));

		    string contig2 = ""; 
		    if (junction_type == "fus")
		      {
			vector<string> contigs;
			tokenize(contig, "-", contigs);
			if (contigs.size() != 2)
			  return false;

			contig = contigs[0];
			contig2 = contigs[1];

			if (junction_strand == "rf" || junction_strand == "rr")
			  orientation = (orientation == '+' ? '-' : '+');
		      }

		    bh = create_hit(name,
				    contig,
				    contig2,
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


int parseCigar(vector<CigarOp>& cigar, const char* cigar_str,
		bool &spliced_alignment)   {
const char* p_cig = cigar_str;
int refspan=0; //alignment span on reference sequence

while (*p_cig)
{
	char* t;
	int op_len = (int)strtol(p_cig, &t, 10);
	if (op_len <= 0)
	{
		fprintf (stderr, "Error: CIGAR op has zero length\n");
		return 0;
	}
	char op_char = toupper(*t);
	CigarOpCode opcode;
	switch (op_char) {
	case '=':
	case 'X':
	case 'M': opcode = MATCH;
	  refspan+=op_len;
	  break;
	case 'I': opcode = INS;
	  break;
	case 'D': opcode = DEL;
	  refspan+=op_len;
	  break;
	case 'N': if (op_len > max_report_intron_length)
	    return 0;
	  opcode = REF_SKIP;
	  spliced_alignment = true;
	  refspan+=op_len;
	  break;
	case 'S': opcode = SOFT_CLIP;
	  break;
	case 'H': opcode = HARD_CLIP;
	  break;
	case 'P': opcode = PAD;
	  break;
	default:  fprintf (stderr, "Error: invalid CIGAR operation\n");
	  return 0;
	  }
	p_cig = t + 1;
	cigar.push_back(CigarOp(opcode, op_len));
} //while cigar codes
 if (*p_cig) {
	  fprintf (stderr, "Error: unmatched CIGAR operation (%s in %s)\n",
		   p_cig, cigar_str);
	  return 0;
    }
return refspan;
}

int getBAMmismatches(const bam1_t* buf, vector<CigarOp>& cigar,
		     vector<bool>& mismatches, int& sam_nm, bool& antisense_splice) {
  int gspan=0;//genomic span of the alignment
  sam_nm=0;
  int num_mismatches=0;
  
  uint8_t* ptr = bam_aux_get(buf, "XS");
  if (ptr) {
    char src_strand_char = bam_aux2A(ptr);
    if (src_strand_char == '-')
      antisense_splice = true;
  }

  ptr = bam_aux_get(buf, "MD");
  if (ptr) {
    const char* p = bam_aux2Z(ptr);
    int bi=0; //base offset position in the read
    while (*p != 0) {
      if (isdigit(*p)) {
	int v=atoi(p);
	do { p++; } while (isdigit(*p));
	bi+=v;
      }
      while (isalpha(*p)) {
	p++;
	num_mismatches++;
	//mismatches.push_back(bi);
	mismatches[bi]=true;
	bi++;
      }
      if (*p=='^') { //reference deletion
	p++;
	while (isalpha(*p)) { //insert read bases
	  p++; bi++;
	}
      }
    }
  }

  /* By convention,the NM field of the SAM record
   *  counts an insertion or deletion. I dont' think
   *  we want the mismatch count in the BowtieHit
   *  record to reflect this. Therefore, subtract out
   *  the mismatches due to in/dels
   */
  for(vector<CigarOp>::const_iterator itr = cigar.begin(); itr != cigar.end(); ++itr){
    switch (itr->opcode)
      {
      case MATCH:
      case REF_SKIP:
      case PAD:
        gspan += itr->length;
        break;
      case DEL:
        gspan += itr->length;
        sam_nm -= itr->length;
        break;
      case INS:
        sam_nm -= itr->length;
        break;
      default:
	break;
      }
  }
  return num_mismatches;
}

int getSAMmismatches(char* &buf, vector<CigarOp>& cigar,
		  vector<bool>& mismatches, int& sam_nm, bool& antisense_splice) {
  int gspan=0;//genomic span of the alignment
  const char* tag_buf = buf;
	sam_nm=0;
	int num_mismatches=0;
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
				sam_nm = atoi(tuple_fields[2].c_str());
			}
			else if (tuple_fields[0] == "NS")
			{
				//ignored for now
			}
			else if (tuple_fields[0] == "MD")
			{
              const char* p=tuple_fields[2].c_str();
              int bi=0; //base offset position in the read
              while (*p != 0) {
            	if (isdigit(*p)) {
            	  int v=atoi(p);
            	  do { p++; } while (isdigit(*p));
            	  bi+=v;
            	  }
            	 while (isalpha(*p)) {
            	  p++;
            	  num_mismatches++;
            	  //mismatches.push_back(bi);
            	  mismatches[bi]=true;
            	  bi++;
            	  }
            	 if (*p=='^') { //reference deletion
                   p++;
                   while (isalpha(*p)) { //insert read bases
                	  p++; bi++;
                    }
            	   }
                }
			}
     //else
			//{
				//fprintf(stderr, "%s attribute not supported\n", tuple_fields[0].c_str());
				//return false;
			//}
		}
	}
	 /* By convention,the NM field of the SAM record
	 *  counts an insertion or deletion. I dont' think
	 *  we want the mismatch count in the BowtieHit
	 *  record to reflect this. Therefore, subtract out
	 *  the mismatches due to in/dels
	 */
	for(vector<CigarOp>::const_iterator itr = cigar.begin(); itr != cigar.end(); ++itr){
    switch (itr->opcode)
    {
      case MATCH:
      case REF_SKIP:
      case PAD:
        gspan += itr->length;
        break;
      case DEL:
        gspan += itr->length;
        sam_nm -= itr->length;
        break;
      case INS:
        sam_nm -= itr->length;
        break;
    default:
      break;
		}
	 }
	return num_mismatches;
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
	string ref_name = text_name, ref_name2 = "";
	int text_offset = atoi(text_offset_str);

	bool end = true;
	unsigned int seg_offset = 0;
	unsigned int seg_num = 0;
	unsigned int num_segs = 0;

	// Copy the tag out of the name field before we might wipe it out
	parseSegReadName(name, name_tags, strip_slash, end, seg_offset, seg_num, num_segs);

	vector<CigarOp> cigar;
	bool spliced_alignment = false;

	int refspan=parseCigar(cigar, cigar_str, spliced_alignment);
	if (refspan==0)
	   return false;
	//vector<string> attributes;
	//tokenize(tag_buf, " \t",attributes);

	bool antisense_splice = false;
	int sam_nm = 0; //the value of the NM tag (edit distance)
	//int mismatches[1024];//array with mismatch positions on the read (0-based from the left aligned end of the read)
	vector<bool> mismatches;
	mismatches.resize(strlen(seq_str), false);
	int num_mismatches=getSAMmismatches(buf, cigar, mismatches, sam_nm, antisense_splice);

	if (spliced_alignment)
	{
		bh = create_hit(name,
				ref_name,
				ref_name2,
				text_offset - 1, 
				cigar,
				sam_flag & 0x0010,
				antisense_splice,
				num_mismatches,
				0,
				end);
	}
	else
	{		
		//assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
		bh = create_hit(name,
				ref_name,
				ref_name2,
				text_offset - 1, // SAM files are 1-indexed 
				cigar,
				sam_flag & 0x0010,
				false,
				num_mismatches,
				0,
				end);
	}
	return true;
}

void cigar_add(vector<CigarOp>& cigar, CigarOp& op) {
 if (op.length<=0) return;
 if (cigar.size()>0 && cigar.back().opcode==op.opcode) {
    cigar.back().length+=op.length;
    }
 cigar.push_back(op);
}

bool spliceCigar(vector<CigarOp>& splcigar, vector<CigarOp>& cigar, vector<bool> mismatches,
		 int &left, int spl_start, int spl_len, CigarOpCode spl_code, int& spl_mismatches) {
  //merge the original 'cigar' with the new insert/gap operation
  //at position spl_start and place the result into splcigar;
  //TODO: ideally this should also get and rebuild the MD string (alignment mismatches)

  //return value: mismatches in the insert region for INS case,
  //or number of mismatches in the anchor region
  //return -1 if somehow the hit seems bad
  
  //these offsets are relative to the beginning of alignment on reference
  int spl_ofs=spl_start-left; //relative position of splice op
  if (spl_code == FUSION_FF || spl_code == FUSION_FR || spl_code == FUSION_RF || spl_code == FUSION_RR)
    spl_ofs = abs(spl_ofs);
  int spl_ofs_end=spl_ofs; //relative position of first ref base AFTER splice op
  CigarOp gapop(spl_code, spl_len); //for DEL, REF_SKIP, FUSIONS
  if (spl_code==INS)
    spl_ofs_end += spl_len;
  
  int ref_ofs=0; //working offset on reference
  int read_ofs=0; //working offset on the read, relative to the leftmost aligned base
  bool xfound=false;
  //if (left<=spl_start+spl_len) {
  if (spl_ofs_end>0) {
    int prev_opcode=0;
    int prev_oplen=0;
    for (size_t c = 0 ; c < cigar.size(); ++c)
      {
        int prev_read_ofs=read_ofs;
        int cur_op_ofs=ref_ofs;
        int cur_opcode=cigar[c].opcode;
        int cur_oplen=cigar[c].length;

        switch (cur_opcode) {
           case MATCH:
             ref_ofs+=cur_oplen;
             read_ofs+=cur_oplen;
             if (spl_code==REF_SKIP || spl_code==DEL ||
		 spl_code==FUSION_FF || spl_code==FUSION_FR ||
		 spl_code==FUSION_RF || spl_code==FUSION_RR) {
	       for (int o=cur_op_ofs;o<ref_ofs;o++) {
		 int rofs=prev_read_ofs+(o-cur_op_ofs);
		 if (abs(spl_ofs-o)<min_anchor_len && mismatches[rofs])
		   spl_mismatches++;
	       }
	     }
	     else if (spl_code==INS) {
	       for (int o=cur_op_ofs;o<ref_ofs;o++) {
		 int rofs=prev_read_ofs+(o-cur_op_ofs);
		 if (o>=spl_ofs && o<spl_ofs_end && mismatches[rofs])
		   spl_mismatches++;
	       }
	     }
             break;
           case DEL:
           case REF_SKIP:
           case PAD:
             ref_ofs+=cur_oplen;
             break;
           case SOFT_CLIP:
           case INS:
             read_ofs+=cur_oplen;
             break;
           //case HARD_CLIP:
           }

        if (cur_op_ofs>=spl_ofs_end || ref_ofs<=spl_ofs) {
	  if (cur_op_ofs==spl_ofs_end) {
	    if (spl_code!=INS) {
	      if (cur_opcode!=INS) {
		xfound=true;
		//we have to insert the gap here first
		cigar_add(splcigar, gapop);
		//also, check
	      }
	    }
	  }

	  CigarOp op(cigar[c]);
	  if (xfound)
	    {
	      if (spl_code == FUSION_FR || spl_code == FUSION_RR)
		{
		  if (op.opcode == MATCH)
		    op.opcode = mATCH;
		  else if (op.opcode == INS)
		    op.opcode = iNS;
		  else if (op.opcode == DEL)
		    op.opcode = dEL;
		  else if (op.opcode == REF_SKIP)
		    op.opcode = rEF_SKIP;
		}
	    }
	  else 
	    {
	      if (spl_code == FUSION_RF || spl_code == FUSION_RR)
		{
		  if (op.opcode == MATCH)
		    op.opcode = mATCH;
		  else if (op.opcode == INS)
		    op.opcode = iNS;
		  else if (op.opcode == DEL)
		    op.opcode = dEL;
		  else if (op.opcode == REF_SKIP)
		    op.opcode = rEF_SKIP;
		}
	    }
	  
	  cigar_add(splcigar, op);
	}
        else //if (ref_ofs>spl_ofs) {
           { //op intersection
           xfound=true;
           if (spl_code==INS) {
                 //we have to shorten cur_opcode
                 // find the overlap between current range
                 int ovl_start = (cur_op_ofs>spl_ofs) ? cur_op_ofs : spl_ofs;
                 int ovl_end = (ref_ofs>spl_ofs_end) ? spl_ofs_end : ref_ofs;
		 
		 CigarOp op(cigar[c]);
		 op.length=spl_ofs-cur_op_ofs;
		 if (spl_ofs>cur_op_ofs)
		   cigar_add(splcigar, op);
		 if (spl_ofs<0)
		   {
		     CigarOp temp = gapop;
		     temp.length += spl_ofs;
		     if (temp.length>0)
		       cigar_add(splcigar, temp);
		   }
		 else
		   cigar_add(splcigar, gapop);
		 op.length=ref_ofs-spl_ofs_end;
		 if (ref_ofs>spl_ofs_end)
		   cigar_add(splcigar,op);
                 }
              else {//DEL or REF_SKIP or FUSION_[FR][FR]
                 //spl_ofs == spl_ofs_end
                 //we have to split cur_opcode
                 //look for mismatches within min_anchor_len distance from splice point
                 CigarOp op(cigar[c]);
		 CigarOpCode opcode = op.opcode;
                 op.length=spl_ofs-cur_op_ofs;
		 if (spl_code == FUSION_RF || spl_code == FUSION_RR)
		   {
		     if (opcode == MATCH)
		       op.opcode = mATCH;
		     else if (opcode == INS)
		       op.opcode = iNS;
		     else if (opcode == DEL)
		       op.opcode = dEL;
		     else if (opcode == REF_SKIP)
		       op.opcode = rEF_SKIP;			 
		   }
                 cigar_add(splcigar, op);

		 cigar_add(splcigar, gapop);

		 op.opcode = opcode;
		 if (spl_code == FUSION_FR || spl_code == FUSION_RR)
		   {
		     if (opcode == MATCH)
		       op.opcode = mATCH;
		     else if (opcode == INS)
		       op.opcode = iNS;
		     else if (opcode == DEL)
		       op.opcode = dEL;
		     else if (opcode == REF_SKIP)
		       op.opcode = rEF_SKIP;			 
		   }
                 op.length=ref_ofs-spl_ofs;
                 cigar_add(splcigar,op);
                 }
            } //op intersection
        prev_opcode=cur_opcode;
        prev_oplen=cur_oplen;
      } //for each cigar opcode
    } //intersection possible

  //if (!xfound) {//no intersection found between splice event and alignment
   if (spl_ofs_end<=0) {
      //alignment starts after the splice event
      if (spl_code==INS) left-=spl_len;
                   else  left+=spl_len;

      splcigar = cigar;
      }
   //else {
     //alignment ends before the splice event
     //nothing to do
    //  }
   //return spl_mismatches;
  // }

   if (splcigar.size() < 3)
     return false;
   else if (splcigar.front().opcode != MATCH && splcigar.front().opcode != mATCH)
     return false;
   else if (splcigar.back().opcode != MATCH && splcigar.back().opcode != mATCH)
     return false;
   else
     return true;
     
}

bool SplicedSAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf,
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
	//int num_mismatches=0;
  //int mismatches[1024]; //list of 0-based mismatch positions in this read
                          //parsed from SAM's MD:Z: tag
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
  text_offset--; //make it 0-based (SAM is 1-based, Bowtie is 0-based)
	bool end = true;
	unsigned int seg_offset = 0;
	unsigned int seg_num = 0;
	unsigned int num_segs = 0;

	// Copy the tag out of the name field before we might wipe it out
	parseSegReadName(name, name_tags, strip_slash, end, seg_offset, seg_num, num_segs);

	vector<CigarOp> samcigar;
	bool spliced_alignment = false;
	int refspan=parseCigar(samcigar, cigar_str, spliced_alignment);

  if (refspan==0)
      return false;
	bool antisense_splice = false;
	int sam_nm = 0;
  vector<bool> mismatches;
  mismatches.resize(strlen(seq_str), false);

  int num_mismatches=getSAMmismatches(buf, samcigar, mismatches, sam_nm, antisense_splice);

	//##############################################

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
			//fprintf(stderr, "\t%s (token: %s)\n", text_name,
			//        toks[num_extra_toks + splice_field].c_str());
			return false;
		}

    string junction_strand = toks[num_extra_toks + strand_field];
    if(junction_strand != "rev" && junction_strand != "fwd"){
      fprintf(stderr, "Malformed insertion record\n");
      return false;
      }

    //
    // check for an insertion hit
    //
    if(toks[num_extra_toks + junction_type_field] == "ins")
      {
	//int8_t spliced_read_len = strlen(seq_str);
	//TODO FIXME: use the CIGAR instead of seq length!
	// The 0-based position of the left edge of the alignment. Note that this
	// value may need to be further corrected to account for the presence of
	// of the insertion.
	int left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
	
	// The 0-based position of the last genomic sequence before the insertion
	int left_splice_pos = atoi(splice_toks[0].c_str());
	
	string insertedSequence = splice_toks[1];
	// The 0-based position of the first genomic sequence after the insertion

	vector<CigarOp> splcigar;
	//this also updates left to the adjusted genomic coordinates
	int spl_num_mismatches=0;
	bool overlapped = spliceCigar(splcigar, samcigar, mismatches,
				      left, left_splice_pos+1, insertedSequence.length(),
				      INS, spl_num_mismatches);
	
	if (!overlapped)
	  return false;

	if (spl_num_mismatches<0) return false;
	num_mismatches-=spl_num_mismatches;
	/*
	  uint32_t right_splice_pos = left_splice_pos + 1;
	  
	  //uint32_t right = left + spliced_read_len - 1;
	  int right = left + refspan - 1;
	  
	  if(left > left_splice_pos){
	  //The genomic position of the left edge of the alignment needs to be corrected
	  //If the alignment does not extend into the insertion, simply subtract the length
	  //of the inserted sequence, otherwise, just set it equal to the right edge
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
	  // Now, right and left should be properly transformed into genomic coordinates
	  // We should be able to deduce how much the alignment matches the insertion
	  // simply based on the length of the read
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
	  
	  //char* pch = strtok( mismatches, ",");
	  //unsigned char num_mismatches = 0;
	  //text_offset holds the left end of the alignment,
	  //RELATIVE TO the start of the contig
	  
	  //The 0-based relative position of the left-most character
	  //before the insertion in the contig
	  int relative_splice_pos = left_splice_pos - atoi(toks[num_extra_toks + left_window_edge_field].c_str());
	  for (size_t i=0;i<mismatches.size();++i) {
	  int mismatch_pos = mismatches[i];
	  // for reversely mapped reads,
	  //find the correct mismatched position.
	  if (sam_flag & 0x0010){
	  mismatch_pos = spliced_read_len - mismatch_pos - 1;
	  }
	  
	  //Only count mismatches outside of the insertion region
	  // If there is a mismatch within the insertion,
	  // disallow this hit
	  if(mismatch_pos + text_offset <= relative_splice_pos ||
	  mismatch_pos + text_offset > relative_splice_pos + insertedSequence.length()){
	  spl_num_mismatches++;
	  }else{
	  return false;
	  }
	  }
	*/
	//vector<CigarOp> splcigar;
	//spliceCigar(splcigar, samcigar, left_match_length, insertion_match_length, right_match_length, INS);
	//splcigar.push_back(CigarOp(MATCH, left_match_length));
	//splcigar.push_back(CigarOp(INS, insertion_match_length));
	//splcigar.push_back(CigarOp(MATCH, right_match_length));
	
	bh = create_hit(name,
			contig,
			"",
			left,
			//splcigar,
			splcigar,
			sam_flag & 0x0010,
			junction_strand == "rev",
			num_mismatches,
			0,
			end);

	return true;
      } //"ins"
    else //"del" or intron
      {
	// The 0-based position of the left edge of the alignment.
	int left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
	
	// The 0-based position of the last genomic sequence before the deletion
	int left_splice_pos = atoi(splice_toks[0].c_str());
	
        int gap_len = atoi(splice_toks[1].c_str()) - left_splice_pos - 1;
        /*
	  if ((sam_flag & 0x0010) == 0) //######
	  {
	  if (left_splice_ofs + seg_offset < _anchor_length)
	  return false;
          }
	  else
          {
	  if (right_splice_ofs + seg_offset < _anchor_length)
	  return false;
          }
	*/
	//uint32_t right = atoi(splice_toks[1].c_str()) + right_splice_pos;
	//atoi(toks[right_window_edge_field].c_str());
	
        /*
        //offset of deletion point, relative to the beginning of the alignment
        int left_splice_ofs = left_splice_pos-left+1;
	
	int mismatches_in_anchor = 0;
	for (size_t i=0;i<mismatches.size();++i) {
	spl_num_mismatches++;
	int mismatch_pos = mismatches[i];
	if (((sam_flag & 0x0010) == 0 && abs(mismatch_pos - left_splice_ofs) < (int)min_anchor_len) ||
	((sam_flag & 0x0010) != 0 &&
	abs(((int)refspan - left_splice_ofs + 1) - mismatch_pos)) < (int)min_anchor_len)
	mismatches_in_anchor++;
	}
	*/
	vector<CigarOp> splcigar;

	CigarOpCode opcode=(toks[num_extra_toks + junction_type_field] == "del")? DEL : REF_SKIP;
	
	int spl_num_mismatches=0;
	bool overlapped = spliceCigar(splcigar, samcigar, mismatches, left,
				      left_splice_pos+1, gap_len, opcode, spl_num_mismatches);

	if (!overlapped)
	  return false;

	if (spl_num_mismatches<0) // || spl_num_mismatches>max_anchor_mismatches)
	  return false;
        /*
	  splcigar.push_back(CigarOp(MATCH, left_splice_pos));
	  if(toks[num_extra_toks + junction_type_field] == "del"){
	  splcigar.push_back(CigarOp(DEL, gap_len));
	  }else{
	  splcigar.push_back(CigarOp(REF_SKIP, gap_len));
	  }
	  splcigar.push_back(CigarOp(MATCH, right_splice_pos));
        */
	bh = create_hit(name,
			contig,
			"",
			left,
			splcigar,
			(sam_flag & 0x0010),
			junction_strand == "rev",
			num_mismatches,
			spl_num_mismatches,
			end);

	return true;
      }
	} //parse splice data
	else
	  {
	    fprintf(stderr, "Warning: found malformed splice record, skipping\n");
	    //fprintf(stderr, "%s\n", orig_bwt_buf);
	    //			continue;
	    return false;
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

  bool antisense_splice=false;
  unsigned char num_mismatches = 0;
  unsigned char num_splice_anchor_mismatches = 0;
  
  uint8_t* ptr = bam_aux_get(hit_buf, "XS");
  if (ptr) {
    char src_strand_char = bam_aux2A(ptr);
    if (src_strand_char == '-')
      antisense_splice = true;
  }
  
  ptr = bam_aux_get(hit_buf, "NM");
  if (ptr) {
    num_mismatches = bam_aux2i(ptr);
  }
  
  ptr = bam_aux_get(hit_buf, "NH");
  if (ptr) {
    num_hits = bam_aux2i(ptr);
  }

  int alignment_score = 0;
  bool has_alignment_score = false;
  ptr = bam_aux_get(hit_buf, "AS");
  if (ptr) {
    alignment_score = bam_aux2i(ptr);
    has_alignment_score = true;
  }
  
  string text_name = _sam_header->target_name[target_id];
  string text_name2 = "";
  
  bool fusion_alignment = false;
  string fusion_cigar_str;
  ptr = bam_aux_get(hit_buf, "XF");
  if (ptr) {
    fusion_alignment = true;
    char* xf = bam_aux2Z(ptr);

    // ignore the second part of a fusion alignment                                                                                                                                                                                                                            
    if (xf[0] == '2')
      return false;

    vector<string> fields;
    tokenize(xf, " ", fields);

    vector<string> contigs;
    tokenize(fields[1], "-", contigs);
    if (contigs.size() >= 2)
      {
	text_name = contigs[0];
	text_name2 = contigs[1];
      }

    text_offset = atoi(fields[2].c_str()) - 1;
    fusion_cigar_str = fields[3].c_str();

    if (seq)
      strcpy(seq, fields[4].c_str());
    if (qual)
      strcpy(qual, fields[5].c_str());
  }

  if (fusion_alignment) {
    const char* p_cig = fusion_cigar_str.c_str();
    while (*p_cig) {
      char* t;
      int length = (int)strtol(p_cig, &t, 10);
      if (length <= 0)
	{
	  //fprintf (stderr, "CIGAR op has zero length\n");
	  return false;
	}
      char op_char = *t;
      CigarOpCode opcode;
      if (op_char == 'M') opcode = MATCH;
      else if(op_char == 'm') opcode = mATCH;
      else if (op_char == 'I') opcode = INS;
      else if (op_char == 'i') opcode = iNS;
      else if (op_char == 'D') opcode = DEL;
      else if (op_char == 'd') opcode = dEL;
      else if (op_char == 'N' || op_char == 'n')
	{
	  if (length > max_report_intron_length)
	    return false;
	  
	  if (op_char == 'N')
	    opcode = REF_SKIP;
	  else
	    opcode = rEF_SKIP;
	  spliced_alignment = true;
	}
      else if (op_char == 'F')
	{
	  opcode = FUSION_FF;
	  length = length - 1;
	}
      else if (op_char == 'S') opcode = SOFT_CLIP;
      else if (op_char == 'H') opcode = HARD_CLIP;
      else if (op_char == 'P') opcode = PAD;
      else
	{
	  fprintf (stderr, "(%d-%d) invalid CIGAR operation\n", length, (int)op_char);
	  return false;
	}
      p_cig = t + 1;
      cigar.push_back(CigarOp(opcode, length));

      if (opcode == INS)
	num_mismatches -= length;
      else if (opcode == DEL)
	num_mismatches -= length;

      if (!has_alignment_score)
	{
	  if (opcode == INS)
	    alignment_score -= (bowtie2_read_gap_open * bowtie2_read_gap_cont * length);
	  else if(opcode == DEL)
	    alignment_score -= (bowtie2_ref_gap_open * bowtie2_ref_gap_cont * length);
	}
      
      /*
       * update fusion direction.
       */
      size_t cigar_size = cigar.size();
      if (cigar_size >= 3 && cigar[cigar_size - 2].opcode == FUSION_FF)
	{
	  CigarOpCode prev = cigar[cigar_size - 3].opcode;
	  CigarOpCode next = cigar[cigar_size - 1].opcode;
	  
	  bool increase1 = false, increase2 = false;
	  if (prev == MATCH || prev == DEL || prev == INS || prev == REF_SKIP)
	    increase1 = true;
	  if (next == MATCH || next == DEL || next == INS || next == REF_SKIP)
	    increase2 = true;
	  
	  if (increase1 && !increase2)
	    cigar[cigar_size - 2].opcode = FUSION_FR;
	  else if (!increase1 && increase2)
	    cigar[cigar_size - 2].opcode = FUSION_RF;
	  else if (!increase1 && !increase2)
	    cigar[cigar_size - 2].opcode = FUSION_RR;
	}
    }
  }
  else {
    for (int i = 0; i < hit_buf->core.n_cigar; ++i) 
      {
	int length = bam1_cigar(hit_buf)[i] >> BAM_CIGAR_SHIFT;
	if (length <= 0)
	  {
	    fprintf (stderr, "insert_id: %s - BAM error: CIGAR op has zero length\n", qname);
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
	
	/*
	 * By convention,the NM field of the SAM record
	 * counts an insertion or deletion. I dont' think
	 * we want the mismatch count in the BowtieHit
	 * record to reflect this. Therefore, subtract out
	 * the mismatches due to in/dels
	 */
	if (opcode == INS)
	  num_mismatches -= length;
	else if (opcode == DEL)
	  num_mismatches -= length;
	
	if (!has_alignment_score)
	  {
	    if (opcode == INS)
	      alignment_score -= (bowtie2_read_gap_open * bowtie2_read_gap_cont * length);
	    else if(opcode == DEL)
	      alignment_score -= (bowtie2_ref_gap_open * bowtie2_ref_gap_cont * length);
	  }
      }
  }

  if (!has_alignment_score)
    alignment_score -= (num_mismatches * (bowtie2_max_penalty + bowtie2_min_penalty) / 2);
  
  string mrnm;
  if (mate_target_id >= 0) {
    if (mate_target_id == target_id) {
      mrnm = _sam_header->target_name[mate_target_id];
    }
    else {
      return false;
    }
  }
  else {
    text_mate_pos = 0;
  }

  if (spliced_alignment) {
    bh = create_hit(qname,
		    text_name,
		    text_name2,
		    text_offset,  // BAM files are 0-indexed
		    cigar,
		    sam_flag & 0x0010,
		    antisense_splice,
		    num_mismatches,
		    num_splice_anchor_mismatches,
		    end);
    
  }
  else {
    bh = create_hit(qname,
		    text_name,
		    text_name2,
		    text_offset,  // BAM files are 0-indexed
		    cigar,
		    sam_flag & 0x0010,
		    false,
		    num_mismatches,
		    0,
		    end);
  }

  bh.alignment_score(alignment_score);
  
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

bool SplicedBAMHitFactory::get_hit_from_buf(const char* orig_bwt_buf, 
				 BowtieHit& bh, bool strip_slash,
				char* name_out, char* name_tags,
				char* seq, char* qual)
{
  if (_sam_header==NULL)
    err_die("Error: no SAM header when BAMHitFactory::get_hit_from_buf()!");
  
  const bam1_t* hit_buf = (const bam1_t*)orig_bwt_buf;
  uint32_t sam_flag = hit_buf->core.flag;
  
  int text_offset = hit_buf->core.pos;
  int text_mate_pos = hit_buf->core.mpos;
  int target_id = hit_buf->core.tid;
  int mate_target_id = hit_buf->core.mtid;
  
  vector<CigarOp> samcigar;
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
  
  bool end = true;
  unsigned int seg_offset = 0;
  unsigned int seg_num = 0;
  unsigned int num_segs = 0;
  // Copy the tag out of the name field before we might wipe it out
  char* name = bam1_qname(hit_buf);
  parseSegReadName(name, name_tags, strip_slash, end, seg_offset, seg_num, num_segs);
  
  if (target_id < 0) {
    //assert(cigar.size() == 1 && cigar[0].opcode == MATCH);
    bh = create_hit(name,
		    "*", //ref_name
		    0, //left coord
		    0, //read_len
		    false, //antisense_aln
		    0, //edit_dist
		    end);
    return true;
  }

  string text_name = _sam_header->target_name[target_id];
  for (int i = 0; i < hit_buf->core.n_cigar; ++i) 
    {
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
	samcigar.push_back(CigarOp(opcode, length));
    }
  
  string mrnm;
  if (mate_target_id >= 0) {
    if (mate_target_id == target_id) {
      mrnm = _sam_header->target_name[mate_target_id];
    }
    else {
      return false;
    }
  }
  else {
    text_mate_pos = 0;
  }
  
  bool antisense_splice = false;
  int sam_nm = 0;
  vector<bool> mismatches;
  mismatches.resize(strlen(seq), false);
  int num_mismatches=getBAMmismatches(hit_buf, samcigar, mismatches, sam_nm, antisense_splice);
  unsigned char num_splice_anchor_mismatches = 0;
  
  //##############################################
  
  // Add this alignment to the table of hits for this half of the
  // Bowtie map
  
  // Parse the text_name field to recover the splice coords
  vector<string> toks;
  tokenize_strict(text_name.c_str(), "|", toks);
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
	  //fprintf(stderr, "\t%s (token: %s)\n", text_name,
	  //        toks[num_extra_toks + splice_field].c_str());
	  return false;
	}
      
      const string& junction_type = toks[num_extra_toks + junction_type_field];
      const string junction_strand = toks[num_extra_toks + strand_field];
      
      //
      // check for an insertion hit
      //
      if(junction_type == "ins")
	{
	  //int8_t spliced_read_len = strlen(seq_str);
	  //TODO FIXME: use the CIGAR instead of seq length!
	  // The 0-based position of the left edge of the alignment. Note that this
	  // value may need to be further corrected to account for the presence of
	  // of the insertion.
	  int left = atoi(toks[num_extra_toks + left_window_edge_field].c_str()) + text_offset;
	  
	  // The 0-based position of the last genomic sequence before the insertion
	  int left_splice_pos = atoi(splice_toks[0].c_str());

	  if (left > left_splice_pos)
	    return false;
	  
	  string insertedSequence = splice_toks[1];
	  // The 0-based position of the first genomic sequence after the insertion
	  
	  vector<CigarOp> splcigar;
	  //this also updates left to the adjusted genomic coordinates
	  int spl_num_mismatches = 0;
	  bool overlapped = spliceCigar(splcigar, samcigar, mismatches,
					left, left_splice_pos+1, insertedSequence.length(), INS, spl_num_mismatches);
	  
	  if (!overlapped)
	    return false;
	  
	  if (spl_num_mismatches < 0)
	    return false;
	  
	  num_mismatches -= spl_num_mismatches;
	  bh = create_hit(name,
			  contig,
			  "",
			  left,
			  splcigar,
			  sam_flag & 0x0010,
			  junction_strand == "rev",
			  num_mismatches,
			  0,
			  end);
	  
	  return true;
	} //"ins"
      else //"del", "intron", or "fusion"
	{
	  char orientation = (sam_flag & 0x0010 ? '-' : '+');
	  if (!(junction_strand == "ff" || junction_strand == "fr" || junction_strand == "rf" || junction_strand == "rr" || junction_strand == "rev" || junction_strand == "fwd")||
	      !(orientation == '-' || orientation == '+'))
	    {
	      fprintf(stderr, "Warning: found malformed splice record, skipping\n");
	      return false;
	    }
	  
	  // The 0-based position of the left edge of the alignment.
	  int left = atoi(toks[num_extra_toks + left_window_edge_field].c_str());
	  if (junction_type != "fus" || (junction_strand != "rf" && junction_strand != "rr"))
	    left += text_offset;
	  else
	    left -= text_offset;
	  
	  vector<CigarOp> splcigar;
	  CigarOpCode opcode;
	  if(junction_type == "del")
	    opcode = DEL;
	  else if(junction_type == "fus")
	    {
	      if (junction_strand == "ff")
		opcode = FUSION_FF;
	      else if (junction_strand == "fr")
		opcode = FUSION_FR;
	      else if (junction_strand == "rf")
		opcode = FUSION_RF;
	      else
		opcode = FUSION_RR;
	    }
	  else
	    opcode = REF_SKIP;
	  
	  int left_splice_pos = atoi(splice_toks[0].c_str());
	  
	  // The 0-based position of the last genomic sequence before the deletion
	  int gap_len = 0;
	  if (junction_type == "fus")
	    gap_len = atoi(splice_toks[1].c_str());
	  else
	    gap_len = atoi(splice_toks[1].c_str()) - left_splice_pos - 1;
	  
	  if (opcode == FUSION_RF || opcode == FUSION_RR)
	    {
	      left_splice_pos -= 1;
	      if (left <= left_splice_pos)
		return false;
	    }
	  else
	    {
	      left_splice_pos += 1;
	      if (left >= left_splice_pos)
		return false;
	    }

	  int spl_num_mismatches = 0;
	  bool overlapped = spliceCigar(splcigar, samcigar, mismatches, left,
				      left_splice_pos, gap_len, opcode, spl_num_mismatches);
	  
	  if (!overlapped)
	    return false;
	  
	  if (spl_num_mismatches < 0) // || spl_num_mismatches>max_anchor_mismatches)
	    return false;

	  // daehwan - remove this
	  if (strcmp(name, "1921") == 0 && false)
	    {
	      cout << text_name << "\t" << left << endl
		   << print_cigar(samcigar) << endl
		   << print_cigar(splcigar) << endl
		   << "splice pos: " << left_splice_pos << endl;
	    }
	  
	  string contig2 = ""; 
	  if (junction_type == "fus")
	    {
	      vector<string> contigs;
	      tokenize(contig, "-", contigs);
	      if (contigs.size() != 2)
		return false;
	      
	      contig = contigs[0];
	      contig2 = contigs[1];
	      
	      if (junction_strand == "rf" || junction_strand == "rr")
		orientation = (orientation == '+' ? '-' : '+');
	    }
	  
	  bh = create_hit(name,
			  contig,
			  contig2,
			  left,
			  splcigar,
			  orientation == '-',
			  junction_strand == "rev",
			  num_mismatches,
			  spl_num_mismatches,
			  end);
	  
	  // daehwan - remove this
	  if (samcigar.size() > 1 && false)
	    {
	      cout << text_name << "\t" << left << endl
		   << print_cigar(samcigar) << endl
		   << print_cigar(splcigar) << endl;
	    }
	  
	  return true;
	}
    } //parse splice data
  else
    {
      fprintf(stderr, "Warning: found malformed splice record, skipping\n");
      //fprintf(stderr, "%s\n", orig_bwt_buf);
      //			continue;
      return false;
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

void print_bamhit(GBamWriter& wbam,
		  const char* read_name,
		  const BowtieHit& bh,
		  const char* ref_name,
		  const char* ref_name2,
		  const char* sequence,
		  const char* qualities,
		  bool from_bowtie,
		  const vector<string>* extra_fields)
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
  else {
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
  CigarOpCode fusion_dir = FUSION_NOTHING;
  for (size_t c = 0; c < bh_cigar.size(); ++c)
    {
      const CigarOp& op = bh_cigar[c];

      char ibuf[64];
      sprintf(ibuf, "%d", op.length);
      switch(op.opcode)
	{
	case MATCH:
	case mATCH:
	  strcat(cigar, ibuf);
	  if (bh_cigar[c].opcode == MATCH)
	    strcat(cigar, "M");
	  else
	    strcat(cigar, "m");
	  break;
	case INS:
	case iNS:
	  strcat(cigar, ibuf);
	  if (bh_cigar[c].opcode == INS)
	    strcat(cigar, "I");
	  else
	    strcat(cigar, "i");
	  indel_distance += bh_cigar[c].length;
	  break;
	case DEL:
	case dEL:
	  strcat(cigar, ibuf);
	  if (bh_cigar[c].opcode == DEL)
	    strcat(cigar, "D");
	  else
	    strcat(cigar, "d");
	  indel_distance += bh_cigar[c].length;
	  break;
	case REF_SKIP:
	case rEF_SKIP:
	  strcat(cigar, ibuf);
	  if (bh_cigar[c].opcode == REF_SKIP)
	    strcat(cigar, "N");
	  else
	    strcat(cigar, "n");
	  break;
	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  fusion_dir = op.opcode;
	  sprintf(ibuf, "%d", bh_cigar[c].length + 1);
	  strcat(cigar, ibuf);
	  strcat(cigar, "F");
	  break;
	default:
	  break;
	}
    }

  char cigar1[256] = {0}, cigar2[256] = {0};
  string left_seq, right_seq, left_qual, right_qual;
  int left1 = -1, left2 = -1;
  extract_partial_hits(bh, seq, quals,
		       cigar1, cigar2, left_seq, right_seq,
		       left_qual, right_qual, left1, left2);
  
  bool containsSplice = false;
  for (vector<CigarOp>::const_iterator itr = bh.cigar().begin(); itr != bh.cigar().end(); ++itr)
    {
      if (itr->opcode == REF_SKIP || itr->opcode == rEF_SKIP)
	{
	  containsSplice = true;
	  break;
	}
    }

  vector<string> auxdata;
  if (extra_fields)
    auxdata.insert(auxdata.end(), extra_fields->begin(), extra_fields->end());

  if (!sam_readgroup_id.empty())
    {
      string nm("RG:Z:");
      nm += sam_readgroup_id;
      auxdata.push_back(nm);
    }
  
  string nm("NM:i:");
  str_appendInt(nm, bh.edit_dist() + indel_distance);
  auxdata.push_back(nm);
  
  if (containsSplice) {
    // do not add more than once
    bool XS_found = false;
    for (size_t i = 0; i < auxdata.size(); ++i)
      {
	if (auxdata[i].substr(0, 2) == "XS")
	  {
	    XS_found = true;
	    break;
	  }
      }

    if (!XS_found)
      {
	nm="XS:A:";
	nm+=(char)(bh.antisense_splice() ? '-' : '+');
	auxdata.push_back(nm);
      }
  }
  
  if (fusion_dir != FUSION_NOTHING)
    {
      char XF[2048] = {0};
      sprintf(XF,
	      "XF:Z:1 %s-%s %u %s %s %s",
	      ref_name,
	      ref_name2,
	      sam_pos,
	      cigar,
	      seq.c_str(),
	      quals.c_str());
      auxdata.push_back(XF);

      GBamRecord *brec = wbam.new_record(read_name, sam_flag, ref_name, left1 + 1, map_quality,
					 cigar1, mate_ref_name.c_str(), mate_pos,
					 insert_size, left_seq.c_str(), left_qual.c_str(), &auxdata);
      
      wbam.write(brec);
      delete brec;

      sprintf(XF,
	      "XF:Z:2 %s-%s %u %s %s %s",
	      ref_name,
	      ref_name2,
	      sam_pos,
	      cigar,
	      seq.c_str(),
	      quals.c_str());
      auxdata.back() = XF;

      brec = wbam.new_record(read_name, sam_flag, ref_name2, left2 + 1, map_quality,
			     cigar2, mate_ref_name.c_str(), mate_pos,
			     insert_size, right_seq.c_str(), right_qual.c_str(), &auxdata);
      
      wbam.write(brec);
      delete brec;
    }
  else
    {
      GBamRecord *brec = wbam.new_record(read_name, sam_flag, ref_name, sam_pos, map_quality,
					 cigar, mate_ref_name.c_str(), mate_pos,
					 insert_size, seq.c_str(), quals.c_str(), &auxdata);
      
      wbam.write(brec);
      delete brec;
    }
}

/**
 * Print a vector of cigar operations to a file.
 * @param bh_cigar A vector of CigarOps
 * @return a string representation of the cigar string
*/
std::string print_cigar(const vector<CigarOp>& bh_cigar){
	char cigar[256];
	cigar[0] = 0;	
	for (size_t c = 0; c < bh_cigar.size(); ++c)
	{
		char ibuf[64];
		sprintf(ibuf, "%d", bh_cigar[c].length);
		strcat(cigar, ibuf);
		switch(bh_cigar[c].opcode)
		{
		case MATCH:
		  strcat(cigar, "M");
		  break;
		case mATCH:
		  strcat(cigar, "m");
		  break;
		case INS:
		  strcat(cigar, "I");
		  break;
		case iNS:
		  strcat(cigar, "i");
		  break;
		case DEL:
		  strcat(cigar, "D");
		  break;
		case dEL:
		  strcat(cigar, "d");
		  break;
		case REF_SKIP:
		  strcat(cigar, "N");
		  break;
		case rEF_SKIP:
		  strcat(cigar, "n");
		  break;
		case FUSION_FF:
		case FUSION_FR:
		case FUSION_RF:
		case FUSION_RR:
		  strcat(cigar, "F");
		  break;
		default:
		  break;
		}
	}
	string result(cigar);
	return result;
}

void extract_partial_hits(const BowtieHit& bh, const string& seq, const string& qual,
			  char* cigar1, char* cigar2, string& seq1, string& seq2,
			  string& qual1, string& qual2, int& left1, int& left2)
{
  const int left = bh.left();
  int right = left;
  int fusion_left = -1, fusion_right = -1;
  
  const vector<CigarOp>& bh_cigar = bh.cigar();
  
  CigarOpCode fusion_dir = FUSION_NOTHING;
  size_t fusion_idx = 0;
  size_t left_part_len = 0;
  for (size_t c = 0; c < bh_cigar.size(); ++c)
    {
      const CigarOp& op = bh_cigar[c];
      switch(op.opcode)
	{
	case MATCH:
	case REF_SKIP:
	case DEL:
	  right += op.length;
	  break;
	case mATCH:
	case rEF_SKIP:
	case dEL:
	  right -= op.length;
	  break;
	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  {
	    fusion_dir = op.opcode;
	    fusion_idx = c;
	    if (op.opcode == FUSION_FF || op.opcode == FUSION_FR)
	      fusion_left = right - 1;
	    else
	      fusion_left = right + 1;
	    fusion_right = right = op.length;
	  }
	  break;
	default:
	  break;
	}

      if (fusion_dir == FUSION_NOTHING)
	{
	  if (op.opcode == MATCH || op.opcode == mATCH || op.opcode == INS || op.opcode == iNS)
	    {
	      left_part_len += op.length;
	    }
	}
    }

  if (fusion_dir == FUSION_FF || fusion_dir == FUSION_FR)
    {
      for (size_t c = 0; c < fusion_idx; ++c)
	{
	  const CigarOp& op = bh_cigar[c];
	  char ibuf[64];
	  sprintf(ibuf, "%d", op.length);
	  strcat(cigar1, ibuf);

	  switch (op.opcode)
	    {
	    case MATCH:
	      strcat(cigar1, "M");
	      break;
	    case INS:
	      strcat(cigar1, "I");
	      break;
	    case DEL:
	      strcat(cigar1, "D");
	      break;
	    case REF_SKIP:
	      strcat(cigar1, "N");
	      break;
	    default:
	      assert (0);
	      break;
	    }
	}
    }
  else if (fusion_dir == FUSION_RF || fusion_dir == FUSION_RR)
    {
      assert (fusion_idx > 0);
      for (int c = fusion_idx - 1; c >=0; --c)
	{
	  const CigarOp& op = bh_cigar[c];
	  char ibuf[64];
	  sprintf(ibuf, "%d", op.length);
	  strcat(cigar1, ibuf);

	  switch (op.opcode)
	    {
	    case mATCH:
	      strcat(cigar1, "M");
	      break;
	    case iNS:
	      strcat(cigar1, "I");
	      break;
	    case dEL:
	      strcat(cigar1, "D");
	      break;
	    case rEF_SKIP:
	      strcat(cigar1, "N");
	      break;
	    default:
	      assert (0);
	      break;
	    }
	}
    }

  if (fusion_dir == FUSION_FF || fusion_dir == FUSION_RF)
    {
      for (size_t c = fusion_idx + 1; c < bh_cigar.size(); ++c)
	{
	  const CigarOp& op = bh_cigar[c];
	  char ibuf[64];
	  sprintf(ibuf, "%d", op.length);
	  strcat(cigar2, ibuf);

	  switch (op.opcode)
	    {
	    case MATCH:
	      strcat(cigar2, "M");
	      break;
	    case INS:
	      strcat(cigar2, "I");
	      break;
	    case DEL:
	      strcat(cigar2, "D");
	      break;
	    case REF_SKIP:
	      strcat(cigar2, "N");
	      break;
	    default:
	      assert (0);
	      break;
	    }
	}
    }
  else if (fusion_dir == FUSION_FR || fusion_dir == FUSION_RR)
    {
      assert (bh_cigar.size() > 0);
      for (size_t c = bh_cigar.size() - 1; c > fusion_idx; --c)
	{
	  const CigarOp& op = bh_cigar[c];
	  char ibuf[64];
	  sprintf(ibuf, "%d", op.length);
	  strcat(cigar2, ibuf);

	  switch (op.opcode)
	    {
	    case mATCH:
	      strcat(cigar2, "M");
	      break;
	    case iNS:
	      strcat(cigar2, "I");
	      break;
	    case dEL:
	      strcat(cigar2, "D");
	      break;
	    case rEF_SKIP:
	      strcat(cigar2, "N");
	      break;
	    default:
	      assert (0);
	      break;
	    }
	}
    }
  
  if (fusion_dir != FUSION_NOTHING)
    {
      seq1 = seq.substr(0, left_part_len);
      qual1 = qual.substr(0, left_part_len);

      if (fusion_dir == FUSION_RF || fusion_dir == FUSION_RR)
	{
	  reverse_complement(seq1);
	  reverse(qual1.begin(), qual1.end());
	}
      
      seq2 = seq.substr(left_part_len);
      qual2 = qual.substr(left_part_len);

      if (fusion_dir == FUSION_FR || fusion_dir == FUSION_RR)
	{
	  reverse_complement(seq2);
	  reverse(qual2.begin(), qual2.end());
	}

      left1 = ((fusion_dir == FUSION_FF || fusion_dir == FUSION_FR) ? left : fusion_left);
      left2 = ((fusion_dir == FUSION_FF || fusion_dir == FUSION_RF) ? fusion_right : right + 1);
    }
}


bool BowtieHit::check_editdist_consistency(const RefSequenceTable& rt, bool bDebug)
{
  RefSequenceTable::Sequence* ref_str1 = rt.get_seq(_ref_id);
  RefSequenceTable::Sequence* ref_str2 = rt.get_seq(_ref_id2);
  
  if (!ref_str1 || !ref_str2)
    return false;

  if (bDebug)
    {
      cout << "check_editdist_consistency" << endl
	   << "insert id: " << _insert_id << endl;
    }
  
  RefSequenceTable::Sequence* ref_str = ref_str1;

  size_t pos_seq = 0;
  size_t pos_ref = _left;
  size_t mismatch = 0;
  size_t N_mismatch = 0;
  bool bSawFusion = false;
  for (size_t i = 0; i < _cigar.size(); ++i)
    {
      CigarOp cigar = _cigar[i];
      switch(cigar.opcode)
	{
	case MATCH:
	  {
	    seqan::Dna5String ref_seq = seqan::infix(*ref_str, pos_ref, pos_ref + cigar.length);
	    for (size_t j = 0; j < cigar.length; ++j)
	      {
		seqan::Dna5 ref_nt = _seq[pos_seq];
		if (ref_nt != ref_seq[j])
		  ++mismatch;

		if (ref_nt == ref_seq[j] && ref_nt == 'N')
		  ++N_mismatch;

		if (bDebug)
		  cout << pos_seq << "\t" << ref_nt << " vs. " << ref_seq[j] << "\t" << mismatch << endl;
	    
		++pos_seq;
	      }

	    pos_ref += cigar.length;
	  }
	  break;
	case mATCH:
	  {
	    seqan::Dna5String ref_seq = seqan::infix(*ref_str, pos_ref - cigar.length + 1, pos_ref + 1);
	    seqan::reverseComplement(ref_seq);

	    for (size_t j = 0; j < cigar.length; ++j)
	      {
		seqan::Dna5 ref_nt = _seq[pos_seq];
		if (ref_nt != ref_seq[j])
		  ++mismatch;

		if (ref_nt == ref_seq[j] && ref_nt == 'N')
		  ++N_mismatch;

		if (bDebug)
		  cout << pos_seq << "\t" << ref_nt << " vs. " << ref_seq[j] << "\t" << mismatch << endl;
	    
		++pos_seq;
	      }

	    pos_ref -= cigar.length;
	  }
	  break;
	case INS:
	case iNS:
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

	case dEL:
	case rEF_SKIP:
	  {
	    pos_ref -= cigar.length;
	  }
	  break;

	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  {
	    // We don't allow a read spans more than two chromosomes.
	    if (bSawFusion)
	      return false;

	    ref_str = ref_str2;  
	    pos_ref = cigar.length;

	    bSawFusion = true;
	  }
	  break;

	default:
	  break;
	}
    }

  if (bDebug)
    cout << "mismatch (real) vs. (calculated):" << mismatch << " vs. " << (int)_edit_dist << endl;

  return mismatch == _edit_dist || mismatch + N_mismatch == _edit_dist;
}

void bowtie_sam_extra(const BowtieHit& bh, const RefSequenceTable& rt, vector<string>& fields)
{
  RefSequenceTable::Sequence* ref_str1 = rt.get_seq(bh.ref_id());
  RefSequenceTable::Sequence* ref_str2 = rt.get_seq(bh.ref_id2());
  
  if (!ref_str1 || !ref_str2)
    return;
  
  RefSequenceTable::Sequence* ref_str = ref_str1;

  size_t pos_seq = 0;
  size_t pos_mismatch = 0;
  size_t pos_ref = bh.left();
  size_t mismatch = 0;
  size_t N_mismatch = 0;
  size_t num_gap_opens = 0;
  size_t num_gap_conts = 0;
  bool bSawFusion = false;

  int AS_score = 0;
  
  const vector<CigarOp>& cigars = bh.cigar();
  const string& seq = bh.seq();
  const string& qual = bh.qual();

  string AS = "AS:i:";
  string MD = "MD:Z:";
  
  for (size_t i = 0; i < cigars.size(); ++i)
    {
      CigarOp cigar = cigars[i];
      switch(cigar.opcode)
	{
	case MATCH:
	case mATCH:
	  {
	    seqan::Dna5String ref_seq;
	    if (cigar.opcode == MATCH)
	      {
		ref_seq = seqan::infix(*ref_str, pos_ref, pos_ref + cigar.length);
		pos_ref += cigar.length;
	      }
	    else
	      {
		ref_seq = seqan::infix(*ref_str, pos_ref - cigar.length + 1, pos_ref + 1);
		seqan::reverseComplement(ref_seq);
		pos_ref -= cigar.length;
	      }
	    
	    for (size_t j = 0; j < cigar.length; ++j)
	      {
		seqan::Dna5 ref_nt = ref_seq[j];
		if (seq[pos_seq] != ref_nt)
		  {
		    ++mismatch;

		    if (pos_seq < qual.length())
		      {
			if (seq[pos_seq] == 'N' || ref_nt == 'N')
			  {
			    AS_score -= (int)bowtie2_penalty_for_N;
			  }
			else
			  {
			    float penalty = bowtie2_min_penalty + (bowtie2_max_penalty - bowtie2_min_penalty) * min((int)(qual[pos_seq] - '!'), 40) / 40.0;
			    AS_score -= (int)penalty;
			  }
		      }

		    str_appendInt(MD, (int)pos_mismatch);
		    MD.push_back((char)ref_nt);
		    pos_mismatch = 0;
		  }
		else
		  {
		    if (ref_nt == 'N')
		      {
			++N_mismatch;
			AS_score -= (int)bowtie2_penalty_for_N;
		      }

		    ++pos_mismatch;
		  }

		++pos_seq;
	      }
	  }
	  break;

	case INS:
	case iNS:
	  {
	    pos_seq += cigar.length;

	    AS_score -= bowtie2_read_gap_open;
	    AS_score -= (int)(bowtie2_read_gap_cont * cigar.length);

	    num_gap_opens += 1;
	    num_gap_conts += cigar.length;
	  }
	  break;
	  
	case DEL:
	case dEL:
	  {
	    AS_score -= bowtie2_ref_gap_open;
	    AS_score -= (int)(bowtie2_ref_gap_cont * cigar.length);

	    num_gap_opens += 1;
	    num_gap_conts += cigar.length;
	      
	    seqan::Dna5String ref_seq;
	    if (cigar.opcode == DEL)
	      {
		ref_seq = seqan::infix(*ref_str, pos_ref, pos_ref + cigar.length);
		pos_ref += cigar.length;
	      }
	    else
	      {
		ref_seq = seqan::infix(*ref_str, pos_ref - cigar.length + 1, pos_ref + 1);
		seqan::reverseComplement(ref_seq);
		pos_ref -= cigar.length;
	      }

	    str_appendInt(MD, (int)pos_mismatch);
	    MD.push_back('^');
	    for (size_t k = 0; k < length(ref_seq); ++k)
	      MD.push_back((char)ref_seq[k]);
	    
	    pos_mismatch = 0;
	  }
	  break;

	case REF_SKIP:
	case rEF_SKIP:
	  {
	    if (cigar.opcode == REF_SKIP)
	      pos_ref += cigar.length;
	    else
	      pos_ref -= cigar.length;
	  }
	  break;

	case FUSION_FF:
	case FUSION_FR:
	case FUSION_RF:
	case FUSION_RR:
	  {
	    // We don't allow a read spans more than two chromosomes.
	    if (bSawFusion)
	      return;

	    ref_str = ref_str2;  
	    pos_ref = cigar.length;

	    bSawFusion = true;
	  }
	  break;

	default:
	  break;
	}
    }

  str_appendInt(AS, AS_score);
  fields.push_back(AS);

  string XM = "XM:i:";
  str_appendInt(XM, (int)mismatch);
  fields.push_back(XM);

  string XO = "XO:i:";
  str_appendInt(XO, (int)num_gap_opens);
  fields.push_back(XO);

  string XG = "XG:i:";
  str_appendInt(XG, (int)num_gap_conts);
  fields.push_back(XG);

  str_appendInt(MD, (int)pos_mismatch);
  fields.push_back(MD);
}
