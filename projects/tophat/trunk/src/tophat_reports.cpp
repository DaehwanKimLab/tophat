/*
 *  tophat_reports.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/20/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
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
#include <getopt.h>

#include "common.h"
#include "bwt_map.h"
#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "align_status.h"
#include "fragments.h"
#include "wiggles.h"
#include "tokenize.h"
#include "reads.h"

#include "inserts.h"

using namespace std;
using namespace seqan;
using std::set;

void fragment_best_alignments(const HitsForRead& hits_for_read,
							  FragmentAlignmentGrade& best_grade,
                              HitsForRead& best_hits)
{
	const vector<BowtieHit>& hits = hits_for_read.hits;
	
	for (size_t i = 0; i < hits.size(); ++i)
	{
		FragmentAlignmentGrade g(hits[i]);
		// Is the new status better than the current best one?
		if (best_grade < g)
		{
			best_hits.hits.clear();
			best_grade = g;
			best_hits.hits.push_back(hits[i]);
		}
		else if (!(g < best_grade)) // is it just as good?
		{
			best_grade.num_alignments++;
			best_hits.hits.push_back(hits[i]);
		}
	}
}

void insert_best_alignments(const HitsForRead& left_hits,
			    const HitsForRead& right_hits,
			    InsertAlignmentGrade& best_grade,
			    HitsForRead& left_best_hits,
			    HitsForRead& right_best_hits)
{
  // max mate inner distance (genomic)
  int min_mate_inner_dist = inner_dist_mean - inner_dist_std_dev;
  if (max_mate_inner_dist == -1)
    {
      max_mate_inner_dist = inner_dist_mean + inner_dist_std_dev;
    }
  
  const vector<BowtieHit>& left = left_hits.hits;
  const vector<BowtieHit>& right = right_hits.hits;
  
  for (size_t i = 0; i < left.size(); ++i)
    {
      const BowtieHit& lh = left[i];
      for (size_t j = 0; j < right.size(); ++j)
	{
	  const BowtieHit& rh = right[j];

	  if (lh.ref_id() != rh.ref_id())
	    continue;
	  
	  InsertAlignmentGrade g(lh, rh, min_mate_inner_dist, max_mate_inner_dist);

	  // Is the new status better than the current best one?
	  if (best_grade < g)
	    {
	      left_best_hits.hits.clear();
	      right_best_hits.hits.clear();
	      
	      best_grade = g;
	      left_best_hits.hits.push_back(lh);
	      right_best_hits.hits.push_back(rh);
	    }
	  else if (!(g < best_grade))
	    {
	      left_best_hits.hits.push_back(lh);
	      right_best_hits.hits.push_back(rh);
	      best_grade.num_alignments++;
	    }
	}
    }
}

enum FragmentType {FRAG_UNPAIRED, FRAG_LEFT, FRAG_RIGHT};

bool rewrite_sam_hit(const RefSequenceTable& rt,
                     const BowtieHit& bh,
                     const char* bwt_buf,
                     char* rebuf,
                     char* read_alt_name,
                     const FragmentAlignmentGrade& grade,
                     FragmentType insert_side,
                     int num_hits,
                     const BowtieHit* next_hit)
{
    // Rewrite this hit, filling in the alt name, mate mapping
    // and setting the pair flag
    vector<string> sam_toks;
    tokenize(bwt_buf, "\t", sam_toks);
    char* slash_pos = strrchr(read_alt_name,'/');
    if (slash_pos)
    {
        *slash_pos = 0;
    }

    *rebuf = 0;

    for (size_t t = 0; t < sam_toks.size(); ++t)
    {
        switch(t)
        {
            case 0: //QNAME
            {
                sam_toks[t] = read_alt_name;
                break;
            }

            case 1:
            {
                // SAM FLAG
                if (insert_side != FRAG_UNPAIRED)
                {
                    int flag = atoi(sam_toks[1].c_str());;
                    // mark this guys as a singleton mate
                    flag |= 0x0001;
                    if (insert_side == FRAG_LEFT)
                      flag |= 0x0040;
                    else if (insert_side == FRAG_RIGHT)
                      flag |= 0x0080;
                    flag |= 0x0008;

                    char flag_buf[64];
                    sprintf(flag_buf, "%d", flag);
                    sam_toks[t] = flag_buf;
                }
                break;
            }

            case 4: //MAPQ
            {
                int mapQ;
                if (grade.num_alignments > 1)
                {
                    double err_prob = 1 - (1.0 / grade.num_alignments);
                    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
                }
                else
                {
                    mapQ = 255;
                }
                char mapq_buf[64];
                sprintf(mapq_buf, "%d", mapQ);
                sam_toks[t] = mapq_buf;
                break;
            }

            default:
                break;
        }
        strcat (rebuf, sam_toks[t].c_str());
        if (t != sam_toks.size() - 1)
            strcat(rebuf, "\t");
    }

    char nh_buf[2048];

    sprintf(nh_buf,
            "\tNH:i:%d",
            num_hits);

    strcat(rebuf, nh_buf);

    if (next_hit)
    {
        const char* nh_ref_name = rt.get_name(next_hit->ref_id());
        assert (nh_ref_name != NULL);
        const char* curr_ref_name = rt.get_name(bh.ref_id());
        assert (curr_ref_name != NULL);

        char mate_buf[2048];
        bool same_contig = !strcmp(curr_ref_name, nh_ref_name);

        sprintf(mate_buf,
                "\tCC:Z:%s\tCP:i:%d",
                same_contig ? "=" : nh_ref_name,
                next_hit->left() + 1);
        strcat(rebuf, mate_buf);
    }

    // FIXME: this code is still a bit brittle, because it contains no
    // consistency check that the mates are on opposite strands (a current protocol
    // requirement, and that the strand indicated by the alignment is consistent
    // with the orientation of the splices (though that should be handled upstream).
    if (bh.contiguous())
    {
        if (library_type == FR_FIRSTSTRAND)
        {
            if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED)
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:+");
                else
                    strcat(rebuf, "\tXS:A:-");
            }
            else
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:-");
                else
                    strcat(rebuf, "\tXS:A:+");
            }
        }

        else if (library_type == FR_SECONDSTRAND)
        {
            if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED)
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:-");
                else
                    strcat(rebuf, "\tXS:A:+");
            }
            else
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:+");
                else
                    strcat(rebuf, "\tXS:A:-");
            }
        }
    }

    strcat(rebuf, "\n");

    return true;
}

bool rewrite_sam_hit(const RefSequenceTable& rt,
                     const BowtieHit& bh,
                     const char* bwt_buf,
                     char* rebuf,
                     char* read_alt_name,
                     const InsertAlignmentGrade& grade,
                     FragmentType insert_side,
                     const BowtieHit* partner,
                     int num_hits,
                     const BowtieHit* next_hit)
{
    // Rewrite this hit, filling in the alt name, mate mapping
    // and setting the pair flag
    vector<string> sam_toks;
    tokenize(bwt_buf, "\t", sam_toks);
    char* slash_pos = strrchr(read_alt_name,'/');
    if (slash_pos)
    {
        *slash_pos = 0;
    }

    *rebuf = 0;

    for (size_t t = 0; t < sam_toks.size(); ++t)
    {
        switch(t)
        {
            case 0: //QNAME
            {
                sam_toks[t] = read_alt_name;
                break;
            }
            case 1: //SAM FLAG
            {
                // 0x0010 (strand of query) is assumed to be set correctly
                // to begin with

                int flag = atoi(sam_toks[1].c_str());
                flag |= 0x0001;
                if (insert_side == FRAG_LEFT)
                    flag |= 0x0040;
                else if (insert_side == FRAG_RIGHT)
                    flag |= 0x0080;

                if (grade.happy() && partner)
                    flag |= 0x0002;

                if (partner)
                {
                    if (partner->antisense_align())
                        flag |= 0x0020;
                }
                else
                {
                    flag |= 0x0008;
                }

                char flag_buf[64];
                sprintf(flag_buf, "%d", flag);
                sam_toks[t] = flag_buf;
                break;
            }

            case 4: //MAPQ
            {
                int mapQ;
                if (grade.num_alignments > 1)
                {
                    double err_prob = 1 - (1.0 / grade.num_alignments);
                    mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
                }
                else
                {
                    mapQ = 255;
                }
                char mapq_buf[64];
                sprintf(mapq_buf, "%d", mapQ);
                sam_toks[t] = mapq_buf;
                break;
            }
            case 6: //MRNM
            {
                if (partner)
                {
                    //FIXME: this won't be true forever.  Someday, we will report
                    //alignments of pairs not on the same contig.
                    sam_toks[t] = "=";
                }
                else
                {
                    sam_toks[t] = "*";
                }
                break;
            }
            case 7: //MPOS
            {
                if (partner)
                {
                    char pos_buf[64];
                    int partner_pos = partner->left() + 1;  // SAM is 1-indexed

                    sprintf(pos_buf, "%d", partner_pos);
                    sam_toks[t] = pos_buf;
                }
                else
                {
                    sam_toks[t] = "0";
                }
		break;
            }
	case 8: //TLEN
	  {
	    if (partner && bh.ref_id() == partner->ref_id())
	      {
		char tlen_buf[64];
		int tlen = 0;

		if (bh.left() < partner->left())
		  tlen = partner->right() - bh.left();
		else
		  tlen = partner->left() - bh.right();
		
		sprintf(tlen_buf, "%d", tlen);
		sam_toks[t] = tlen_buf;
	      }
	    else
	      {
		sam_toks[t] = "0";
	      }
	  }
	  
	default:
	  break;
        }
        strcat (rebuf, sam_toks[t].c_str());
        if (t != sam_toks.size() - 1)
            strcat(rebuf, "\t");
    }

    char nh_buf[2048];
    sprintf(nh_buf, "\tNH:i:%d", num_hits);
    strcat(rebuf, nh_buf);

    if (next_hit)
    {
        const char* nh_ref_name = rt.get_name(next_hit->ref_id());
        assert (nh_ref_name != NULL);
        const char* curr_ref_name = rt.get_name(bh.ref_id());
        assert (curr_ref_name != NULL);

        char mate_buf[2048];
        bool same_contig = !strcmp(curr_ref_name, nh_ref_name);

        assert (num_hits > 1);

        sprintf(mate_buf,
                "\tCC:Z:%s\tCP:i:%d",
                same_contig ? "=" : nh_ref_name,
                next_hit->left() + 1);

        strcat(rebuf, mate_buf);
    }

    // FIXME: this code is still a bit brittle, because it contains no
    // consistency check that the mates are on opposite strands (a current protocol
    // requirement, and that the strand indicated by the alignment is consistent
    // with the orientation of the splices (though that should be handled upstream).
    if (bh.contiguous() && grade.opposite_strands)
    {
        if (library_type == FR_FIRSTSTRAND)
        {
            if (insert_side == FRAG_LEFT)
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:+");
                else
                    strcat(rebuf, "\tXS:A:-");
            }
            else
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:-");
                else
                    strcat(rebuf, "\tXS:A:+");
            }
        }

        else if (library_type == FR_SECONDSTRAND)
        {
            if (insert_side == FRAG_LEFT)
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:-");
                else
                    strcat(rebuf, "\tXS:A:+");
            }
            else
            {
                if (bh.antisense_align())
                    strcat(rebuf, "\tXS:A:+");
                else
                    strcat(rebuf, "\tXS:A:-");
            }
        }
    }
    strcat(rebuf, "\n");

    return true;
}


void add_auxData(vector<string>& auxdata, vector<string>& sam_toks,
    const RefSequenceTable& rt,const BowtieHit& bh, FragmentType insert_side,
     int num_hits, const BowtieHit* next_hit) {

  if (sam_toks.size()>11) {
    for (size_t i=11;i<sam_toks.size();++i) {
       if (sam_toks[i].find("NH:i:")==string::npos)
            auxdata.push_back(sam_toks[i]);
       }
    }
  string aux("NH:i:");
  str_appendInt(aux, num_hits);
  auxdata.push_back(aux);
  if (next_hit) {
      const char* nh_ref_name = "=";
      nh_ref_name = rt.get_name(next_hit->ref_id());
      assert (nh_ref_name != NULL);
      bool same_contig=(next_hit->ref_id()==bh.ref_id());
      aux="CC:Z:";
      aux+= (same_contig ? "=" : nh_ref_name);
      auxdata.push_back(aux);
      aux="CP:i:";
      int nh_gpos=next_hit->left() + 1;
      str_appendInt(aux, nh_gpos);
      auxdata.push_back(aux);
      } //has next_hit
  // FIXME: this code is still a bit brittle, because it contains no
  // consistency check that the mates are on opposite strands (a current protocol
  // requirement, and that the strand indicated by the alignment is consistent
  // with the orientation of the splices (though that should be handled upstream).
  const string xs_f("XS:A:+");
  const string xs_r("XS:A:-");
  if (bh.contiguous())  {
      if (library_type == FR_FIRSTSTRAND) {
          if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED) {
              if (bh.antisense_align())
                  auxdata.push_back(xs_f);
              else
                  auxdata.push_back(xs_r);
          }
          else {
              if (bh.antisense_align())
                auxdata.push_back(xs_r);
              else
                auxdata.push_back(xs_f);
          }
      }
      else if (library_type == FR_SECONDSTRAND)   {
          if (insert_side == FRAG_LEFT || insert_side == FRAG_UNPAIRED){
              if (bh.antisense_align())
                auxdata.push_back(xs_r);
              else
                auxdata.push_back(xs_f);
          }
          else
          {
              if (bh.antisense_align())
                auxdata.push_back(xs_f);
              else
                auxdata.push_back(xs_r);
          }
      }
   } //bh.contiguous()
}

bool rewrite_sam_record(GBamWriter& bam_writer, const RefSequenceTable& rt,
                     const BowtieHit& bh,
                     const char* bwt_buf,
                     //char* rebuf,
                     char* read_alt_name,
                     const FragmentAlignmentGrade& grade,
                     FragmentType insert_side,
                     int num_hits,
                     const BowtieHit* next_hit)
{
	// Rewrite this hit, filling in the alt name, mate mapping
	// and setting the pair flag
	vector<string> sam_toks;
	tokenize(bwt_buf, "\t", sam_toks);
	char* slash_pos = strrchr(read_alt_name,'/');
	if (slash_pos) *slash_pos = 0;
	//*rebuf = 0;
	//read_alt_name as QNAME
    int flag=atoi(sam_toks[1].c_str()); //FLAG
    if (insert_side != FRAG_UNPAIRED) {
        //flag = atoi(sam_toks[1].c_str());
        // mark this as a singleton mate
        flag |= 0x0001;
        if (insert_side == FRAG_LEFT)
          flag |= 0x0040;
        else if (insert_side == FRAG_RIGHT)
          flag |= 0x0080;
        flag |= 0x0008;
        //char flag_buf[64];
        //sprintf(flag_buf, "%d", flag);
        //sam_toks[t] = flag_buf;
        }
    int gpos=isdigit(sam_toks[3][0]) ? atoi(sam_toks[3].c_str()) : 0;
    int mapQ=255;
    if (grade.num_alignments > 1)  {
      double err_prob = 1 - (1.0 / grade.num_alignments);
      mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
      }
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit);
    int insert_size=atoi(sam_toks[8].c_str());
    int mate_pos=atoi(sam_toks[7].c_str());
    GBamRecord* bamrec=bam_writer.new_record(read_alt_name, flag, sam_toks[2].c_str(), gpos, mapQ,
        sam_toks[5].c_str(), sam_toks[6].c_str(), mate_pos,
        insert_size, sam_toks[9].c_str(), sam_toks[10].c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
    return true;
}

bool rewrite_sam_record(GBamWriter& bam_writer, const RefSequenceTable& rt,
                     const BowtieHit& bh,
                     const char* bwt_buf,
                     //char* rebuf,
                     char* read_alt_name,
                     const InsertAlignmentGrade& grade,
                     FragmentType insert_side,
                     const BowtieHit* partner,
                     int num_hits,
                     const BowtieHit* next_hit)
{
    // Rewrite this hit, filling in the alt name, mate mapping
    // and setting the pair flag
    vector<string> sam_toks;
    tokenize(bwt_buf, "\t", sam_toks);
    char* slash_pos = strrchr(read_alt_name,'/');
    if (slash_pos) *slash_pos = 0;
    //read_alt_name as QNAME
    int flag = atoi(sam_toks[1].c_str());
    // 0x0010 (strand of query) is assumed to be set correctly
    // to begin with
    flag |= 0x0001; //it must be paired
    if (insert_side == FRAG_LEFT)
        flag |= 0x0040;
    else if (insert_side == FRAG_RIGHT)
        flag |= 0x0080;
    int gpos=isdigit(sam_toks[3][0]) ? atoi(sam_toks[3].c_str()) : 0;
    int mapQ=255;
    if (grade.num_alignments > 1) {
        double err_prob = 1 - (1.0 / grade.num_alignments);
        mapQ = (int)(-10.0 * log(err_prob) / log(10.0));
        }
    int mate_pos=atoi(sam_toks[7].c_str());
    if (partner) {
        if (partner->ref_id()==bh.ref_id())
             sam_toks[6] = "=";
        else {
          sam_toks[6] = rt.get_name(partner->ref_id());
          if (sam_toks[6].empty()) {
             //FIXME
             sam_toks[6] = "=";
             fprintf(stderr, "Warning: partner ref_id %d has no entry in ref table?\n", partner->ref_id());
             }
          }
        mate_pos = partner->left() + 1;
        if (grade.happy())
             flag |= 0x0002;
        if (partner->antisense_align())
             flag |=  0x0020;
        }
      else {
        sam_toks[6] = "*";
        mate_pos = 0;
        flag |= 0x0008;
        }
    vector<string> auxdata;
    add_auxData(auxdata, sam_toks, rt, bh, insert_side, num_hits, next_hit);
    int insert_size=atoi(sam_toks[8].c_str());
    //strcat(rebuf, "\n");
    GBamRecord* bamrec=bam_writer.new_record(read_alt_name, flag, sam_toks[2].c_str(), gpos, mapQ,
        sam_toks[5].c_str(), sam_toks[6].c_str(), mate_pos,
        insert_size, sam_toks[9].c_str(), sam_toks[10].c_str(), &auxdata);
    bam_writer.write(bamrec);
    delete bamrec;
    return true;
}

struct lex_hit_sort
{
  lex_hit_sort(const RefSequenceTable& rt, const HitsForRead& hits)
    : _rt(rt), _hits(hits)
  {}
  
  bool operator()(const uint32_t& l, const uint32_t& r) const
  {
    const BowtieHit& lhs = _hits.hits[l];
    const BowtieHit& rhs = _hits.hits[r];
      
    uint32_t l_id = lhs.ref_id();
    uint32_t r_id = rhs.ref_id();
    if (l_id != r_id)
      {
	return (strcmp(_rt.get_name(lhs.ref_id()), _rt.get_name(rhs.ref_id())) < 0);
      }
    return lhs.left() < rhs.left();
  }
  
  const RefSequenceTable& _rt;
  const HitsForRead& _hits;
};

void print_sam_for_hits(const RefSequenceTable& rt,
                        const HitsForRead& hits,
			const FragmentAlignmentGrade& grade,
			FragmentType frag_type,
			FILE* reads_file,
			//FILE* fout)
			GBamWriter& bam_writer)
{
  static const int buf_size = 2048;
  char read_name[buf_size];
  char read_seq[buf_size];
  char read_alt_name[buf_size];
  char read_quals[buf_size];
  
  //char rebuf[buf_size];

  lex_hit_sort s(rt, hits);
  vector<uint32_t> index_vector;
  for (size_t i = 0; i < hits.hits.size(); ++i)
      index_vector.push_back(i);
  
  sort(index_vector.begin(), index_vector.end(), s);
  
  bool got_read = get_read_from_stream(hits.insert_id, 
				       reads_file,
				       reads_format,
				       false,
				       read_name, 
				       read_seq,
				       read_alt_name,
				       read_quals);
  
  assert (got_read);
  
  for (size_t i = 0; i < hits.hits.size(); ++i)
    {
      size_t index = index_vector[i];
      const BowtieHit& bh = hits.hits[index];
      /*if (rewrite_sam_hit(rt,
              bh,
              bh.hitfile_rec().c_str(),
              rebuf,
              read_alt_name,
              grade,
              frag_type,
              hits.hits.size(),
              (i < hits.hits.size() - 1) ? &(hits.hits[index_vector[i+1]]) : NULL))
        {
      fprintf(fout, "%s", rebuf);
        }
    */
    rewrite_sam_record(bam_writer, rt,
         bh,
         bh.hitfile_rec().c_str(),
         read_alt_name,
         grade,
         frag_type,
         hits.hits.size(),
         (i < hits.hits.size()-1) ? &(hits.hits[index_vector[i+1]]) : NULL);
    }
}

void print_sam_for_hits(const RefSequenceTable& rt,
                        const HitsForRead& left_hits,
			const HitsForRead& right_hits,
			const InsertAlignmentGrade& grade,
			FILE* left_reads_file,
			FILE* right_reads_file,
			//FILE* fout)
			GBamWriter& bam_writer)
{
  assert (left_hits.insert_id == right_hits.insert_id);

  static const int buf_size = 2048;
  char left_read_name[buf_size];
  char left_read_seq[buf_size];
  char left_read_alt_name[buf_size];
  char left_read_quals[buf_size];
  //char left_rebuf[buf_size];
  
  char right_read_name[buf_size];
  char right_read_seq[buf_size];
  char right_read_alt_name[buf_size];
  char right_read_quals[buf_size];
  //char right_rebuf[buf_size];

  assert (left_hits.hits.size() == right_hits.hits.size() || 
	  (left_hits.hits.empty() || right_hits.hits.empty()));
    
  vector<uint32_t> index_vector;
  if(right_hits.hits.size() > 0)
    {
        lex_hit_sort s(rt, right_hits);
	for (size_t i = 0; i < right_hits.hits.size(); ++i)
	  index_vector.push_back(i);

	sort(index_vector.begin(), index_vector.end(), s);
    }
  else if (left_hits.hits.size() > 0)
    {
      lex_hit_sort s(rt, left_hits);
      for (size_t i = 0; i < left_hits.hits.size(); ++i)
	  index_vector.push_back(i);
  
      sort(index_vector.begin(), index_vector.end(), s);
    }
  
  bool got_left_read = get_read_from_stream(left_hits.insert_id, 
					    left_reads_file,
					    reads_format,
					    false,
					    left_read_name, 
					    left_read_seq,
					    left_read_alt_name,
					    left_read_quals);
  
  bool got_right_read = get_read_from_stream(right_hits.insert_id, 
					     right_reads_file,
					     reads_format,
					     false,
					     right_read_name, 
					     right_read_seq,
					     right_read_alt_name,
					     right_read_quals);

  assert (got_left_read && got_right_read);
   
  if (left_hits.hits.size() == right_hits.hits.size())
    {
      for (size_t i = 0; i < right_hits.hits.size(); ++i)
	{
	  size_t index = index_vector[i];
	  const BowtieHit& right_bh = right_hits.hits[index];
	  const BowtieHit& left_bh = left_hits.hits[index];
      /*if (rewrite_sam_hit(rt,
                  right_bh,
                  right_bh.hitfile_rec().c_str(),
                  right_rebuf,
                  right_read_alt_name,
                  grade,
                  FRAG_RIGHT,
                  &left_bh,
                  right_hits.hits.size(),
                  (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL))
            {
          fprintf(fout, "%s", right_rebuf);
        }

      if (rewrite_sam_hit(rt,
                  left_bh,
                  left_bh.hitfile_rec().c_str(),
                  left_rebuf,
                  left_read_alt_name,
                  grade,
                  FRAG_LEFT,
                  &right_bh,
                  left_hits.hits.size(),
                  (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL))
        {
          fprintf(fout, "%s", left_rebuf);
        }
      */
      rewrite_sam_record(bam_writer, rt,
                  right_bh,
                  right_bh.hitfile_rec().c_str(),
                  right_read_alt_name,
                  grade,
                  FRAG_RIGHT,
                  &left_bh,
                  right_hits.hits.size(),
                  (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL);
      rewrite_sam_record(bam_writer, rt,
                 left_bh,
                 left_bh.hitfile_rec().c_str(),
                 left_read_alt_name,
                 grade,
                 FRAG_LEFT,
                 &right_bh,
                 left_hits.hits.size(),
                 (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL);
	 }
    }
  else if (left_hits.hits.empty())
    {
      for (size_t i = 0; i < right_hits.hits.size(); ++i)
      {
	  size_t index = index_vector[i];
	  const BowtieHit& bh = right_hits.hits[index];
      /*if (rewrite_sam_hit(rt,
                  bh,
                  bh.hitfile_rec().c_str(),
                  right_rebuf,
                  right_read_alt_name,
                  grade,
                  FRAG_RIGHT,
                  NULL,
                  right_hits.hits.size(),
                  (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL))
        fprintf(fout, "%s", right_rebuf);
      */
      rewrite_sam_record(bam_writer, rt,
              bh,
              bh.hitfile_rec().c_str(),
              right_read_alt_name,
              grade,
              FRAG_RIGHT,
              NULL,
              right_hits.hits.size(),
              (i < right_hits.hits.size() - 1) ? &(right_hits.hits[index_vector[i+1]]) : NULL);

	  }
    }
  else if (right_hits.hits.empty())
    {
      for (size_t i = 0; i < left_hits.hits.size(); ++i)
      {
	  size_t index = index_vector[i];
	  const BowtieHit& bh = left_hits.hits[index];
      /*if (rewrite_sam_hit(rt,
                  bh,
                  bh.hitfile_rec().c_str(),
                  left_rebuf,
                  left_read_alt_name,
                  grade,
                  FRAG_LEFT,
                  NULL,
                  left_hits.hits.size(),
                  (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL))

        fprintf(fout, "%s", left_rebuf);
      */
      rewrite_sam_record(bam_writer, rt,
                  bh,
                  bh.hitfile_rec().c_str(),
                  left_read_alt_name,
                  grade,
                  FRAG_LEFT,
                  NULL,
                  left_hits.hits.size(),
                  (i < left_hits.hits.size() - 1) ? &(left_hits.hits[index_vector[i+1]]) : NULL);

	  }
    }
  else
    {
      assert (false);
    }
}

/**
 * Given all of the hits for a particular read, update the read counts for insertions and deletions.
 * @param hits hits The alignments for a particular read
 * @param insertions Maps from an insertion to the number of supporting reads for that insertion
 * @param deletions Maps from a deletion to the number of supporting reads for that deletion
 */ 
void update_insertions_and_deletions(const HitsForRead& hits,
					InsertionSet& insertions,
					DeletionSet& deletions)
{
	for (size_t i = 0; i < hits.hits.size(); ++i)
	{
		const BowtieHit& bh = hits.hits[i];
		insertions_from_alignment(bh, insertions);
		deletions_from_alignment(bh, deletions);
	}
}					



void update_junctions(const HitsForRead& hits,
					  JunctionSet& junctions)
{
	for (size_t i = 0; i < hits.hits.size(); ++i)
	{
		const BowtieHit& bh = hits.hits[i];
		junctions_from_alignment(bh, junctions);
	}
}

// Extracts junctions from all the SAM hits (based on REF_SKIPs) in the hit file
// resets the stream when finished.
void exclude_hits_on_filtered_junctions(const JunctionSet& junctions,
										HitsForRead& hits)
{
	HitsForRead remaining;
	remaining.insert_id = hits.insert_id;

	for (size_t i = 0; i < hits.hits.size(); ++i)
	{
		BowtieHit& bh = hits.hits[i];
		bool filter_hit = false;
		if (!bh.contiguous())
		{
			JunctionSet bh_juncs;
			junctions_from_alignment(bh, bh_juncs);
			for (JunctionSet::iterator itr = bh_juncs.begin();
				 itr != bh_juncs.end();
				 itr++)
			{
				const Junction& j = itr->first;
				JunctionSet::const_iterator target = junctions.find(j);
				if (target == junctions.end() || !target->second.accepted)
				{
					filter_hit = true;
					break;
				}
			}
		}
		if (!filter_hit)
		{
			remaining.hits.push_back(bh);
		}
		else
		{
			//fprintf(stderr, "Filtering hit\n");
		}
	}
	hits = remaining;
}

void get_junctions_from_best_hits(HitStream& left_hs,
				  HitStream& right_hs,
				  ReadTable& it, 
				  JunctionSet& junctions)
{
	HitsForRead curr_left_hit_group;
	HitsForRead curr_right_hit_group;
	
	left_hs.next_read_hits(curr_left_hit_group);
	right_hs.next_read_hits(curr_right_hit_group);
	
	uint32_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	uint32_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	
	// While we still have unreported hits...
	while(curr_left_obs_order != 0xFFFFFFFF || 
		  curr_right_obs_order != 0xFFFFFFFF)
	{		
		// Chew up left singletons
		while (curr_left_obs_order < curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_left_obs_order;
			FragmentAlignmentGrade grade;
			
			// Process hits for left singleton, select best alignments
			fragment_best_alignments(curr_left_hit_group, grade, best_hits);
			
			update_junctions(best_hits, junctions);
			
			// Get next hit group
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
		}
		
		// Chew up right singletons
		while (curr_left_obs_order > curr_right_obs_order &&
			   curr_right_obs_order != 0xFFFFFFFF)
		{
			HitsForRead best_hits;
			best_hits.insert_id = curr_right_obs_order;
			FragmentAlignmentGrade grade;
			
			// Process hit for right singleton, select best alignments
			fragment_best_alignments(curr_right_hit_group,grade, best_hits);
			
			update_junctions(best_hits, junctions);
			
			// Get next hit group
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
		
		// Since we have both left hits and right hits for this insert,
		// Find the best pairing and print both
		while (curr_left_obs_order == curr_right_obs_order &&
			   curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
		{			 
			if (curr_left_hit_group.hits.empty())
			{
				HitsForRead right_best_hits;
				right_best_hits.insert_id = curr_right_obs_order;
				
				FragmentAlignmentGrade grade;
				fragment_best_alignments(curr_right_hit_group, grade, right_best_hits);
				
				update_junctions(right_best_hits, junctions);
			}
			else if (curr_right_hit_group.hits.empty())
			{
				HitsForRead left_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				
				FragmentAlignmentGrade grade;
				// Process hits for left singleton, select best alignments
				fragment_best_alignments(curr_left_hit_group, grade, left_best_hits);
				
				update_junctions(left_best_hits, junctions);
			}
			else
			{		
				HitsForRead left_best_hits;
				HitsForRead right_best_hits;
				left_best_hits.insert_id = curr_left_obs_order;
				right_best_hits.insert_id = curr_right_obs_order;
				
				InsertAlignmentGrade grade;
				insert_best_alignments(curr_left_hit_group, 
						       curr_right_hit_group, 
						       grade,
						       left_best_hits,
						       right_best_hits);
				
				update_junctions(left_best_hits, junctions);
				update_junctions(right_best_hits, junctions);
			}
			
			left_hs.next_read_hits(curr_left_hit_group);
			curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
			
			right_hs.next_read_hits(curr_right_hit_group);
			curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
		}
	}
	left_hs.reset();
	right_hs.reset();
}


void driver(GBamWriter& bam_writer,
//void driver(FILE* accepted_hits_out,
        string& left_map_fname,
	    FILE* left_reads,
        string& right_map_fname,
	    FILE* right_reads,
        FILE* junctions_out,
	    FILE* insertions_out,
	    FILE* deletions_out)
{	
    ReadTable it;
    RefSequenceTable rt(sam_header, true);
    
    //SAMHitFactory hit_factory(it,rt);
    BAMHitFactory hit_factory(it,rt);
    
    JunctionSet junctions;
    {
      HitStream l_hs(left_map_fname, &hit_factory, false, true, true, true);
      HitStream r_hs(right_map_fname, &hit_factory, false, true, true, true);
      get_junctions_from_best_hits(l_hs, r_hs, it, junctions);
      //this resets the streams
    }
    HitStream left_hs(left_map_fname, &hit_factory, false, true, true, true);
    HitStream right_hs(right_map_fname, &hit_factory, false, true, true, true);

    size_t num_unfiltered_juncs = junctions.size();
    
    fprintf(stderr, "Loaded %lu junctions\n", (long unsigned int) num_unfiltered_juncs);
    
    HitsForRead curr_left_hit_group;
    HitsForRead curr_right_hit_group;
    
    left_hs.next_read_hits(curr_left_hit_group);
    right_hs.next_read_hits(curr_right_hit_group);
    
    uint32_t curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
    uint32_t curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
    
    // Read hits, extract junctions, and toss the ones that arent strongly enough supported.
    
    filter_junctions(junctions);
    //size_t num_juncs_after_filter = junctions.size();
    //fprintf(stderr, "Filtered %lu junctions\n", num_unfiltered_juncs - num_juncs_after_filter);
    
    size_t small_overhangs = 0;
    for (JunctionSet::iterator i = junctions.begin(); i != junctions.end(); ++i)
      {
	if (i->second.accepted && 
	    (i->second.left_extent < min_anchor_len || i->second.right_extent < min_anchor_len))
	  {
	    small_overhangs++;
	  }
      }
    
    if (small_overhangs >0)
      fprintf(stderr, "Warning: %lu small overhang junctions!\n", (long unsigned int)small_overhangs);
    
    JunctionSet final_junctions; // the junctions formed from best hits
    InsertionSet final_insertions;
    DeletionSet final_deletions;
    
    fprintf (stderr, "Reporting final accepted alignments...");
    
    // While we still have unreported hits...
    while(curr_left_obs_order != 0xFFFFFFFF || 
	  curr_right_obs_order != 0xFFFFFFFF)
      {		
	// Chew up left singletons
	while (curr_left_obs_order < curr_right_obs_order &&
	       curr_left_obs_order != 0xFFFFFFFF)
	  {
	    HitsForRead best_hits;
	    best_hits.insert_id = curr_left_obs_order;
	    FragmentAlignmentGrade grade;
	    
	    exclude_hits_on_filtered_junctions(junctions, curr_left_hit_group);
	    
	    // Process hits for left singleton, select best alignments
	    fragment_best_alignments(curr_left_hit_group, grade, best_hits);

	    if (best_hits.hits.size() <= max_multihits)
	      {
		update_junctions(best_hits, final_junctions);
		update_insertions_and_deletions(best_hits, final_insertions, final_deletions);
	    
		print_sam_for_hits(rt,
				   best_hits, 
				   grade,
				   (right_map_fname.empty() ? FRAG_UNPAIRED : FRAG_LEFT),
				   left_reads,
				   //accepted_hits_out);
				   bam_writer);
	      }
	    
	    // Get next hit group
	    left_hs.next_read_hits(curr_left_hit_group);
	    curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	  }

	// Chew up right singletons
	while (curr_left_obs_order > curr_right_obs_order &&
	       curr_right_obs_order != 0xFFFFFFFF)
	  {
	    HitsForRead best_hits;
	    best_hits.insert_id = curr_right_obs_order;
	    FragmentAlignmentGrade grade;

	    exclude_hits_on_filtered_junctions(junctions, curr_right_hit_group);

	    // Process hit for right singleton, select best alignments
	    fragment_best_alignments(curr_right_hit_group, grade, best_hits);

	    if (best_hits.hits.size() <= max_multihits)
	      {
		update_junctions(best_hits, final_junctions);
		update_insertions_and_deletions(best_hits, final_insertions, final_deletions);
	    
		print_sam_for_hits(rt,
				   best_hits, 
				   grade, 
				   FRAG_RIGHT,
				   right_reads, 
				   //accepted_hits_out);
				   bam_writer);
	      }
	    
	    // Get next hit group
	    right_hs.next_read_hits(curr_right_hit_group);
	    curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	  }

	// Since we have both left hits and right hits for this insert,
	// Find the best pairing and print both
	while (curr_left_obs_order == curr_right_obs_order &&
	       curr_left_obs_order != 0xFFFFFFFF && curr_right_obs_order != 0xFFFFFFFF)
	  {			 
	    exclude_hits_on_filtered_junctions(junctions, curr_left_hit_group);
	    exclude_hits_on_filtered_junctions(junctions, curr_right_hit_group);
	    
	    if (curr_left_hit_group.hits.empty())
	      {
		HitsForRead right_best_hits;
		right_best_hits.insert_id = curr_right_obs_order;
		
		FragmentAlignmentGrade grade;
		fragment_best_alignments(curr_right_hit_group, grade, right_best_hits);

		if (right_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(right_best_hits, final_junctions);
		    update_insertions_and_deletions(right_best_hits, final_insertions, final_deletions);
				
		    print_sam_for_hits(rt,
				       right_best_hits, 
				       grade, 
				       FRAG_RIGHT,
				       right_reads, 
				       //accepted_hits_out);
				       bam_writer);
		  }
	      }
	    else if (curr_right_hit_group.hits.empty())
	      {
		HitsForRead left_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
		
		FragmentAlignmentGrade grade;
		// Process hits for left singleton, select best alignments
		fragment_best_alignments(curr_left_hit_group, grade, left_best_hits);

		if (left_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(left_best_hits, final_junctions);
		    update_insertions_and_deletions(left_best_hits, final_insertions, final_deletions);			
		    
		    print_sam_for_hits(rt,
				       left_best_hits, 
				       grade,
				       FRAG_LEFT,
				       left_reads, 
				       //accepted_hits_out);
				       bam_writer);
		  }
	      }
	    else
	      {		
		HitsForRead left_best_hits;
		HitsForRead right_best_hits;
		left_best_hits.insert_id = curr_left_obs_order;
		right_best_hits.insert_id = curr_right_obs_order;
		
		InsertAlignmentGrade grade;
		insert_best_alignments(curr_left_hit_group, 
				       curr_right_hit_group, 
				       grade,
				       left_best_hits,
				       right_best_hits);

		if (left_best_hits.hits.size() <= max_multihits && right_best_hits.hits.size() <= max_multihits)
		  {
		    update_junctions(left_best_hits, final_junctions);
		    update_junctions(right_best_hits, final_junctions);
		    update_insertions_and_deletions(left_best_hits, final_insertions, final_deletions);
		    update_insertions_and_deletions(right_best_hits, final_insertions, final_deletions);
		
		    print_sam_for_hits(rt,
				       left_best_hits,
				       right_best_hits,
				       grade,
				       left_reads, 
				       right_reads, 
				       //accepted_hits_out);
				       bam_writer);
		  }
	      }
	    
	    left_hs.next_read_hits(curr_left_hit_group);
	    curr_left_obs_order = it.observation_order(curr_left_hit_group.insert_id);
	    
	    right_hs.next_read_hits(curr_right_hit_group);
	    curr_right_obs_order = it.observation_order(curr_right_hit_group.insert_id);
	  }
      }

    fprintf (stderr, "done\n");
    
    small_overhangs = 0;
    for (JunctionSet::iterator i = final_junctions.begin(); i != final_junctions.end();)
      {
	const JunctionStats& stats = i->second;
	if (i->second.supporting_hits == 0 || 
	    i->second.left_extent < 8 || 
	    i->second.right_extent < 8)
	  {
	    final_junctions.erase(i++);
	  }
	else
	  {
	    ++i;
	  }
      }
    
    
    
//	if (small_overhangs > 0)
//		fprintf(stderr, "Warning: %lu small overhang junctions!\n", small_overhangs);
	
    fprintf (stderr, "Printing junction BED track...");
    print_junctions(junctions_out, final_junctions, rt);
    fprintf (stderr, "done\n");
    
    fprintf (stderr, "Printing insertions...");
    print_insertions(insertions_out, final_insertions,rt);
    fclose(insertions_out);
    fprintf (stderr, "done\n");
    
    fprintf (stderr, "Printing deletions...");
    print_deletions(deletions_out, final_deletions, rt);
    fclose(deletions_out);
    fprintf (stderr, "done\n");
    
    fprintf(stderr, "Found %lu junctions from happy spliced reads\n", (long unsigned int)final_junctions.size());
}

void print_usage()
{
    fprintf(stderr, "Usage:   tophat_reports <junctions.bed> <insertions.vcf> <deletions.vcf> <accepted_hits.sam> <left_map1,...,left_mapN> <left_reads.fq>  [right_map1,...,right_mapN] [right_reads.fq]\n");
	
	//    fprintf(stderr, "Usage:   tophat_reports <coverage.wig> <junctions.bed> <accepted_hits.sam> <map1.bwtout> [splice_map1.sbwtout]\n");
}

int main(int argc, char** argv)
{
  fprintf(stderr, "tophat_reports v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION); 
  fprintf(stderr, "---------------------------------------\n");
  
  reads_format = FASTQ;
  
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
	
    string accepted_hits_file_name = argv[optind++];
	
    if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string left_map_filename = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
        return 1;
    }
  //FZPipe left_map_file;
  //string unbamcmd=getBam2SamCmd(left_map_filename);
  //left_map_file.openRead(left_map_filename, unbamcmd);
  string left_reads_filename = argv[optind++];
  string unzcmd=getUnpackCmd(left_reads_filename, false);
  
  string right_map_filename;
  string right_reads_filename;
  FZPipe right_reads_file;
	
  if (optind < argc)
    {
      right_map_filename = argv[optind++];      
      if(optind >= argc) {
	     print_usage();
	     return 1;
	     }
      //right_map_file.openRead(right_map_filename, unbamcmd);
  //if (optind<argc) {
     right_reads_filename=argv[optind++];
     right_reads_file.openRead(right_reads_filename,unzcmd);
   //  }
    }
  FILE* junctions_file = fopen(junctions_file_name.c_str(), "w");
  if (junctions_file == NULL)
    {
      fprintf(stderr, "Error: cannot open BED file %s for writing\n",
	      junctions_file_name.c_str());
      exit(1);
    }

    FILE* insertions_file = fopen(insertions_file_name.c_str(), "w");
    if (insertions_file == NULL)
    {
	fprintf(stderr, "Error: cannot open VCF file %s for writing\n",
		insertions_file_name.c_str());
	exit(1);
    }


    FILE* deletions_file = fopen(deletions_file_name.c_str(), "w");
    if (deletions_file == NULL)
    {
	fprintf(stderr, "Error: cannot open VCF file %s for writing\n",
		deletions_file_name.c_str());
	exit(1);
    }
    /*
    FILE* accepted_hits_file=NULL;
    if (accepted_hits_file_name == "-") {
         accepted_hits_file=stdout;
         }
       else {
         //accepted_hits_file = fopen(accepted_hits_file_name.c_str(), "a");
         accepted_hits_file = fopen(accepted_hits_file_name.c_str(), "w");
         if (accepted_hits_file == NULL) {
           fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
                           accepted_hits_file_name.c_str());
           exit(1);
           }
         }
   writeSamHeader(accepted_hits_file);
  */
  bool uncompressed_bam=(accepted_hits_file_name=="-");
  GBamWriter bam_writer(accepted_hits_file_name.c_str(), sam_header.c_str(), uncompressed_bam);
  FZPipe left_reads_file(left_reads_filename, unzcmd);
  if (left_reads_file.file==NULL)
    {
      fprintf(stderr, "Error: cannot open reads file %s for reading\n",
        left_reads_filename.c_str());
        exit(1);
    }

  driver(bam_writer, left_map_filename,
  //driver(accepted_hits_file, left_map_filename,
	 left_reads_file.file,
     right_map_filename,
	 right_reads_file.file,
	 junctions_file,
	 insertions_file,
	 deletions_file);
  //if (accepted_hits_file!=stdout)
  //  fclose(accepted_hits_file);
  return 0;
}


