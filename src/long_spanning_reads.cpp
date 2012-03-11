#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  long_spanning_reads.cpp
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
#include <cstdlib>
#include <cstring>
#include <bitset>
//#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <getopt.h>

#include <boost/thread.hpp>

#include "common.h"
#include "utils.h"
#include "bwt_map.h"
#include "tokenize.h"
#include "segments.h"
#include "reads.h"

#include "junctions.h"
#include "insertions.h"
#include "deletions.h"
#include "fusions.h"

using namespace seqan;
using namespace std;

// daehwan
bool bDebug = false;

void print_usage()
{
    fprintf(stderr, "Usage:   long_spanning_reads <reads.fa/.fq> <possible_juncs1,...,possible_juncsN> <possible_insertions1,...,possible_insertionsN> <possible_deletions1,...,possible_deletionsN> <seg1.bwtout,...,segN.bwtout> [spliced_seg1.bwtout,...,spliced_segN.bwtout]\n");
}

bool key_lt(const pair<uint32_t, HitsForRead>& lhs,
      const pair<uint32_t, HitsForRead>& rhs)
{
  return lhs.first < rhs.first;
}

void get_seqs(istream& ref_stream,
        RefSequenceTable& rt,
        bool keep_seqs = true)
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
      seqan::read(ref_stream, *ref_str, Fasta());

      rt.get_id(name, keep_seqs ? ref_str : NULL, 0);
    }
}


void look_right_for_hit_group(ReadTable& unmapped_reads,
            vector<HitStream>& contig_hits,
            size_t curr_file,
            vector<HitStream>& spliced_hits,
            const HitsForRead& targets,
            vector<HitsForRead>& seg_hits_for_read)
{
  int right_file = curr_file + 1;

  HitStream& right = contig_hits[right_file];
  uint64_t curr_next_group_id = targets.insert_id;
  int curr_order = unmapped_reads.observation_order(curr_next_group_id);

  assert (curr_order != -1);
  while(true)
    {
      HitsForRead hit_group;
      uint64_t right_next_group_id = right.next_group_id();
      int right_order = unmapped_reads.observation_order(right_next_group_id);

      // If we would have seen the hits by now, bail out.
      if (curr_order < right_order || right_order == -1)
  {
    break;
  }
      if (right.next_read_hits(hit_group))
  {
    if (hit_group.insert_id == targets.insert_id)
      {
        // Some of the targets may be missing, we need to
        // process them individually
        seg_hits_for_read[right_file] = hit_group;
        break;
      }
  }
    }

  HitsForRead& curr_seg_hits = seg_hits_for_read[right_file];

  if (right_file < (int)spliced_hits.size() && right_file >= 0)
    {    
      // Scan forward in the spliced hits file for this hit group
      HitsForRead spliced_group;
      HitsForRead curr_spliced_group;
      while (spliced_hits[right_file].next_group_id() > 0 &&
       spliced_hits[right_file].next_group_id() <= (uint32_t)curr_order)
  {
    spliced_hits[right_file].next_read_hits(curr_spliced_group);

    if (curr_spliced_group.insert_id == (uint32_t)curr_order)
      {
        spliced_group = curr_spliced_group;
        break;
      }
  }

      if (!spliced_group.hits.empty())
  {
    curr_seg_hits.insert_id = spliced_group.insert_id;
    curr_seg_hits.hits.insert(curr_seg_hits.hits.end(),
            spliced_group.hits.begin(),
            spliced_group.hits.end());
  }
    }

  if (curr_seg_hits.hits.empty())
    return;
  else if (right_file + 1 < (int)contig_hits.size())
    {
      look_right_for_hit_group(unmapped_reads,
            contig_hits,
            curr_file + 1,
            spliced_hits,
            curr_seg_hits,
            seg_hits_for_read);
    }
}

BowtieHit merge_chain_color(RefSequenceTable& rt,
			    const string& read_seq,
			    const string& read_quals,
			    std::set<Junction>& possible_juncs,
			    std::set<Insertion>& possible_insertions,
			    list<BowtieHit>& hit_chain)
{
  bool antisense = hit_chain.front().antisense_align();
  uint32_t reference_id = hit_chain.front().ref_id();
  uint64_t insert_id = hit_chain.front().insert_id();

  int left = hit_chain.front().left();

  list<BowtieHit>::iterator prev_hit = hit_chain.begin();
  list<BowtieHit>::iterator curr_hit = ++(hit_chain.begin());

  string seq;
  string qual;
  int old_read_length = 0;
  int first_seg_length = hit_chain.front().seq().length();
  for (list<BowtieHit>::iterator i = hit_chain.begin(); i != hit_chain.end(); ++i)
    {
      seq += i->seq();
      qual += i->qual();
      old_read_length += i->read_len();
    }

  string rev_read_seq, rev_read_quals;
  if (color && antisense)
    {
      rev_read_seq = read_seq;
      reverse(rev_read_seq.begin() + 1, rev_read_seq.end());

      rev_read_quals = read_quals;
      reverse(rev_read_quals.begin(), rev_read_quals.end());
    }

  while (curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
    {
      /*
       * Note that the gap size may be negative, since left() and right() return
       * signed integers, this will be OK.
       */
      int gap  = curr_hit->left() - prev_hit->right();
      if (gap < -(int)max_insertion_length ||
    (gap > (int)max_deletion_length &&
     (gap < min_report_intron_length || gap > max_report_intron_length)))
  {
    return BowtieHit();
  }

      ++prev_hit;
      ++curr_hit;
    }

  prev_hit = hit_chain.begin();
  curr_hit = ++(hit_chain.begin());

  RefSequenceTable::Sequence* ref_str = rt.get_seq(prev_hit->ref_id());
  if (!ref_str)
    return BowtieHit();

  int curr_seg_index = 1;
  while (curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
    {
      /*
       * This code is assuming that the cigar strings end and start with a match
       * While segment alignments can actually end with a junction, insertion or deletion, the hope
       * is that in those cases, the right and left ends of the alignments will correctly
       * line up, so we won't get to this bit of code
       */
      if (prev_hit->cigar().back().opcode != MATCH || curr_hit->cigar().front().opcode != MATCH)
  {
    return BowtieHit();
  }

      if (prev_hit->is_spliced() && curr_hit->is_spliced() && prev_hit->antisense_splice() != curr_hit->antisense_splice())
  {
    /*
     * There is no way that we can splice these together into a valid
     * alignment
     */
    return BowtieHit();
  }

      bool found_closure = false;
      /*
       * Note that antisense_splice is the same for prev and curr
       */
      bool antisense_closure = prev_hit->is_spliced() ? prev_hit->antisense_splice() : curr_hit->antisense_splice();
      vector<CigarOp> new_cigar;
      int new_left = -1;
      int mismatch = 0;

      /*
       * Take the length of matched bases in each segment into consideration for closures,
       * this can be a problem for reads of variable lengths.
       */
      int prev_right_end_match_length = prev_hit->cigar().back().length;
      int curr_left_end_match_length = curr_hit->cigar().front().length;

      if (prev_hit->right() > curr_hit->left())
  {
    std::set<Insertion>::iterator lb, ub;
    /*
     * Note, this offset is determined by the min-anchor length supplied to
     * juncs_db, which is currently hard-coded at 3 in tophat.py
     * as this value determines what sort of read segments
     * should be able to align directly to the splice sequences
     */
    int left_boundary = prev_hit->right() - 4;
    int right_boundary = curr_hit->left() + 4;

    /*
     * Create a dummy sequence to represent the maximum possible insertion
     */
    std::string maxInsertedSequence = "";
    maxInsertedSequence.resize(max_insertion_length,'A');

    lb = possible_insertions.upper_bound(Insertion(reference_id, left_boundary, ""));
    ub = possible_insertions.upper_bound(Insertion(reference_id, right_boundary, maxInsertedSequence));

    int reference_mismatch = 0;

    while (lb != ub && lb != possible_insertions.end())
      {
        /*
         * In the following code, we will check to make sure that the segments have the proper
         * separation and sequence for the insertions, and generate the appropriate merged bowtie hit
         * In general, reads with insertions must match the inserted sequence exactly.
         */
        if (((int)lb->sequence.size()) == (prev_hit->right() - curr_hit->left()))
    {
      /*
       * Check we have enough matched bases on prev or curr segment.
       */
      int insert_to_prev_right = prev_hit->right() - lb->left - 1;
      int curr_left_to_insert = lb->left - curr_hit->left() + 1;
      if (insert_to_prev_right > prev_right_end_match_length || curr_left_to_insert > curr_left_end_match_length)
        {
          ++lb;
          continue;
        }

      /*
       * Keep track of how many mismatches were made to the genome in the region
       * where we should actually be matching against the insertion
       */
      int this_reference_mismatch = 0;
      int insertion_mismatch = 0;
      int insertion_len = lb->sequence.length();
      const seqan::Dna5String insertionSequence = seqan::Dna5String(lb->sequence);

      /*
       * First check to see if we need to adjust number of observed errors for the left (prev)
       * hit. This is only the case if this segment runs into the insertion. To be consistent
       * with bwt_map.cpp, we will not allow a segment to have errors in the insertion region
       */
      string colorSegmentSequence_prev;
      if (insert_to_prev_right > 0)
        {
          const seqan::Dna5String referenceSequence = seqan::infix(*ref_str, lb->left + 1, prev_hit->right());
          const seqan::Dna5String oldSegmentSequence = seqan::Dna5String(prev_hit->seq().substr(prev_hit->seq().length() - insert_to_prev_right));

          if (color)
      {
        string color;

        if (antisense)
            color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - insert_to_prev_right - 2, insert_to_prev_right + 1);
        else
            color = read_seq.substr(curr_seg_index * segment_length - insert_to_prev_right, insert_to_prev_right + 1);

        color[0] = prev_hit->seq()[segment_length - insert_to_prev_right - 1];
        colorSegmentSequence_prev = convert_color_to_bp(color);
      }

          const seqan::Dna5String newSegmentSequence = color ? colorSegmentSequence_prev : oldSegmentSequence;

          /*
           * Scan right in the read until we run out of read
           */
          for (int read_index = 0; read_index < insert_to_prev_right; ++read_index)
      {
        /*
         * Any mismatch to the insertion is a failure
         */
        if (referenceSequence[read_index] == 'N' || referenceSequence[read_index] != oldSegmentSequence[read_index])
          {
            ++this_reference_mismatch;
          }

        if (read_index < insertion_len)
          {
            if (insertionSequence[read_index] == 'N' || insertionSequence[read_index] != newSegmentSequence[read_index])
        {
          ++insertion_mismatch;
          break;
        }
          }
        else
          {
            if (referenceSequence[read_index - insertion_len] == 'N' ||
          referenceSequence[read_index - insertion_len] != newSegmentSequence[read_index])
        {
          --this_reference_mismatch;
        }
          }
      }
        }

      string colorSegmentSequence_curr;
      if (curr_left_to_insert > 0)
        {
          const seqan::Dna5String referenceSequence = seqan::infix(*ref_str, curr_hit->left(), lb->left + 1);
          const seqan::Dna5String oldSegmentSequence = seqan::Dna5String(curr_hit->seq().substr(0, curr_left_to_insert));

          if (color)
      {
        string color;
        if (antisense)
            color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length), curr_left_to_insert);
        else
            color = read_seq.substr(curr_seg_index * segment_length + 2, curr_left_to_insert);

        color.push_back(curr_hit->seq()[curr_left_to_insert]);
        reverse(color.begin(), color.end());
        string bp = convert_color_to_bp(color);
        reverse(bp.begin(), bp.end());
        colorSegmentSequence_curr = bp;
      }

          const seqan::Dna5String newSegmentSequence = color ? colorSegmentSequence_curr : oldSegmentSequence;

          /*
           * Scan left in the read until
           * We ran out of read sequence (insertion extends past segment)
           */
          for (int read_index = 0; read_index < curr_left_to_insert; ++read_index)
      {
        int segmentPosition = curr_left_to_insert - read_index - 1;
        int insertionPosition = insertion_len - read_index - 1;

        if (referenceSequence[segmentPosition] == 'N' || (referenceSequence[segmentPosition] != oldSegmentSequence[segmentPosition]))
          {
            ++this_reference_mismatch;
          }

        if (read_index < insertion_len)
          {
            if (insertionSequence[insertionPosition] == 'N' || (insertionSequence[insertionPosition] != newSegmentSequence[segmentPosition]))
        {
          ++insertion_mismatch;
          break;
        }
          }
        else
          {
            if (referenceSequence[segmentPosition + insertion_len] == 'N' ||
          (referenceSequence[segmentPosition + insertion_len] != newSegmentSequence[segmentPosition]))
        {
          --this_reference_mismatch;
        }
          }
      }
        }

      if (found_closure)
        {
    fprintf(stderr, "Warning: multiple closures found for insertion read # %d\n", (int)insert_id);
    return BowtieHit();
        }

      if (insertion_mismatch == 0)
        {
    reference_mismatch = this_reference_mismatch;
    mismatch = -reference_mismatch;
    found_closure = true;
    new_left = prev_hit->left();
    new_cigar = prev_hit->cigar();

    /*
     * Need to make a new insert operation between the two match character that begin
     * and end the intersection of these two junction. Note that we necessarily assume
     * that this insertion can't span beyond the boundaries of these reads. That should
     * probably be better enforced somewhere
     */

    new_cigar.back().length -= insert_to_prev_right;
    if (new_cigar.back().length <= 0)
      new_cigar.pop_back();

    new_cigar.push_back(CigarOp(INS, lb->sequence.size()));
    vector<CigarOp> new_right_cigar = curr_hit->cigar();
    new_right_cigar.front().length += (insert_to_prev_right - lb->sequence.size());

    /*
     * Finish stitching together the new cigar string
     */
    size_t c = new_right_cigar.front().length > 0 ? 0 : 1;
    for (; c < new_right_cigar.size(); ++c)
      {
        new_cigar.push_back(new_right_cigar[c]);
      }

    if (color)
      {
        if (insert_to_prev_right > 0)
          seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length - insert_to_prev_right, insert_to_prev_right, colorSegmentSequence_prev);

        if (curr_left_to_insert > 0)
          seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length, curr_left_to_insert, colorSegmentSequence_curr);
      }
        }
    }
        ++lb;
  }

  if (!found_closure)
  {
    return BowtieHit();
  }
      }

      /*
       * Stitch segments together using juctions or deletions if necessary.
       */
      else if (prev_hit->right() < curr_hit->left())
  {
    std::set<Junction>::iterator lb, ub;

    int left_boundary = prev_hit->right() - 4;
    int right_boundary = curr_hit->left() + 4;

    lb = possible_juncs.upper_bound(Junction(reference_id, left_boundary, right_boundary - 8, true));
    ub = possible_juncs.lower_bound(Junction(reference_id, left_boundary + 8, right_boundary, false));

    int new_diff_mismatches = 0xff;
    while (lb != ub && lb != possible_juncs.end())
      {
        int dist_to_left = lb->left - prev_hit->right() + 1;
        int dist_to_right = lb->right - curr_hit->left();

        if (abs(dist_to_left) <= 4 && abs(dist_to_right) <= 4 && dist_to_left == dist_to_right)
    {
      /*
       * Check we have enough matched bases on prev or curr segment.
       */
      if (dist_to_left > curr_left_end_match_length || -dist_to_left > prev_right_end_match_length )
        {
          ++lb;
          continue;
        }

      Dna5String new_cmp_str, old_cmp_str;
      int new_mismatch = 0, old_mismatch = 0;
      string new_patch_str; // this is for colorspace reads
      if (dist_to_left > 0)
        {
          new_cmp_str = seqan::infix(*ref_str, prev_hit->right(), lb->left + 1);
          old_cmp_str = seqan::infix(*ref_str, curr_hit->left(), lb->right);

          string new_seq;
          if (color)
      {
        string ref = DnaString_to_string(seqan::infix(*ref_str, prev_hit->right() - 1, lb->left + 1));

        string color, qual;
        if (antisense)
          {
            color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - 1, dist_to_left);
            qual = rev_read_quals.substr(rev_read_quals.length() - (curr_seg_index * segment_length) - 1, dist_to_left);
          }
          else
          {
            color = read_seq.substr(1 + curr_seg_index * segment_length, dist_to_left);
            qual = read_quals.substr(curr_seg_index * segment_length, dist_to_left);
          }

        BWA_decode(color, qual, ref, new_seq);
        new_seq = new_seq.substr(1);
      }

          const string& curr_old_seq = curr_hit->seq();
          const string& curr_seq = color ? new_seq : curr_hit->seq();
          for (int i = 0; i < dist_to_left; ++i)
      {
        if (curr_seq[i] != new_cmp_str[i])
          ++new_mismatch;

        if (curr_old_seq[i] != old_cmp_str[i])
          ++old_mismatch;
      }

          if (color)
      new_patch_str = curr_seq.substr(0, dist_to_left);
        }
      else if (dist_to_left < 0)
        {
          new_cmp_str = seqan::infix(*ref_str, lb->right, curr_hit->left());
          old_cmp_str = seqan::infix(*ref_str, lb->left + 1, prev_hit->right());

          size_t abs_dist = -dist_to_left;
          string new_seq;
          if (color)
      {
        string ref = DnaString_to_string(seqan::infix(*ref_str, lb->left, lb->left + 1));
        ref += DnaString_to_string(seqan::infix(*ref_str, lb->right, curr_hit->left()));

        string color, qual;
        if (antisense)
          {
            color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - 1 - abs_dist, abs_dist);
            qual = rev_read_quals.substr(rev_read_quals.length() - (curr_seg_index * segment_length) - 1 - abs_dist, abs_dist);
          }
          else
          {
            color = read_seq.substr(1 + curr_seg_index * segment_length - abs_dist, abs_dist);
            qual = read_quals.substr(curr_seg_index * segment_length - abs_dist, abs_dist);
          }

        BWA_decode(color, qual, ref, new_seq);
        new_seq = new_seq.substr(1);
      }

          const string& prev_old_seq = prev_hit->seq();
          size_t prev_old_seq_len = prev_old_seq.length();
          const string& prev_seq = color ? new_seq : prev_hit->seq();
          size_t prev_seq_len = prev_seq.length();
          for (size_t i = 0; i < abs_dist; ++i)
      {
        if (prev_seq[prev_seq_len - (abs_dist - i)] != new_cmp_str[i])
          ++new_mismatch;
        if (prev_old_seq[prev_old_seq_len - (abs_dist - i)] != old_cmp_str[i])
          ++old_mismatch;
      }

          if (color)
      new_patch_str = prev_seq.substr(prev_seq_len - abs_dist, abs_dist);
        }

      int temp_diff_mismatches = new_mismatch - old_mismatch;
      if (temp_diff_mismatches >= new_diff_mismatches || new_mismatch >= 2)
        {
          ++lb;
          continue;
        }

      if (color)
        {
          /*
           * We need to recover the origianl sequence.
           */
          if (found_closure)
      {
        seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length - 4, 8,
              prev_hit->seq().substr(prev_hit->seq().length() - 4) + curr_hit->seq().substr(0, 4));
      }

          if (dist_to_left > 0)
      seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length, dist_to_left, new_patch_str);
          else if (dist_to_left < 0)
      seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length + dist_to_left, -dist_to_left, new_patch_str);
        }

      new_diff_mismatches = temp_diff_mismatches;

      new_left = prev_hit->left();
      new_cigar = prev_hit->cigar();

      int new_left_back_len = new_cigar.back().length;
      new_left_back_len += dist_to_left;

      vector<CigarOp> new_right_cig = curr_hit->cigar();
      int new_right_front_len = new_right_cig.front().length;
      new_right_front_len -= dist_to_right;
      if (new_left_back_len > 0)
        new_cigar.back().length = new_left_back_len;
      else
        new_cigar.pop_back();

      /*
       * FIXME, currently just differentiating between a deletion and a
       * reference skip based on length. However, would probably be better
       * to denote the difference explicitly, this would allow the user
       * to supply their own (very large) deletions
       */
      if ((lb->right - lb->left - 1) <= max_deletion_length)
        {
          new_cigar.push_back(CigarOp(DEL, lb->right - lb->left - 1));
          antisense_closure = prev_hit->is_spliced() ? prev_hit->antisense_splice() : curr_hit->antisense_splice();
        }
      else
        {
          new_cigar.push_back(CigarOp(REF_SKIP, lb->right - lb->left - 1));
          antisense_closure = lb->antisense;
        }

      new_right_cig.front().length = new_right_front_len;
      size_t c = new_right_front_len > 0 ? 0 : 1;
      for (; c < new_right_cig.size(); ++c)
        new_cigar.push_back(new_right_cig[c]);

      mismatch = new_diff_mismatches;
      found_closure = true;
    }
        ++lb;
      }

    if (!found_closure)
      {
        return BowtieHit();
      }
  }

      if (found_closure)
  {
    bool end = false;
    BowtieHit merged_hit(reference_id,
			 reference_id,
             insert_id,
             new_left,
             new_cigar,
             antisense,
             antisense_closure,
             prev_hit->edit_dist() + curr_hit->edit_dist() + mismatch,
             prev_hit->splice_mms() + curr_hit->splice_mms(),
             end);

    if (curr_seg_index > 1)
      merged_hit.seq(seq.substr(first_seg_length + (curr_seg_index - 1) * segment_length, 2 * segment_length));
    else
      merged_hit.seq(seq.substr(0, first_seg_length + segment_length));

    prev_hit = hit_chain.erase(prev_hit, ++curr_hit);
    /*
     * prev_hit now points PAST the last element removed
     */
    prev_hit = hit_chain.insert(prev_hit, merged_hit);
    /*
     * merged_hit has been inserted before the old position of
     * prev_hit. New location of prev_hit is merged_hit
     */
    curr_hit = prev_hit;
    ++curr_hit;
    ++curr_seg_index;
    continue;
  }

      ++prev_hit;
      ++curr_hit;
      ++curr_seg_index;
    }

  bool saw_antisense_splice = false;
  bool saw_sense_splice = false;
  vector<CigarOp> long_cigar;
  int num_mismatches = 0;
  int num_splice_mms = 0;
  for (list<BowtieHit>::iterator s = hit_chain.begin(); s != hit_chain.end(); ++s)
    {
      num_mismatches += s->edit_dist();
      num_splice_mms += s->splice_mms();

      /*
       * Check whether the sequence contains any reference skips. Previously,
       * this was just a check to see whether the sequence was contiguous; however
       * we don't want to count an indel event as a splice
       */
      bool containsSplice = s->is_spliced();
      if (containsSplice)
  {
    if (s->antisense_splice())
      {
        if (saw_sense_splice)
    return BowtieHit();
        saw_antisense_splice = true;
      }
    else
      {
        if (saw_antisense_splice)
    return BowtieHit();
        saw_sense_splice = true;
      }
  }
      const vector<CigarOp>& cigar = s->cigar();
      if (long_cigar.empty())
  {
    long_cigar = cigar;
  }
      else
  {
    CigarOp& last = long_cigar.back();
    /*
     * If necessary, merge the back and front
     * cigar operations
     */
    if(last.opcode == cigar[0].opcode){
      last.length += cigar[0].length;
      for (size_t b = 1; b < cigar.size(); ++b)
        {
    long_cigar.push_back(cigar[b]);
        }
    }else{
      for(size_t b = 0; b < cigar.size(); ++b)
        {
    long_cigar.push_back(cigar[b]);
        }
    }
  }
    }

  bool end = false;
  BowtieHit new_hit(reference_id,
		    reference_id,
        insert_id,
        left,
        long_cigar,
        antisense,
        saw_antisense_splice,
        num_mismatches,
        num_splice_mms,
        end);

  new_hit.seq(seq);
  new_hit.qual(qual);

  int new_read_len = new_hit.read_len();
  if (new_read_len != old_read_length || !new_hit.check_editdist_consistency(rt))
    {
      fprintf(stderr, "Warning: malformed closure\n");
      return BowtieHit();
    }

  return new_hit;
}

BowtieHit merge_chain(RefSequenceTable& rt,
		      const string& read_seq,
		      const string& read_quals,
		      std::set<Junction>& possible_juncs,
		      std::set<Insertion>& possible_insertions,
		      std::set<Fusion>& possible_fusions,
		      list<BowtieHit>& hit_chain,
		      int fusion_dir = FUSION_NOTHING)
{
  bool antisense = hit_chain.front().antisense_align();
  uint64_t insert_id = hit_chain.front().insert_id();
  
  const int left = hit_chain.front().left();	
  
  list<BowtieHit>::iterator prev_hit = hit_chain.begin();
  list<BowtieHit>::iterator curr_hit = ++(hit_chain.begin());
  
  string seq;
  string qual;
  int old_read_length = 0;
  int first_seg_length = hit_chain.front().seq().length();
  for (list<BowtieHit>::iterator i = hit_chain.begin(); i != hit_chain.end(); ++i)
    {
      seq += i->seq();
      qual += i->qual();
      old_read_length += i->read_len();
    }

  string rev_read_seq, rev_read_quals;
  if (color && antisense)
    {
      rev_read_seq = read_seq;
      reverse(rev_read_seq.begin() + 1, rev_read_seq.end());
      
      rev_read_quals = read_quals;
      reverse(rev_read_quals.begin(), rev_read_quals.end());
    }

  size_t num_fusions = prev_hit->fusion_opcode() == FUSION_NOTHING ? 0 : 1;
  bool fusion_passed = false;
  while (curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
    {
      if (prev_hit->ref_id() != prev_hit->ref_id2() || prev_hit->ref_id2() != curr_hit->ref_id())
	  fusion_passed = true;

      if (prev_hit->ref_id2() != curr_hit->ref_id())
	++num_fusions;
      
      if (curr_hit->fusion_opcode() != FUSION_NOTHING)
	++num_fusions;

      if (prev_hit->ref_id2() == curr_hit->ref_id())
	{
	  bool reversed = false;
	  if ((fusion_dir == FUSION_FR && fusion_passed) || (fusion_dir == FUSION_RF && !fusion_passed))
	    reversed = true;
	  
	  /*
	   * Note that the gap size may be negative, since left() and right() return
	   * signed integers, this will be OK.
	   */
	  int gap;
	  if (reversed)
	    gap = prev_hit->right() - curr_hit->left();
	  else
	    gap = curr_hit->left() - prev_hit->right();

	  // daehwan
	  if (bDebug)
	    {
	      cout << "prev: " << prev_hit->ref_id() << ":" << prev_hit->left() << ":" << (prev_hit->antisense_align() ? "-" : "+")
		   << "\t" << prev_hit->ref_id2() << ":" << prev_hit->right() << ":" << (prev_hit->antisense_align2() ? "-" : "+") << endl
		   << "curr: " << curr_hit->ref_id() << ":" << curr_hit->left() << ":" << (curr_hit->antisense_align() ? "-" : "+")
		   << "\t" << curr_hit->ref_id2() << ":" << curr_hit->right() << ":" << (curr_hit->antisense_align2() ? "-" : "+") << endl
		   << "gap: " << gap << endl;
	    }
	  
	  if (gap < -(int)max_insertion_length ||
	      (gap > (int)max_deletion_length && 
	       (gap < min_report_intron_length || gap > min(max_report_intron_length, (int)fusion_min_dist))))
	    {
	      fusion_passed = true;
	      ++num_fusions;
	    }
	}

      if (num_fusions >= 2)
	  return BowtieHit();
      
      ++prev_hit;
      ++curr_hit;
    }
  
  prev_hit = hit_chain.begin();
  curr_hit = ++(hit_chain.begin());

  // daehwan
  if (bDebug)
    {
      cout << "daehwan - test" << endl;
    }
  
  int curr_seg_index = 1;
  fusion_passed = false;
  while (curr_hit != hit_chain.end() && prev_hit != hit_chain.end())
    {
      antisense = prev_hit->antisense_align();

      // daehwan
      if (bDebug)
	{
	  cout << "daehwan - start - stitch" << endl;
	  cout << "prev right: " << prev_hit->right() << endl;
	  cout << "curr left: " << curr_hit->left() << endl;
	  cout << "prev back: " << prev_hit->cigar().back().opcode << endl;
	  cout << "curr front: " << curr_hit->cigar().front().opcode << endl;
	  cout << "prev refs: " << prev_hit->ref_id() << "-" << prev_hit->ref_id2() << endl;
	  cout << "curr refs: " << curr_hit->ref_id() << "-" << curr_hit->ref_id2() << endl;
	}

      if (prev_hit->fusion_opcode() != FUSION_NOTHING || prev_hit->ref_id2() != curr_hit->ref_id())
	fusion_passed = true;
      
      /*
       * This code is assuming that the cigar strings end and start with a match
       * While segment alignments can actually end with a junction, insertion or deletion, the hope
       * is that in those cases, the right and left ends of the alignments will correctly
       * line up, so we won't get to this bit of code
       */
      if (!(prev_hit->cigar().back().opcode == MATCH || curr_hit->cigar().front().opcode == MATCH ||
	    prev_hit->cigar().back().opcode == mATCH || curr_hit->cigar().front().opcode == mATCH))
	{
	  return BowtieHit();
	}

      // daehwan
      if (bDebug)
	{
	  cout << "daehwan - pass - enough matched bases" << endl;
	}
      
      if (prev_hit->is_spliced() && curr_hit->is_spliced() && prev_hit->antisense_splice() != curr_hit->antisense_splice())
	{
	  /*
	   * There is no way that we can splice these together into a valid
	   * alignment
	   */
	  return BowtieHit();
	}

      bool found_closure = false;
      /*
       * Note that antisense_splice is the same for prev and curr
       */
      bool antisense_closure = prev_hit->is_spliced() ? prev_hit->antisense_splice() : curr_hit->antisense_splice();
      vector<CigarOp> new_cigar;
      int new_left = -1;
      int mismatch = 0;
      
      /*
       * Take the length of matched bases in each segment into consideration for closures,
       * this can be a problem for reads of variable lengths.
       */
      int prev_right_end_match_length = prev_hit->cigar().back().length;
      int curr_left_end_match_length = curr_hit->cigar().front().length;

      bool check_fusion = prev_hit->ref_id2() != curr_hit->ref_id();

      if (prev_hit->ref_id2() == curr_hit->ref_id())
	{
	  // daehwan
	  if (bDebug)
	    {
	      cout << "daehwan - start - junction or insertion" << endl;
	      cout << "prev right: " << prev_hit->right() << endl;
	      cout << "curr left: " << curr_hit->left() << endl;
	      cout << "prev refs: " << prev_hit->ref_id() << "-" << prev_hit->ref_id2() << endl;
	      cout << "curr refs: " << curr_hit->ref_id() << "-" << curr_hit->ref_id2() << endl;
	    }

	  bool reversed = false;
	  if ((fusion_dir == FUSION_FR && fusion_passed) || (fusion_dir == FUSION_RF && !fusion_passed))
	    reversed = true;

	  uint32_t reference_id = prev_hit->ref_id2();
	  RefSequenceTable::Sequence* ref_str = rt.get_seq(reference_id);
	  
	  int left_boundary, right_boundary;
	  if (reversed)
	    {
	      left_boundary = curr_hit->left() - 4;
	      right_boundary = prev_hit->right() + 4;
	    }
	  else
	    {
	      left_boundary = prev_hit->right() - 4;
	      right_boundary = curr_hit->left() + 4;
	    }

	  int dist_btw_two;
	  if (reversed)
	    dist_btw_two = prev_hit->right() - curr_hit->left();
	  else
	    dist_btw_two = curr_hit->left() - prev_hit->right();

	  if (dist_btw_two < 0 && dist_btw_two >= -max_insertion_length && prev_hit->antisense_align2() == curr_hit->antisense_align())
	    {
	      std::set<Insertion>::iterator lb, ub;
	      
	      /*
	       * Create a dummy sequence to represent the maximum possible insertion
	       */
	      std::string maxInsertedSequence = "";
	      maxInsertedSequence.resize(max_insertion_length,'A');

	      lb = possible_insertions.upper_bound(Insertion(reference_id, left_boundary, ""));
	      ub = possible_insertions.upper_bound(Insertion(reference_id, right_boundary, maxInsertedSequence));	
	      
	      int reference_mismatch = 0;
	      
	      while (lb != ub && lb != possible_insertions.end())
		{
		  /*
		   * In the following code, we will check to make sure that the segments have the proper
		   * separation and sequence for the insertions, and generate the appropriate merged bowtie hit
		   * In general, reads with insertions must match the inserted sequence exactly.
		   */
		  if (((int)lb->sequence.size()) == (reversed ? curr_hit->left() - prev_hit->right() : prev_hit->right() - curr_hit->left()))
		    {
		      /*
		       * Check we have enough matched bases on prev or curr segment.
		       */
		      int insert_to_prev_right, curr_left_to_insert;
		      if (reversed)
			{
			  insert_to_prev_right = lb->left - prev_hit->right();
			  curr_left_to_insert = curr_hit->left() - lb->left;
			}
		      else
			{
			  insert_to_prev_right = prev_hit->right() - lb->left - 1;
			  curr_left_to_insert = lb->left - curr_hit->left() + 1;
			}

		      if (insert_to_prev_right > prev_right_end_match_length || curr_left_to_insert > curr_left_end_match_length)
			{
			  ++lb;
			  continue;
			}

		      // daehwan
		      if (bDebug)
			{
			  cout << "insert_to_prev_right: " << insert_to_prev_right << endl;
			  cout << "curr_left_to_insert: " << curr_left_to_insert << endl;
			  cout << "curr_seg_index: " << curr_seg_index << endl;
			}
		      
		      /*
		       * Keep track of how many mismatches were made to the genome in the region
		       * where we should actually be matching against the insertion
		       */
		      int this_reference_mismatch = 0;
		      int insertion_mismatch = 0;
		      int insertion_len = lb->sequence.length();
		      seqan::Dna5String insertionSequence = seqan::Dna5String(lb->sequence);
		      if (reversed)
			{
			  seqan::reverseComplement(insertionSequence);
			}
		      
		      /*
		       * First check to see if we need to adjust number of observed errors for the left (prev)
		       * hit. This is only the case if this segment runs into the insertion. To be consistent
		       * with bwt_map.cpp, we will not allow a segment to have errors in the insertion region
		       */
		      string colorSegmentSequence_prev;
		      if (insert_to_prev_right > 0)
			{
			  seqan::Dna5String referenceSequence, oldSegmentSequence;

			  if (reversed)
			    {
			      referenceSequence = seqan::infix(*ref_str, prev_hit->right() + 1, lb->left + 1);
			      seqan::reverseComplement(referenceSequence);

			      string temp;

			      // daehwan
			      if (bDebug)
				{
				  cout << "reversed: " << read_seq.length() << " " << read_seq << endl;
				}
			      
			      temp = read_seq.substr(curr_seg_index * segment_length - insert_to_prev_right, insert_to_prev_right);
			      oldSegmentSequence = seqan::Dna5String(temp);
			    }
			  else
			    {
			      referenceSequence = seqan::infix(*ref_str, lb->left + 1, prev_hit->right());

			      // daehwan
			      if (bDebug)
				{
				  cout << "non-reversed: " << prev_hit->seq() << endl;
				}

			      oldSegmentSequence = seqan::Dna5String(prev_hit->seq().substr(prev_hit->seq().length() - insert_to_prev_right));
			    }
			  
			  if (color)
			    {
			      string color;
			      
			      if (antisense)
				color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - insert_to_prev_right - 2, insert_to_prev_right + 1);
			      else
				color = read_seq.substr(curr_seg_index * segment_length - insert_to_prev_right, insert_to_prev_right + 1);
			      
			      color[0] = prev_hit->seq()[segment_length - insert_to_prev_right - 1];
			      colorSegmentSequence_prev = convert_color_to_bp(color);
			    }
			  
			  const seqan::Dna5String newSegmentSequence = color ? colorSegmentSequence_prev : oldSegmentSequence;

			  // daehwan
			  if (bDebug)
			    {
			      cout << "ref: " << referenceSequence << endl;
			      cout << "old: " << oldSegmentSequence << endl;
			      cout << "ins: " << insertionSequence << endl;
			    }
			  
			  /*
			   * Scan right in the read until we run out of read
			   */
			  for (int read_index = 0; read_index < insert_to_prev_right; ++read_index)
			    { 
			      /*  
			       * Any mismatch to the insertion is a failure
			       */
			      if (referenceSequence[read_index] == 'N' || referenceSequence[read_index] != oldSegmentSequence[read_index])
				{
				  ++this_reference_mismatch;
				}
			      
			      if (read_index < insertion_len)
				{
				  if (insertionSequence[read_index] == 'N' || insertionSequence[read_index] != newSegmentSequence[read_index])
				    {
				      ++insertion_mismatch;
				      break;
				    }
				}
			      else
				{
				  if (referenceSequence[read_index - insertion_len] == 'N' ||
				      referenceSequence[read_index - insertion_len] != newSegmentSequence[read_index])
				    {
				      --this_reference_mismatch;
				    }
				}
			    }
			}
		      
		      string colorSegmentSequence_curr;
		      if (curr_left_to_insert > 0)
			{
			  seqan::Dna5String referenceSequence, oldSegmentSequence;
			  if (reversed)
			    {
			      referenceSequence = seqan::infix(*ref_str, lb->left + 1, curr_hit->left() + 1);
			      seqan::reverseComplement(referenceSequence);

			      string temp = read_seq.substr(curr_seg_index * segment_length, curr_left_to_insert);
			      oldSegmentSequence = seqan::Dna5String(temp);
			    }
			  else
			    {
			      referenceSequence = seqan::infix(*ref_str, curr_hit->left(), lb->left + 1);
			      oldSegmentSequence = seqan::Dna5String(curr_hit->seq().substr(0, curr_left_to_insert));
			    }
			  
			  if (color)
			    {
			      string color;
			      if (antisense)
				color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length), curr_left_to_insert);
			      else
				color = read_seq.substr(curr_seg_index * segment_length + 2, curr_left_to_insert);
			      
			      color.push_back(curr_hit->seq()[curr_left_to_insert]);
			      reverse(color.begin(), color.end());
			      string bp = convert_color_to_bp(color);
			      reverse(bp.begin(), bp.end());
			      colorSegmentSequence_curr = bp;
			    }
			  
			  const seqan::Dna5String newSegmentSequence = color ? colorSegmentSequence_curr : oldSegmentSequence;

			  // daehwan
			  if (bDebug)
			    {
			      cout << "ref: " << referenceSequence << endl;
			      cout << "old: " << oldSegmentSequence << endl;
			      cout << "ins: " << insertionSequence << endl;
			    }
		      
			  /*
			   * Scan left in the read until
			   * We ran out of read sequence (insertion extends past segment)
			   */
			  for (int read_index = 0; read_index < curr_left_to_insert; ++read_index)
			    {
			      int segmentPosition = curr_left_to_insert - read_index - 1;
			      int insertionPosition = insertion_len - read_index - 1;
			      
			      if (referenceSequence[segmentPosition] == 'N' || (referenceSequence[segmentPosition] != oldSegmentSequence[segmentPosition]))
				{
				  ++this_reference_mismatch;
				}
			      
			      if (read_index < insertion_len)
				{
				  if (insertionSequence[insertionPosition] == 'N' || (insertionSequence[insertionPosition] != newSegmentSequence[segmentPosition]))
				    {
				      ++insertion_mismatch;
				      break;
				    }
				}
			      else
				{
				  if (referenceSequence[segmentPosition + insertion_len] == 'N' ||
				      (referenceSequence[segmentPosition + insertion_len] != newSegmentSequence[segmentPosition]))
				    {
				      --this_reference_mismatch;
				    }
				}
			    }
			}
		      
		      if (found_closure)
			{
			  // fprintf(stderr, "Warning: multiple closures found for insertion read # %d\n", (int)insert_id);
			  return BowtieHit();
			}

		      if (insertion_mismatch == 0)
			{
			  reference_mismatch = this_reference_mismatch;
			  mismatch = -reference_mismatch;
			  found_closure = true;
			  new_left = prev_hit->left();
			  new_cigar = prev_hit->cigar();
			  
			  /*
			   * Need to make a new insert operation between the two match character that begin
			   * and end the intersection of these two junction. Note that we necessarily assume
			   * that this insertion can't span beyond the boundaries of these reads. That should
			   * probably be better enforced somewhere
			   */

			  new_cigar.back().length -= insert_to_prev_right;
			  
			  if (new_cigar.back().length <= 0)
			    new_cigar.pop_back();

			  if (reversed)
			    new_cigar.push_back(CigarOp(iNS, lb->sequence.size()));
			  else
			    new_cigar.push_back(CigarOp(INS, lb->sequence.size()));
			  
			  vector<CigarOp> new_right_cigar = curr_hit->cigar();
			  new_right_cigar.front().length += (insert_to_prev_right - lb->sequence.size());
			  
			  /*
			   * Finish stitching together the new cigar string
			   */
			  size_t c = new_right_cigar.front().length > 0 ? 0 : 1;
			  for (; c < new_right_cigar.size(); ++c)
			    {
			      new_cigar.push_back(new_right_cigar[c]);
			    }
			  
			  if (color)
			    {
			      if (insert_to_prev_right > 0)
				seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length - insert_to_prev_right, insert_to_prev_right, colorSegmentSequence_prev);
			      
			      if (curr_left_to_insert > 0)
				seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length, curr_left_to_insert, colorSegmentSequence_curr);
			    }
			}
		    }
		  ++lb;
		}
	      
	      if (!found_closure)
		{
		  return BowtieHit();
		}
	    }

	  /*
	   * Stitch segments together using juctions or deletions if necessary.
	   */
	  else if (dist_btw_two > 0 && dist_btw_two <= max_report_intron_length && prev_hit->antisense_align2() == curr_hit->antisense_align())
	    {
	      std::set<Junction>::iterator lb, ub;

	      // daehwan
	      if (bDebug)
		{
		  cout << "junction" << endl;
		  cout << "min: " << left_boundary << "-" << right_boundary - 8 << endl;
		  cout << "max: " << left_boundary + 8 << "-" << right_boundary << endl;
		}

	      lb = possible_juncs.upper_bound(Junction(reference_id, left_boundary, right_boundary - 8, true));
	      ub = possible_juncs.lower_bound(Junction(reference_id, left_boundary + 8, right_boundary, false));
	      
	      int new_diff_mismatches = 0xff;
	      while (lb != ub && lb != possible_juncs.end())
		{
		  int dist_to_left, dist_to_right;
		  if (reversed)
		    {
		      dist_to_left = lb->left - curr_hit->left();
		      dist_to_right = lb->right - prev_hit->right() - 1;
		    }
		  else
		    {
		      dist_to_left = lb->left - prev_hit->right() + 1;
		      dist_to_right = lb->right - curr_hit->left();
		    }
		  
		  if (abs(dist_to_left) <= 4 && abs(dist_to_right) <= 4 && dist_to_left == dist_to_right)
		    {
		      /*
		       * Check we have enough matched bases on prev or curr segment.
		       */
		      if ((reversed && (dist_to_left > prev_right_end_match_length || -dist_to_left > curr_left_end_match_length)) ||
			  (!reversed && (dist_to_left > curr_left_end_match_length || -dist_to_left > prev_right_end_match_length)))
			{
			  ++lb;
			  continue;
			}

		      // daehwan
		      if (bDebug)
			{
			  cout << "candidate junction: " << endl;
			  cout << "coords: " << lb->left << "-" << lb->right << endl;
			  cout << "dist to left: " << dist_to_left << endl;
			}
		      
		      Dna5String new_cmp_str, old_cmp_str;
		      int new_mismatch = 0, old_mismatch = 0;
		      string new_patch_str; // this is for colorspace reads
		      if (dist_to_left > 0)
			{
			  if (reversed)
			    {
			      new_cmp_str = seqan::infix(*ref_str, curr_hit->left() + 1, lb->left + 1);
			      seqan::reverseComplement(new_cmp_str);
			      
			      old_cmp_str = seqan::infix(*ref_str, prev_hit->right() + 1, lb->right);
			      seqan::reverseComplement(old_cmp_str);
			    }
			  else
			    {
			      new_cmp_str = seqan::infix(*ref_str, prev_hit->right(), lb->left + 1);
			      old_cmp_str = seqan::infix(*ref_str, curr_hit->left(), lb->right);
			    }
			  
			  string new_seq;
			  if (color)
			    {
			      string ref = DnaString_to_string(seqan::infix(*ref_str, prev_hit->right() - 1, lb->left + 1));
			      
			      string color, qual;
			      if (antisense)
				{
				  color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - 1, dist_to_left);
				  qual = rev_read_quals.substr(rev_read_quals.length() - (curr_seg_index * segment_length) - 1, dist_to_left);
				}
			      else
				{
				  color = read_seq.substr(1 + curr_seg_index * segment_length, dist_to_left);
				  qual = read_quals.substr(curr_seg_index * segment_length, dist_to_left);
				}
			      
			      BWA_decode(color, qual, ref, new_seq);
			      new_seq = new_seq.substr(1);
			    }

			  string curr_hit_seq;
			  if (reversed)
			    curr_hit_seq = read_seq.substr(curr_seg_index * segment_length - dist_to_left, dist_to_left);
			  else
			    curr_hit_seq = curr_hit->seq();
			  
			  const string& curr_old_seq = curr_hit_seq;
			  const string& curr_seq = color ? new_seq : curr_hit_seq;
			  for (int i = 0; i < dist_to_left; ++i)
			    {
			      if (curr_seq[i] != new_cmp_str[i])
				++new_mismatch;
			      
			      if (curr_old_seq[i] != old_cmp_str[i])
				++old_mismatch;
			    }
			  
			  if (color)
			    new_patch_str = curr_seq.substr(0, dist_to_left);
			}
		      else if (dist_to_left < 0)
			{
			  if (reversed)
			    {
			      new_cmp_str = seqan::infix(*ref_str, lb->right, prev_hit->right() + 1);
			      seqan::reverseComplement(new_cmp_str);

			      old_cmp_str = seqan::infix(*ref_str, lb->left + 1, curr_hit->left() + 1);
			      seqan::reverseComplement(old_cmp_str);
			    }
			  else
			    {
			      new_cmp_str = seqan::infix(*ref_str, lb->right, curr_hit->left());
			      old_cmp_str = seqan::infix(*ref_str, lb->left + 1, prev_hit->right());
			    }

			  size_t abs_dist = -dist_to_left;
			  string new_seq;
			  if (color)
			    {
			      string ref = DnaString_to_string(seqan::infix(*ref_str, lb->left, lb->left + 1));
			      ref += DnaString_to_string(seqan::infix(*ref_str, lb->right, curr_hit->left()));
			      
			      string color, qual;
			      if (antisense)
				{
				  color = rev_read_seq.substr(rev_read_seq.length() - (curr_seg_index * segment_length) - 1 - abs_dist, abs_dist);
				  qual = rev_read_quals.substr(rev_read_quals.length() - (curr_seg_index * segment_length) - 1 - abs_dist, abs_dist);
				}
			      else
				{
				  color = read_seq.substr(1 + curr_seg_index * segment_length - abs_dist, abs_dist);
				  qual = read_quals.substr(curr_seg_index * segment_length - abs_dist, abs_dist);
				}
			      
			      BWA_decode(color, qual, ref, new_seq);
			      new_seq = new_seq.substr(1);
			    }

			  string prev_hit_seq;
			  if (reversed)
			    prev_hit_seq = read_seq.substr(curr_seg_index * segment_length, abs_dist);
			  else
			    prev_hit_seq = prev_hit->seq();

			  // daehwan
			  if (bDebug)
			    {
			      cout << "reverse: " << (int)reversed << endl;
			      cout << "new cmp str: " << new_cmp_str << endl;
			      cout << "old cmp str: " << old_cmp_str << endl;
			      cout << "hit seq: " << prev_hit_seq << endl;
			      cout << "curr seq: " << curr_hit->seq() << endl;

			      cout << read_seq
				   << endl;
			      cout << read_seq.substr(first_seg_length + (curr_seg_index - 1) * segment_length, segment_length)
				   << endl;
			    }
			  
			  const string& prev_old_seq = prev_hit_seq;
			  size_t prev_old_seq_len = prev_old_seq.length();
			  const string& prev_seq = color ? new_seq : prev_hit_seq;
			  size_t prev_seq_len = prev_seq.length();
			  for (size_t i = 0; i < abs_dist; ++i)
			    {
			      if (prev_seq[prev_seq_len - (abs_dist - i)] != new_cmp_str[i])
				++new_mismatch;
			      if (prev_old_seq[prev_old_seq_len - (abs_dist - i)] != old_cmp_str[i])
				++old_mismatch;
			    }

			  if (color)
			    new_patch_str = prev_seq.substr(prev_seq_len - abs_dist, abs_dist);
			}
		      
		      int temp_diff_mismatches = new_mismatch - old_mismatch;

		      // daehwan
		      if (bDebug)
			{
			  cout << "new mismatch: " << new_mismatch << endl;
			  cout << "old mismatch: " << old_mismatch << endl;
			  cout << "new_diff_mismatch: " << new_diff_mismatches << endl;
			  cout << "temp mismatch: " << temp_diff_mismatches << endl;
			}
		      
		      if (temp_diff_mismatches >= new_diff_mismatches || new_mismatch >= 2)
			{
			  ++lb;
			  continue;
			}

		      if (color)
			{
			  /*
			   * We need to recover the origianl sequence.
			   */
			  if (found_closure)
			    {
			      seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length - 4, 8,
					  prev_hit->seq().substr(prev_hit->seq().length() - 4) + curr_hit->seq().substr(0, 4));
			    }
			  
			  if (dist_to_left > 0)
			    seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length, dist_to_left, new_patch_str);
			  else if (dist_to_left < 0)
			    seq.replace(first_seg_length + (curr_seg_index - 1) * segment_length + dist_to_left, -dist_to_left, new_patch_str);
			}
		      
		      new_diff_mismatches = temp_diff_mismatches;
		      
		      new_left = prev_hit->left();
		      new_cigar = prev_hit->cigar();
		      
		      int new_left_back_len = new_cigar.back().length;
		      if (reversed)
			new_left_back_len -= dist_to_left;
		      else
			new_left_back_len += dist_to_left;
		      
		      vector<CigarOp> new_right_cig = curr_hit->cigar();
		      int new_right_front_len = new_right_cig.front().length;
		      if (reversed)
			new_right_front_len += dist_to_right;
		      else
			new_right_front_len -= dist_to_right;
		      if (new_left_back_len > 0)
			new_cigar.back().length = new_left_back_len;					  
		      else
			new_cigar.pop_back();
		      
		      /*
		       * FIXME, currently just differentiating between a deletion and a 
		       * reference skip based on length. However, would probably be better
		       * to denote the difference explicitly, this would allow the user
		       * to supply their own (very large) deletions
		       */
		      if ((lb->right - lb->left - 1) <= max_deletion_length)
			{
			  if (reversed)
			    new_cigar.push_back(CigarOp(dEL, lb->right - lb->left - 1));
			  else
			    new_cigar.push_back(CigarOp(DEL, lb->right - lb->left - 1));
			  antisense_closure = prev_hit->is_spliced() ? prev_hit->antisense_splice() : curr_hit->antisense_splice();
			}
		      else
			{
			  if (reversed)
			    new_cigar.push_back(CigarOp(rEF_SKIP, lb->right - lb->left - 1));
			  else
			    new_cigar.push_back(CigarOp(REF_SKIP, lb->right - lb->left - 1));
			  antisense_closure = lb->antisense;
			}
		      
		      new_right_cig.front().length = new_right_front_len;
		      size_t c = new_right_front_len > 0 ? 0 : 1;
		      for (; c < new_right_cig.size(); ++c)
			new_cigar.push_back(new_right_cig[c]);
		      
		      mismatch = new_diff_mismatches;
		      found_closure = true;
		    }
		  ++lb;
		}

	      if (!found_closure)
		{
		  return BowtieHit();
		}
	    }
	  else if (!(dist_btw_two == 0 && prev_hit->antisense_align2() == curr_hit->antisense_align()))
	     check_fusion = true;
	}

      if (check_fusion)
	{
	  std::set<Fusion>::iterator lb, ub;

	  uint32_t ref_id1 = prev_hit->ref_id2();
	  uint32_t ref_id2 = curr_hit->ref_id();
	  uint32_t left = prev_hit->right() - 4;
	  uint32_t right = curr_hit->left() - 4;

	  // daehwan
	  if (bDebug)
	    {
	      cout << "daehwan - start - fusion" << endl
		   << "ref_id1: " << ref_id1 << endl
		   << "ref_id2: " << ref_id2 << endl
		   << "left: " << left << endl
		   << "right: " << right << endl
		   << "dir: " << fusion_dir << endl;
	    }

	  bool reversed = false;
	  if (fusion_dir != FUSION_FF &&
	      (ref_id2 < ref_id1 || (ref_id1 == ref_id2 && left > right)))
	    {
	      reversed = true;
	      
	      uint32_t temp = ref_id1;
	      ref_id1 = ref_id2;
	      ref_id2 = temp;
	      
	      temp = left;
	      left = right;
	      right = temp;
	    }
	      
	  lb = possible_fusions.upper_bound(Fusion(ref_id1, ref_id2, left, right));
	  ub = possible_fusions.lower_bound(Fusion(ref_id1, ref_id2, left + 8, right + 8));

	  RefSequenceTable::Sequence* ref_str = rt.get_seq(prev_hit->ref_id2());
	  RefSequenceTable::Sequence* ref_str2 = rt.get_seq(curr_hit->ref_id());
	  
	  int new_diff_mismatches = 0xff;
	  while (lb != ub && lb != possible_fusions.end())
	    {
	      int lb_left = lb->left;
	      int lb_right = lb->right;
	      
	      if (reversed)
		{
		  lb_left = lb->right;
		  lb_right = lb->left;
		}

	      int dist_to_left, dist_to_right;
	      if (fusion_dir == FUSION_RF)
		dist_to_left = prev_hit->right() - lb_left + 1;
	      else
		dist_to_left = lb_left - prev_hit->right() + 1;
	      
	      if (fusion_dir == FUSION_FR)
		dist_to_right = curr_hit->left() - lb_right;
	      else
		dist_to_right = lb_right - curr_hit->left();

	      // daehwan
	      if (bDebug)
		{
		  cout << "daehwan - fusion gap" << endl;
		  cout << "dist left: " << dist_to_left << endl;
		  cout << "dist right: " << dist_to_right << endl;
		}
	      
	      if (abs(dist_to_left) <= 4 && abs(dist_to_right) <= 4 && dist_to_left == dist_to_right)
		{
		  /*
		   * Check we have enough matched bases on prev or curr segment.
		   */
		  if (dist_to_left > curr_left_end_match_length || -dist_to_left > prev_right_end_match_length)
		    {
		      ++lb;
		      continue;
		    }

		  Dna5String new_cmp_str, old_cmp_str;
		  int new_mismatch = 0, old_mismatch = 0;
		  string new_patch_str; // this is for colorspace reads
		  if (dist_to_left > 0)
		    {
		      if (fusion_dir == FUSION_RF)
			{
			  new_cmp_str = seqan::infix(*ref_str, lb_left, prev_hit->right() + 1);
			  seqan::reverseComplement(new_cmp_str);
			}
		      else
			new_cmp_str = seqan::infix(*ref_str, prev_hit->right(), lb_left + 1);

		      if (fusion_dir == FUSION_FR)
			{
			  old_cmp_str = seqan::infix(*ref_str2, lb_right + 1, curr_hit->left() + 1);
			  seqan::reverseComplement(old_cmp_str);
			}
		      else
			old_cmp_str = seqan::infix(*ref_str2, curr_hit->left(), lb_right);

		      // daehwan
		      if (bDebug)
			{
			  cout << "new str: " << new_cmp_str << endl;
			  cout << "old str: " << old_cmp_str << endl;
			  cout << "curr seq: " << curr_hit->seq() << endl;
			}

		      string curr_hit_seq;
		      if (fusion_dir == FUSION_FF || fusion_dir == FUSION_RR)
			curr_hit_seq = curr_hit->seq();
		      else
			curr_hit_seq = read_seq.substr(curr_seg_index * segment_length, segment_length);
		      
		      string new_seq;
		      const string& curr_old_seq = curr_hit_seq;
		      const string& curr_seq = color ? new_seq : curr_hit_seq;
		      for (int i = 0; i < dist_to_left; ++i)
			{
			  if (curr_seq[i] != new_cmp_str[i])
			    ++new_mismatch;
			  
			  if (curr_old_seq[i] != old_cmp_str[i])
			    ++old_mismatch;
			}
		    }
		  else if (dist_to_left < 0)
		    {
		      if (fusion_dir == FUSION_FR)
			{
			  new_cmp_str = seqan::infix(*ref_str2, curr_hit->left() + 1, lb_right + 1);
			  seqan::reverseComplement(new_cmp_str);
			}
		      else
			new_cmp_str = seqan::infix(*ref_str2, lb_right, curr_hit->left());

		      if (fusion_dir == FUSION_RF)
			{
			  old_cmp_str = seqan::infix(*ref_str, prev_hit->right() + 1, lb_left);
			  seqan::reverseComplement(old_cmp_str);
			}
		      else
			old_cmp_str = seqan::infix(*ref_str, lb_left + 1, prev_hit->right());

		      string prev_hit_seq;
		      if (fusion_dir == FUSION_FF || fusion_dir == FUSION_RR)
			prev_hit_seq = prev_hit->seq();
		      else
			prev_hit_seq = read_seq.substr((curr_seg_index - 1) * segment_length, segment_length);

		      size_t abs_dist = -dist_to_left;
		      string new_seq;
			  
		      const string& prev_old_seq = prev_hit_seq;
		      size_t prev_old_seq_len = prev_old_seq.length();
		      const string& prev_seq = color ? new_seq : prev_hit_seq;
		      size_t prev_seq_len = prev_seq.length();
		      for (size_t i = 0; i < abs_dist; ++i)
			{
			  if (prev_seq[prev_seq_len - (abs_dist - i)] != new_cmp_str[i])
			    ++new_mismatch;
			  if (prev_old_seq[prev_old_seq_len - (abs_dist - i)] != old_cmp_str[i])
			    ++old_mismatch;
			}
		    }

		  int temp_diff_mismatches = new_mismatch - old_mismatch;
		  if (temp_diff_mismatches >= new_diff_mismatches || new_mismatch >= 2)
		    {
		      ++lb;
		      continue;
		    }
		  new_diff_mismatches = temp_diff_mismatches;

		  new_left = prev_hit->left();
		  new_cigar = prev_hit->cigar();
		      
		  int new_left_back_len = new_cigar.back().length;
		  new_left_back_len += dist_to_left;
		  
		  vector<CigarOp> new_right_cig = curr_hit->cigar();
		  int new_right_front_len = new_right_cig.front().length;
		  new_right_front_len -= dist_to_right;
		  if (new_left_back_len > 0)
		    new_cigar.back().length = new_left_back_len;					  
		  else
		    new_cigar.pop_back();

		  new_cigar.push_back(CigarOp((CigarOpCode)fusion_dir, lb_right));
		  antisense_closure = prev_hit->is_spliced() ? prev_hit->antisense_splice() : curr_hit->antisense_splice();
		  
		  new_right_cig.front().length = new_right_front_len;
		  size_t c = new_right_front_len > 0 ? 0 : 1;
		  for (; c < new_right_cig.size(); ++c)
		    new_cigar.push_back(new_right_cig[c]);
		  
		  mismatch = new_diff_mismatches;
		  found_closure = true;
		  ++num_fusions;

		  // daehwan
		  if (bDebug)
		    {
		      cout << "daehwan - fusion gap - found" << endl;
		    }

		}
	      ++lb;
	    }

	  // daehwan
	  if (bDebug)
	    {
	      cout << "daehwan2 - end - fusion: " << (found_closure ? "found" : "not found") << endl;
	    }
	  
	  if (!found_closure)
	    {
	      return BowtieHit();
	    }
	}

      if (found_closure)
	{
	  bool end = false;
	  BowtieHit merged_hit(prev_hit->ref_id(),
			       curr_hit->ref_id2(),
			       insert_id,
			       new_left,
			       new_cigar,
			       antisense,
			       antisense_closure,
			       prev_hit->edit_dist() + curr_hit->edit_dist() + mismatch,
			       prev_hit->splice_mms() + curr_hit->splice_mms(),
			       end);

	  // daehwan - should fix this for SOLiD dataset
	  merged_hit.seq(prev_hit->seq() + curr_hit->seq());

	  // daehwan
	  if (bDebug)
	    {
	      cout << "fusing of " << merged_hit.left() << " and " << merged_hit.right() << endl;
	      cout << print_cigar(merged_hit.cigar()) << endl;
	      if (!merged_hit.check_editdist_consistency(rt, bDebug))
		{
		  prev_hit->check_editdist_consistency(rt, bDebug);
		  curr_hit->check_editdist_consistency(rt, bDebug);
		  cout << "btw " << print_cigar(prev_hit->cigar()) << " and " << print_cigar(curr_hit->cigar()) << endl;
		  cout << "this is a malformed hit" << endl;
		  exit(1);
		}
	    }

	  prev_hit = hit_chain.erase(prev_hit, ++curr_hit);
	  /*
	   * prev_hit now points PAST the last element removed
	   */
	  prev_hit = hit_chain.insert(prev_hit, merged_hit);
	  /*
	   * merged_hit has been inserted before the old position of
	   * prev_hit. New location of prev_hit is merged_hit
	   */
	  curr_hit = prev_hit;
	  ++curr_hit;
	  ++curr_seg_index;
	  continue;
	}

      // daehwan
      if (bDebug)
	{
	  cout << "daehwan - test 0.3" << endl;
	}
      
      ++prev_hit;
      ++curr_hit;
      ++curr_seg_index;
    }

  // daehwan
  if (bDebug)
    {
      cout << "daehwan - test2" << endl;
    }
  
  bool saw_antisense_splice = false;
  bool saw_sense_splice = false;
  vector<CigarOp> long_cigar;
  int num_mismatches = 0;
  int num_splice_mms = 0;
  for (list<BowtieHit>::iterator s = hit_chain.begin(); s != hit_chain.end(); ++s)
    {
      num_mismatches += s->edit_dist();
      num_splice_mms += s->splice_mms();

      /*
       * Check whether the sequence contains any reference skips. Previously,
       * this was just a check to see whether the sequence was contiguous; however
       * we don't want to count an indel event as a splice
       */
      bool containsSplice = s->is_spliced();
      if (containsSplice)
	{
	  if (s->antisense_splice())
	    {
	      if (saw_sense_splice)
		return BowtieHit();
	      saw_antisense_splice = true;
	    }
	  else
	    {
	      if (saw_antisense_splice)
		return BowtieHit();
	      saw_sense_splice = true;
	    }
	}
      const vector<CigarOp>& cigar = s->cigar();
      if (long_cigar.empty())
	{
	  long_cigar = cigar;
	}
      else
	{
	  CigarOp& last = long_cigar.back();
	  /*
	   * If necessary, merge the back and front
	   * cigar operations
	   */
	  if(last.opcode == cigar[0].opcode){
	    last.length += cigar[0].length;
	    for (size_t b = 1; b < cigar.size(); ++b)
	      {
		long_cigar.push_back(cigar[b]);
	      }
	  }else{
	    for(size_t b = 0; b < cigar.size(); ++b)
	      {
		long_cigar.push_back(cigar[b]);
	      }
	  }
	}
    }

  bool end = false;
  BowtieHit new_hit(hit_chain.front().ref_id(),
		    hit_chain.back().ref_id2(),
		    insert_id, 
		    left, 
		    long_cigar, 
		    antisense,
		    saw_antisense_splice,
		    num_mismatches,
		    num_splice_mms,
		    end);

  if (fusion_dir == FUSION_NOTHING || fusion_dir == FUSION_FF || fusion_dir == FUSION_RR)
    {
      new_hit.seq(seq);
      if (bowtie2)
	{
	  // for the time being, let's compare "seq" and "read_seq"
	  if (seq != read_seq)
	    {
	      string temp_qual = read_quals;
	      reverse(temp_qual.begin(), temp_qual.end());
	      new_hit.qual(temp_qual);
	    }
	  else
	    new_hit.qual(read_quals);
	}
      else
	new_hit.qual(qual);
    }
  else
    {
      new_hit.seq(read_seq);
      new_hit.qual(read_quals);
    }

  bool do_reverse = new_hit.ref_id() > new_hit.ref_id2();
  if (new_hit.ref_id() == new_hit.ref_id2())
    {
      vector<Fusion> fusions;
      bool auto_sort = false;
      fusions_from_spliced_hit(new_hit, fusions, auto_sort);
      if (fusions.size() > 0)
	{
	  const Fusion& fusion = fusions[0];
	  do_reverse = fusion.left > fusion.right;
	}
    }

  if (do_reverse)
    {
      new_hit = new_hit.reverse();
    }

  if (fusion_dir == FUSION_RF || fusion_dir == FUSION_RR)
    {
      new_hit.antisense_align(!new_hit.antisense_align());
    }

  // daehwan
  if (bDebug)
    {
      cout << "daehwan - test3" << endl;
      cout << new_hit.left() << " " << print_cigar(new_hit.cigar()) << endl;
      cout << new_hit.ref_id() << "-" << new_hit.ref_id2() << ": " << new_hit.fusion_opcode() << endl;
    }

  int new_read_len = new_hit.read_len();
  if (new_read_len != old_read_length || !new_hit.check_editdist_consistency(rt, bDebug))
    {
      // daehwan
      if (bDebug)
	{
	  cout << "Warning: " << new_hit.insert_id() << " malformed closure: " << print_cigar(new_hit.cigar()) << endl;
	  exit(1);
	}
      
      fprintf(stderr, "Warning: %d malformed closure\n", new_hit.insert_id());
      return BowtieHit();
    }

  return new_hit;
}


int multi_closure = 0;
int anchor_too_short = 0;
int gap_too_short = 0;

bool valid_hit(const BowtieHit& bh)
{
  if (bh.insert_id())
    {	
      /*	
       * validate the cigar chain - no gaps shorter than an intron, etc.
       * also,
       * 	-Don't start or end with an indel or refskip
       *	-Only a match operation is allowed is allowed
       *	 adjacent to an indel or refskip
       *      -Indels should confrom to length restrictions
       */
      const CigarOp* prevCig = &(bh.cigar()[0]); 
      const CigarOp* currCig = &(bh.cigar()[1]);
      for (size_t i = 1; i < bh.cigar().size(); ++i){
	currCig = &(bh.cigar()[i]);
	if(!(currCig->opcode == MATCH || currCig->opcode == mATCH) &&
	   !(prevCig->opcode == MATCH || prevCig->opcode == mATCH)){
	  return false;
	}
	if(currCig->opcode == INS || currCig->opcode == iNS){
	  if(currCig->length > max_insertion_length){
	    return false;
	  }
	}
	if(currCig->opcode == DEL || currCig->opcode == dEL){
	  if(currCig->length > max_deletion_length){
	    return false;
	  }
	}
	if(currCig->opcode == REF_SKIP || currCig->opcode == rEF_SKIP){
	  if(currCig->length < (uint64_t)min_report_intron_length){
	    gap_too_short++;
	    return false;
	  }
	}
	prevCig = currCig;
      }
      if (!(bh.cigar().front().opcode == MATCH || bh.cigar().front().opcode == mATCH) || 
	  !(bh.cigar().back().opcode == MATCH || bh.cigar().back().opcode == mATCH)/* ||
										      (int)bh.cigar().front().length < min_anchor_len||
										      (int)bh.cigar().back().length < min_anchor_len*/ )
	{
	  anchor_too_short++;
	  return false;
	}	
    }
  else
    {
      multi_closure++;
      return false;
    }
  
  return true;
}

void merge_segment_chain(RefSequenceTable& rt,
			 const string& read_seq,
			 const string& read_quals,
			 std::set<Junction>& possible_juncs,
			 std::set<Insertion>& possible_insertions,
			 std::set<Fusion>& possible_fusions,
			 vector<BowtieHit>& hits,
			 vector<BowtieHit>& merged_hits,
			 int fusion_dir = FUSION_NOTHING)
{
  if (hits.size() == 0)
    return;

  BowtieHit bh;
  if (hits.size() > 1)
    {
      list<BowtieHit> hit_chain;
      if (fusion_dir == FUSION_NOTHING || fusion_dir == FUSION_FF || fusion_dir == FUSION_RR)
	{
	  if (hits.front().antisense_align())
	    copy(hits.rbegin(), hits.rend(), back_inserter(hit_chain));
	  else
	    copy(hits.begin(), hits.end(), back_inserter(hit_chain));
	}
      else
	{
	  bool bSawFusion = false;
	  for (size_t i = 0; i < hits.size(); ++i)
	    {
	      bool pushed = false;
	      if (!bSawFusion)
		{
		  if (i > 0)
		    {
		      if (hits[i-1].ref_id() != hits[i].ref_id())
			bSawFusion = true;
		      else if(hits[i-1].antisense_align() != hits[i].antisense_align())
			bSawFusion = true;
		      else
			{
			  int dist = 0;
			  if (hits[i].antisense_align())
			    dist = hits[i-1].left() - hits[i].right();
			  else
			    dist = hits[i].left() - hits[i-1].right();

			  if (dist >= max_report_intron_length || dist < -(int)max_insertion_length)
			    bSawFusion = true;
			}
		    }
		}
	      
	      if (hits[i].fusion_opcode() == FUSION_NOTHING &&
		  ((fusion_dir == FUSION_FR && bSawFusion) || (fusion_dir == FUSION_RF && !bSawFusion)) &&
		  hits[i].left() < hits[i].right())
		{
		  hit_chain.push_back(hits[i].reverse());
		  pushed = true;
		}

	      if (i > 0 &&
		  hits[i].fusion_opcode() != FUSION_NOTHING &&
		  hits[i].ref_id() != hits[i-1].ref_id())
		{
		  hit_chain.push_back(hits[i].reverse());
		  pushed = true;
		}

	      if (!bSawFusion)
		{
		  if (hits[i].fusion_opcode() != FUSION_NOTHING)
		    bSawFusion = true;
		}

	      if (!pushed)
		hit_chain.push_back(hits[i]);
	    }
	}

      // todo: merge_chain_color needs to be merged into merge_chain fuction.
      if (color)
	bh = merge_chain_color(rt,
			 read_seq,
			 read_quals,
			 possible_juncs,
			 possible_insertions,
			       hit_chain);
      else
	bh = merge_chain(rt,
			 read_seq,
			 read_quals,
			 possible_juncs,
			 possible_insertions,
			 possible_fusions,
			 hit_chain,
			 fusion_dir);
    }
  else
    {
      bh = hits[0];
      bool do_reverse = bh.ref_id() > bh.ref_id2();
      if (bh.ref_id() == bh.ref_id2())
	{
	  vector<Fusion> fusions;
	  bool auto_sort = false;
	  fusions_from_spliced_hit(bh, fusions, auto_sort);
	  if (fusions.size() > 0)
	    {
	      const Fusion& fusion = fusions[0];
	      do_reverse = fusion.left > fusion.right;
	    }
	}
      
      if (do_reverse)
	bh = bh.reverse();
    }

  if (valid_hit(bh))
      merged_hits.push_back(bh);
}

bool dfs_seg_hits(RefSequenceTable& rt,
		  const string& read_seq,
		  const string& read_quals,
		  std::set<Junction>& possible_juncs,
		  std::set<Insertion>& possible_insertions,
		  std::set<Fusion>& possible_fusions, 
		  vector<HitsForRead>& seg_hits_for_read,
		  size_t curr,
		  vector<BowtieHit>& seg_hit_stack,
		  vector<BowtieHit>& joined_hits,
		  int fusion_dir = FUSION_NOTHING)
{
  assert (!seg_hit_stack.empty());
  bool join_success = false;

  if (curr < seg_hits_for_read.size())
    {
      for (size_t i = 0; i < seg_hits_for_read[curr].hits.size(); ++i)
	{
	  /*
	   * As we reverse segments depending on directions like FR or RF,
	   * it's necessary to recover the original segments.
	   */
	  BowtieHit bh = seg_hits_for_read[curr].hits[i];
	  BowtieHit bh_prev = seg_hit_stack.back();
	  
	  BowtieHit* prevHit = &bh_prev;
	  BowtieHit* currHit = &bh;

	  // daehwan
	  // if (bh.insert_id() == 11921733)
	  //  bDebug = true;

	  /*
	   * Each segment has at most one fusion by assumption,
	   */
	  bool prevHit_fused = prevHit->fusion_opcode() != FUSION_NOTHING;
	  bool currHit_fused = currHit->fusion_opcode() != FUSION_NOTHING;

	  /*
	   * Count the number of fusions on prev and curr segments,
	   * but this doesn't take into account the gap (if exists) between the two segments.
	   */
	  size_t num_fusions = prevHit_fused ? 1 : 0;
	  num_fusions += currHit_fused ? 1 : 0;

	  int dir = prevHit_fused ? prevHit->fusion_opcode() : currHit->fusion_opcode();

	  /*
	   * We don't allow reads that span more than two fusion points.
	   */
	  if (num_fusions >= 2)
	    continue;

	  if (fusion_dir != FUSION_NOTHING && currHit_fused)
	    continue;

	  if (fusion_dir == FUSION_FF || fusion_dir == FUSION_RR)
	    {
	      if ((currHit->antisense_align() && currHit->ref_id() != prevHit->ref_id()) ||
		  (!currHit->antisense_align() && currHit->ref_id() != prevHit->ref_id2()))
		continue;
	    }

	  if (bDebug)
	    {
	      cout << "daehwan - prev ref: " << prevHit->ref_id() << "-" << prevHit->ref_id2() << ": "
		   << print_cigar(prevHit->cigar()) << endl;
	      cout << "daehwan - prev sense: "
		   << (prevHit->antisense_align() ? "-" : "+")
		   << "\t"
		   << (prevHit->antisense_align2() ? "-" : "+")
		   << endl;
	      cout << "daehwan - prev coords: "
		   << prevHit->left()
		   << "\t"
		   << prevHit->right()
		   << endl;
	      
	      cout << "daehwan - curr ref: " << currHit->ref_id() << "-" << currHit->ref_id2() << ": "
		   << print_cigar(currHit->cigar()) << endl;
	      cout << "daehwan - curr sense: "
		   << (currHit->antisense_align() ? "-" : "+")
		   << "\t"
		   << (currHit->antisense_align2() ? "-" : "+")
		   << endl;
	      cout << "daehwan - curr coords: "
		   << currHit->left()
		   << "\t"
		   << currHit->right()
		   << endl;
	    }

	  if ((fusion_dir == FUSION_FR || fusion_dir == FUSION_RF) &&
	      prevHit->ref_id2() != currHit->ref_id())
	    continue;

	  if ((fusion_dir == FUSION_FR && !currHit->antisense_align()) ||
	      (fusion_dir == FUSION_RF && currHit->antisense_align()))
	    continue;

	  if (currHit_fused && dir == FUSION_RR)
	    *currHit = currHit->reverse();

	  if (fusion_dir == FUSION_FR || fusion_dir == FUSION_RF ||
	      (currHit_fused && currHit->ref_id() == currHit->ref_id2() && (dir == FUSION_FR || dir == FUSION_RF)))
	    {
	      if (currHit_fused)
		{
		  if ((dir == FUSION_FR && currHit->antisense_align()) ||
		      (dir == FUSION_RF && !currHit->antisense_align()))
		    *currHit = currHit->reverse();
		}
	      else
		{
		  if (fusion_dir == FUSION_FR && currHit->antisense_align())
		    *currHit = currHit->reverse();
		}
	    }
	  /*
	   * Switch prevHit and currHit in FUSION_NOTHING and FUSION_FF cases
	   * to make it easier to check the distance in the gap between the two segments.
	   */
	  else if ((num_fusions == 0 && prevHit->antisense_align() && currHit->antisense_align() && prevHit->ref_id() == currHit->ref_id() && prevHit->left() <= currHit->right() + max_report_intron_length && prevHit->left() + max_insertion_length >= currHit->right()) ||
		   (num_fusions == 1 && (dir == FUSION_FF || dir == FUSION_RR)&&
		    ((!prevHit_fused && prevHit->antisense_align()) || (!currHit_fused && currHit->antisense_align())))
		   )
	    {
	      BowtieHit* tempHit = prevHit;
	      prevHit = currHit;
	      currHit = tempHit;
	    }
	  else if (num_fusions == 0)
	    {
	      if (prevHit->ref_id2() == currHit->ref_id() && prevHit->antisense_align() == currHit->antisense_align())
		{
		  int dist = 0;
		  if (prevHit->antisense_align())
		    dist = prevHit->left() - currHit->right();
		  else
		    dist = currHit->left() - prevHit->right();
		  
		  if (dist > max_report_intron_length || dist < -(int)max_insertion_length)
		    {
		      if ((prevHit->antisense_align() && prevHit->left() > currHit->left()) ||
			  (!prevHit->antisense_align() && prevHit->left() < currHit->left()))
			dir = FUSION_FF;
		      else
			dir = FUSION_RR;
		    }
		}
	      else
		{
		  if (prevHit->antisense_align() == currHit->antisense_align())
		    {
		      if ((prevHit->antisense_align() && prevHit->ref_id() > currHit->ref_id()) ||
			  (!prevHit->antisense_align() && prevHit->ref_id() < currHit->ref_id()))
			dir = FUSION_FF;
		      else
			dir = FUSION_RR;
		    }
		  else if (!prevHit->antisense_align())
		    dir = FUSION_FR;
		  else
		    dir = FUSION_RF;

		  if (dir == FUSION_FR)
		    *currHit = currHit->reverse();
		  else if(dir == FUSION_RF)
		    *prevHit = prevHit->reverse();
		}
	    }

	  // daehwan - test
	  if (bDebug)
	    {
	      cout << "insert id: " << prevHit->insert_id() << endl;
	      cout << "(" << curr - 1 << ") prev: " << prevHit->seq() << " : " << (prevHit->fusion_opcode() != FUSION_NOTHING ? "fused" : "no") << endl;
	      cout << "(" << curr << ") curr: " << currHit->seq() << " : " << (currHit->fusion_opcode() != FUSION_NOTHING ? "fused" : "no") << endl;
	      cout << "prev ref: " << prevHit->ref_id() << "-" << prevHit->ref_id2() << ": " << print_cigar(prevHit->cigar()) << endl;
	      cout << "curr ref: " << currHit->ref_id() << "-" << currHit->ref_id2() << ": " << print_cigar(currHit->cigar()) << endl;
	      cout << "prev coords: "
		   << prevHit->left()
		   << "\t"
		   << prevHit->right()
		   << endl;
	      cout << "curr corrds: "
		   << currHit->left()
		   << "\t"
		   << currHit->right()
		   << endl;
	      cout << "prev sense: "
		   << (prevHit->antisense_align() ? "-" : "+")
		   << "\t"
		   << (prevHit->antisense_align2() ? "-" : "+")
		   << endl;
	      cout << "curr sense: "
		   << (currHit->antisense_align() ? "-" : "+")
		   << "\t"
		   << (currHit->antisense_align2() ? "-" : "+")
		   << endl;
	    }

	  if (num_fusions == 1)
	    {
	      // daehwan
	      if (bDebug)
		{
		  cout << "direction: " << (int)dir << endl;
		}

	      /*
	       * orient the fused segment, which depends on a fusion direction.
	       */
	      if (dir != FUSION_FF && dir != FUSION_RR)
		{
		  bool prevHit_rep = false;
		  bool currHit_rep = false;

		  if (prevHit_fused)
		    {
		      if ((dir == FUSION_FR && !currHit->antisense_align()) ||
			  (dir == FUSION_RF && currHit->antisense_align()))
			continue;

		      if (prevHit->ref_id2() != currHit->ref_id())
			prevHit_rep = true;
		      else if ((dir == FUSION_FR && prevHit->antisense_align()) ||
			       (dir == FUSION_RF && !prevHit->antisense_align()))
			prevHit_rep = true;
		    }
		  
		  if (currHit_fused)
		    {
		      if ((dir == FUSION_FR && prevHit->antisense_align()) ||
			  (dir == FUSION_RF && !prevHit->antisense_align()))
			continue;
		      
		      if (currHit->ref_id() != prevHit->ref_id2())
			currHit_rep = true;
		    }
		  
		  if (bDebug)
		    {
		      if (prevHit_rep) cout << "1. reversed in prev" << endl;
		      if (currHit_rep) cout << "1. reversed in curr" << endl;
		    }

		  if (prevHit_rep)
		    *prevHit = prevHit->reverse();
		  if (currHit_rep)
		    *currHit = currHit->reverse();

		  prevHit_rep = false;
		  currHit_rep = false;

		  if (prevHit_fused)
		    {
		      if (prevHit->is_forwarding_right() != currHit->is_forwarding_left())
			currHit_rep = true;
		    }
		  else
		    {
		      if (prevHit->is_forwarding_right() != currHit->is_forwarding_left())
			prevHit_rep = true;
		    }

		  if (prevHit_rep)
		    *prevHit = prevHit->reverse();
		  if (currHit_rep)
		    *currHit = currHit->reverse();

		  // daehwan
		  if (bDebug)
		    {
		      if (prevHit_rep) cout << "2. reversed in prev" << endl;
		      if (currHit_rep) cout << "2. reversed in curr" << endl;
		    }
		}	      
	    }

	  bool same_contig = prevHit->ref_id2() == currHit->ref_id();
	  if (!same_contig && num_fusions > 0)
	    continue;

	  if (same_contig && num_fusions >= 1)
	    {
	      if (prevHit->antisense_align2() != currHit->antisense_align())
		continue;
	    }

	  int bh_l = 0, back_right = 0, dist = 0;
	  if (same_contig)
	    {
	      if ((fusion_dir == FUSION_FR || fusion_dir == FUSION_RF || dir == FUSION_FR || dir == FUSION_RF) &&
		  prevHit->antisense_align2())
		{
		  bh_l = prevHit->right() + 1;
		  back_right = currHit->left() + 1;
		}
	      else
		{
		  bh_l = currHit->left();
		  back_right = prevHit->right();
		}

	      dist = bh_l - back_right;
	    }

	  // daehwan - pass
	  if (bDebug)
	    {
	      cout << "daehwan - pass" << endl;
	      cout << "prev coords: " << prevHit->left() << "\t" << prevHit->right() << endl;
	      cout << "curr coords: " << currHit->left() << "\t" << currHit->right() << endl;
	    }

	  if (!same_contig ||
	      (same_contig && num_fusions == 0 && dir != FUSION_NOTHING && fusion_dir == FUSION_NOTHING) ||
	      (same_contig && dist <= max_report_intron_length && dist >= -(int)max_insertion_length && prevHit->is_forwarding_right() == currHit->is_forwarding_left()))
	    {
	      // daehwan
	      if (bDebug)
		{
		  cout << "daehwan - really passed!!" << endl;
		}

	      BowtieHit tempHit = seg_hit_stack.back();
	      seg_hit_stack.back() = bh_prev;
	      
	      // these hits are compatible, so push bh onto the 
	      // stack, recurse, and pop it when done.
	      seg_hit_stack.push_back(bh);
	      bool success = dfs_seg_hits(rt,
					  read_seq,
					  read_quals,
					  possible_juncs,
					  possible_insertions,
					  possible_fusions,
					  seg_hits_for_read, 
					  curr + 1,
					  seg_hit_stack,
					  joined_hits,
					  dir == FUSION_NOTHING ? fusion_dir : dir);
	      
	      if (success)
		join_success = true;
	      
	      seg_hit_stack.pop_back();
	      seg_hit_stack.back() = tempHit;
	    }
	}
    }
  else
    {
      merge_segment_chain(rt,
			  read_seq,
			  read_quals,
			  possible_juncs,
			  possible_insertions,
			  possible_fusions,
			  seg_hit_stack,
			  joined_hits,
			  fusion_dir);

      return join_success = true;
    }
  return join_success;
}

bool join_segments_for_read(RefSequenceTable& rt,
			    const string& read_seq,
			    const string& read_quals,
			    std::set<Junction>& possible_juncs,
			    std::set<Insertion>& possible_insertions,
			    std::set<Fusion>& possible_fusions,
			    vector<HitsForRead>& seg_hits_for_read,
			    vector<BowtieHit>& joined_hits)
{	
  vector<BowtieHit> seg_hit_stack;
  bool join_success = false;
  
  // ignore segments that map to more than this many places.
  if (bowtie2)
    {
      for (size_t s = 0; s < seg_hits_for_read.size(); ++s)
	{
	  if (seg_hits_for_read[s].hits.size() > max_seg_multihits)
	    return join_success;
	}
    }

  for (size_t i = 0; i < seg_hits_for_read[0].hits.size(); ++i)
    {
      BowtieHit& bh = seg_hits_for_read[0].hits[i];

      // daehwan - remove this
      //if (bh.insert_id() == 16487)
      //	bDebug = true;
      
      if (bh.fusion_opcode() == FUSION_RR)
	seg_hit_stack.push_back(bh.reverse());
      else
	seg_hit_stack.push_back(bh);
      
      bool success = dfs_seg_hits(rt,
				  read_seq,
				  read_quals,
				  possible_juncs,
				  possible_insertions,
				  possible_fusions,
				  seg_hits_for_read, 
				  1, 
				  seg_hit_stack,
				  joined_hits);
      if (success)
	join_success = true;
      seg_hit_stack.pop_back();
    }
  
  return join_success;
}

struct JoinSegmentsWorker
{
  void operator()()
  {
    ReadTable it;
    GBamWriter bam_writer(bam_output_fname.c_str(), sam_header_fname.c_str());
    ReadStream readstream(reads_fname);
    if (readstream.file() == NULL)
      err_die("Error: cannot open %s for reading\n", reads_fname.c_str());

    if (read_offset > 0)
      readstream.seek(read_offset);

    bool need_seq = true;
    bool need_qual = true;

    vector<HitStream> contig_hits;
    vector<HitStream> spliced_hits;
    vector<HitFactory*> factories;
    for (size_t i = 0; i < segmap_fnames.size(); ++i)
      {
	HitFactory* fac = new BAMHitFactory(it, *rt);
	factories.push_back(fac);
	HitStream hs(segmap_fnames[i], fac, false, false, false, need_seq, need_qual);
	
	if (seg_offsets[i] > 0)
	  hs.seek(seg_offsets[i]);

	contig_hits.push_back(hs);
      }
    
    for (size_t i = 0; i < spliced_segmap_fnames.size(); ++i)
      {
	int anchor_length = 0;
	HitFactory* fac = new SplicedBAMHitFactory(it, *rt, anchor_length);
	factories.push_back(fac);
	
	HitStream hs(spliced_segmap_fnames[i], fac, true, false, false, need_seq, need_qual);
	if (spliced_seg_offsets[i] > 0)
	  hs.seek(spliced_seg_offsets[i]);

	spliced_hits.push_back(hs);
      }

    uint32_t curr_contig_obs_order = VMAXINT32;
    HitStream* first_seg_contig_stream = NULL;
    uint64_t next_contig_id = 0;
    
    if (contig_hits.size() > 0)
      {
	first_seg_contig_stream = &(contig_hits.front());
	next_contig_id = first_seg_contig_stream->next_group_id();
	curr_contig_obs_order = it.observation_order(next_contig_id);
      }

    HitsForRead curr_hit_group;
    
    uint32_t curr_spliced_obs_order = VMAXINT32;
    HitStream* first_seg_spliced_stream = NULL;
    uint64_t next_spliced_id = 0;
    
    if (spliced_hits.size() > 0)
      {
	first_seg_spliced_stream = &(spliced_hits.front());
	next_spliced_id = first_seg_spliced_stream->next_group_id();
	curr_spliced_obs_order = it.observation_order(next_spliced_id);
      }

    while((curr_contig_obs_order != VMAXINT32 || curr_spliced_obs_order != VMAXINT32) &&
	  (curr_contig_obs_order < end_id || curr_spliced_obs_order < end_id))
      {
	uint32_t read_in_process;
	vector<HitsForRead> seg_hits_for_read;
	seg_hits_for_read.resize(contig_hits.size());

	if (curr_contig_obs_order < curr_spliced_obs_order)
	  {
	    first_seg_contig_stream->next_read_hits(curr_hit_group);
	    seg_hits_for_read.front() = curr_hit_group;
	    
	    next_contig_id = first_seg_contig_stream->next_group_id();
	    uint32_t next_order = it.observation_order(next_contig_id);
	    
	    read_in_process = curr_contig_obs_order;
	    curr_contig_obs_order = next_order;
	  }
	else if  (curr_spliced_obs_order < curr_contig_obs_order)
	  {
	    first_seg_spliced_stream->next_read_hits(curr_hit_group);
	    seg_hits_for_read.front() = curr_hit_group;
		
	    next_spliced_id = first_seg_spliced_stream->next_group_id();
	    uint32_t next_order = it.observation_order(next_spliced_id);

	    read_in_process = curr_spliced_obs_order;
	    curr_spliced_obs_order = next_order;

	    if (read_in_process < begin_id)
	      continue;
	  }
	else if (curr_contig_obs_order == curr_spliced_obs_order &&
		 curr_contig_obs_order != VMAXINT32 &&
		 curr_spliced_obs_order != VMAXINT32)
	  {
	    first_seg_contig_stream->next_read_hits(curr_hit_group);
	    
	    HitsForRead curr_spliced_group;
	    first_seg_spliced_stream->next_read_hits(curr_spliced_group);
	    
	    curr_hit_group.hits.insert(curr_hit_group.hits.end(),
				       curr_spliced_group.hits.begin(),
				       curr_spliced_group.hits.end());
	    seg_hits_for_read.front() = curr_hit_group;
	    read_in_process = curr_spliced_obs_order;
	    
	    next_contig_id = first_seg_contig_stream->next_group_id();
	    uint32_t next_order = it.observation_order(next_contig_id);
	    
	    next_spliced_id = first_seg_spliced_stream->next_group_id();
	    uint32_t next_spliced_order = it.observation_order(next_spliced_id);
	    
	    curr_spliced_obs_order = next_spliced_order;
	    curr_contig_obs_order = next_order;
	  }
	else
	  {
	    break;
	  }
	
	if (contig_hits.size() > 1)
	  {
	    look_right_for_hit_group(it,
				     contig_hits,
				     0,
				     spliced_hits,
				     curr_hit_group,
				     seg_hits_for_read);
	  }

	size_t last_non_empty = seg_hits_for_read.size() - 1;
	while(last_non_empty >= 0 && seg_hits_for_read[last_non_empty].hits.empty())
	  {
	    --last_non_empty;
	  }
	
	seg_hits_for_read.resize(last_non_empty + 1);
	if (!seg_hits_for_read[last_non_empty].hits[0].end())
	  continue;
	
	if (!seg_hits_for_read.empty() && !seg_hits_for_read[0].hits.empty())
	  {
	    uint64_t insert_id = seg_hits_for_read[0].hits[0].insert_id();
	    if (insert_id >= begin_id && insert_id < end_id)
	      {
		Read read;
		if (readstream.getRead(insert_id, read))
		  {
		    vector<BowtieHit> joined_hits;
		    join_segments_for_read(*rt,
					   read.seq.c_str(),
					   read.qual.c_str(),
					   *possible_juncs,
					   *possible_insertions,
					   *possible_fusions,
					   seg_hits_for_read,
					   joined_hits);
		    
		    sort(joined_hits.begin(), joined_hits.end());
		    vector<BowtieHit>::iterator new_end = unique(joined_hits.begin(), joined_hits.end());
		    joined_hits.erase(new_end, joined_hits.end());
		    for (size_t i = 0; i < joined_hits.size(); i++)
		      {
			const char* ref_name = rt->get_name(joined_hits[i].ref_id());
			const char* ref_name2 = "";
			if (joined_hits[i].fusion_opcode() != FUSION_NOTHING)
			  ref_name2 = rt->get_name(joined_hits[i].ref_id2());

			vector<string> extra_fields;

			if (!color)
			  bowtie_sam_extra(joined_hits[i], *rt, extra_fields);

			if (color)
			  print_bamhit(bam_writer, read.name.c_str(), joined_hits[i], ref_name, ref_name2,
				       joined_hits[i].seq().c_str(), joined_hits[i].qual().c_str(), true, &extra_fields);
			else
			  print_bamhit(bam_writer, read.name.c_str(), joined_hits[i], ref_name, ref_name2,
				       read.seq.c_str(), read.qual.c_str(), false,  &extra_fields);
		      }
		  }
		else
		  {
		    err_die("Error: could not get read # %d from stream\n",
			    read_in_process);
		  }
	      }
	  }
	else
	  {
	    //fprintf(stderr, "Warning: couldn't join segments for read # %d\n", read_in_process);
	  }
      }

    for (size_t fac = 0; fac < factories.size(); ++fac)
      {
	delete factories[fac];
      }
    factories.clear();
  }

  string bam_output_fname;
  string sam_header_fname;
  string reads_fname;

  vector<string> segmap_fnames;
  vector<string> spliced_segmap_fnames;
  
  std::set<Junction>* possible_juncs;
  std::set<Insertion>* possible_insertions;
  std::set<Fusion>* possible_fusions;
  
  RefSequenceTable* rt;

  uint64_t begin_id;
  uint64_t end_id;
  int64_t read_offset;
  vector<int64_t> seg_offsets;
  vector<int64_t> spliced_seg_offsets;
};

void driver(const string& bam_output_fname,
	    istream& ref_stream,
	    vector<FILE*>& possible_juncs_files,
	    vector<FILE*>& possible_insertions_files,
	    vector<FILE*>& possible_deletions_files,
	    vector<FILE*>& possible_fusions_files,
	    vector<string>& spliced_segmap_fnames, //.bam files
	    vector<string>& segmap_fnames, //.bam files
	    const string& reads_fname)
{
  if (!parallel)
    num_threads = 1;

  if (segmap_fnames.size() == 0)
  {
    fprintf(stderr, "No hits to process, exiting\n");
    exit(0);
  }

  RefSequenceTable rt(sam_header, true);
  fprintf (stderr, "Loading reference sequences...\n");
  get_seqs(ref_stream, rt, true);
    fprintf (stderr, "        reference sequences loaded.\n");

  fprintf(stderr, "Loading junctions...");
  std::set<Junction> possible_juncs;

  for (size_t i = 0; i < possible_juncs_files.size(); ++i)
  {
    char buf[2048];
    while(!feof(possible_juncs_files[i]) &&
        fgets(buf, sizeof(buf), possible_juncs_files[i]))
    {
      char junc_ref_name[256];
      int left;
      int right;
      char orientation;
      int ret = sscanf(buf, "%s %d %d %c", junc_ref_name, &left, &right, &orientation);
      if (ret != 4)
        continue;
      uint32_t ref_id = rt.get_id(junc_ref_name, NULL, 0);
      possible_juncs.insert(Junction(ref_id, left, right, orientation == '-'));
    }
  }
  fprintf(stderr, "done\n");

  fprintf(stderr, "Loading deletions...");

  for (size_t i = 0; i < possible_deletions_files.size(); ++i)
  {
    char splice_buf[2048];
    FILE* deletions_file = possible_deletions_files[i];
    if(!deletions_file){
      continue;
    }
    while(fgets(splice_buf, 2048, deletions_file)){
      char* nl = strrchr(splice_buf, '\n');
      char* buf = splice_buf;
      if (nl) *nl = 0;

			char* ref_name = get_token((char**)&buf, "\t");
			char* scan_left_coord = get_token((char**)&buf, "\t");
			char* scan_right_coord = get_token((char**)&buf, "\t");

      if (!scan_left_coord || !scan_right_coord)
      {
        err_die("Error: malformed deletion coordinate record\n");
      }

      uint32_t ref_id = rt.get_id(ref_name,NULL,0);
      uint32_t left_coord = atoi(scan_left_coord);
      uint32_t right_coord = atoi(scan_right_coord);
      possible_juncs.insert((Junction)Deletion(ref_id, left_coord - 1,right_coord, false));
    }
  }
  fprintf(stderr, "done\n");

  /*
   * Read the insertions from the list of insertion
   * files into a set
   */
  fprintf(stderr, "Loading insertions...");
  std::set<Insertion> possible_insertions;
  for (size_t i=0; i < possible_insertions_files.size(); ++i)
  {
    char splice_buf[2048];
    FILE* insertions_file = possible_insertions_files[i];
    if(!insertions_file){
      continue;
    }
    while(fgets(splice_buf, 2048, insertions_file)){
      char* nl = strrchr(splice_buf, '\n');
      char* buf = splice_buf;
      if (nl) *nl = 0;

			char* ref_name = get_token((char**)&buf, "\t");
			char* scan_left_coord = get_token((char**)&buf, "\t");
			char* scan_right_coord = get_token((char**)&buf, "\t");
			char* scan_sequence = get_token((char**)&buf, "\t");

      if (!scan_left_coord || !scan_sequence || !scan_right_coord)
      {
        err_die("Error: malformed insertion coordinate record\n");
      }

      uint32_t ref_id = rt.get_id(ref_name,NULL,0);
      uint32_t left_coord = atoi(scan_left_coord);
      std::string sequence(scan_sequence);
      possible_insertions.insert(Insertion(ref_id, left_coord, sequence));
    }
  }
  fprintf(stderr, "done\n");

  vector<uint64_t> read_ids;
  vector<vector<int64_t> > offsets;
  if (num_threads > 1)
    {
      vector<string> fnames;
      fnames.push_back(reads_fname);
      fnames.insert(fnames.end(), spliced_segmap_fnames.rbegin(), spliced_segmap_fnames.rend());
      fnames.insert(fnames.end(), segmap_fnames.rbegin(), segmap_fnames.rend());
      bool enough_data = calculate_offsets(fnames, read_ids, offsets);
      if (!enough_data)
	num_threads = 1;
    }

  std::set<Fusion> possible_fusions;
  if (fusion_search)
    {
      fprintf(stderr, "Loading fusions...");
      for (size_t i=0; i < possible_fusions_files.size(); ++i)
	{
	  char splice_buf[2048];	
	  FILE* fusions_file = possible_fusions_files[i];
	  if(!fusions_file){
	    continue;
	  } 
	  while(fgets(splice_buf, 2048, fusions_file)){
	    char* nl = strrchr(splice_buf, '\n');
	    char* buf = splice_buf;
	    if (nl) *nl = 0;
	    
	    char* ref_name1 = strsep((char**)&buf, "\t");
	    char* scan_left_coord = strsep((char**)&buf, "\t");
	    char* ref_name2 = strsep((char**)&buf, "\t");
	    char* scan_right_coord = strsep((char**)&buf, "\t");
	    char* scan_dir = strsep((char**)&buf, "\t");
	    
	    if (!ref_name1 || !scan_left_coord || !ref_name2 || !scan_right_coord || !scan_dir)
	      {
		fprintf(stderr,"Error: malformed insertion coordinate record\n");
		exit(1);
	      }
	    
	    uint32_t ref_id1 = rt.get_id(ref_name1,NULL,0);
	    uint32_t ref_id2 = rt.get_id(ref_name2,NULL,0);
	    uint32_t left_coord = atoi(scan_left_coord);
	    uint32_t right_coord = atoi(scan_right_coord);
	    uint32_t dir = FUSION_FF;
	    if (strcmp(scan_dir, "fr") == 0)
	      dir = FUSION_FR;
	    else if(strcmp(scan_dir, "rf") == 0)
	      dir = FUSION_RF;
	    else if (strcmp(scan_dir, "rr") == 0)
	      dir = FUSION_RR;
	    
	    possible_fusions.insert(Fusion(ref_id1, ref_id2, left_coord, right_coord, dir));
	  }
	}
      fprintf(stderr, "done\n");
    }
      
  vector<boost::thread*> threads;
  for (int i = 0; i < num_threads; ++i)
    {
      JoinSegmentsWorker worker;

      if (num_threads == 1)
	worker.bam_output_fname = bam_output_fname;
      else
	{
	  string filename_base = bam_output_fname.substr(0, bam_output_fname.length() - 4);
	  char filename[1024] = {0};
	  sprintf(filename, "%s%d.bam", filename_base.c_str(), i);
	  worker.bam_output_fname = filename;
	}
      
      worker.sam_header_fname = sam_header;
      worker.reads_fname = reads_fname;
      worker.segmap_fnames = segmap_fnames;
      worker.spliced_segmap_fnames = spliced_segmap_fnames;
      worker.possible_juncs = &possible_juncs;
      worker.possible_insertions = &possible_insertions;
      worker.possible_fusions = &possible_fusions;
      worker.rt = &rt;
      if (i == 0)
	{
	  worker.begin_id = 0;
	  worker.seg_offsets = vector<int64_t>(segmap_fnames.size(), 0);
	  worker.spliced_seg_offsets = vector<int64_t>(spliced_segmap_fnames.size(), 0);
	  worker.read_offset = 0;
	}
      else
	{
	  worker.begin_id = read_ids[i-1];
	  worker.seg_offsets.insert(worker.seg_offsets.end(),
				    offsets[i-1].rbegin(),
				    offsets[i-1].rbegin() + segmap_fnames.size());
	  worker.spliced_seg_offsets.insert(worker.spliced_seg_offsets.end(),
					    offsets[i-1].rbegin() + segmap_fnames.size(),
					    offsets[i-1].rend() - 1);
	  worker.read_offset = offsets[i-1][0];
	}
      
      worker.end_id = (i+1 < num_threads) ? read_ids[i] : std::numeric_limits<uint64_t>::max();

      if (num_threads > 1)
	threads.push_back(new boost::thread(worker));
      else
	worker();
    }

  for (size_t i = 0; i < threads.size(); ++i)
    {
      threads[i]->join();
      delete threads[i];
      threads[i] = NULL;
    }
  threads.clear();
  
} //driver

int main(int argc, char** argv)
{
  fprintf(stderr, "long_spanning_reads v%s (%s)\n", PACKAGE_VERSION, SVN_REVISION);
  fprintf(stderr, "--------------------------------------------\n");

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


  string reads_file_name = argv[optind++];

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string juncs_file_list = argv[optind++];

  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string insertions_file_list = argv[optind++];
  if(optind >= argc)
      {
	print_usage();
	return 1;
      }
  
  string deletions_file_list = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }

  string fusions_file_list = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string bam_output_fname = argv[optind++];
  
  if(optind >= argc)
    {
      print_usage();
      return 1;
    }
  
  string segmap_file_list = argv[optind++];
  string spliced_segmap_flist;
  if(optind < argc)
    {
      spliced_segmap_flist = argv[optind++];
    }

  ifstream ref_stream(ref_file_name.c_str(), ifstream::in);
  if (!ref_stream.good())
      err_die("Error: cannot open %s for reading\n",ref_file_name.c_str());

  checkSamHeader();

  //FILE* reads_file = fopen(reads_file_name.c_str(), "r");
  vector<string> segmap_file_names;
  tokenize(segmap_file_list, ",",segmap_file_names);

  vector<string> juncs_file_names;
  vector<FILE*> juncs_files;
  tokenize(juncs_file_list, ",",juncs_file_names);
  for (size_t i = 0; i < juncs_file_names.size(); ++i) {
      //fprintf(stderr, "Opening %s for reading\n",
      //  juncs_file_names[i].c_str());
      FILE* juncs_file = fopen(juncs_file_names[i].c_str(), "r");
      if (juncs_file == NULL) {
        fprintf(stderr, "Warning: cannot open %s for reading\n",
        juncs_file_names[i].c_str());
        continue;
        }
      juncs_files.push_back(juncs_file);
    }

    /*
     * Read in the deletion file names
     */
    vector<string> deletions_file_names;
    vector<FILE*> deletions_files;
    tokenize(deletions_file_list, ",",deletions_file_names);
    for (size_t i = 0; i < deletions_file_names.size(); ++i)
    {
      //fprintf(stderr, "Opening %s for reading\n",
      //     deletions_file_names[i].c_str());
      FILE* deletions_file = fopen(deletions_file_names[i].c_str(), "r");
      if (deletions_file == NULL) {
            fprintf(stderr, "Warning: cannot open %s for reading\n",
                    deletions_file_names[i].c_str());
            continue;
            }
      deletions_files.push_back(deletions_file);
    }

    /*
     * Read in the list of filenames that contain
     * insertion coordinates
     */

    vector<string> insertions_file_names;
    vector<FILE*> insertions_files;
    tokenize(insertions_file_list, ",",insertions_file_names);
    for (size_t i = 0; i < insertions_file_names.size(); ++i)
    {
     //fprintf(stderr, "Opening %s for reading\n",
     //        insertions_file_names[i].c_str());
     FILE* insertions_file = fopen(insertions_file_names[i].c_str(), "r");
     if (insertions_file == NULL)
        {
            fprintf(stderr, "Warning: cannot open %s for reading\n",
                    insertions_file_names[i].c_str());
            continue;
        }
        insertions_files.push_back(insertions_file);
    }

  vector<string> spliced_segmap_file_names;
  vector<FZPipe> spliced_segmap_files;
  string unzcmd;
  tokenize(spliced_segmap_flist, ",",spliced_segmap_file_names);

  vector<string> fusions_file_names;
  vector<FILE*> fusions_files;
  tokenize(fusions_file_list, ",",fusions_file_names);
  for (size_t i = 0; i < fusions_file_names.size(); ++i)
    {
      fprintf(stderr, "Opening %s for reading\n",
	      fusions_file_names[i].c_str());
      FILE* fusions_file = fopen(fusions_file_names[i].c_str(), "r");
      if (fusions_file == NULL)
        {
	  fprintf(stderr, "Warning: cannot open %s for reading\n",
		  fusions_file_names[i].c_str());
	  continue;
        }
      fusions_files.push_back(fusions_file);
    }
  
  
  driver(bam_output_fname,
	 ref_stream,
	 juncs_files,
	 insertions_files,
	 deletions_files,
	 fusions_files,
	 spliced_segmap_file_names,
	 segmap_file_names,
	 reads_file_name);

  return 0;
}
