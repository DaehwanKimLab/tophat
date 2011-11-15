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

BowtieHit merge_chain(RefSequenceTable& rt,
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
     *   -Don't start or end with an indel or refskip
     *  -Only a match operation is allowed is allowed
     *   adjacent to an indel or refskip
     *      -Indels should confrom to length restrictions
     */
    const CigarOp* prevCig = &(bh.cigar()[0]);
    const CigarOp* currCig = &(bh.cigar()[1]);
    for (size_t i = 1; i < bh.cigar().size(); ++i){
      currCig = &(bh.cigar()[i]);
      if(currCig->opcode != MATCH && prevCig->opcode != MATCH){
        return false;
      }
      if(currCig->opcode == INS){
        if(currCig->length > max_insertion_length){
          return false;
        }
      }
      if(currCig->opcode == DEL){
        if(currCig->length > max_deletion_length){
          return false;
        }
      }
      if(currCig->opcode == REF_SKIP){
        if(currCig->length < (uint64_t)min_report_intron_length){
          gap_too_short++;
          return false;
        }
      }
      prevCig = currCig;
    }
    if (bh.cigar().front().opcode != MATCH ||
      bh.cigar().back().opcode != MATCH /* ||
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
       vector<BowtieHit>& hits,
       vector<BowtieHit>& merged_hits)
{
  if (hits.size() == 0)
    return;

  BowtieHit bh;
  if (hits.size() > 1)
    {
      list<BowtieHit> hit_chain;

      if (hits.front().antisense_align())
  copy(hits.rbegin(), hits.rend(), back_inserter(hit_chain));
      else
  copy(hits.begin(), hits.end(), back_inserter(hit_chain));

      bh = merge_chain(rt, read_seq, read_quals, possible_juncs, possible_insertions, hit_chain);
    }
  else
    {
      bh = hits[0];
    }

  if (valid_hit(bh))
      merged_hits.push_back(bh);
}

bool dfs_seg_hits(RefSequenceTable& rt,
      const string& read_seq,
      const string& read_quals,
      std::set<Junction>& possible_juncs,
      std::set<Insertion>& possible_insertions,
      vector<HitsForRead>& seg_hits_for_read,
      size_t curr,
      vector<BowtieHit>& seg_hit_stack,
      vector<BowtieHit>& joined_hits)
{
  assert (!seg_hit_stack.empty());
  bool join_success = false;

  if (curr < seg_hits_for_read.size())
  {
    for (size_t i = 0; i < seg_hits_for_read[curr].hits.size(); ++i)
    {
      BowtieHit& bh = seg_hits_for_read[curr].hits[i];
      BowtieHit& back = seg_hit_stack.back();

      bool consistent_sense = bh.antisense_align() == back.antisense_align();
      bool same_contig = bh.ref_id() == back.ref_id();

      // FIXME: when we have stranded reads, we need to fix this condition
      //bool consistent_strand = (bh.contiguous() || back.contiguous() ||
      //              (bh.antisense_splice() == back.antisense_splice()));
      if (consistent_sense && same_contig /*&& consistent_strand*/)
      {
        if (bh.antisense_align())
        {
          unsigned int bh_r = bh.right();
          unsigned int back_left = seg_hit_stack.back().left();
          if ((bh_r + max_report_intron_length >= back_left &&
            back_left >= bh_r) || (back_left + max_insertion_length >= bh_r && back_left < bh_r))
          {
            // these hits are compatible, so push bh onto the
            // stack, recurse, and pop it when done.
            seg_hit_stack.push_back(bh);
            bool success = dfs_seg_hits(rt,
                      read_seq,
                      read_quals,
                      possible_juncs,
                      possible_insertions,
                      seg_hits_for_read,
                      curr + 1,
                      seg_hit_stack,
                      joined_hits);
            if (success)
              join_success = true;

            seg_hit_stack.pop_back();
          }
        }
        else
        {
          unsigned int bh_l = bh.left();

          unsigned int back_right = seg_hit_stack.back().right();
          if ((back_right + max_report_intron_length >= bh_l  &&
            bh_l >= back_right) || (bh_l + max_insertion_length >= back_right && bh_l < back_right))
          {
            // these hits are compatible, so push bh onto the
            // stack, recurse, and pop it when done.
            seg_hit_stack.push_back(bh);
            bool success = dfs_seg_hits(rt,
                      read_seq,
                      read_quals,
                      possible_juncs,
                      possible_insertions,
                      seg_hits_for_read,
                      curr + 1,
                      seg_hit_stack,
                      joined_hits);

            if (success)
              join_success = true;

            seg_hit_stack.pop_back();
          }
        }
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
            seg_hit_stack,
            joined_hits);
    return join_success = true;
  }
  return join_success;
}

bool join_segments_for_read(RefSequenceTable& rt,
          const string& read_seq,
          const string& read_quals,
          std::set<Junction>& possible_juncs,
          std::set<Insertion>& possible_insertions,
          vector<HitsForRead>& seg_hits_for_read,
          vector<BowtieHit>& joined_hits)
{
  vector<BowtieHit> seg_hit_stack;
  bool join_success = false;

  for (size_t i = 0; i < seg_hits_for_read[0].hits.size(); ++i)
    {
      BowtieHit& bh = seg_hits_for_read[0].hits[i];
      seg_hit_stack.push_back(bh);
      bool success = dfs_seg_hits(rt,
          read_seq,
          read_quals,
          possible_juncs,
          possible_insertions,
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

void join_segment_hits(GBamWriter& bam_writer, std::set<Junction>& possible_juncs,
           std::set<Insertion>& possible_insertions,
           ReadTable& unmapped_reads,
           RefSequenceTable& rt,
           //FILE* reads_file,
           ReadStream& readstream,
           vector<HitStream>& contig_hits,
           vector<HitStream>& spliced_hits)
{
  uint32_t curr_contig_obs_order = VMAXINT32;
  HitStream* first_seg_contig_stream = NULL;
  uint64_t next_contig_id = 0;

  if (contig_hits.size())
    {
      first_seg_contig_stream = &(contig_hits.front());
      next_contig_id = first_seg_contig_stream->next_group_id();
      curr_contig_obs_order = unmapped_reads.observation_order(next_contig_id);
    }

  HitsForRead curr_hit_group;

  uint32_t curr_spliced_obs_order = VMAXINT32;
  HitStream* first_seg_spliced_stream = NULL;
  uint64_t next_spliced_id = 0;

  if (spliced_hits.size())
    {
      first_seg_spliced_stream = &(spliced_hits.front());
      next_spliced_id = first_seg_spliced_stream->next_group_id();
      curr_spliced_obs_order = unmapped_reads.observation_order(next_spliced_id);
    }

  while(curr_contig_obs_order != VMAXINT32 ||
  curr_spliced_obs_order != VMAXINT32)
    {
      uint32_t read_in_process;
      vector<HitsForRead> seg_hits_for_read;
      seg_hits_for_read.resize(contig_hits.size());

      if (curr_contig_obs_order < curr_spliced_obs_order)
  {
    first_seg_contig_stream->next_read_hits(curr_hit_group);
    seg_hits_for_read.front() = curr_hit_group;

    next_contig_id = first_seg_contig_stream->next_group_id();
    uint32_t next_order = unmapped_reads.observation_order(next_contig_id);

    read_in_process = curr_contig_obs_order;
    curr_contig_obs_order = next_order;
  }
      else if  (curr_spliced_obs_order < curr_contig_obs_order)
  {
    first_seg_spliced_stream->next_read_hits(curr_hit_group);
    seg_hits_for_read.front() = curr_hit_group;

    next_spliced_id = first_seg_spliced_stream->next_group_id();
    uint32_t next_order = unmapped_reads.observation_order(next_spliced_id);

    read_in_process = curr_spliced_obs_order;
    curr_spliced_obs_order = next_order;
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
    uint32_t next_order = unmapped_reads.observation_order(next_contig_id);

    next_spliced_id = first_seg_spliced_stream->next_group_id();
    uint32_t next_spliced_order = unmapped_reads.observation_order(next_spliced_id);

    curr_spliced_obs_order = next_spliced_order;
    curr_contig_obs_order = next_order;
  }
      else
  {
    break;
  }

      if (contig_hits.size() > 1)
      {
        look_right_for_hit_group(unmapped_reads,
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
    Read read;
    /* if (get_read_from_stream(insert_id,
           reads_file,
           reads_format,
           false,
           read)) */
    if (readstream.getRead(insert_id, read))
      {
        vector<BowtieHit> joined_hits;
        join_segments_for_read(rt,
             read.seq.c_str(),
             read.qual.c_str(),
             possible_juncs,
             possible_insertions,
             seg_hits_for_read,
             joined_hits);

        sort(joined_hits.begin(), joined_hits.end());
        vector<BowtieHit>::iterator new_end = unique(joined_hits.begin(), joined_hits.end());
        joined_hits.erase(new_end, joined_hits.end());

        for (size_t i = 0; i < joined_hits.size(); i++)
        {
          const char* ref_name = rt.get_name(joined_hits[i].ref_id());
          if (color && !color_out)
            //print_hit(stdout, read_name, joined_hits[i], ref_name, joined_hits[i].seq().c_str(), joined_hits[i].qual().c_str(), true);
            print_bamhit(bam_writer, read.name.c_str(), joined_hits[i], ref_name, joined_hits[i].seq().c_str(),
                                         joined_hits[i].qual().c_str(), true);
          else
            print_bamhit(bam_writer, read.name.c_str(), joined_hits[i], ref_name,
                                          read.seq.c_str(), read.qual.c_str(), false);
            //print_hit(stdout, read_name, joined_hits[i], ref_name, read_seq, read_quals, false);
        }
      }
    else
      {
        err_die("Error: could not get read # %d from stream\n",
               read_in_process);
      }
  }
      else
  {
    //fprintf(stderr, "Warning: couldn't join segments for read # %d\n", read_in_process);
  }
 }
}

void driver(GBamWriter& bam_writer, istream& ref_stream,
      vector<FILE*>& possible_juncs_files,
      vector<FILE*>& possible_insertions_files,
      vector<FILE*>& possible_deletions_files,
      vector<FZPipe>& spliced_seg_files,
      vector<FZPipe>& seg_files,
      ReadStream& readstream)
{
  if (seg_files.size() == 0)
  {
    fprintf(stderr, "No hits to process, exiting\n");
    exit(0);
  }

  RefSequenceTable rt(true, true);
  fprintf (stderr, "Loading reference sequences...\n");
  get_seqs(ref_stream, rt, true, false);
    fprintf (stderr, "        reference sequences loaded.\n");
  ReadTable it;

  bool need_seq = true;
  bool need_qual = color;
  //rewind(reads_file);
  readstream.rewind();

  //vector<HitTable> seg_hits;
  vector<HitStream> contig_hits;
  vector<HitStream> spliced_hits;

  vector<HitFactory*> factories;
  for (size_t i = 0; i < seg_files.size(); ++i)
  {
    HitFactory* fac = new BowtieHitFactory(it, rt);
    factories.push_back(fac);
    HitStream hs(seg_files[i].file, fac, false, false, false, need_seq, need_qual);
    contig_hits.push_back(hs);
  }

  fprintf(stderr, "Loading spliced hits...");

  //vector<vector<pair<uint32_t, HitsForRead> > > spliced_hits;
  for (size_t i = 0; i < spliced_seg_files.size(); ++i)
  {
    int anchor_length = 0;

    // we allow read alignments even 1bp overlapping with juctions.
    // if (i == 0 || i == spliced_seg_files.size() - 1)
    // anchor_length = min_anchor_len;

    HitFactory* fac = new SplicedBowtieHitFactory(it,
                    rt,
                    anchor_length);
    factories.push_back(fac);

    HitStream hs(spliced_seg_files[i].file, fac, true, false, false, need_seq, need_qual);
    spliced_hits.push_back(hs);

  }
  fprintf(stderr, "done\n");

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

  join_segment_hits(bam_writer, possible_juncs, possible_insertions, it, rt, readstream, contig_hits, spliced_hits);
  //join_segment_hits(possible_juncs, possible_insertions, it, rt, reads_file.file, contig_hits, spliced_hits);

  for (size_t i = 0; i < seg_files.size(); ++i)
      seg_files[i].close();
    for (size_t i = 0; i < spliced_seg_files.size(); ++i)
        spliced_seg_files[i].close();
    readstream.close();

  for (size_t fac = 0; fac < factories.size(); ++fac)
  {
    delete factories[fac];
  }
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


  string segment_file_list = argv[optind++];
  string spliced_segment_file_list;
  if(optind < argc)
    {
      spliced_segment_file_list = argv[optind++];
    }

  ifstream ref_stream(ref_file_name.c_str(), ifstream::in);
  if (!ref_stream.good())
      err_die("Error: cannot open %s for reading\n",ref_file_name.c_str());

  checkSamHeader();

  vector<string> segment_file_names;
  tokenize(segment_file_list, ",",segment_file_names);

  // string unzcmd=getUnpackCmd(reads_file_name, spliced_segment_file_list.empty() &&
  //              segment_file_names.size()<4);

  //FZPipe reads_file(reads_file_name, unzcmd);
  ReadStream readstream(reads_file_name);
  if (readstream.file()==NULL)
     err_die("Error: cannot open %s for reading\n",
        reads_file_name.c_str());

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

    //vector<FILE*> segment_files;
    vector<FZPipe> segment_files;
    string unzcmd;
    for (size_t i = 0; i < segment_file_names.size(); ++i)
    {
      if (unzcmd.empty())
          unzcmd=getUnpackCmd(segment_file_names[i],false);
      fprintf(stderr, "Opening %s for reading\n",
        segment_file_names[i].c_str());
      //FILE* seg_file = fopen(segment_file_names[i].c_str(), "r");
      FZPipe seg_file(segment_file_names[i], unzcmd);
      if (seg_file.file == NULL)
        {
        err_die("Error: cannot open %s for reading\n",
              segment_file_names[i].c_str());
        }
      segment_files.push_back(seg_file);
    }

  vector<string> spliced_segment_file_names;
  //vector<FILE*> spliced_segment_files;
  vector<FZPipe> spliced_segment_files;
  unzcmd.clear();
  tokenize(spliced_segment_file_list, ",",spliced_segment_file_names);
  for (size_t i = 0; i < spliced_segment_file_names.size(); ++i)
    {
      fprintf(stderr, "Opening %s for reading\n",
        spliced_segment_file_names[i].c_str());
      if (unzcmd.empty())
          unzcmd=getUnpackCmd(spliced_segment_file_names[i],false);
      //FILE* spliced_seg_file = fopen(spliced_segment_file_names[i].c_str(), "r");
      FZPipe spliced_seg_file(spliced_segment_file_names[i], unzcmd);
      if (spliced_seg_file.file == NULL)
        {
        err_die("Error: cannot open %s for reading\n",
             spliced_segment_file_names[i].c_str());
        }
        spliced_segment_files.push_back(spliced_seg_file);
    }

  if (spliced_segment_files.size())
    {
      if (spliced_segment_files.size() != segment_files.size())
        {
        err_die("Error: each segment file must have a corresponding spliced segment file\n");
        }
    }
  GBamWriter bam_writer("-", sam_header.c_str());
  driver(bam_writer, ref_stream, juncs_files, insertions_files, deletions_files, spliced_segment_files, segment_files, readstream);
  //driver(ref_stream, juncs_files, insertions_files, deletions_files, spliced_segment_files, segment_files, reads_file);
  return 0;
}
