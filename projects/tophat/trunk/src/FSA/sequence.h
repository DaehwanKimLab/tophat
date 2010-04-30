
/**
 * \file sequence.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef SEQ_SEQUENCE_INCLUDED
#define SEQ_SEQUENCE_INCLUDED

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "misc.h"
#include "fsa_alphabet.h"

namespace fsa {

  /**
   * \brief Represent a single named sequence.
   *
   * Coordinates are 0-based.
   * Interavls are always fully-closed, [start, end].
   */
  struct Sequence {

  public:

    std::string name;   ///< Sequence name.
    std::string seq;    ///< Sequence itself.

    /**
     * \brief Constructor.
     */
    Sequence (std::string name, std::string seq);

    /**
     * \brief Detect whether a file seems to be in FASTA format.
     */
    static bool detect_fasta (const std::string& filename);

    /**
     * \brief Write sequence in FASTA format.
     */
    void write_fasta (std::ostream& o) const;

    /**
     * \brief Get sequence length.
     */
    inline size_t length() const { return seq.length(); }

    /**
     * \brief Append sequence data.
     */
    inline void append (const std::string& s) { seq += s; }

    /**
     * \brief Get subsequence.
     *
     */
    Sequence subsequence (const unsigned& start, const unsigned& end) const;

    /**
     * \brief Reverse-complement sequence under the passed alphabet.
     * \param alphabet Alphabet to reverse-complement under
     */
    void revcomp (const Alphabet& alphabet);

    /**
     * \brief Complement a strand.
     */
    inline static void complement_strand (char& strand);



    static const std::string fasta_seq_start;   ///< designate sequence name in FASTA files


    /**
     * \brief Initialize hardmasking by finding hardmasked intervals in original sequence.
     *
     * \param min_hardmask_length minimum length of hardmasked sequence to strip out
     * \param is_hardmask_char function pointer defining a hardmasked character
     */
     void init_hardmasking (size_t min_hardmask_length, bool (*is_hardmask_char) (char));

    /**
     * \brief Is the sequence at a position hardmasked?
     */
    bool is_pos_hardmasked (const unsigned& orig_pos) const;

    /**
     * \brief Map coordinates from stripped sequence back to original sequence.
     *
     * Uses a linear-time search over the masked intervals.
     * NB: The coordinate mappings could be made faster with a binary search.
     */
    unsigned map_stripped_to_orig (const unsigned& pos) const;

    /**
     * \brief Map coordinates from original sequence to stripped sequence.
     *
     * Uses a linear-time search over the masked intervals.
     * Returns sequence length (out of bounds)
     * if the coordinate falls in a masked interval.
     */
    unsigned map_orig_to_stripped (const unsigned& pos) const;

    /**
     * \brief Strip out hardmasked sequence.
     */
    Sequence get_stripped_sequence() const;


  private:

    bool __is_hardmasked;  ///< is this sequence hardmasked?

    /**
     * \brief coordinates of hardmasked intervals
     *
     * 0-based, fully-closed coordinates [start, end].
     */
    std::vector<Interval> __masked_intervals;


  };


  /**
   * \brief Represent a translated nucleotide sequence.
   *
   * Assumes and enforces a non-degenerate DNA alphabet:
   * Degenerate characters are randomized and a DNA alphabet is imposed.
   */
  struct Translated_sequence {

  public:

    /**
     * \brief constructor
     */
    Translated_sequence (const Sequence& orig_seq);

    /**
     * \brief Get translated sequence on forward strand.
     * \param frame reading frame 0, 1, 2
     */
    inline Sequence get_forward (const unsigned& frame) const {
      assert (frame < 3);
      return Sequence (orig_seq.name, __forward[frame]);
    }

    /**
     * \brief Get translated sequence on reverse strand.
     */
    inline Sequence get_reverse (const unsigned& frame) const {
      assert (frame < 3);
      return Sequence (orig_seq.name, __reverse[frame]);
    }

    /**
     * \brief Map a coordinate in a translated sequence back to the position in the original sequence.
     *
     * All coordinates are assumed to be 0-based.
     * Maps back to the first nucleotide of the corresponding codon.
     * \param forward true for forward strand, false for reverse
     * \param frame reading frame (0, 1 or 2)
     * \param pos position in translated sequence
     * \return position in original sequence
     */
    inline unsigned map_coords_to_orig (const bool forward, const unsigned& frame, const unsigned& pos) const;

    /**
     * \brief Map an interval in a translated sequence back to the interval in the original sequence.
     *
     * All coordinates are assumed to be 0-based.
     * Maps back to the first nucleotide of first corresponding codon and last nucleotide of the last codon.
     * \param forward true for forward strand, false for reverse
     * \param frame reading frame (0, 1 or 2)
     * \param interval position in translated sequence
     * \return interval in original sequence
     */
    inline Interval map_interval_to_orig (const bool forward, const unsigned& frame, const Interval& interval) const;

  private:

    static const char __codon_map[4][4][4];   ///< map from codons to amino acids

    const Sequence& orig_seq;                 ///< original sequence
    const Alphabet __alphabet;                ///< DNA alphabet

    std::vector<std::string> __forward;       ///< 3 reading frames of forward strand
    std::vector<std::string> __reverse;       ///< 3 reading frames of reverse strand

  };

  /**
   * \brief Represent a set of sequences.
   *
   * This container holds all sequence data,
   * and as such should be instantiated at the "top" of a program.
   */
  struct Sequence_database {

  public:

    /**
     * \brief Constructor.
     */
    Sequence_database() { }

    /**
     * \brief Read a FASTA-format file.
     *
     * Strips out whitespace. Strips out gaps if requested with, e.g.,
     * seq_db.read_fasta ("rob.fasta", Alignment::is_gap_char)
     * Does NOT clear sequence data; this method can therefore be
     * used to read a series of different files into the same Sequence_database.
     * \param is_gap_char function pointer defining a gap character
     */
    void read_fasta (const std::string& filename, bool (*is_gap_char) (char) = NULL);

    /**
     * \brief Write in FASTA format.
     */
    void write_fasta (std::ostream& o) const;

    /**
     * \brief Add sequence to database.
     */
    void add_seq (const Sequence& sequence);

    /**
     * \brief Append sequence data.
     * \param name sequence name
     * \param seq sequence data to append
     */
    inline void append_seq (const std::string& name, const std::string seq);

    /**
     * \brief Does a sequence exist in the database?
     */
    inline bool exists_seq (const std::string& name) const;

    /**
     * \brief Get sequence from database.
     */
    inline const Sequence& get_seq (const std::string& name) const;

    /**
     * \brief Get sequence from database.
     */
    inline Sequence& get_seq (const std::string& name);

    /**
     * \brief Get sequence from database.
     * \param i index of sequence in database
     */
    inline const Sequence& get_seq (const size_t& i) const;

    /**
     * \brief Get sequence from database.
     * \param i index of sequence in database
     */
    inline Sequence& get_seq (const size_t& i);

    /**
     * \brief Get index of sequence in database.
     * 
     * Sequences are indexed according to the order in which they were first stored.
     */
    inline size_t get_seq_index (const std::string& name) const;

    /**
     * \brief Reverse-complement sequences under the passed alphabet.
     * \param alphabet Alphabet to reverse-complement under
     */
    void revcomp (const Alphabet& alphabet);

    /**
     * \brief Do the sequences appear to match a particular alphabet?
     * \param alphabet alphabet to match against
     * \param threshold fraction of sequence characters required to match alphabet
     */
    bool matches_alphabet (const Alphabet& alphabet, const double threshold = 0.95) const;

    /**
     * \brief Number of sequences in database.
     */
    size_t size() const { return __sequences.size(); }

    /**
     * \brief Average length of sequences in database.
     */
    size_t meanlength() const;

    /**
     * \brief Clear all sequence data.
     */
    void clear();

    /**
     * \brief Translate all sequences.
     *
     * "Overhangs" (incomplete codons at the end of the sequence) are dropped.
     * \return the translated database
     */
    Sequence_database translate() const;

    /**
     * \brief Get iterator to start of __sequences.
     */
    std::vector<Sequence>::iterator begin() {
      return __sequences.begin();
    }

    /**
     * \brief Get iterator to start of __sequences.
     */
    std::vector<Sequence>::const_iterator begin() const {
      return __sequences.begin();
    }

    /**
     * \brief Get iterator to end of __sequences.
     */
    std::vector<Sequence>::iterator end() {
      return __sequences.end();
    }

    /**
     * \brief Get iterator to end of __sequences.
     */
    std::vector<Sequence>::const_iterator end() const {
      return __sequences.end();
    }


  private:

    std::vector<Sequence> __sequences;                ///< sequences
    std::map<std::string, size_t> __sequence_index;   ///< lookup index in *this from sequence name

  };




  /****************************************
   * Function definitions.
   ****************************************/

  inline void Sequence::complement_strand (char& strand) {

    if (strand == '+')
      strand = '-';
    else if (strand == '-')
      strand = '+';

  }

  inline unsigned Translated_sequence::map_coords_to_orig (const bool forward, const unsigned& frame, const unsigned& pos) const {

    unsigned orig_pos;
    if (forward)
      orig_pos = (3 * pos + frame);
    else
      orig_pos = ((int)orig_seq.length() - 1) - (3 * pos + frame); // - 1 because 0-based coordinates

    assert (orig_pos < orig_seq.length());
    return orig_pos;
  }

  inline Interval Translated_sequence::map_interval_to_orig (const bool forward, const unsigned& frame, const Interval& interval) const {

    assert (interval.first <= interval.second);

    Interval orig_interval;
    if (forward) {
      orig_interval.first = map_coords_to_orig (forward, frame, interval.first);        // first nt of first codon
      orig_interval.second = map_coords_to_orig (forward, frame, interval.second) + 2;  // last nt of last codon
    }
    else {
      // for reverse strand:
      // first nt of first codon on forward strand is last nt of last codon in revcomp'ed sequence
      orig_interval.first = map_coords_to_orig (forward, frame, interval.second) - 2;
      orig_interval.second = orig_interval.first + 3 * (interval.second - interval.first + 1) - 1; // remember 0-based coordinates
    }

    assert (orig_interval.first <= orig_interval.second);

    return orig_interval;
  }


  inline bool Sequence_database::exists_seq (const std::string& name) const {

    return (__sequence_index.find (name) != __sequence_index.end());

  }

  inline size_t Sequence_database::get_seq_index (const std::string& name) const {
    assert (exists_seq (name));
    return __sequence_index.find (name)->second;
  }

  inline const Sequence& Sequence_database::get_seq (const std::string& name) const {

    assert (exists_seq (name));
    return __sequences[__sequence_index.find (name)->second];

  }

  inline Sequence& Sequence_database::get_seq (const std::string& name) {

    assert (exists_seq (name));
    return __sequences[__sequence_index.find (name)->second];

  }

  inline const Sequence& Sequence_database::get_seq (const size_t& i) const {

    assert (i < __sequences.size());
    return __sequences[i];

  }

  inline Sequence& Sequence_database::get_seq (const size_t& i) {

    assert (i < __sequences.size());
    return __sequences[i];

  }

  inline void Sequence_database::append_seq (const std::string& name, const std::string seq) {
    if (!exists_seq (name))
      add_seq (Sequence (name, seq));
    else
      get_seq (name).append (seq);
    return;
  }

}

#endif /* SEQ_SEQUENCE_INCLUDED */
