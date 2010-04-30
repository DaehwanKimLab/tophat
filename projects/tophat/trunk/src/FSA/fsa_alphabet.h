
/**
 * \file alphabet.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 * The alphabet representation code is loosely based on Ian Holmes's
 * Alphabet class.
 */

#ifndef SEQ_ALPHABET_INCLUDED
#define SEQ_ALPHABET_INCLUDED

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdint.h>

#include "misc.h"

namespace fsa {

  /**
   * \brief Represent a sequence alphabet.
   */
  struct Alphabet {

  public:

    /**
     * \brief Constructor.
     *
     * Unless the alphabet is designated as case-sensitive,
     * alphabet information is stored as lower-case internally.
     * \param name alphabet name
     * \param size alphabet size
     * \param case_sensitive is the alphabet case-sensitive?
     */
    Alphabet (std::string name, size_t size, bool case_sensitive = false);

    /**
     * \brief Make a sequence nondegenerate.
     */
    void make_nondegen (std::string& sequence) const;

    /**
     * \brief Get nondegenerate sequence.
     */
    std::string get_nondegen (const std::string& sequence) const;

    /**
     * \brief Reverse-complement sequence under alphabet.
     *
     * Degenerate characters are set to __unknown_char
     * after reverse-complementing.
     * \param seq sequence to be reverse-complemented
     */
    void revcomp (std::string& seq) const;

    /**
     * \brief Reverse-complement sequence under alphabet.
     * \see revcomp
     */
    std::string revcomp (const std::string& seq) const;

    /**
     * \brief Get alphabet size.
     */
    inline size_t size() const { return __size; };

    /**
     * \brief Is this a non-degenerate character in the alphabet?
     */
    inline bool is_nondegen_char (const char c) const;

    /**
     * \brief Is this a degenerate character in the alphabet?
     */
    inline bool is_degen_char (const char c) const;

    /**
     * \brief Does the alphabet contain a possibly-degenerate character?
     * \see is_nondegen_char
     * \see is_degen_char
     */
    inline bool contains_char (const char c) const;

    /**
     * \brief Is this an unknown character?
     */
    inline bool is_unknown_char (const char c) const;

    /**
     * \brief Convert a possibly-degenerate character to a non-degenerate character.
     *
     * Randomizes unknown characters across the entire alphabet.
     * Preserves character case.
     */
    inline char get_nondegen_char (const char c) const;

    /**
     * \brief Get the numerical index for a character.
     *
     * If the character is degenerate or unknown,
     * randomizes the character prior to getting the index.
     */
    size_t get_char_index (const char c) const;

  protected:

    /**
     * \brief Initialize characters in alphabet.
     */
    void init_chars (const std::string& chars, const std::string complement = "");

    /**
     * \brief Add a degenerate character to the alphabet.
     */
    void add_degen_char (const char c, const std::string& nondegen);

    /**
     * \brief Set the unknown character for the alphabet.
     */
    void set_unknown_char (const char c);

  private:

    std::string __name;          ///< alphabet name
    size_t __size;               ///< alphabet size (# of non-degenerate characters)
    bool __case_sensitive;       ///< is the alphabet case-sensitive?
    bool __has_complement;       ///< does the alphabet have a complement?

    std::vector<char> __char_list;                          ///< ordered list of characters
    std::vector<char> __char_complement_list;               ///< ordered list of complementary characters
    std::map<char, size_t> __char_index;                    ///< map from characters to index in __char_list
    char __unknown_char;                                    ///< unknown degenerate character
    std::map<char, std::vector<char> > __degen_char_map;    ///< map from degenerate to non-degenerate characters

  };

  /**
   * \brief Represent a DNA alphabet.
   */
  struct DNA_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    DNA_alphabet();

    /**
     * \brief Define the hardmasked character (ignores case!).
     */
    static inline bool is_hardmask_char (const char c) {
      return toupper (c) == 'N';
    }

  private:

    static const std::string DNA_alphabet_name;   ///< alphabet name

  };

  /**
   * \brief Represent a RNA alphabet.
   */
  struct RNA_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    RNA_alphabet();

  private:

    static const std::string RNA_alphabet_name;   ///< alphabet name

  };

  /**
   * \brief Represent a protein alphabet.
   *
   * Characters are in IUPAC alphabetical order, 'acdefghiklmnpqrstvwy'.
   */
  struct Protein_alphabet : public Alphabet {

  public:

    /**
     * \brief Constructor.
     */
    Protein_alphabet();

  private:

    static const std::string Protein_alphabet_name;   ///< alphabet name

  };

  inline bool Alphabet::is_nondegen_char (const char c) const {
    return (__char_index.find (__case_sensitive ? c : tolower (c)) != __char_index.end());
  }

  inline bool Alphabet::is_degen_char (const char c) const {
    return (__degen_char_map.find (__case_sensitive ? c : tolower (c)) != __degen_char_map.end());
  }

  inline bool Alphabet::contains_char (const char c) const {
    return (is_nondegen_char (c) || is_degen_char (c));
  }

  inline bool Alphabet::is_unknown_char (const char c) const {
    return (!is_nondegen_char (c) && !is_degen_char (c));
  }

  inline char Alphabet::get_nondegen_char (const char c) const {

    // if upper-case
    if (isupper (c)) {

      // if nondegenerate
      if (is_nondegen_char (c))
	return c;

      // else if degenerate
      else if (is_degen_char (c)) {
	const std::vector<char>& nondegen = __degen_char_map.find (__case_sensitive ? c : tolower (c))->second;
	return toupper (nondegen[Util::rand ((uint32_t)(nondegen.size() - 1))]);
      }
      // else we don't know anything about it, so treat as unknown
      else
	return toupper (__char_list[Util::rand ((uint32_t)(__size - 1))]);

    }

    // if lower-case
    else {

      // if nondegenerate
      if (is_nondegen_char (c))
	return c;

      // else if degenerate
      else if (is_degen_char (c)) {
	const std::vector<char>& nondegen = __degen_char_map.find (c)->second;
	return tolower (nondegen[Util::rand ((uint32_t)(nondegen.size() - 1))]);
      }
      // else we don't know anything about it, so treat as unknown
      else {
	return tolower (__char_list[Util::rand ((uint32_t)(__size - 1))]);
      }

    }

  }

}

#endif /* SEQ_ALPHABET_INCLUDED */
