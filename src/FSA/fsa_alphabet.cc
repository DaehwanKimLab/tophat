
/**
 * \file sequence.cc
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#include <algorithm>

#include "fsa_alphabet.h"

using namespace fsa;

const std::string DNA_alphabet::DNA_alphabet_name = "DNA";
const std::string RNA_alphabet::RNA_alphabet_name = "RNA";
const std::string Protein_alphabet::Protein_alphabet_name = "Protein";



Alphabet::Alphabet (std::string name, size_t size, bool case_sensitive /* = false */)
  : __name (name), __size (size), __case_sensitive (case_sensitive) {

}

void Alphabet::init_chars (const std::string& chars, const std::string complement /* = "" */) {

  __has_complement = (complement.length() > 0);

  // check sane
  assert (__size == chars.length());
  if (__has_complement) {
    assert (__size == complement.length());
    __char_complement_list.resize (__size);
  }

  // store the alphabet
  __char_list.resize (__size);  
  for (size_t i = 0; i < chars.length(); ++i) {
    const char& c = chars[i];
    __char_list[i] = __case_sensitive ? c : tolower (c);
    if (__has_complement)
      __char_complement_list[i] = complement[i];
    __char_index.insert (std::make_pair (__case_sensitive ? c : tolower (c), i));
  }
  
}

void Alphabet::add_degen_char (const char c, const std::string& nondegen) {

  std::vector<char> nondegen_chars (nondegen.length());
  for (size_t i = 0; i < nondegen.length(); ++i)
    nondegen_chars[i] = nondegen[i];
  __degen_char_map.insert (std::make_pair (c, nondegen_chars));
  
}

void Alphabet::set_unknown_char (const char c) {

  __unknown_char = c;

}

void Alphabet::make_nondegen (std::string& seq) const {

  for (std::string::iterator c = seq.begin(); c != seq.end(); ++c)
    *c = get_nondegen_char (*c);

}

std::string Alphabet::get_nondegen (const std::string& seq) const {

  std::string seq_nondegen;
  seq_nondegen.resize (seq.length());
  for (size_t i = 0; i < seq.length(); ++i)
    seq_nondegen[i] = get_nondegen_char (seq[i]);

  return seq_nondegen;

}

size_t Alphabet::get_char_index (const char c) const {

  return __char_index.find (__case_sensitive ? get_nondegen_char (c) : tolower (get_nondegen_char (c)))->second;

}

void Alphabet::revcomp (std::string& seq) const {

  if (!__has_complement) {
    cerr << "Tried to reverse-complement under an alphabet without a complementary alphabet." << endl;
    exit (1);
  }

  // first reverse
  std::reverse (seq.begin(), seq.end());

  // then complement
  for (std::string::iterator c = seq.begin(); c != seq.end(); ++c) {
    // complement nondegen chars
    if (is_nondegen_char (*c)) {
      if (isupper (*c))
	*c = toupper (__char_complement_list[__char_index.find (tolower (__case_sensitive ? *c : tolower (*c)))->second]);
      else
	*c = tolower (__char_complement_list[__char_index.find (*c)->second]);
    }
    // and set degen & unknown chars to __unknown_char
    else {
      *c = isupper (*c) ? toupper (__unknown_char) : __unknown_char;
    }
  }

}

std::string Alphabet::revcomp (const std::string& seq) const {
  
  std::string seq_revcomp = seq;
  revcomp (seq_revcomp);
  return seq_revcomp;

}

DNA_alphabet::DNA_alphabet()
  : Alphabet (DNA_alphabet_name, 4, false) {

  init_chars ("acgt", "tgca");
  set_unknown_char ('n');
  add_degen_char ('u', "t");
  add_degen_char ('r', "ag");
  add_degen_char ('y', "ct");
  add_degen_char ('m', "ac");
  add_degen_char ('k', "gt");
  add_degen_char ('s', "cg");
  add_degen_char ('w', "at");
  add_degen_char ('h', "act");
  add_degen_char ('b', "cgt");
  add_degen_char ('v', "acg");
  add_degen_char ('d', "agt");

}

RNA_alphabet::RNA_alphabet()
  : Alphabet (RNA_alphabet_name, 4, false) {

  init_chars ("acgu", "ugca");
  set_unknown_char ('n');
  add_degen_char ('t', "u");
  add_degen_char ('r', "ag");
  add_degen_char ('y', "ct");
  add_degen_char ('m', "ac");
  add_degen_char ('k', "gt");
  add_degen_char ('s', "cg");
  add_degen_char ('w', "at");
  add_degen_char ('h', "act");
  add_degen_char ('b', "cgt");
  add_degen_char ('v', "acg");
  add_degen_char ('d', "agt");

}

Protein_alphabet::Protein_alphabet()
  : Alphabet (Protein_alphabet_name, 20, false) {

  init_chars ("acdefghiklmnpqrstvwy");
  set_unknown_char ('x');
  add_degen_char ('b', "nd");
  add_degen_char ('z', "qe");

}
