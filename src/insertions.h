#ifndef INSERTIONS_H
#define INSERTIONS_H
/*
 *  insertions.h
 *  TopHat
 * 
 *  Adapted from junctions.h
 */

#include <cstdio>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <cstring>
#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>

#include "bwt_map.h"

using namespace std;


/**
 * Data structure to represent an insertion.
 * Need to keep track of the actual position of the insertion
 * as well the actual inserted sequence.
 */
struct Insertion 
{
  Insertion (uint32_t ref, uint32_t l, const std::string& seq)
  : refid(ref), left(l), sequence(seq){}
Insertion() : refid(0), left(0), sequence("") {}
  /**
   * The ID of the assoicated reference sequence (eg, chr22). In order to actually turn this into
   * a string, need to use the map associated with the reference table
   */
  uint32_t refid;
  /**
   * The position of the insertion.
   * This is the 0-based position of the last genomic nucleotide before the insertion
   */
  uint32_t left;
  
  /**
   * The actual inserted sequence.
   */
  std::string sequence;
  
  bool operator<(const Insertion& rhs) const
  {
    if (refid < rhs.refid)
      return true;
    else if (refid > rhs.refid)
      return false;
    
    if (left < rhs.left)
      return true;
    else if (left > rhs.left)
      return false;
    
    if (sequence.length() < rhs.sequence.length())
      return true;
    return false;
  }
  
  
  bool operator==(const Insertion& rhs) const
  {
    return  (refid == rhs.refid && left == rhs.left && sequence == rhs.sequence);
  }
  
};

/**
 * A function used to compare Insertions, specifically for use
 * with std::sets. My C++ is a little weak, but I imagine
 * there should be someway to directly address the overloaded
 * operator function and ditch this code.
 */ 
struct insertion_comparison 
{
  bool operator()(const Insertion& lhs, const Insertion& rhs)
  {
    return lhs < rhs;
  }
};

typedef std::map<Insertion, uint32_t> InsertionSet;

void insertions_from_alignment(const BowtieHit& bh, InsertionSet& insertions);
void print_insertions(FILE* insertions_out, const InsertionSet& insertions, RefSequenceTable& ref_sequences);
void insertions_from_spliced_hit(const BowtieHit& bh, vector<Insertion>& insertions);
void merge_with(InsertionSet& insertions, const InsertionSet& other);

#endif
