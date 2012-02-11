#ifndef DELETIONS_H
#define DELETIONS_H
/*
 *  deletions.h
 *  TopHat
 *
 *  Created by Ryan Kelley on 11/09/2010 
 *
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
#include "junctions.h"
using namespace std;

typedef Junction Deletion;
typedef std::map<Junction, uint32_t> DeletionSet;

void deletions_from_alignment(const BowtieHit& spliced_alignment, DeletionSet& junctions);
void deletions_from_spliced_hit(const BowtieHit& bh, vector<Deletion>& deletions);
void print_deletions(FILE* deletions_out, const DeletionSet& deletions, RefSequenceTable& ref_sequences);
void merge_with(DeletionSet& deletions, const DeletionSet& other);

#endif
