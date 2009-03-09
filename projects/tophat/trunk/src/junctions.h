#ifndef JUNCTIONS_H
#define JUNCTIONS_H
/*
 *  junctions.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/22/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include <cstdio>
#include <vector>
#include <string>
//#include <map>
//#include <algorithm>
#include <set>
//#include <stdexcept>
#include <iostream>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>

#include "bwt_map.h"

using namespace std;

struct Junction
{
	
	Junction (uint32_t ref, uint32_t l, uint32_t r, bool a)
	: refid(ref), left(l), right(r), antisense(a) {}
	Junction() {}
	uint32_t refid;
	uint32_t left;
	uint32_t right;
	bool antisense;
	
	bool operator<(const Junction& rhs) const
	{
		if (refid < rhs.refid)
			return true;
		else if (refid > rhs.refid)
			return false;
		
		if (left < rhs.left)
			return true;
		else if (left > rhs.left)
			return false;
		
		if (right < rhs.right)
			return true;
		else if (right > rhs.right)
			return false;
		
		return antisense < rhs.antisense;
	}
	
#if !NDEBUG
	bool valid() const
	{
		return refid != 0xFFFFFFFF && left < right && (left != right);
	}
#endif
};

struct JunctionStats
{
	uint8_t left_extent;
	uint8_t right_extent;
	set<BowtieHit*> supporting_hits;
	bool accepted;
};

typedef std::map<Junction, JunctionStats> JunctionSet;

// This routine DOES NOT set the real refid!  
pair<Junction, JunctionStats> junction_from_spliced_hit(const BowtieHit& h);

void print_junction(FILE* junctions_out, 
					const string& name, 
					const Junction& j, 
					const JunctionStats& s, 
					uint32_t junc_id);




void junctions_from_alignment(BowtieHit& spliced_alignment,

							 JunctionSet& junctions);

void junctions_from_alignments(HitTable& hits,
							   JunctionSet& junctions);

void accept_valid_junctions(JunctionSet& junctions,
							const uint32_t refid,
							const vector<short>& DoC,
							double min_isoform_fraction);

void accept_all_junctions(JunctionSet& junctions,
							const uint32_t refid);

void print_junctions(FILE* junctions_out, 
					 const JunctionSet& junctions,
					 RefSequenceTable& ref_sequences);

#endif

