#ifndef ISLANDS_H
#define ISLANDS_H
/*
 *  islands.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 10/17/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

// Keeps metadata about each island in the mapping for later reporting in
// the gff file
struct Island
{
	Island(const char* r, int s, int l, int d, unsigned int i) 
	: ref_src(r), start_pos(s), length(l), total_depth(d), id(i), d_stat(0.0) {}
	const char* ref_src;
	int start_pos;
	int length; // This is the *unextended* length of the island
	int total_depth;
	unsigned int id;
	long double d_stat;
};

#endif
