#ifndef WIGGLES_H
#define WIGGLES_H

/*
 *  wiggles.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 12/12/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#include "common.h"
#include "bwt_map.h"

void print_wiggle_header(FILE* coverage_out);
void print_wiggle_for_ref(FILE* coverage_out,
						  const string& ref_name,
						  const vector<unsigned int>& DoC);


#endif
