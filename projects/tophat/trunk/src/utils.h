#ifndef UTILS_H
#define UTILS_H
/*
 *  utils.h
 *  TopHat
 *
 *  Created by Daehwan Kim on 12/28/11.
 *  Copyright 2011 Daehwan Kim. All rights reserved.
 *
 */

#include <vector>
#include <string>
using namespace std;

#include "common.h"

// this is for parallelization purposes in segment_juncs, long_spanning_reads, and tophat_reports.
// given "index" files, it calculates "read ids" in increasing order and
// their corresponding file offsets.
bool calculate_offsets(const vector<string>& fnames,
		       vector<uint64_t>& ids,
		       vector<vector<int64_t> >& offsets);

// given "read ids" as reference read ids,
// it finds the closest read ids (with file offsets) not greater than them.
void calculate_offsets_from_ids(const string& fname,
				const vector<uint64_t>& ids,
				vector<int64_t>& offsets);

#endif
