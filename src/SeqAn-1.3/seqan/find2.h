// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Umbrella header for the find module.
// ==========================================================================

#ifndef SEQAN_SEQAN_FIND_H_
#define SEQAN_SEQAN_FIND_H_

// Prerequisites.
#include <cmath>

#include <deque>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/score.h>
//#include <seqan/modifier.h>
//#include <seqan/graph_types.h>
//#include <seqan/graph_algorithms.h>
//#include <seqan/map.h>
//#include <seqan/find.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/find2/find2_generated_forwards.h>
#endif

#include <seqan/find2/find_base.h>
#include <seqan/find2/find_finder_default.h>

// Exact pattern matching.
#include <seqan/find2/find_exact_shiftand.h>
#include <seqan/find2/find_exact_simple.h>

// Complex pattern matching.
#include <seqan/find2/find_pattern_wild_shiftand.h>

// Multiple exact pattern search.
#include <seqan/find2/find_multiple_exact_simple.h>

// Approximate pattern matching with Hamming distance.
#include <seqan/find2/find_hamming_simple.h>

// Approximate matching with linear/affine gap costs, edit distance etc.
#include <seqan/find2/find_approx_find_begin.h>  // findBegin() support
#include <seqan/find2/find_approx_dpsearch.h>

#endif  // SEQAN_SEQAN_FIND_H_
