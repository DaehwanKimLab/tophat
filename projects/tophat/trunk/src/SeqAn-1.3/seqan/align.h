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

#ifndef SEQAN_ALIGN_H_
#define SEQAN_ALIGN_H_

// ----------------------------------------------------------------------
// Prerequisites.

#include <cmath>
#include <stack>
#include <ostream>

#include <seqan/sequence.h>
#include <seqan/score.h>

#include <seqan/map.h>
#include <seqan/graph_align.h>
#include <seqan/modifier.h>

// ----------------------------------------------------------------------
// The module's headers.

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/align/align_generated_forwards.h>
#endif

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>

#include <seqan/align/gaps_array.h>
#include <seqan/align/gaps_sumlist.h>

#include <seqan/align/align_base.h>
#include <seqan/align/align_cols_base.h>
#include <seqan/align/align_iterator_base.h>

#include <seqan/align/matrix_base.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#include <seqan/align/align_algorithms.h>
#include <seqan/align/align_dynprog.h>
#include <seqan/align/align_local_dynprog.h>
#include <seqan/align/align_local_dynprog_banded.h>
#include <seqan/align/hirschberg_set.h>
#include <seqan/align/align_myers.h>
#include <seqan/align/align_hirschberg.h>

#endif  // SEQAN_ALIGN_H_
