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

#ifndef SEQAN_HEADER_INDEX_H
#define SEQAN_HEADER_INDEX_H

//____________________________________________________________________________
// prerequisites

#include <seqan/sequence.h>
#include <seqan/pipe.h>
#include <seqan/modifier.h>

#include <seqan/find.h>
#include <seqan/misc/misc_set.h>

#include <climits>
#include <functional>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>
#include <string.h> // memset


//////////////////////////////////////////////////////////////////////////////
// INDEX CONSTRUCTION
//////////////////////////////////////////////////////////////////////////////


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/index/index_manual_forwards.h>
#include <seqan/index/index_generated_forwards.h>
#endif

#include <seqan/index/index_base.h>

//____________________________________________________________________________
// suffix array creators

#include <seqan/index/radix.h>
#include <seqan/index/index_sa_btree.h>
#include <seqan/index/index_sa_lss.h>
#include <seqan/index/index_sa_mm.h>
#include <seqan/index/index_sa_qsort.h>
#include <seqan/index/index_sa_bwtwalk.h>

#include <seqan/index/pump_extender3.h>
#include <seqan/index/pipe_merger3.h>
#include <seqan/index/index_skew3.h>

#include <seqan/index/pump_extender7.h>
#include <seqan/index/pipe_merger7.h>
#include <seqan/index/index_skew7.h>

#include <seqan/index/pump_separator7.h>
#include <seqan/index/index_skew7_multi.h>

//____________________________________________________________________________
// enhanced table creators

#include <seqan/index/pump_lcp_core.h>
#include <seqan/index/index_lcp.h>
#include <seqan/index/index_lcp_tree.h>

#include <seqan/index/index_childtab.h>
#include <seqan/index/index_bwt.h>

//____________________________________________________________________________
// q-gram index creator

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>
#include <seqan/index/shape_onegapped.h>
#include <seqan/index/shape_predefined.h>
#include <seqan/index/shape_threshold.h>
#include <seqan/index/index_qgram.h>
#include <seqan/index/index_qgram_openaddressing.h>
//#include <seqan/index/index_qgram_nested.h>


//////////////////////////////////////////////////////////////////////////////
// INDEX USAGE
//////////////////////////////////////////////////////////////////////////////

#include <seqan/index/index_shims.h>

//____________________________________________________________________________
// (virtual) suffix trees

#include <seqan/index/index_esa_base.h>
#include <seqan/index/index_esa_stree.h>
#include <seqan/index/index_wotd.h>
#include <seqan/index/index_dfi.h>

//____________________________________________________________________________
// suffix tree algorithms

#include <seqan/index/index_esa_algs.h>
#include <seqan/index/index_esa_algs_multi.h>
#include <seqan/index/index_esa_drawing.h>
#include <seqan/index/repeat_base.h>

//____________________________________________________________________________
// Pizza & Chili interface (compressed indices)

#include <seqan/index/index_pizzachili.h>
#include <seqan/index/index_pizzachili_find.h>

//____________________________________________________________________________
// Shawarma interface (suffix array creators)

#include <seqan/index/index_shawarma.h>


//////////////////////////////////////////////////////////////////////////////
// FINDER INTERFACE
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
// index based finders

#include <seqan/index/find_index_esa.h>
#include <seqan/index/find_index_approx.h>
#include <seqan/index/find_index_qgram.h>
#include <seqan/index/find_index.h>
#include <seqan/index/find_quasar.h>
#include <seqan/index/find_swift.h>

#endif //#ifndef SEQAN_HEADER_...
