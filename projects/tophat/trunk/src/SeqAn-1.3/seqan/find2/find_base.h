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
// Basic definitions and types for the find module.
// ==========================================================================

#ifndef SEQAN_FIND2_FIND_BASE_H_
#define SEQAN_FIND2_FIND_BASE_H_

namespace seqan {

// Contains the state for finder and patterns.
struct FindState_ {
    enum TState {
        STATE_EMPTY,           // Finder/pattern is empty.
        STATE_INITIAL,         // Finer/pattern has just been initialized.
        STATE_FOUND,           // Found the end position of a hit.
        STATE_NOTFOUND,        // No hit found, no more hits possible.
        STATE_BEGIN_FOUND,     // Found begin position.
        STATE_BEGIN_NOTFOUND,  // Found end but not begin, should not happen.
        STATE_NO_HIT           // Set manually to non-hit.
    };
};


template <typename TNeedle, typename TSpec>
struct Pattern;


template <typename THaystack, typename TSpec = Default>
struct Finder;


template <typename TPattern>
struct Needle;


template <typename TNeedle, typename TSpec>
struct Needle<Pattern<TNeedle, TSpec> > {
    typedef TNeedle Type;
};


// Metafunction for retrieving scoring schemes.
template <typename TPattern>
struct ScoringScheme;


// TODO(holtgrew): Implement HasFeature<TPattern, SetEndPosition> and HasFeature<TPattern, SetBeginPosition>, HasFeature<TPattern, ScoreLimit>, HasFeature<TPattern, Score>?

}  // namespace seqan

#endif  // SEQAN_FIND2_FIND_BASE_H_
