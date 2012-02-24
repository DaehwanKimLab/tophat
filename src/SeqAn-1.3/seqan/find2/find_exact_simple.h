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
// Exact pattern matching using a naive implementation.
// ==========================================================================

#ifndef SEQAN_FIND2_FIND_EXACT_SIMPLE_H_
#define SEQAN_FIND2_FIND_EXACT_SIMPLE_H_

namespace seqan {

template <typename TNeedle_>
struct Pattern<TNeedle_, Simple> : FindState_ {
    typedef TNeedle_ TNeedle;

    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;

    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    explicit
    Pattern(TNeedle & ndl)
        : _state(STATE_INITIAL),
          _host(ndl) {
        SEQAN_CHECKPOINT;
    }
};


template <typename TNeedle>
struct Needle<Pattern<TNeedle, Simple> > {
    typedef TNeedle Type;
};


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, Simple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, Simple> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, Simple> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, Simple> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, Simple> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._beginPosition = 0u;
        finder._endPosition = length(needle(pattern));
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder)))
            return false;
        finder._beginPosition += 1;
    } else {
        finder._beginPosition += 1;
    }

    // Search the needle in the haystack naively.
    for (TPosition i = 0u; i < length(needle(pattern));) {
        // Break out of loop if no more match is possible.
        if (finder._beginPosition >= length(haystack(finder)) - length(needle(pattern))) {
            finder._state = TFinder::STATE_NOTFOUND;
            pattern._state = TPattern::STATE_NOTFOUND;
            return false;
        }
        // Otherwise, go on searching.
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._beginPosition += 1;
            i = 0u;
            continue;
        }
        i += 1;
    }
    finder._endPosition = finder._beginPosition + length(needle(pattern));
    finder._state = TFinder::STATE_FOUND;
    pattern._state = TPattern::STATE_FOUND;
    return true;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, Simple> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    // State of finder and pattern should be in sync and in "found" or
    // "found begin" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    if (pattern._state == TPattern::STATE_FOUND) {
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    }
    finder._state = TFinder::STATE_BEGIN_NOTFOUND;
    pattern._state = TPattern::STATE_BEGIN_NOTFOUND;
    return false;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, Simple> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // End position must not be right of the end of the haystack.
    SEQAN_ASSERT_LEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(haystack(finder)));
    // Begin position must not be left of the beginning of the haystack.
    SEQAN_ASSERT_GEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(needle(pattern)));

    // Set the end position.
    finder._endPosition = pos;
    finder._beginPosition = pos - length(needle(pattern));

    // Check whether there is a hit at this position and update the
    // state accordingly.
    typedef typename Position<THaystack>::Type THaystackPos;
    for (THaystackPos i = 0u; i < length(needle(pattern)); ++i) {
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._state = TPattern::STATE_NO_HIT;
            pattern._state = TFinder::STATE_NO_HIT;
            return false;
        }
    }
    finder._state = TPattern::STATE_FOUND;
    pattern._state = TPattern::STATE_FOUND;
    return true;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setBeginPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, Simple> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    return setEndPosition(finder, pattern, pos + length(needle(pattern)));
}


/*
  Build the alignment resulting from the search result as specified by the
  finder and the pattern.  If the state is not "begin found" then no alignment
  is built and false is returned.
*/
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, Simple> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    typedef Align<TAlignSeq, TAlignSpec> TAlign;
    typedef typename Row<TAlign>::Type TRow;

    // Both finder and pattern must be in the "found begin position" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // Can only build alignment if the state is "begin found".
    if (finder._state != TFinder::STATE_BEGIN_FOUND)
        return false;

    // Initialize alignment with the two sequences.
    resize(rows(outAlignment), 2);
    assignSource(row(outAlignment, 0), haystack(finder));
    assignSource(row(outAlignment, 1), needle(pattern));
    // Insert gap into the needle.
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_EXACT_SIMPLE_H_
