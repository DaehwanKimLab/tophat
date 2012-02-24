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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Approximate string matching with linear and affine gap costs based
// on Sellers' (Sellers, 1980) algorithm based on the Needleman-Wunsch
// (Needleman and Wunsch, 1970) algorithm.  Described in the SeqAn book
// Section 9.3.1, p.154ff.
//
//   Needleman, S. B. and C. D. Wunsch (1970).  A general method
//     applicable to the search for similarities in the amino acid
//     sequence of two proteins. J. Molecular Biol. 48, 443-453.
//   Sellers, P. H. (1980).  The theory and computations of evolutionary
//     distances: Pattern recognition. Journal of Algorithms 1, 359-373.
// ==========================================================================

#ifndef SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_
#define SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_

// TODO(holtgrew): Document accessor for scoring scheme on whitepaper?

namespace seqan {

struct FindInfix_;
typedef Tag<FindInfix_> FindInfix;


struct FindPrefix_;
typedef Tag<FindPrefix_> FindPrefix;


template <typename TScore, typename TSpec = FindInfix, typename TFindBeginPatternSpec = True>
struct DPSearch;


// The simple dynamic programming implementation following Sellers'
// algorithm.
//
// TODO(holtgrew): Rename TScore to TScoringScheme?  Also: Rename Score class to ScoringScheme?
// TODO(holtgrew): Support affine gap costs, documented in book, code does not allow this.
// TODO(holtgrew): Add FindPrefix support.
template <typename TNeedle_, typename TScore, typename TSpec, typename TSupportFindBegin>
struct Pattern<TNeedle_, DPSearch<TScore, TSpec, TSupportFindBegin> >
        : public ApproxFindBegin_<TNeedle_, TScore, TSupportFindBegin> {
    typedef Pattern<TNeedle_, DPSearch<TScore, TSpec, TSupportFindBegin> > TPattern;
    typedef TNeedle_ TNeedle;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TMatrixColumn;

    // The pattern's state.
    FindState_::TState _state;

    // The needle we work on.
    Holder<TNeedle> _host;

    // The minimal score of a match.
    TScoreValue _scoreLimit;

    // The scoring scheme to use.
    TScore _scoringScheme;

    // The current score of a match, in book: v.
    TScoreValue _currentScore;

    // The matrix column from the DP search algorithm.
    TMatrixColumn _matrixColumn;

    Pattern() : _state(FindState_::STATE_EMPTY) { SEQAN_CHECKPOINT; }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, int scoreLimit)
        : _state(FindState_::STATE_INITIAL), _host(needle), _scoreLimit(scoreLimit) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_EQ_MSG(scoreGapOpen(_scoringScheme), scoreGapExtend(_scoringScheme), "The DPSearch Pattern only works with linear gap costs.");
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 & needle, int scoreLimit, TScore scoringScheme)
        : _state(FindState_::STATE_INITIAL), _host(needle), _scoreLimit(scoreLimit),
          _scoringScheme(scoringScheme) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_EQ_MSG(scoreGapOpen(_scoringScheme), scoreGapExtend(_scoringScheme), "The DPSearch Pattern only works with linear gap costs.");
    }
};


template <typename TNeedle, typename TScore, typename TSpec, typename TSupportFindBegin>
struct ScoringScheme<Pattern<TNeedle, DPSearch<TScore, TSpec, TSupportFindBegin> > > {
	typedef TScore Type;
};


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
TScore const & scoringScheme(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoringScheme;
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Value<TScore>::Type const & scoreLimit(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    return pattern._scoreLimit;
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Value<TScore>::Type const & score(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    // State of pattern should be in "found", "found begin" or "found
    // no begin" state.
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_NOTFOUND);
    return pattern._currentScore;
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
TNeedle const & host(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
TNeedle const & host(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
TNeedle const & needle(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    return host(pattern);
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Position<TNeedle>::Type length(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND || pattern._state == TPattern::STATE_FOUND);
    return length(needle(pattern));
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND || pattern._state == TPattern::STATE_FOUND);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND || pattern._state == TPattern::STATE_FOUND);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND || pattern._state == TPattern::STATE_FOUND);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND || pattern._state == TPattern::STATE_FOUND);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
inline void _initPattern(Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
	typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TPattern>::Type TSize;
	typedef typename TPattern::TMatrixColumn TMatrixColumn;
	typedef typename Iterator<TMatrixColumn, Standard>::Type TIterator;
	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;

    // Get some shortcuts.
	TScore const & myScoringScheme = scoringScheme(pattern);
	TScoreValue const gapScore = scoreGap(myScoringScheme);
	TMatrixColumn & matrixColumn = pattern._matrixColumn;

	// Allocate memory for the matrix column and fill it with the
    // scores for aligning with an end position of 0 and only deletes.
	resize(matrixColumn, length(needle(pattern)));
    TScoreValue s = gapScore;
    for (TIterator it = begin(matrixColumn, Standard()); it != end(matrixColumn, Standard()); ++it) {
        *it = s;
        s += gapScore;
    }
    SEQAN_ASSERT_EQ(static_cast<TScoreValue>(length(needle(pattern)) + 1) * gapScore, s);
}


template <typename THaystack, typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    typedef typename Position<TNeedle>::Type TPosition;
	typedef typename Value<TScore>::Type TScoreValue;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder and pattern if state is "initial".  Otherwise
    // advance at least by one (if not set to of haystack with
    // setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._beginPosition = 0u;
        finder._endPosition = 0u;
        _initPattern(pattern);
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        if (finder._endPosition == length(haystack(finder)))
            return false;
    }

    // The actual find implementation follows.
    //
    // Get some shortcuts.
    TScore const & myScore = scoringScheme(pattern);

    // Search for the next match.
    for (TPosition j = finder._endPosition; j < length(haystack(finder)); ++j) {
        // Compute best score if the end of needle would align at j.
        // pattern._currentScore is the "v" from the book.  This
        // initialization depends on whether we are doing a prefix
        // search or an infix search.
        TScoreValue d = IsSameType<TDPSearchSpec, FindInfix>::VALUE ? 0 : (j * scoreGap(myScore));
        pattern._currentScore = IsSameType<TDPSearchSpec, FindInfix>::VALUE ? 0 : ((j + 1) * scoreGap(myScore));
        for (TPosition i = 0; i < length(pattern._matrixColumn); ++i) {
            TScoreValue h = pattern._matrixColumn[i];
            pattern._currentScore = _max(d + score(myScore, haystack(finder)[j], needle(pattern)[i]),
                                         _max(pattern._currentScore, h) + scoreGap(myScore));
            pattern._matrixColumn[i] = pattern._currentScore;
            d = h;
        }
        if (pattern._currentScore >= scoreLimit(pattern)) {
            // Found a match, update the state and report the match.
            finder._state = TFinder::STATE_FOUND;
            pattern._state = TPattern::STATE_FOUND;
            // Update finder._endPosition to end instead of "last" position.
            finder._endPosition = j + 1;
            return true;
        }
    }
    // No match found, update state accordingly.
    finder._state = TFinder::STATE_NOTFOUND;
    pattern._state = TPattern::STATE_NOTFOUND;
	return false;
}


// findBeginScore() only works if the Pattern object has the necessary
// state information via ApproxFindBegin_<..., True>.
template <typename TNeedle, typename TScore>
typename Value<typename ScoringScheme<Pattern<TNeedle, DPSearch<TScore, FindInfix, True> > >::Type>::Type
findBeginScore(Pattern<TNeedle, DPSearch<TScore, FindInfix, True> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, FindInfix, True> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND);
    return pattern._findBeginCurrentScore;
}


// findBeginScore() also works if we are doing a prefix search.  In
// this case, findBegin() is superflous.
template <typename TNeedle, typename TScore, typename TFindBeginSpec>
typename Value<typename ScoringScheme<Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginSpec> > >::Type>::Type
findBeginScore(Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginSpec> > const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginSpec> > TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_BEGIN_FOUND);
    return score(pattern);
}


// findBegin() only works if the Pattern object has the necessary
// state information via ApproxFindBegin_<..., True>
template <typename THaystack, typename TNeedle, typename TScore, typename TFindBeginScoreLimitTag>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, DPSearch<TScore, FindInfix, True> > & pattern,
               TFindBeginScoreLimitTag const & tag) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, DPSearch<TScore, FindInfix, True> > TPattern;
    // State of finder and pattern should be in sync and in "found" or
    // "found begin" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    // Initialize findBegin() data if necessary.
    if (pattern._state == TPattern::STATE_FOUND)
        _initFindBegin(pattern, length(needle(pattern)), scoringScheme(pattern), scoreLimit(pattern));
    bool res = _findBeginImpl(pattern, scoringScheme(pattern), score(pattern), endPosition(finder), haystack(finder), needle(pattern), tag);
    if (res) {
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        finder._beginPosition = finder._endPosition - pattern._findBeginMatchLength;
    }
    return res;
}


// findBegin() also works if we are doing a prefix search only.  In
// this case, an actual search is superflous and we can always return
// true after updating the state.
template <typename THaystack, typename TNeedle, typename TScore, typename TFindBeginSpec, typename TFindBeginScoreLimitTag>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginSpec> > & pattern,
               TFindBeginScoreLimitTag const &) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginSpec> > TPattern;
    // State of finder and pattern should be in sync and in "found" or
    // "found begin" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    finder._beginPosition = 0;
    if (pattern._state == TPattern::STATE_FOUND) {
        finder._state = TFinder::STATE_BEGIN_FOUND;
        pattern._state = TPattern::STATE_BEGIN_FOUND;
        return true;
    } else {
        finder._state = TFinder::STATE_BEGIN_NOTFOUND;
        pattern._state = TPattern::STATE_BEGIN_NOTFOUND;
        return false;
    }
}


// If not given,t he find begin score limit tag is FindBeginScore and
// the score from the find() search is used for the score limit of the
// findBegin() search.
template <typename THaystack, typename TNeedle, typename TScore, typename TDPSearchSpec>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, True> > & pattern) {
    SEQAN_CHECKPOINT;
    return findBegin(finder, pattern, UseScore());
}


// The execution of this function is costly, O(length(needle)) for most realistic scoring schemes.
template <typename THaystack, typename TNeedle, typename TScore, typename TDPSearchSpec, typename TSupportFindBegin>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > & pattern,
                    typename Position<THaystack>::Type const & pos) {
    SEQAN_CHECKPOINT;
    typedef typename Position<THaystack>::Type TPosition;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, TSupportFindBegin> > TPattern;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Iterator<TNeedle, Standard>::Type TIterator;
    SEQAN_ASSERT_LEQ_MSG(pos, length(haystack(finder)), "End position must be valid.");
    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(pattern._scoringScheme), scoreGapExtend(pattern._scoringScheme), "The DPSearch Pattern only works with linear gap costs.");

    // Get some shortcuts.
    TScore const & myScore = scoringScheme(pattern);

    // In order to correctly initialize the matrix column, we have to
    // start far enough left and perform some DP steps.  The exact
    // number of steps depends on the scoring scheme.
    //
    // We make the assumption that the highest score is achieved if
    // all entries of the needle are matched with themselves.  This
    // score is taken as the base number for the computation.  We then
    // compute how many indels X are required to make the score fall
    // below the score limit.  We then start the DP search
    // length(needle) + X positions left of position and search until
    // the pattern's end aligns at pos.
    //
    // Start with computing the "base score" and getting the gap score.
    TScoreValue baseScore = -pattern._scoreLimit;
    for (TIterator it = begin(needle(pattern), Standard()); it != end(needle(pattern), Standard()); ++it)
        baseScore += score(myScore, *it, *it);
    TScoreValue gapScore = scoreGap(myScore);
    // In the worst case, the search start position is 0.  This is the
    // case if the gap penalty is non-negative or we are doing a
    // prefix search.  Otherwise, compute the search start position.
    TPosition searchStartPosition = 0;
    if (!IsSameType<TDPSearchSpec, FindPrefix>::VALUE) {
        double frac = (gapScore == 0) ? 1.0 : (1.0 * baseScore / gapScore);
        if (frac <= 0.0) {
            TPosition delta = length(needle(pattern)) - static_cast<TPosition>(frac);
            if (delta < pos)
                searchStartPosition = pos - delta;
        }
    }

    // Re-initialize the pattern.
    _initPattern(pattern);

    // Perform the DP search until the pattern aligns at pos.
    //
    // Search for the next match.
    TPosition j;
    for (j = searchStartPosition; j < pos; ++j) {
        // Compute best score if the end of needle would align at j.
        // pattern._currentScore is the "v" from the book.
        TScoreValue d = IsSameType<TDPSearchSpec, FindInfix>::VALUE ? 0 : (j * scoreGap(myScore));
        pattern._currentScore = IsSameType<TDPSearchSpec, FindInfix>::VALUE ? 0 : ((j + 1) * scoreGap(myScore));
        for (TPosition i = 0; i < length(pattern._matrixColumn); ++i) {
            TScoreValue h = pattern._matrixColumn[i];
            pattern._currentScore = _max(d + score(myScore, haystack(finder)[j], needle(pattern)[i]),
                                         _max(pattern._currentScore, h) + scoreGap(myScore));
            pattern._matrixColumn[i] = pattern._currentScore;
            d = h;
        }
    }
    // Update finder._endPosition to end instead of "last" position.
    finder._endPosition = j;
    // Update the state and return value depending on whether pos is a hit.
    if (pattern._currentScore >= scoreLimit(pattern)) {
        // Found a match, update the state and report the match.
        finder._state = TFinder::STATE_FOUND;
        pattern._state = TPattern::STATE_FOUND;
        return true;
    } else {
        // Found no match.
        finder._state = TFinder::STATE_NO_HIT;
        pattern._state = TPattern::STATE_NO_HIT;
        return false;
    }
}


// Build the alignment resulting from the search result as specified by the
// finder and the pattern.  If the state is not "begin found" then no alignment
// is built and false is returned.
//
// The alignment can only be built if the Pattern has the necessary state
// information.  This function is pretty costly since we internally use the
// Needleman-Wunsch alignment algorithm for generating the alignment between
// the computed begin and end positions in pattern and haystack.
template <typename THaystack, typename TNeedle, typename TScore, typename TDPSearchSpec, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, DPSearch<TScore, TDPSearchSpec, True> > &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, Simple> TPattern;
    typedef Align<TAlignSeq, TAlignSpec> TAlign;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Position<THaystack>::Type TPosition;

    // Both finder and pattern must be in the "found begin position" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // Can only build alignment if the state is "begin found".
    if (finder._state != TFinder::STATE_BEGIN_FOUND)
        return false;

    // Initialize alignment with the two sequences.
    resize(rows(outAlignment), 2);
    assignSource(row(outAlignment, 0), haystack(finder));
    assignSource(row(outAlignment, 1), needle(pattern));
    // Now shift the needle's beginning to the begin position.
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    // Build the alignment with the matching infixes of finder and pattern.
    typedef Align<Segment<THaystack const, InfixSegment> > TTempAlign;
    TTempAlign tempAlignment;
    resize(rows(tempAlignment), 2);
    assignSource(row(tempAlignment, 0), infix(finder));
    assignSource(row(tempAlignment, 1), infix(pattern));
    globalAlignment(tempAlignment, scoringScheme(pattern));
//     std::cout << tempAlignment;

    // Copy over the gaps from the alignment.
    typedef typename Cols<TTempAlign>::Type TColumns;
    typedef typename Iterator<TColumns>::Type TIterator;
    TPosition i = beginPosition(finder);
//     std::cout << "length(cols(tempAlignment)) = " << length(cols(tempAlignment)) << std::endl;
    // TODO(holtgrew): This does not work as expected.
    for (TPosition j = 0; j < length(cols(tempAlignment)); ++j) {
        if (value(col(tempAlignment, j), 0) == '-') {
//             std::cout << "j = " << j << ", insertGaps(row(outAlignment, 0), " << i << ", 1)" << std::endl;
            insertGaps(row(outAlignment, 0), i, 1);
            i += 1;
        } else if (value(col(tempAlignment, j), 1) == '-') {
//             std::cout << "j = " << j << ", insertGaps(row(outAlignment, 1), " << i << ", 1)" << std::endl;
            insertGaps(row(outAlignment, 1), i, 1);
            i += 1;
        }
        i += 1;
    }
//     std::cout << tempAlignment;

    return true;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_APPROX_DPSEARCH_H_
