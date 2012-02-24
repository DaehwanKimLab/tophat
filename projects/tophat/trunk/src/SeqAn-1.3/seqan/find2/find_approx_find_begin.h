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
// A simple implementation of the algorithms and datastructures required for
// findBegin().  It is a "backwards" implementation of Sellers' algorithm (see
// find_approx_dpsearch.h for references) for approximate prefix search.
//
// In order to support findBegin(), the approximate DP-based algorithms need
// to store more information than only for the forward search.  If findBegin()
// does not need to be supported, this information does not have to be stored.
// For this reason, the functionality is centralized in the struct template
// ApproxFindBegin_ that is to be used as the base class for the approximate
// search patterns.
//
// If the parameter HasFindBeginSupport is set to True then the class gets the
// information and functionality for findBegin(), if it is set to False, it
// does not get this information.
// ============================================================================

#ifndef SEQAN_FIND2_FIND_APPROX_DPSEARCH_
#define SEQAN_FIND2_FIND_APPROX_DPSEARCH_

namespace seqan {

// TODO(holtgrew): Add support for affine gap costs.
template <typename TNeedle, typename TScore, typename HasFindBeginSupport>
struct ApproxFindBegin_;


// Implementation for no findBegin() support.
template <typename TNeedle, typename TScore>
struct ApproxFindBegin_<TNeedle, TScore, False> {};


// Implementation for findBegin() support.
template <typename TNeedle, typename TScore>
struct ApproxFindBegin_<TNeedle, TScore, True> : public FindState_ {
    // We use the position of the needle for the finder, too.  This is
    // the best we can do with the split-up interface.
    typedef typename Position<TNeedle>::Type TPosition;
    // We will use a matrix column of score values for the dynamic
    // programming and types for this.
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TMatrixColumn;

    // The length of the match for findBegin().
    TPosition _findBeginMatchLength;

    // The score limit when searching backwards.
    int _findBeginScoreLimit;

    // The current score of the findBegin() call.
    int _findBeginCurrentScore;

    // The matrix column for the findBegin computation.
    TMatrixColumn _findBeginMatrixColumn;

    // The state of the begin search.
    TState _findBeginState;
};


// Tag for findBegin().  Use score of forward search result when doing
// find begin.
struct UseScore_;
typedef Tag<UseScore_> UseScore;


// Tag for findBegin().  Use score limit of forward search when
// doing find begin.
struct UseScoreLimit_;
typedef Tag<UseScoreLimit_> UseScoreLimit;


// Called to initialize a ApproxFindBegin_ datastructure.
template <typename TNeedle, typename TScore>
void _initFindBegin(ApproxFindBegin_<TNeedle, TScore, True> & findBeginStruct,
                    typename Position<TNeedle>::Type const & needleLength, TScore const & scoringScheme, int scoreLimit) {
    SEQAN_CHECKPOINT;
    typedef ApproxFindBegin_<TNeedle, TScore, True> TApproxFindBegin;
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TMatrixColumn;
    typedef typename Iterator<TMatrixColumn>::Type TIterator;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapOpen(scoringScheme), "findBegin() only supports linear gap costs at the moment.");

    // Assumption: The end position must have been found already!  The
    // caller has to ensure this.
    findBeginStruct._findBeginState = TApproxFindBegin::STATE_FOUND;
    findBeginStruct._findBeginMatchLength = 0;
    findBeginStruct._findBeginScoreLimit = scoreLimit;
    // Re-initialize the matrix column.
    resize(findBeginStruct._findBeginMatrixColumn, needleLength);
    TScoreValue gapValue = scoreGap(scoringScheme);
    for (TIterator it = begin(findBeginStruct._findBeginMatrixColumn, Standard()); it != end(findBeginStruct._findBeginMatrixColumn, Standard()); ++it) {
        *it = gapValue;
        gapValue += scoreGap(scoringScheme);
    }
}


template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type _findBeginScoreLimit(ApproxFindBegin_<TNeedle, TScore, True> const &,
                                                  typename Value<TScore>::Type const & findScore,
                                                  typename Value<TScore>::Type const &,
                                                  UseScore const &) {
    SEQAN_CHECKPOINT;
    return findScore;
}


template <typename TNeedle, typename TScore>
inline typename Value<TScore>::Type _findBeginScoreLimit(ApproxFindBegin_<TNeedle, TScore, True> const &,
                                                  typename Value<TScore>::Type const &,
                                                  typename Value<TScore>::Type const & findScoreLimit,
                                                  UseScoreLimit const &) {
    SEQAN_CHECKPOINT;
    return findScoreLimit;
}


// Actual implementation of findBegin() through ApproxFindBegin_.
template <typename TNeedle, typename TScore, typename THaystack, typename TTag>
bool _findBeginImpl(ApproxFindBegin_<TNeedle, TScore, True> & findBeginStruct,
                    TScore const & scoringScheme,
                    typename Value<TScore>::Type const & findScore,
                    typename Position<TNeedle>::Type const & endPosition,
                    THaystack const & haystack, TNeedle const & needle,
                    Tag<TTag> const & tag) {
    SEQAN_CHECKPOINT;
    typedef ApproxFindBegin_<TNeedle, TScore, True> TApproxFindBegin;
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TScore>::Type TScoreValue;
    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme), "findBegin() only supports linear gap costs at the moment.");

    // Finding the start of the match follows Seller's algorithm but
    // works from the right to the left.

    // TODO(holtgrew): Remove all commented out debug output code from this function.
//     std::cout << "== start of findBegin() == " << std::endl;
//     for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//         std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//     std::cout << std::endl;
//     std::cout << "=====" << std::endl;

    TScoreValue findBeginScoreLimit = _findBeginScoreLimit(findBeginStruct, findScore, findBeginStruct._findBeginScoreLimit, tag);

    // Search for the next match.
    for (TPosition & j = findBeginStruct._findBeginMatchLength; j < endPosition; ++j) {
        // Compute best score if the begin of the match was j position
        // left of the end.
        findBeginStruct._findBeginCurrentScore = (j + 1) * scoreGap(scoringScheme);  // is v in SeqAn book
        TScoreValue d = (j) * scoreGap(scoringScheme);
        for (TPosition i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i) {
            TScoreValue h = findBeginStruct._findBeginMatrixColumn[i];
            findBeginStruct._findBeginCurrentScore = _max(d + score(scoringScheme, haystack[endPosition - j - 1], needle[length(needle) - 1 - i]),
                                                          _max(findBeginStruct._findBeginCurrentScore, h) + scoreGap(scoringScheme));
            findBeginStruct._findBeginMatrixColumn[i] = findBeginStruct._findBeginCurrentScore;
            d = h;
        }
//         std::cout << "== within of findBegin() == " << std::endl;
//         for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//             std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//         std::cout << std::endl;
//         std::cout << "=====" << std::endl;
        if (findBeginStruct._findBeginCurrentScore >= findBeginScoreLimit) {
//             std::cout << "== end of findBegin() == " << std::endl;
//             for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//                 std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//             std::cout << std::endl;
//             std::cout << "=====" << std::endl;
            findBeginStruct._findBeginMatchLength += 1;
            // Found a match, update the state and report the match.
            findBeginStruct._findBeginState = TApproxFindBegin::STATE_BEGIN_FOUND;
            return true;
        }
    }

//     std::cout << "== end of findBegin() == " << std::endl;
//     for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//         std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//     std::cout << std::endl;
//     std::cout << "=====" << std::endl;

    // No match found, update state accordingly.
    findBeginStruct._findBeginState = TApproxFindBegin::STATE_BEGIN_NOTFOUND;
	return false;
}

}

#endif  // SEQAN_FIND2_FIND_APPROX_DPSEARCH_
