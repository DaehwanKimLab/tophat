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
// Implementation of the default finder.
// ==========================================================================

#ifndef SEQAN_FIND2_FINDER_DEFAULT_H_
#define SEQAN_FIND2_FINDER_DEFAULT_H_

namespace seqan {

template <typename THaystack_>
struct Finder<THaystack_, Default> : FindState_ {
    typedef THaystack_ THaystack;
    typedef typename Position<THaystack>::Type TPosition;
    typedef typename Iterator<THaystack>::Type TIterator;

    TState _state;

    Holder<THaystack> _holder;

    // TODO(holtgrew): Maybe switch to an iterator based implementation?
    TPosition _beginPosition;
    TPosition _endPosition;

    Finder(THaystack & haystack)
        : _state(STATE_INITIAL),
          _holder(haystack) {
        SEQAN_CHECKPOINT;
    }
};


template <typename THaystack>
THaystack const & haystack(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    return value(finder._holder);
}


template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type begin(Finder<THaystack, Default> const & finder,
                                               Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return begin(haystack(finder), spec) + finder._beginPosition;
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type begin(Finder<THaystack, Default> & finder,
                                               Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return begin(const_cast<TFinder const &>(finder), spec);
}


template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type end(Finder<THaystack, Default> const & finder,
                                             Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return begin(haystack(finder), spec) + finder._endPosition;
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename THaystack, typename TTag>
typename Iterator<THaystack const>::Type end(Finder<THaystack, Default> & finder,
                                             Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return end(const_cast<TFinder const &>(finder), spec);
}


template <typename THaystack>
typename Position<THaystack const>::Type beginPosition(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return finder._beginPosition;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename THaystack>
typename Position<THaystack>::Type beginPosition(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return beginPosition(const_cast<TFinder const &>(finder));
}


template <typename THaystack>
typename Position<THaystack>::Type endPosition(Finder<THaystack, Default> const & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND || finder._state == TFinder::STATE_FOUND);
    return finder._endPosition;
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename THaystack>
typename Position<THaystack>::Type endPosition(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    return endPosition(const_cast<TFinder const &>(finder));
}


template <typename THaystack>
Segment<THaystack const, InfixSegment> infix(Finder<THaystack, Default> & finder) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    SEQAN_ASSERT_TRUE(finder._state == TFinder::STATE_BEGIN_FOUND);
    return infix(haystack(finder), beginPosition(finder), endPosition(finder));
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FINDER_DEFAULT_H_
