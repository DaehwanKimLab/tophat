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
// Exact pattern matching using the Shift-And algorithm.
// ==========================================================================

#ifndef SEQAN_FIND2_FIND_EXACT_SHIFTAND_H_
#define SEQAN_FIND2_FIND_EXACT_SHIFTAND_H_

// TODO(holtgrew): Try to achieve a speedup using the LongWord class.

namespace seqan {

struct ShiftAnd_;
typedef Tag<ShiftAnd_> ShiftAnd;


template <typename TNeedle_>
struct Pattern<TNeedle_, ShiftAnd> : FindState_ {
    typedef TNeedle_ TNeedle;
    typedef unsigned TWord;

    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;
    
    String<TWord> _table;         // Lookup table for each character in the alphabet (called B in "Navarro")
    String<TWord> _prefSufMatch;  // Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
    TWord _blockCount;            // number of unsigneds required to store needle  


    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    explicit
    Pattern(TNeedle & ndl)
        : _state(STATE_INITIAL),
          _host(ndl) {
        SEQAN_CHECKPOINT;
        _initPattern(*this);
    }
};


template <typename TNeedle>
struct Needle<Pattern<TNeedle, ShiftAnd> > {
    typedef TNeedle Type;
};


template <typename TNeedle>
void _initPattern(Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
	typedef typename TPattern::TWord TWord;
	typedef typename Position<TNeedle>::Type TPosition;
	typedef typename Value<TNeedle>::Type TAlphabet;

    // Get some shortcuts.
    TNeedle const & ndl = needle(pattern);
    TPosition needleLength = length(ndl);
    // Compute number of blocks required for the bitmasks.
    pattern._blockCount = needleLength / BitsPerValue<TWord>::VALUE;
    if (needleLength % BitsPerValue<TWord>::VALUE > 0)
        pattern._blockCount += 1;
    // Resize and initialize the "match/mismatch" word.
	resize(pattern._prefSufMatch, pattern._blockCount, 0u, Exact());
    // Resize and initialize the bitmask table.
	resize(pattern._table, pattern._blockCount * ValueSize<TAlphabet>::VALUE, 0u, Exact());
	for (TWord j = 0; j < needleLength; ++j) {
		// Determine character position in array table.
		TWord pos = ordValue(getValue(ndl, j));
		pattern._table[pattern._blockCount * pos + j / BitsPerValue<TWord>::VALUE] |= (1 << (j % BitsPerValue<TWord>::VALUE));
	}
}


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, ShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, ShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, ShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, ShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, ShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    SEQAN_ASSERT_TRUE(pattern._state == TPattern::STATE_FOUND ||
                      pattern._state == TPattern::STATE_BEGIN_FOUND);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle>
inline bool _findShiftAndSmallNeedle(Finder<THaystack, Default> & finder,
                                       Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
	typedef typename Value<TNeedle>::Type TValue;
	typedef unsigned int TWord;
    typedef typename Position<TNeedle>::Type TPosition;
    THaystack const & hstck = haystack(finder);
    TPosition needleLength = length(needle(pattern));
    TPosition haystackLength = length(hstck);
	TWord compare = 1 << (needleLength - 1);
	while (finder._endPosition < haystackLength) {
		TWord pos = ordValue(hstck[finder._endPosition]);
		pattern._prefSufMatch[0] = ((pattern._prefSufMatch[0] << 1) | 1) & pattern._table[pattern._blockCount * pos];
		if ((pattern._prefSufMatch[0] & compare) != 0)
			return true;
        finder._endPosition += 1;
	}
	return false;
}


template <typename THaystack, typename TNeedle>
inline bool _findShiftAndLargeNeedle(Finder<THaystack, Default> & finder,
                                      Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
	typedef typename Value<TNeedle>::Type TValue;
	typedef unsigned int TWord;
    typedef typename Position<TNeedle>::Type TPosition;
    THaystack const & hstck = haystack(finder);
    TPosition needleLength = length(needle(pattern));
    TPosition haystackLength = length(hstck);
	TWord compare = 1 << ((needleLength - 1) % BitsPerValue<TWord>::VALUE);
	while (finder._endPosition < haystackLength) {
		TWord pos = ordValue(hstck[finder._endPosition]);
		TWord carry = 1;
		for (TWord block = 0; block < pattern._blockCount; ++block) {
			bool newCarry = (pattern._prefSufMatch[block] & (1 << (BitsPerValue<TWord>::VALUE - 1))) != 0;
			pattern._prefSufMatch[block] <<= 1;
			pattern._prefSufMatch[block] |= carry;
			carry = newCarry;
		}
		for (TWord block = 0; block < pattern._blockCount; ++block)
            pattern._prefSufMatch[block] &= pattern._table[pattern._blockCount * pos + block];
		if ((pattern._prefSufMatch[pattern._blockCount - 1] & compare) != 0)
			return true; 
        finder._endPosition += 1;
	}
	return false;
}


template <typename THaystack, typename TNeedle>
inline bool find(Finder<THaystack, Default> & finder,
                 Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._endPosition = 0;
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder))) {
            finder._state = TPattern::STATE_NOTFOUND;
            return false;
        }
    }

    bool res;
    if (pattern._blockCount == 1)
        res = _findShiftAndSmallNeedle(finder, pattern);
    else
        res = _findShiftAndLargeNeedle(finder, pattern);
    // Advance end position to make an "end" position from a "last" position.
    finder._endPosition += 1;
    if (res) {
        finder._beginPosition = finder._endPosition - length(needle(pattern));
        finder._state = TFinder::STATE_FOUND;
        pattern._state = TPattern::STATE_FOUND;
        return true;
    } else {
        finder._state = TFinder::STATE_NOTFOUND;
        pattern._state = TPattern::STATE_NOTFOUND;
        return false;
    }
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, ShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
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



template <typename THaystack, typename TNeedle>
inline bool _setEndPositionShiftAndSmallNeedle(Finder<THaystack, Default> & finder,
                                                 Pattern<TNeedle, ShiftAnd> & pattern,
                                                 typename Position<THaystack>::Type const & pos) {
    SEQAN_CHECKPOINT;
	typedef typename Value<TNeedle>::Type TValue;
    typedef typename Position<THaystack>::Type TPosition;
	typedef unsigned int TWord;
    THaystack const & hstck = haystack(finder);
    TPosition needleLength = length(needle(pattern));
	TWord compare = 1 << (needleLength - 1);
	while (finder._endPosition < pos) {
		TWord pos = ordValue(hstck[finder._endPosition]);
		pattern._prefSufMatch[0] = ((pattern._prefSufMatch[0] << 1) | 1) & pattern._table[pattern._blockCount * pos];
        finder._endPosition += 1;
	}
	return (pattern._prefSufMatch[0] & compare) != 0;
}


template <typename THaystack, typename TNeedle>
inline bool _setEndPositionShiftAndLargeNeedle(Finder<THaystack, Default> & finder,
                                                Pattern<TNeedle, ShiftAnd> & pattern,
                                                typename Position<THaystack>::Type const & pos) {
    SEQAN_CHECKPOINT;
	typedef typename Value<TNeedle>::Type TValue;
	typedef unsigned int TWord;
    typedef typename Position<TNeedle>::Type TPosition;
    THaystack const & hstck = haystack(finder);
    TPosition needleLength = length(needle(pattern));
	TWord compare = 1 << ((needleLength - 1) % BitsPerValue<TWord>::VALUE);
	while (finder._endPosition < pos) {
		TWord pos = ordValue(hstck[finder._endPosition]);
		TWord carry = 1;
		for (TWord block = 0; block < pattern._blockCount; ++block) {
			bool newCarry = (pattern._prefSufMatch[block] & (1 << (BitsPerValue<TWord>::VALUE - 1))) != 0;
			pattern._prefSufMatch[block] <<= 1;
			pattern._prefSufMatch[block] |= carry;
			carry = newCarry;
		}
		for (TWord block = 0; block < pattern._blockCount; ++block)
            pattern._prefSufMatch[block] &= pattern._table[pattern._blockCount * pos + block];
        finder._endPosition += 1;
	}
	return ((pattern._prefSufMatch[pattern._blockCount - 1] & compare) != 0);
}


// Sets end position and adjusts the bit pattern state.
template <typename THaystack, typename TNeedle>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, ShiftAnd> & pattern,
                    typename Position<THaystack>::Type const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    typedef typename Position<THaystack>::Type TPosition;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // End position must not be right of the end of the haystack.
    SEQAN_ASSERT_LEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(haystack(finder)));
    // Begin position must not be left of the beginning of the haystack.
    SEQAN_ASSERT_GEQ(static_cast<typename MakeUnsigned_<TPosition>::Type>(pos), length(needle(pattern)));

    // Set finder's end position if it is in the initial state.
    if (finder._state == TPattern::STATE_INITIAL)
        finder._endPosition = 0;

    // Set the end position to the start position and search from
    // there until pos is reached.
    if (finder._endPosition > pos)
        finder._endPosition = finder._endPosition - pos;
    else
        finder._endPosition = 0u;

    bool ret;
    if (pattern._blockCount == 1)
        ret = _setEndPositionShiftAndSmallNeedle(finder, pattern, pos);
    else
        ret = _setEndPositionShiftAndLargeNeedle(finder, pattern, pos);

    // Adjust state.
    if (ret) {
        finder._beginPosition = finder._endPosition - length(needle(pattern));
        finder._state = TPattern::STATE_FOUND;
        pattern._state = TPattern::STATE_FOUND;
        return true;
    } else {
        finder._state = TPattern::STATE_NO_HIT;
        pattern._state = TFinder::STATE_NO_HIT;
        return false;
    }
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setBeginPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, ShiftAnd> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
    return setEndPosition(finder, pattern, pos + length(needle(pattern)));
}

/*
  Build the alignment resulting from the search result as specified by the
  finder and the pattern.  If the state is not "begin found" then no alignment
  is built and false is returned.
*/
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, ShiftAnd> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, ShiftAnd> TPattern;
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
          
#endif  // SEQAN_FIND2_FIND_EXACT_SHIFTAND_H_
