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
// Author: Stefan Aiche <aiche@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Wildcard pattern matching using a modification of the Shift-And Algorithm.
// ==========================================================================

#ifndef SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
#define SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_

// TODO(holtgrew): This algorithm yields end positions. What about the start position?

// Deactivate debugging for WildShiftAnd by default.
#ifndef SEQAN_WILD_SHIFTAND_DEBUG
#define SEQAN_WILD_SHIFTAND_DEBUG 0
#endif  // SEQAN_WILD_SHIFTAND_DEBUG

namespace seqan {

struct WildShiftAnd_;
typedef Tag<WildShiftAnd_> WildShiftAnd;


template <typename TNeedle_>
struct Pattern<TNeedle_, WildShiftAnd> : FindState_ {
    typedef TNeedle_ TNeedle;
	typedef unsigned TWord;

    // The pattern's state.
    TState _state;

    // Holder with the needle's data.
	Holder<TNeedle> data_host;

    // TODO(holtgrew): Adjust naming of members to the style guide.
	String<TWord> table;			// Look up table for each character in the alphabet (called B in "Navarro")
	
	String<TWord> s_table;			// marks all positions, that can remain active, after reading a specific character (called S in "Navarro")
	String<TWord> a_table;			// marks all positions of optional characters in the pattern (called A in "Navarro")
	String<TWord> i_table;			// marks all positions in the pattern, that preceed a block of optional characters (called I in "Navarro")
	String<TWord> f_table;			// marks all end-positions of blocks of optional characters (called F in "Navarro")
	
	String<TWord> prefSufMatch;		// Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
	String<TWord> df;				// additional bit mask to enable flooding of bits
	
	TWord needleLength;				// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	TWord character_count;			// number of normal characters in the needle
	TWord blockCount;				// number of unsigned ints required to store needle	

    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    explicit
    Pattern(TNeedle & ndl)
        : _state(STATE_INITIAL),
          data_host(ndl) {
        SEQAN_CHECKPOINT;
        _initializePattern(*this);
    }
};


template <typename TNeedle>
struct Needle<Pattern<TNeedle, WildShiftAnd> > {
    typedef typename Value<TNeedle>::Type Value;
};


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern.data_host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern.data_host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, WildShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, WildShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, WildShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, WildShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


// template <typename THaystack, typename TNeedle>
// bool findBegin(Finder<THaystack, Default> & finder,
//                Pattern<TNeedle, WildShiftAnd> & pattern) {
//     SEQAN_CHECKPOINT;
//     // State of finder and pattern should be in sync.
//     SEQAN_ASSERT_EQ(finder._state, pattern._state);
//     typedef Pattern<TNeedle, WildShiftAnd> TPattern;
//     return finder._state == TPattern::STATE_BEGIN_FOUND;
// }


template <typename THaystack, typename TNeedle, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, WildShiftAnd> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
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
    finder._state = TPattern::STATE_BEGIN_FOUND;
    pattern._state = TPattern::STATE_BEGIN_FOUND;
    return true;
}


// Build the alignment resulting from the search result as specified by the
// finder and the pattern.  If the state is not "begin found" then no alignment
// is built and false is returned.
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, WildShiftAnd> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
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
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}


// TODO(holtgrew): Should probably to into some utility header.
inline bool _isUnsigned(CharString const & number) {
    SEQAN_CHECKPOINT;
    if (length(number) == 0u)
        return false;

    typedef Position<CharString>::Type TPosition;
	for (TPosition i = 0; i < length(number); ++i) {
        if (number[i] > '9' || number[i] < '0')
            return false;
	}
	return true;
}


// The states used by the finite state machine that validates the regular
// expressions of the WildShiftAnd algorithm.
//
// These regular expressions consist of sequences of characters, character
// classes or the placeholder ".".  Each can be followed by a quantifier,
// where quantifiers are characters from the set {*, +, ?} or explicit numbers
// as in a{i, j}.
enum FindWildShiftAndParserStates_ {
    STATE_NO_CHAR,                    // Expecting character, character class or ".".
    STATE_CHAR,                       // Read character, character class or ".".
    STATE_ESCAPED,                    // Just read a backslash, next character is escaped.
    STATE_CHAR_CLASS,                 // In a character class, read at least one character.
    STATE_CHAR_CLASS_ESCAPED,         // Escape character in character class, a character must follow.
    STATE_CHAR_CLASS_DASH,            // Just read a dash in a character class.
    STATE_CHAR_CLASS_NO_CHAR,         // In a character class, expecting a char, no "-".
    STATE_QUANTIFIER_BEGIN,           // Just read "{" from a quantifier.
    STATE_QUANTIFIER_NUM1,            // Read at least one character from the first number.
    STATE_QUANTIFIER_COMMA,           // Read comma, expecting first character of second number in quantifier.
    STATE_QUANTIFIER_NUM2             // Read at least one character from second number.
};


// Check whether the pattern is valid.
//
// We use a finite state machine for the validation with the states of
// the enum FindWildShiftAndParserStates_.
inline bool _findWildShiftAndIsValid(CharString const & needle) {
    SEQAN_CHECKPOINT;

    if (length(needle) == 0)
        return false;

    typedef Iterator<CharString, Standard>::Type TIterator;

    FindWildShiftAndParserStates_ state = STATE_NO_CHAR;

    // TODO(holtgrew): Check i <= j in a{i, j}?

    for (TIterator it = begin(needle, Standard()); it != end(needle, Standard()); ++it) {
        switch (state) {
            case STATE_NO_CHAR:
                if (*it == '\\') {
                    state = STATE_ESCAPED;
                } else if (*it == '[') {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                } else if (*it != '+' && *it != '*' && *it != '?' && *it != '{' && *it != '}' && *it != '-') {
                    // No character, char class or "." here, must not be quantifier.
                    state = STATE_CHAR;
                } else {
                    return false;
                }
                break;
            case STATE_CHAR:
                if (*it == '\\') {
                    state = STATE_ESCAPED;
                } else if (*it == '{') {
                    state = STATE_QUANTIFIER_BEGIN;
                } else if (*it == '[') {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                } else if (*it == '+' || *it == '*' || *it == '?') {
                    state = STATE_NO_CHAR;
                } else if (*it == ']' || *it == '}' || *it == '-') {
                    // Other special character have to be escaped.
                    return false;
                }
                break;
            case STATE_ESCAPED:
                state = STATE_CHAR;
                break;
            case STATE_CHAR_CLASS:
                if (*it == ']') {
                    state = STATE_CHAR;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else if (*it == '-') {
                    state = STATE_CHAR_CLASS_DASH;
                } else if (*it == '+' || *it == '*' || *it == '?' || *it == '{' || *it == '}' || *it == '.' || *it == '[') {
                    // Other special character have to be escaped.
                    return false;
                }
                break;
            case STATE_CHAR_CLASS_ESCAPED:
                state = STATE_CHAR_CLASS;
                break;
            case STATE_CHAR_CLASS_DASH:
                if (*it == '+' || *it == '*' || *it == '?' || *it == '{' ||
                    *it == '}' || *it == '.' || *it == '[' || *it == ']' ||
                    *it == '-') {
                    return false;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                }
                break;
            case STATE_CHAR_CLASS_NO_CHAR:
                if (*it == ']') {
                    state = STATE_CHAR;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else if (*it == '+' || *it == '*' || *it == '?' || *it == '{' || *it == '}' || *it == '.' || *it == '-' || *it == '[') {
                    return false;
                } else {
                    state = STATE_CHAR_CLASS;
                }
                break;
            case STATE_QUANTIFIER_BEGIN:
                if (*it >= '0' && *it <= '9') {
                    state = STATE_QUANTIFIER_NUM1;
                } else {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_NUM1:
                if (*it == ',') {
                    state = STATE_QUANTIFIER_COMMA;
                } else if (*it < '0' || *it > '9') {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_COMMA:
                if (*it >= '0' && *it <= '9') {
                    state = STATE_QUANTIFIER_NUM2;
                } else {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_NUM2:
                if (*it == '}') {
                    state = STATE_NO_CHAR;
                } else if (*it < '0' || *it > '9') {
                    return false;
                }
                break;
        }
    }

    if (state == STATE_CHAR || state == STATE_NO_CHAR)
        return true;
    return false;
}


// Determine the pattern length without wildcard characters.  needle
// must be a valid pattern.
inline Position<CharString>::Type _findWildShiftAndLengthWithoutWildcards(CharString const & needle) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(needle));

    typedef Position<CharString>::Type TPosition;
    TPosition result = 0u;

    for (TPosition i = 0u; i < length(needle); ++i) {
        char c = needle[i];

        // Used in the switch statement below, must be defined before.
        CharString number;
        unsigned n, m;

        switch (c) {
            case '+':
            case '*':
            case '?':
                // Skip quantifier characters.
                continue;
            case '[':
                // Skip over character set, means one character.
                while (needle[i] != ']')
                    ++i;
                result += 1;
                break;
            case '{':
                // Skip over explicit quantifier, extract number of quantified
                // characters.
                ++i;
                // Read first number.
                while (needle[i] != ',') {
                    appendValue(number, needle[i]);
                    i += 1;
                }
                n = atoi(toCString(number));
                i += 1;  // Skip comma.
                // Read second number.
                clear(number);
                while (needle[i] != '}') {
                    appendValue(number, needle[i]);
                    i += 1;
                }
                m = atoi(toCString(number));

                SEQAN_ASSERT_GEQ(m, n);
                result += m - 1;
                break;
            case '\\':
                // Escape-sequences count as one character.
                result += 1;
                i += 1;
            default:
                result += 1;
                break;
        }
    }

    return result;
}


// Build a string with all characters in the range [begin, end) in host into
// result.
//
// TODO(holtgrew): We could simply use segments in the caller and a template argument for host instead of begin/end here.
inline void _findWildShiftAndGetCharacterClass(
        CharString & result, CharString const & host,
        Position<CharString>::Type begin, Position<CharString>::Type end) {
    SEQAN_CHECKPOINT;
    typedef Position<CharString>::Type TPosition;
    clear(result);

    for (TPosition pos = begin; pos < end; ++pos) {
        if (host[pos] == '\\') {
            SEQAN_ASSERT_LT(pos + 1, end);
            pos += 1;
            appendValue(result, host[pos]);
        } else if (host[pos] != '-') {
            appendValue(result, host[pos]);
        } else {
            // Read end of range character.
            SEQAN_ASSERT_LT(pos + 1, end);
            char last;
            if (host[pos + 1] == '\\') {
                SEQAN_ASSERT_LT(pos + 2, end);
                pos += 2;
                last = host[pos];
            } else {
                pos += 1;
                last = host[pos];
            }
            SEQAN_ASSERT_LEQ(back(result), last);
            for (char c = back(result) + 1; c <= last; ++c)
                appendValue(result, c);
        }
    }
}


// Called after the host has been set, does the precomputation.
template <typename TNeedle>
void _initializePattern(Pattern<TNeedle, WildShiftAnd> & me) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_INITIAL, me._state);
    
    TNeedle const & needle = value(me.data_host);

	SEQAN_ASSERT_TRUE(_findWildShiftAndIsValid(needle));

	typedef unsigned TWord;
    // TODO(holtgrew): TValue will always be char?!?!
	typedef typename Value<TNeedle>::Type TValue;
	
	me.needleLength = length(needle);
	me.character_count = _findWildShiftAndLengthWithoutWildcards(needle);

	if (me.character_count<1) me.blockCount=1;
	else me.blockCount=((me.character_count-1) / BitsPerValue<TWord>::VALUE)+1;
	
	clear(me.table);
	resize(me.table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

	clear(me.s_table);
	resize(me.s_table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

	clear(me.a_table);
	resize(me.a_table,me.blockCount,0,Exact());

	clear(me.prefSufMatch);
	resize(me.prefSufMatch, me.blockCount, 0, Exact());

	clear(me.df);
	resize(me.df, me.blockCount, 0, Exact());

	int i = -1;
	String <char> last_char; // stores the character (or characters) that were read in the last step
	TWord j=0;
	while(j < me.needleLength){
		if (convert<char>(getValue(needle,j)) == '+'){
			TWord len = length(last_char);
			for (unsigned int k = 0; k < len; ++k)
				me.s_table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
		} 
		else if (convert<char>(getValue(needle,j)) == '?'){
			me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
		}
		else if (convert<char>(getValue(needle,j)) == '*'){
			TWord len = length(last_char);
			for (unsigned int k = 0; k < len; ++k)
				me.s_table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
			me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
		}
		else if(convert<char>(getValue(needle,j)) == '['){
			/* find characters in class */
			TWord e = j;
			while(convert<char>(getValue(needle,e)) != ']') ++e;
			/* get character codes of class */
			_findWildShiftAndGetCharacterClass(last_char, needle, j+1, e);
			TWord len = length(last_char);			
			
			/* add class to the mask */
			++i;
			for (unsigned int k = 0; k < len; ++k){
				me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
			}
			j = e;
		}
		else if(convert<char>(getValue(needle,j)) == '.'){ // matches all characters in the current alphabet
			clear(last_char);
			++i;
			for(unsigned int l = 0;l < ValueSize<TValue>::VALUE;++l){
				append(last_char,l);
				me.table[me.blockCount*l + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
			}
			
		}
		else if(convert<char>(getValue(needle,j)) == '\\'){ // handle escape characters
			/* goto next character use this for the bit mask */
			++i;++j;
			clear(last_char);
			append(last_char, convert<TWord>(convert<TValue>(getValue(needle,j))));
			me.table[me.blockCount*last_char[0] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));			
		}
		else if(convert<char>(getValue(needle,j)) == '{'){ // handle bounded character repeats
			String <char> number;
			TWord n,m,r;
			TWord len = length(last_char);			
			n = m = 0;
			++j;
			while(convert<char>(getValue(needle,j)) != '}' && convert<char>(getValue(needle,j)) != ',') {			
				append(number,convert<char>(getValue(needle,j)));
				++j;
			}
			n = atoi(toCString(number));
			if (convert<char>(getValue(needle,j)) == ','){
				++j;
				clear(number);
				while(convert<char>(getValue(needle,j)) != '}') {			
					append(number,convert<char>(getValue(needle,j)));
					++j;
				}
				m = atoi(toCString(number));
			}
			// we already have seen one required occurence of the character (last_char)
			n -= 1;
			r = 0;
			while(r < n){ // add n normal characters
				++i;
				for (unsigned int k = 0; k < len; ++k){
					me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
				}
				++r;
			}
			++r; // correct the -1 of n to get in the correct relation to m
			while (r < m){ // if there was no m specified this won't be used
				// add m - n charaters and make them optional
				++i;
				for (unsigned int k = 0; k < len; ++k){
					me.table[me.blockCount*last_char[k] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
				}
				me.a_table[i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
				++r;
			}			
		}

		else // we have a character here
		{
			// determine character position in array table
			clear(last_char);
			append(last_char, convert<TWord>(convert<TValue>(getValue(needle,j))));
			++i;
			me.table[me.blockCount*last_char[0] + i / BitsPerValue<TWord>::VALUE] |= (1<<(i%BitsPerValue<TWord>::VALUE));
		}
		++j;
	}

	clear(me.i_table);
	resize(me.i_table,me.blockCount,0,Exact());

	clear(me.f_table);
	resize(me.f_table,me.blockCount,0,Exact());

	for (unsigned int i = 0; i < me.character_count; ++i){
		if ((me.a_table[i / BitsPerValue<TWord>::VALUE] & (1 << (i % BitsPerValue<TWord>::VALUE))) != 0){
			if ((me.f_table[i / BitsPerValue<TWord>::VALUE] & (1 << ((i-1) % BitsPerValue<TWord>::VALUE))) == 0){
				if(i > 0)
					me.i_table[(i-1) / BitsPerValue<TWord>::VALUE] |= 1 << ((i-1) % BitsPerValue<TWord>::VALUE);
				me.f_table[i / BitsPerValue<TWord>::VALUE] |= 1 << (i % BitsPerValue<TWord>::VALUE);
#if SEQAN_WILD_SHIFTAND_DEBUG
				std::cout << "Update F and I" << std::endl;
				_printMask(me.f_table,0,"F ");
				_printMask(me.i_table,0,"I ");
				_printMask(me.a_table,0,"A ");
				std::cout << std::endl;
#endif
			}
			else{
				TWord curBlock = i / BitsPerValue<TWord>::VALUE;
				for (unsigned int k = 0; k < me.blockCount; ++k){
					if(k != curBlock)
						me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~0;
					else
						me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~(1 << ((i-1) % BitsPerValue<TWord>::VALUE));

				}
				//me.f_table[i / BitsPerValue<TWord>::VALUE] &= ~(1 << ((i-1) % BitsPerValue<TWord>::VALUE));
				me.f_table[i / BitsPerValue<TWord>::VALUE] |= 1 << (i % BitsPerValue<TWord>::VALUE);
#if SEQAN_WILD_SHIFTAND_DEBUG
				std::cout << "Update F" << std::endl;
				_printMask(me.f_table,0,"F ");
				std::cout << std::endl;
#endif
			}
		}
	}

#if SEQAN_WILD_SHIFTAND_DEBUG	
	// Debug code
	std::cout << "Alphabet size: " << ValueSize<TValue>::VALUE << ::std::endl;
	std::cout << "Needle length (with wildcards): " << me.needleLength << ::std::endl;
	std::cout << "Needle length (wo wildcards): " << me.character_count << ::std::endl;
	std::cout << "Block count: " << me.blockCount << ::std::endl;

	std::cout << "Needle:" << needle << ::std::endl;

	_printMask(me.f_table,0,"F ");
	_printMask(me.i_table,0,"I ");
	_printMask(me.a_table,0,"A ");
	std::cout << std::endl << std::endl;

	for(unsigned i=0;i<ValueSize<TValue>::VALUE;++i) {
		if (((i<97) && (4 < i) ) || (i>122)) continue;
		std::cout << static_cast<TValue>(i) << ": ";
		for(unsigned int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned>::VALUE;++bit_pos) {
				std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
	}
#endif
}


template <typename THaystack, typename TNeedle>
inline bool _findWildShiftAndSmallNeedle(Finder<THaystack, Default> & finder,
                                           Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    THaystack const & hstck = haystack(finder);
	typedef unsigned TWord;
	TWord compare = (1 << (pattern.character_count - 1));
    for (; finder._endPosition < length(hstck); ++finder._endPosition) {
		TWord pos = convert<TWord>(hstck[finder._endPosition]);
		// added  | (pattern.prefSufMatch[0] & pattern.s_table[pattern.blockCount*pos]) at the end of the line
		pattern.prefSufMatch[0] = (((pattern.prefSufMatch[0] << 1) | 1) & pattern.table[pattern.blockCount*pos]) | (pattern.prefSufMatch[0] & pattern.s_table[pattern.blockCount*pos]);

		// additional bit operations
		pattern.df[0] = pattern.prefSufMatch[0] | pattern.f_table[0];
		pattern.prefSufMatch[0] |= ((pattern.a_table[0] & (~(pattern.df[0] - pattern.i_table[0]))) ^ pattern.df[0]);
		if ((pattern.prefSufMatch[0] & compare) != 0) {
            // Found a match: Update states, positions and report the match.
            finder._endPosition += 1;
            finder._state = TFinder::STATE_FOUND;
            pattern._state = TFinder::STATE_FOUND;
			return true; 
		}
	}
    finder._state = TFinder::STATE_NOTFOUND;
    pattern._state = TFinder::STATE_NOTFOUND;
	return false;
}


template <typename THaystack, typename TNeedle>
inline bool _findWildShiftAndLargeNeedle(Finder<THaystack, Default> & finder,
                                           Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    THaystack const & hstck = haystack(finder);
	typedef unsigned TWord;
	const TWord all1 = ~0;	
	TWord compare = (1 << ((pattern.character_count-1) % BitsPerValue<TWord>::VALUE));

    for (; finder._endPosition < length(hstck); ++finder._endPosition) {
		TWord pos = convert<TWord>(hstck[finder._endPosition]);
		TWord carry = 1;
		TWord wc_carry = 0;
		// shift of blocks with carry
		for(TWord block=0;block<pattern.blockCount;++block) {
			bool newCarry = ((pattern.prefSufMatch[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			pattern.prefSufMatch[block] = (((pattern.prefSufMatch[block] << 1) | carry) & pattern.table[pattern.blockCount*pos+block]) | (pattern.prefSufMatch[block] & pattern.s_table[pattern.blockCount*pos+block]) ;
			carry = newCarry;
			
			pattern.df[block] = pattern.prefSufMatch[block] | pattern.f_table[block];
			TWord Z = pattern.df[block] - pattern.i_table[block] - wc_carry;
			wc_carry = ((pattern.df[block] < Z) || (pattern.i_table[block]==all1 && wc_carry)) ? 1 : 0;
			pattern.prefSufMatch[block] |= (pattern.a_table[block] & (~Z ^ pattern.df[block]));
		}

#if SEQAN_WILD_SHIFTAND_DEBUG
		std::cout << "reading " << *finder << std::endl;
		_printMask(pattern.prefSufMatch,position(finder),"D ");
		_printMask(pattern.df,position(finder),"Df");
		std::cout << std::endl;
#endif
		
		// check for match
		if ((pattern.prefSufMatch[pattern.blockCount-1] & compare) != 0)
			return true; 
	}
	return false;
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        _initializePattern(pattern);
        finder._beginPosition = 0u;
        finder._endPosition = 0u;
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder)))
            return false;
        finder._beginPosition += 1;
    } else {
        finder._beginPosition += 1;
    }

	// Use fast algorithm for needles < machine word if possible.
	if (pattern.blockCount == 1)
		return _findWildShiftAndSmallNeedle(finder, pattern);
	else
		return _findWildShiftAndLargeNeedle(finder, pattern);
    return false;
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
