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
// Journaled String implementation.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes
// ============================================================================

/**
.Spec.Journaled String
..summary:Journaled versions of arbitrary underlying string.
..signature:String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >
..include:seqan/sequence_journaled.h
..cat:Sequences
 */

template <typename THostSpec, typename TJournalSpec = SortedArray, typename TBufferSpec = Alloc<void> >
struct Journaled {};

template <typename TValue_, typename THostSpec_, typename TJournalSpec_, typename TBufferSpec_>
class String<TValue_, Journaled<THostSpec_, TJournalSpec_, TBufferSpec_> >
{
public:
    typedef TValue_ TValue;
    typedef THostSpec_ THostSpec;
    typedef TJournalSpec_ TJournalSpec;
    typedef TBufferSpec_ TBufferSpec;

    typedef String<TValue, THostSpec> THost;
    typedef typename Size<THost>::Type TSize;
    typedef typename Position<THost>::Type TPosition;
    typedef String<TValue, TBufferSpec> TInsertionBuffer;
    typedef JournalEntry<TSize, TPosition> TJournalEntry;
    typedef JournalEntries<TJournalEntry, TJournalSpec> TJournalEntries;

    // The underlying host string.
    Holder<THost> _holder;
    // A buffer for inserted strings.
    TInsertionBuffer _insertionBuffer;
    // The journal is a sorted set of TJournalEntry objects, the exact types
    // depends on TJournalSpec.  Note that the entries resemble a partial
    // sum datastructure.
    TJournalEntries _journalEntries;
    // The journaled string's size.
    TSize _length;

    String() {}

    String(THost & host)
    {
        SEQAN_CHECKPOINT;
        setHost(*this, host);
    }

    // TODO(holtgrew): Actually, we want to have a proxy for non-const.
    TValue operator[](TPosition const & pos) const
    {
        SEQAN_CHECKPOINT;
        return getValue(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

///.Metafunction.Host.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef String<TValue, THostSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef String<TValue, THostSpec> const Type;
};

/**
.Metafunction.InsertionBuffer
..summary:Return type of insertion buffer string for a journaled string.
..param.T:Spec.Journaled String
..include:sequan/sequence_journal.h
*/
template <typename T>
struct InsertionBuffer;

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct InsertionBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef String<TValue, TBufferSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct InsertionBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef String<TValue, TBufferSpec> const Type;
};

///.Metafunction.Size.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef typename String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >::TSize Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
    : Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > > {};

///.Metafunction.Position.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef typename String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >::TPosition Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
    : Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > > {};

///.Metafunction.Value.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Value<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef TValue Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Value<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
  typedef TValue Type;
};

///.Metafunction.GetValue.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef TValue Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
  typedef TValue Type;
};

///.Metafunction.Reference.param.T:Spec.Journaled String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
  typedef TValue & Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Reference<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
  typedef TValue const & Type;
};

/**
.Metafunction.JournalType
..signature:JournalType<T>::Type
..summary:Metafunction for returning the type of the journal of a Journaled String.
..param.T:Spec.Journaled String
..include:seqan/string_journaled.h
 */
template <typename T>
struct JournalType;

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >
{
    typedef typename Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type TSize_;
    typedef typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type TPosition_;
    typedef JournalEntry<TSize_, TPosition_> TJournalEntry_;

    typedef JournalEntries<TJournalEntry_, TJournalSpec> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>
{
    typedef typename JournalType<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
TStream &
operator<<(TStream & stream, String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & s)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > TString;
    typedef typename TString::TJournalEntries TJournalEntries;
    typedef typename Iterator<TJournalEntries const, Standard>::Type TIterator;

    for (TIterator it = begin(s._journalEntries), itend = end(s._journalEntries); it != itend; ++it) {
        if (value(it).segmentSource == SOURCE_ORIGINAL) {
            stream << infix(value(s._holder), value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        } else {
            SEQAN_ASSERT_EQ(value(it).segmentSource, SOURCE_PATCH);
            stream << infix(s._insertionBuffer, value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        }
    }
    return stream;
}

/**
.Function.setHost:
..param.object.type:Spec.Journaled String
..param.host.type:Class.String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TSequence2>
inline
void
setHost(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString, TSequence2 & str)
{
    SEQAN_CHECKPOINT;
    setValue(journaledString._holder, str);
    journaledString._length = length(str);
    reinit(journaledString._journalEntries, length(str));
}

/**
.Function.host:
..param.object.type:Spec.Journaled String
..include:seqan/sequence_journaled.h
*/
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type &
host(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    return value(journaledString._holder);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Host<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type const &
host(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return value(journaledString._holder);
}

/**
.Function.clear:
..param.object.type:Spec.Journaled String
..include:seqan/sequence_journaled.h
 */
// TODO(holtgrew): Behaviour is to clear the journal, not the string!
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline void
clear(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    reinit(journaledString._journalEntries, length(host(journaledString)));
}

/**
.Function.flatten:
..summary:Apply the journal to the underlying string, destructively on the underlying string.
..signature:flatten(journaledString)
..param.journaledString:The journaled string to flatten.
...type:Spec.Journaled String
..include:seqan/sequence_journaled.h
 */
// TODO(holtgrew): Write me! What about non-destructive version that creates a new copy and sets holder to it?

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline void
erase(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
      TPos const & pos,
      TPos const & posEnd)
{
    SEQAN_CHECKPOINT;
	SEQAN_ASSERT_GEQ(static_cast<TPos>(journaledString._length), pos);
	SEQAN_ASSERT_GEQ(static_cast<TPos>(journaledString._length), posEnd);
    SEQAN_ASSERT_GEQ(static_cast<TPos>(journaledString._length), posEnd - pos);
    journaledString._length -= posEnd - pos;
    recordErase(journaledString._journalEntries, pos, posEnd);
    if (length(journaledString._journalEntries) == 0)
        clear(journaledString._insertionBuffer);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline void
erase(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
      TPos const & pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(journaledString._length, 1u);
    erase(journaledString, pos, pos + 1);
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TString, typename TPos>
inline void
insert(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
       TPos const & pos,
       TString const & seq)
{
    SEQAN_CHECKPOINT;
    journaledString._length += length(seq);
    TPos beginPos = length(journaledString._insertionBuffer);
    append(journaledString._insertionBuffer, seq);
    recordInsertion(journaledString._journalEntries, pos, beginPos, length(seq));
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos, typename TValue2>
inline void
insertValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TPos const & pos,
            TValue2 const & value)
{
    SEQAN_CHECKPOINT;
    journaledString._length += 1;
    TPos beginPos = length(journaledString._insertionBuffer);
    appendValue(journaledString._insertionBuffer, value);
    recordInsertion(journaledString._journalEntries, pos, beginPos, 1u);
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos, typename TSequence2>
inline void
assignInfix(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TPos const & beginPos,
            TPos const & endPos,
            TSequence2 const & valueString)
{
    SEQAN_CHECKPOINT;
    erase(journaledString, beginPos, endPos);
    insert(journaledString, beginPos, valueString);
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos, typename TValue2>
inline void
assignValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString,
            TPos const & pos,
            TValue2 const & value)
{
    SEQAN_CHECKPOINT;
    erase(journaledString, pos);
    insertValue(journaledString, pos, value);
}


// TODO(holtgrew): Batch-Assignment of values through segments?

// TODO(holtgrew): begin
// TODO(holtgrew): empty
// TODO(holtgrew): end
// TODO(holtgrew): flatten
// TODO(holtgrew): fill
// TODO(holtgrew): getValue
// TODO(holtgrew): infix
// TODO(holtgrew): infixWithLength
// TODO(holtgrew): iter

// TODO(holtgrew): Unused, remove?
/*
template <typename TSequence, typename TJournalSpec>
inline
typename Value<TSequence>::Type const &
front(String...<TSequence, TJournalSpec> const & journaledString)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode frontNode = front(journaledString._journalEntries);
    if (frontNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), frontNode->virtualPosition + frontNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(frontNode->segmentSource, SOURCE_PATCH);
        return getValue(journaledString._insertionBuffer, frontNode->virtualPosition + frontNode->length - 1);
    }
}

// front/back clash with general sequence definitions.
template <typename TSequence, typename TJournalSpec>
inline
TValue const &
back(SequenceJournal<TSequence, TJournalSpec> const & journaledString)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode backNode = back(journaledString._journalEntries);
    if (backNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), backNode->virtualPosition + backNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(backNode->segmentSource, SOURCE_PATCH);
        return getValue(journaledString._insertionBuffer, backNode->virtualPosition + backNode->length - 1);
    }
}
*/

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Size<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
length(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return journaledString._length;
}


// TODO(holtgrew): toCString
// TODO(holtgrew): value

// TOOD(holtgrew): operator<
// TOOD(holtgrew): operator>
// TOOD(holtgrew): operator<=
// TOOD(holtgrew): operator>=
// TOOD(holtgrew): operator==
// TOOD(holtgrew): operator!=

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename GetValue<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const>::Type
getValue(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
         typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type const & pos)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const TJournaledString;
    typedef typename TJournaledString::TJournalEntry TJournalEntry;
    typedef typename Position<TJournaledString>::Type TPos;

    TJournalEntry entry = findJournalEntry(journaledString._journalEntries, pos);
    TPos relativePos = pos - entry.virtualPosition;

    if (entry.segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journaledString._holder), entry.physicalPosition + relativePos);
    } else {
        return getValue(journaledString._insertionBuffer, entry.physicalPosition + relativePos);
    }
}


// Note that if pos is in a gap, we return the position of the entry
// after the gap in the host.
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
typename Position<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > >::Type
virtualToHostPosition(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
                      TPos const & pos)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): With a better journal entries datastructure, we could solve the main problem here. At the moment, we delegate completely.
    return virtualToHostPosition(journaledString._journalEntries, pos);
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec, typename TPos>
inline
bool
isGapInHost(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString,
            TPos const & pos)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): With a better journal entries datastructure, we could solve the main problem here. At the moment, we delegate completely.
    return isGapInHost(journaledString._journalEntries, pos);
}


template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
const void *
id(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journaledString)
{
    SEQAN_CHECKPOINT;
    return id(value(journaledString._holder));
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
const void *
id(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journaledString)
{
    SEQAN_CHECKPOINT;
    return id(value(journaledString._holder));
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
