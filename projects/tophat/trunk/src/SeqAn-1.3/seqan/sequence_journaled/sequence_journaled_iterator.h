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
// Code for the Journaled string iterator.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

template <typename TJournaledStringSpec>
struct JournaledStringIterSpec;

template <typename TJournaledString, typename TJournalSpec>
class Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
{
public:
    typedef Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > TIterator;
    typedef typename TJournaledString::TValue TValue;
    typedef typename JournalType<TJournaledString>::Type TJournalEntries;
    // We need a rooted iterator for iterating the journal tree since we need atEnd().
    typedef typename Iterator<TJournalEntries, Rooted>::Type TJournalEntriesIterator;
    typedef typename Host<TJournaledString>::Type THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;
    typedef typename InsertionBuffer<TJournaledString>::Type TInsertionBuffer;
    typedef typename Iterator<TInsertionBuffer, Standard>::Type TInsertionBufferIterator;

    // The journal string we iterate over.
    TJournaledString * _journalStringPtr;
    // Iterator over the segments in the journal tree.
    TJournalEntriesIterator _journalEntriesIterator;
    // Begin and end iterator in the host string of the journal string.
    THostIterator _hostSegmentBegin;
    THostIterator _hostSegmentEnd;
    // Current iterator in the host segment.
    THostIterator _currentHostIt;
    // Begin and end iterator in the insertion buffer of the journal string.
    TInsertionBufferIterator _insertionBufferSegmentBegin;
    TInsertionBufferIterator _insertionBufferSegmentEnd;
    // Current iterator in the insertion buffer.
    TInsertionBufferIterator _currentInsertionBufferIt;

    Iter() : _journalStringPtr(0) { SEQAN_CHECKPOINT; }

    Iter(TIterator const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename IterComplementConst<TIterator>::Type const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    explicit
    Iter(TJournaledString & journalString)
    {
        SEQAN_CHECKPOINT;
        _initJournaledStringIterator(*this, journalString);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// For String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >

///.Metafunction.Iterator.param.T:Spec.Journal String
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>
{
    typedef Iter<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, JournaledStringIterSpec<TJournalSpec> > Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
struct Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>
{
    typedef Iter<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, JournaledStringIterSpec<TJournalSpec> > Type;
};

// For Iter<TJournaledString, TJournaledStringIterSpec>
template <typename TJournaledString, typename TJournaledStringIterSpec>
struct GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > >
{
    typedef typename GetValue<TJournaledString>::Type Type;
};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > const>
        : GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > > {};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct Value<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > >
{
    typedef typename Value<TJournaledString>::Type Type;
};

template <typename TJournaledString, typename TJournaledStringIterSpec>
struct Value<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > const>
        : Value<Iter<TJournaledString, JournaledStringIterSpec<TJournaledStringIterSpec> > > {};

// ============================================================================
// Functions
// ============================================================================

// For String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type
begin(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type TResult;
    return TResult(journalString);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type
begin(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type TResult;
    return TResult(journalString);
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type
end(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const & journalString, Standard)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > const, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIteratorEnd(result, journalString);
    return result;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBufferSpec>
inline
typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type
end(String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> > & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >, Standard>::Type TResult;
    TResult result;
    _initJournaledStringIteratorEnd(result, journalString);
    return result;
}

// For Iter<TJournaledString, JournaledStringIterSpec>

template <typename TJournaledString, typename TJournalSpec>
inline
void
_initJournaledStringIterator(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
                           TJournaledString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = begin(journalString._journalEntries);
    // Update iterators on the segment.
    _updateSegmentIterators(iterator);
}

template <typename TJournaledString, typename TJournalSpec>
inline
void
_initJournaledStringIteratorEnd(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
                              TJournaledString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = end(journalString._journalEntries);
}

template <typename TJournaledString, typename TJournalSpec>
inline
void
_updateSegmentIterators(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (atEnd(iterator._journalEntriesIterator))
        return;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
//             static_cast<int>(begin(host(*iterator._journalStringPtr), Standard()));
//             static_cast<int>(host(*iterator._journalStringPtr));
//             static_cast<int>(*iterator._journalStringPtr);
//             static_cast<int>(iterator._journalStringPtr);
//             static_cast<int>(iterator._hostSegmentBegin);
            iterator._hostSegmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._hostSegmentEnd = iterator._hostSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentHostIt = iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            iterator._insertionBufferSegmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._insertionBufferSegmentEnd = iterator._insertionBufferSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentInsertionBufferIt = iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Value<TJournaledString>::Type
value(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Value<TJournaledString>::Type
value(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return getValue(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return getValue(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return getValue(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return getValue(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator++(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            ++iterator._currentHostIt;
            if (iterator._currentHostIt == iterator._hostSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        case SOURCE_PATCH:
            ++iterator._currentInsertionBufferIt;
            if (iterator._currentInsertionBufferIt == iterator._insertionBufferSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
    return iterator;
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator++(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    ++iterator;
    return temp;    
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Value<TJournaledString>::Type
operator*(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournaledString, typename TJournalSpec>
inline
typename Value<TJournaledString>::Type
operator*(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > &
operator+=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > & iterator,
          typename Size<TJournaledString>::Type len)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<TJournaledString>::Type TSize;
    while (len > 0) {
        TSize remaining;
        switch (value(iterator._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                remaining = iterator._hostSegmentEnd - iterator._currentHostIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentHostIt += len;
                    len = 0;
                }
                break;
            case SOURCE_PATCH:
                remaining = iterator._insertionBufferSegmentEnd - iterator._currentInsertionBufferIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentInsertionBufferIt += len;
                    len = 0;
                }
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid segment source!");
        }
    }
    return iterator;
}

template <typename TJournaledString, typename TJournalSpec>
inline
Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> >
operator+(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & iterator,
          typename Size<TJournaledString>::Type const & len)
{
    SEQAN_CHECKPOINT;
    Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > temp(iterator);
    temp += len;
    return temp;
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    if (atEnd(a._journalEntriesIterator) && atEnd(b._journalEntriesIterator))
        return true;
    if (a._journalEntriesIterator != b._journalEntriesIterator)
        return false;
    if (value(a._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        if (a._currentHostIt != b._currentHostIt)
            return false;
    } else {
        SEQAN_ASSERT_EQ(value(a._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        if (a._currentInsertionBufferIt != b._currentInsertionBufferIt)
            return false;
    }
    return true;
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    typedef typename IterMakeConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

template <typename TJournaledString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournaledString, JournaledStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
