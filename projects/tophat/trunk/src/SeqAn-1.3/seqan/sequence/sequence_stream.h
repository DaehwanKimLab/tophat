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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Code that allows strings to be treated like streams.
//
// TODO(holtgrew): This could be called ${module}_adapt_sequence.h, where $module is the name of a future streaming/I/O module.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_STREAM_H
#define SEQAN_HEADER_SEQUENCE_STREAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline bool
_streamEOF(Iter<TContainer, TSpec> const & iter)
{
SEQAN_CHECKPOINT
    return atEnd(iter);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TContainer, typename TSpec>
inline ::std::streamsize
_streamRead(TValue * target,
            Iter<TContainer, TSpec> & source,
            ::std::streamsize limit)
{
SEQAN_CHECKPOINT
    if (position(target) + limit > length(container(target)))
        limit = length(container(target)) - position(target);
    Iter<TContainer, TSpec> sourceEnd = source + limit;
    for (; source != sourceEnd; ++source, ++target)
        *target = *source;
    return limit;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Value<Iter<TContainer, TSpec> >::Type
_streamGet(Iter<TContainer, TSpec> & source)
{
SEQAN_CHECKPOINT
    typename Value<Iter<TContainer, TSpec> >::Type _val = getValue(source);
    goNext(source);
    return _val;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Value<Iter<TContainer, TSpec> >::Type
_streamPeek(Iter<TContainer, TSpec> & source)
{
SEQAN_CHECKPOINT
    return getValue(source);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec, typename TChar>
inline void
_streamPut(Iter<TContainer, TSpec> & target,
           TChar character)
{
SEQAN_CHECKPOINT
    if (atEnd(target))
    {
        typename Container<Iter<TContainer, TSpec> >::Type & container_ = container(target);
        appendValue(container_, character);
        target = begin(container_) + (length(container_) - 1);
    } else
        *target = character;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Position<Iter<TContainer, TSpec> >::Type
_streamTellG(Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
    return position(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline typename Position<Iter<TContainer, TSpec> >::Type
_streamTellP(Iter<TContainer, TSpec> & me)
{
SEQAN_CHECKPOINT
    return position(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeekG(Iter<TContainer, TSpec> & me,
     typename Position<Iter<TContainer, TSpec> >::Type pos)
{
SEQAN_CHECKPOINT
    me = begin(container(me)) + pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeekP(Iter<TContainer, TSpec> & me,
     typename Position<Iter<TContainer, TSpec> >::Type pos)
{
SEQAN_CHECKPOINT
    me = begin(container(me)) + pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TSpec>
inline void
_streamSeek2G(Iter<TContainer, TSpec> & me,
     int off)
{
SEQAN_CHECKPOINT
    me = begin(container(me)) + (position(me) + off);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
