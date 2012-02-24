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

#ifndef SEQAN_HEADER_STORE_BASE_H
#define SEQAN_HEADER_STORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Base structs
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Class.GapAnchor
..summary:Stores the position of an alignment character in sequence-space and in gap-space.
..cat:Alignments
..signature:GapAnchor<TPos>
..param.TPos:Type to store gapped/ungapped positions.
..remarks:Value types of the $gaps$ strings in @Class.ReadStoreElement@ and @Class.ContigStoreElement@.

.Memfunc.GapAnchor#GapAnchor
..summary:Constructor
..signature:GapAnchor<TPos> ()
..signature:GapAnchor<TPos> (TPos seqPos, TPos gapPos)
..param.seqPos:Sequence character position in the ungapped sequence.
..param.gapPos:Sequence character position in the gapped sequence.
..remarks:Default constructor sets both positions to $0$.
..class:Class.GapAnchor
.Memvar.GapAnchor#seqPos
..summary:Sequence character position in the ungapped sequence.
..class:Class.GapAnchor
.Memvar.GapAnchor#gapPos
..summary:Sequence character position in the gapped sequence.
..class:Class.GapAnchor
..include:seqan/store.h
*/

// We store gap anchors only for the first text character behind a gap or a clipped sequence character

template <typename TPos>
struct GapAnchor {
	TPos	seqPos;			// sequence character position in the ungapped sequence
	TPos	gapPos;			// sequence character position in the gapped sequence

	GapAnchor() : seqPos(0), gapPos(0) {}
	GapAnchor(TPos sP, TPos gP) : seqPos(sP), gapPos(gP) {}
	
	template <typename TPos_>
	GapAnchor(GapAnchor<TPos_> const &other)
	{
		seqPos = other.seqPos;
		gapPos = other.gapPos;
	} 

	template <typename TPos_>
	inline GapAnchor const &
	operator = (GapAnchor<TPos_> const &other)
	{
		seqPos = other.seqPos;
		gapPos = other.gapPos;
		return *this;
	} 

	template <typename TOther>
	inline bool
	operator == (TOther const &other) const
	{
		return seqPos == other.seqPos && gapPos == other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator != (TOther const &other) const
	{
		return !(*this == other);
	} 

	template <typename TOther>
	inline bool
	operator < (TOther const &other) const
	{
		return seqPos < other.seqPos || gapPos < other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator > (TOther const &other) const
	{
		return seqPos > other.seqPos || gapPos > other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator <= (TOther const &other) const
	{
		return seqPos < other.seqPos || gapPos <= other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator >= (TOther const &other) const
	{
		return seqPos > other.seqPos || gapPos >= other.gapPos;
	}
};

template <typename TPos>
struct Size< GapAnchor<TPos> >
{
	typedef TPos Type;
};

template <typename TPos>
struct Size< GapAnchor<TPos> const>:
	public Size< GapAnchor<TPos> > {};

template <typename TPos>
struct Position< GapAnchor<TPos> >
{
	typedef TPos Type;
};

template <typename TPos>
struct Position< GapAnchor<TPos> const>:
	public Position< GapAnchor<TPos> > {};


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Sorting tags (just for lower_bound and upper_bound, positions are always sorted)
//////////////////////////////////////////////////////////////////////////////


struct SortSeqPos_;
typedef Tag<SortSeqPos_> const SortSeqPos;

struct SortGapPos_;
typedef Tag<SortGapPos_> const SortGapPos;

//////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor, typename TTag>
struct _LessGapAnchor;

template <typename TGapAnchor>
struct _LessGapAnchor<TGapAnchor, SortSeqPos> :
	public ::std::binary_function<TGapAnchor, TGapAnchor, bool>
{
	inline bool 
	operator() (TGapAnchor const& a1, TGapAnchor const& a2) const {
		return (a1.seqPos) < (a2.seqPos);
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor>
struct _LessGapAnchor<TGapAnchor, SortGapPos> :
	public ::std::binary_function<TGapAnchor, TGapAnchor, bool>
{
	inline bool 
	operator() (TGapAnchor const& a1, TGapAnchor const& a2) const {
		return (a1.gapPos) < (a2.gapPos);
	}
};

//////////////////////////////////////////////////////////////////////////////
// Lower and upper bound search
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
lowerBoundGapAnchor(TGapAnchor const& gaps, 
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortSeqPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
upperBoundGapAnchor(TGapAnchor const& gaps, 
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortSeqPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
lowerBoundGapAnchor(TGapAnchor const& gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortGapPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
upperBoundGapAnchor(TGapAnchor const& gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortGapPos const>() );
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
