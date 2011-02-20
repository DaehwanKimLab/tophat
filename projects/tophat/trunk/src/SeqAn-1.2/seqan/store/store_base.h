 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_STORE_BASE_H
#define SEQAN_HEADER_STORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Base structs
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

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
		return seqPos <= other.seqPos || gapPos <= other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator >= (TOther const &other) const
	{
		return seqPos >= other.seqPos || gapPos >= other.gapPos;
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


struct _SortSeqPos;
typedef Tag<_SortSeqPos> const SortSeqPos;

struct _SortGapPos;
typedef Tag<_SortGapPos> const SortGapPos;

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
