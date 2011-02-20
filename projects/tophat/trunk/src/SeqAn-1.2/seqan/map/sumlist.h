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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SUMLIST_H
#define SEQAN_HEADER_SUMLIST_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Default Specialization Tag

template <typename TSpec = Default>
struct SkipSumList;

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec = SkipSumList<> >
class SumList;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct DIMENSION
{
	enum { VALUE = 1}; //dummy implementation to make VC++ happy
};

template <unsigned int _DIM, typename TValue, typename TSpec>
struct DIMENSION< SumList<_DIM, TValue, TSpec> >
{
	enum { VALUE = _DIM};
};
template <unsigned int _DIM, typename TValue, typename TSpec>
struct DIMENSION< SumList<_DIM, TValue, TSpec> const>
{
	enum { VALUE = _DIM};
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
struct Value< SumList<DIM, TValue, TSpec> >
{
	typedef TValue Type;
};
template <unsigned int DIM, typename TValue, typename TSpec>
struct Value< SumList<DIM, TValue, TSpec> const >
{
	typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////
// SumListValues: 
// Stores one DIM tupel of Values
//////////////////////////////////////////////////////////////////////////////


template <unsigned int DIM, typename TValue>
struct SumListValues
{
	TValue values[DIM];

	SumListValues(MinimalCtor)
	{}
	SumListValues()
	{
		clear(*this);
	}
	SumListValues(SumListValues const & other)
	{
		*this = other;
	}
	SumListValues(TValue const * arr)
	{
		*this = arr;
	}
	~SumListValues()
	{
	}
	SumListValues const & operator = (SumListValues const & other)
	{
		arrayCopy(other.values, other.values + DIM, values);
		return *this;
	}
	SumListValues const & operator = (TValue const * arr)
	{
		arrayCopy(arr, arr + DIM, values);
		return *this;
	}

	TValue & operator [] (unsigned int pos)
	{
		return values[pos];
	}
	TValue const & operator [] (unsigned int pos) const
	{
		return values[pos];
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Values;

template <unsigned int DIM, typename TValue, typename TSpec>
struct Values< SumList<DIM, TValue, TSpec> >
{
	typedef SumListValues<DIM, TValue> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline SumListValues<DIM, TValue> const &
operator += (SumListValues<DIM, TValue> & left,
			 SumListValues<DIM, TValue> const & right)
{
	for (unsigned int i = 0; i < DIM; ++i)
	{
		left.values[i] += right.values[i];
	}
	return left;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline SumListValues<DIM, TValue> const &
operator -= (SumListValues<DIM, TValue> & left,
			 SumListValues<DIM, TValue> const & right)
{
	for (unsigned int i = 0; i < DIM; ++i)
	{
		left.values[i] -= right.values[i];
	}
	return left;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline bool
operator == (SumListValues<DIM, TValue> const & left,
			 SumListValues<DIM, TValue> const & right)
{
	for (unsigned int i = 0; i < DIM; ++i)
	{
		if (left.values[i] != right.values[i]) return false;
	}
	return true;
}
template <unsigned int DIM, typename TValue>
inline bool
operator != (SumListValues<DIM, TValue> const & left,
			 SumListValues<DIM, TValue> const & right)
{
	for (unsigned int i = 0; i < DIM; ++i)
	{
		if (left.values[i] != right.values[i]) return true;
	}
	return false;
}
//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline void
clear(SumListValues<DIM, TValue> & me)
{
	arrayFill(me.values, me.values + DIM, 0);
}


//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
