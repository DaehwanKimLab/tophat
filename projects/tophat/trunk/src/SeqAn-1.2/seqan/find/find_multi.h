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
  $Id: find_multi.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MULTI_H
#define SEQAN_HEADER_FIND_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct _MultipatternFinder;
typedef Tag<_MultipatternFinder> MultipatternFinder;
	
//____________________________________________________________________________

template <typename THaystack>
class Finder<THaystack, MultipatternFinder>
{
//____________________________________________________________________________
private:
	unsigned int data_pattern;

public:
	Finder():
		data_pattern(0)
	{
SEQAN_CHECKPOINT
	}

	Finder(Finder const & other_):
		data_pattern(other_.data_pattern)
	{
SEQAN_CHECKPOINT
	}

	~Finder()
	{
SEQAN_CHECKPOINT
	}
//____________________________________________________________________________

	Finder & 
	operator = (Finder const & other_)
	{
SEQAN_CHECKPOINT
		data_pattern = other_.data_pattern;
		return *this;
	}
//____________________________________________________________________________



//////////////////////////////////////////////////////////////////////////////
/*
template <typename THaystack, typename TNeedle>
bool
findNext(Finder & me,
		 THaystack & hstk,
		 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	++hstk;
	return find(me, hstk, ndl);
}*/

//////////////////////////////////////////////////////////////////////////////

};

template <typename THaystack>
inline unsigned int &
needle(Finder<THaystack, MultipatternFinder> & me)
{
SEQAN_CHECKPOINT
	return me.data_pattern;
}
template <typename THaystack>
inline unsigned int const &
needle(Finder<THaystack, MultipatternFinder> const & me)
{
SEQAN_CHECKPOINT
	return me.data_pattern;
}
//____________________________________________________________________________

template <typename THaystack>
inline void
setNeedle(Finder<THaystack, MultipatternFinder> & me, unsigned int const needleIndex_)
{
SEQAN_CHECKPOINT
	me.data_pattern = needleIndex_;
}

//____________________________________________________________________________

template <typename THaystack>
inline void
init(Finder<THaystack, MultipatternFinder> & me)
{
SEQAN_CHECKPOINT
	me.data_pattern = 0;
}
//____________________________________________________________________________


template <typename THaystack, typename THaystack2, typename TNeedle>
inline bool
find(Finder<THaystack, MultipatternFinder> & me,
	 THaystack2 & hstk,
	 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	while ( needle(me) < length(ndl) )
	{
		Finder<THaystack2, Horspool> horspool(ndl[needle(me)]);
		bool found = find(horspool, hstk, ndl[needle(me)]);
		if (found)
		{
			return true;
		}
		setPosition(hstk, 0);
		++needle(me);
	}
	return false;
}
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
