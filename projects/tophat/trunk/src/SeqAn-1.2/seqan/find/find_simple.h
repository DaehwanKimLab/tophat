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
  $Id: find_bndm.h 2416 2008-06-17 15:30:38Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SIMPLE_H
#define SEQAN_HEADER_FIND_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Simple Finder
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Simple Finder:
..summary:A brute force online searching algorithm.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, Simple>
..param.TNeedle:The needle type.
...type:Class.String
..remarks:This specialization should only be used if no other is applicable.
*/

///.Class.Pattern.param.TSpec.type:Spec.Simple Finder

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Simple> {
//____________________________________________________________________________
public:
	Holder<TNeedle> data_host;

//____________________________________________________________________________

	Pattern() {}

	Pattern(Pattern const & other):
		data_host(other.data_host)
	{
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

	~Pattern(){}

	Pattern const & 
	operator = (Pattern const & other)
	{
		data_host = other.data_host;
		return *this;
	}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, Simple> & me, 
			  TNeedle2 & needle) 
{
	setValue(me.data_host, needle);
}
template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, Simple> & me, 
			  TNeedle2 const & needle) 
{
	setValue(me.data_host, needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, 
				 Pattern<TNeedle, Simple> & me) 
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;
	typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;

	if (empty(finder))
	{
		_setFinderLength(finder, length(needle(me)));
		_finderSetNonEmpty(finder);
	}
	else ++finder;

	THaystack const & hstk = haystack(finder);
	TNeedle const & ndl = needle(me);

	THaystackIterator res = ::std::search(begin(hstk, Standard())+position(finder), end(hstk, Standard()), begin(ndl, Standard()), end(ndl, Standard()));

	if (res == end(hstk)) return false;

	_setFinderEnd(finder, (res - begin(hstk, Standard())) + length(ndl));
	setPosition(finder, beginPosition(finder));
	return true; 

/*
	TSize n = length(hstk);
	TSize m = length(ndl);
	while (position(finder)+m <= n)
	{
		if (ndl == infix(hstk, position(finder), position(finder)+m))
		{
			_setFinderEnd(finder, position(finder)+m);
			return true; 
		}
		++finder;
	}
	return false;
*/
}

//____________________________________________________________________________

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
