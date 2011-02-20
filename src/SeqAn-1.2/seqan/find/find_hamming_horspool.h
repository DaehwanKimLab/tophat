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
  $Id:$
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_HAMMING_HORSPOOL_H
#define SEQAN_HEADER_FIND_HAMMING_HORSPOOL_H



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// approximate Boyer-Moore-Horspool
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.HammingHorspool:
..summary:Hamming distance string matching using approximate Boyer-Moore-Horspool algorithm
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, HammingHorspool>
..param.TNeedle:The needle type.
...type:Class.String
*/

///.Class.Pattern.param.TSpec.type:Spec.HammingHorspool

struct _HammingHorspool;
typedef Tag<_HammingHorspool> HammingHorspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, HammingHorspool>
{
//____________________________________________________________________________
public:
	typedef typename Size<TNeedle>::Type TSize;

	Holder<TNeedle>	data_host; 
	String<TSize>	shift_table; 	// c x (k+1) shift table (c = alphabet size)   
	TSize			shift;			// current shift
	unsigned int	k;				// maximal number of allowed mismatches
	
//____________________________________________________________________________

public:
	Pattern() {}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int const k)
	{
		setHost(*this, ndl, k);
	}

//____________________________________________________________________________
};

template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingHorspool> & me, 
		TNeedle2 const & ndl, 
		int const k)
{
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TValue;
	typename typename Iterator<TNeedle2 const, Standard>::Type TNeedle2Iterator;

	SEQAN_ASSERT(!empty(ndl));

	me.k = k;
	me.shift = 0;
	TSize m = length(ndl);
	TSize s = m - me.k;

	String<TSize> ready; // keep last updates

	TSize alphabet_size = ValueSize<TValue>::VALUE;
	TSize shift_size = alphabet_size * (me.k + 1);

	// get Space for ready-table and shift table
	resize(ready, alphabet_size);
	resize(me.shift_table, shift_size); // c x (k+1)

	// initialize shift table and ready table
	arrayFill(begin(ready), begin(ready) + alphabet_size, length(ndl)+1);
    arrayFill(begin(me.shift_table), begin(me.shift_table) + shift_size, length(ndl));

	// fill shift table by scanning the needle
	TNeedle2Iterator it = begin(ndl, Standard())+ length(ndl)-1;
	
	for (TValue i = m; i >= 1; --i)
	{
		--it;  
		TValue pos = ordValue(*it); //conversion value type to unsigned int
		TValue v = (i > s) ? i : s;

		for (TValue j = ready[pos]-1; j >= v ; --j)
		{
			me.shift_table[pos*(me.k+1)+j-s]= j-i+1;
		}
		ready[pos]= v;

	}

	me.data_host = ndl;
}



template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingHorspool> & horsp, TNeedle2 & ndl, unsigned int k)
{
	setHost(horsp, reinterpret_cast<TNeedle2 const &>(ndl),k);
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _finderInit (Pattern<TNeedle, HammingHorspool> & me) 
{
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, HammingHorspool>const>::Type & 
host(Pattern<TNeedle, HammingHorspool> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, HammingHorspool>const>::Type & 
host(Pattern<TNeedle, HammingHorspool> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

//____________________________________________________________________________

template <typename TFinder, typename TNeedle>
bool
find(TFinder & finder, Pattern<TNeedle, HammingHorspool> & me)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & haystack = container(finder);
	typename Position<TFinder>::Type pos = position(finder);
	typename Size<THaystack>::Type hstk_size = length(haystack);
	typename Size<TNeedle>::Type ndl_size = length(host(me));

	if (ndl_size > hstk_size)
	{//needle larger than haystack: nothing to find
		return false;
	}

	typename Position<THaystack>::Type pos_max = hstk_size - ndl_size; 
	SEQAN_ASSERT2(pos <= pos_max, "invalid search position");

	//helper variables
	typename Position<TFinder>::Type new_pos;
	TFinder h = finder + ndl_size -1;

	typename Iterator<TNeedle const, Standard>::Type it;
	it = begin(host(me), Standard())+ ndl_size -1;
	unsigned int s = ndl_size - me.k;

	if (empty(finder))
	{
		_finderInit(me);
		_setFinderLength(finder, length(needle(me)));
		_finderSetNonEmpty(finder);
	}
	
	while(true){
	
		//move to next position: shift to the right in text
		finder += me.shift;
		pos +=  me.shift;

		if (pos > pos_max)
		{//pos out of text: found nothing
			return false;
		}
		// h scans the text, it scans the pattern
		h = finder + ndl_size - 1;		
		it = begin(host(me), Standard())+ ndl_size -1 ; 

		unsigned int neq = 0;			// # mismatches
		unsigned int i = ndl_size;
		me.shift = s;				    // initialize shift

		while(i > 0 && neq <= me.k){
			// alignment of pattern and text
			// right-to-left-scan over the pattern 
			unsigned int char_t = *h;   // char in text
			unsigned int char_p = *it;  // char in pattern
			if(i >= s)
			{	// minimize over all component shifts
				unsigned int d_k = me.shift_table[char_t*(me.k+1)+i-s];
				me.shift = (d_k < me.shift) ? d_k : me.shift;
			} 

			if(char_t != char_p)
			{   // count mismatches
				neq+=1; 
			}
			--i;
			--it;  // proceed to left in pattern 
			--h;   // proceed to left in text
		
		}
		if(neq <= me.k)	{
			//found a hit
			return true;
		}
	
	}
	


}



//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, HammingHorspool> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, HammingHorspool> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
