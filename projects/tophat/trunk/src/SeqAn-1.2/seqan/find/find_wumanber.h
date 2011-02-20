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
  $Id: find_wumanber.h 4763 2009-09-03 14:11:58Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_WUMANBER_H
#define SEQAN_HEADER_FIND_WUMANBER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// WuManber
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.WuManber:
..general:Class.Pattern
..cat:Searching
..summary:Online-algorithm for multi-pattern search.
..signature:Pattern<TNeedle, WuManber>
..param.TNeedle:The needle type.
...type:Class.String
*/

///.Class.Pattern.param.TSpec.type:Spec.WuManber

struct _WuManber;
typedef Tag<_WuManber> WuManber;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, WuManber> 
{
//____________________________________________________________________________
public:
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Position<TNeedle>::Type TNeedlePosition;
	typedef typename Size<TKeyword>::Type TSize;
	typedef typename Value<TKeyword>::Type TValue;

	//searching data: these members are initialized in _patternInit
	TNeedlePosition position; //last found keyword
	TNeedlePosition * to_verify_begin; //next entry in verify
	TNeedlePosition * to_verify_end; //end of list in verify
			//note: if to_verify_begin == to_verify_end then searching in Haystack must go on

	//preprocessed data: these members are initialized in setHost
	Holder<TNeedle> needle;
	String<TNeedlePosition> verify_tab; //table of keywords to verify depending on the last value (HASH)
	String<TNeedlePosition *> verify; //directory into verify_tab
	String<TSize> shift; //table of skip widths (SHIFT)

	TSize lmin;	//min length of keyword
	unsigned char q; //length of hashed q-grams
//____________________________________________________________________________

	Pattern():
		lmin(0)
	{
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		SEQAN_CHECKPOINT
		setHost(*this, ndl);
	}

	~Pattern() 
	{
	}
//____________________________________________________________________________

private:
	Pattern(Pattern const& other);
	Pattern const & operator=(Pattern const & other);

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

//forward
template <typename TNeedle, int Q>
struct _WuManber_Hash;


//////////////////////////////////////////////////////////////////////////////
//implementation kernel of WuManber 

template <typename TNeedle, int Q>
struct _WuManber_Imp
{
//____________________________________________________________________________

	typedef Pattern<TNeedle, WuManber> TPattern;

	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Position<TNeedle>::Type TNeedlePosition;
	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;

	typedef typename Size<TKeyword>::Type TSize;
	typedef typename Value<TKeyword>::Type TValue;
	typedef typename Iterator<TKeyword, Standard>::Type TIterator;

//____________________________________________________________________________

	enum 
	{
		C = BitsPerValue<TValue>::VALUE, //bits per value
		W = (C*Q <= 16) ? C*Q : 16,	//bits per hash values

		DIR_SIZE = 1 << W,			//size of verify_dir and shift_dir

		//shift width for Q = 2:
		SHIFT = (W > 2*C) ? C : W-C,			

		//shift widths for Q = 3:
		SHIFT_LEFT = (W > 3*C) ? 2*C : W-C,
		SHIFT_MIDDLE = SHIFT_LEFT / 2
	};
//____________________________________________________________________________

	inline static unsigned short
	hash(TValue * vals)
	{
		return _WuManber_Hash<TNeedle, Q>::hash(vals);
	}
//____________________________________________________________________________

	inline static void
	initialize(TPattern & me)
	{
		//note: me.needle, me.q, and me.lmin were already set in setHost before initialize was called

		//some variables
		TNeedleIterator pit;
		TNeedleIterator pit_end = end(needle(me));

		TSize k = length(needle(me));

		//resize and init tables
		resize(me.verify_tab, (unsigned) k);
		resize(me.verify, (unsigned) DIR_SIZE+1);
		resize(me.shift, (unsigned) DIR_SIZE);

		arrayFill(begin(me.shift), begin(me.shift) + (unsigned) DIR_SIZE, me.lmin-Q+1); //maximal shift width is me.lmin-B+1

		//init counters
		unsigned int verify_count[DIR_SIZE];
		arrayFill(verify_count, verify_count + DIR_SIZE, 0);

		//1: first scan: fill me.shift and count verify
		for (pit = begin(needle(me)); pit != pit_end; ++pit)
		{
			unsigned short hash;
			TIterator kit = begin(*pit);
			for (unsigned int i = 0; i <= me.lmin-Q; ++i)
			{
				hash = _WuManber_Hash<TNeedle, Q>::hash(kit + i);
				if (me.shift[hash] > me.lmin-Q - i)
				{
					me.shift[hash] = me.lmin-Q - i;
				}
			}
			++(verify_count[hash]);
		}

		//2: add up and convert to pointers
		TNeedlePosition * verify_ptr = begin(me.verify_tab);
		
		me.verify[0] = verify_ptr;
		unsigned int sum = 0;
		for (unsigned int i = 0; i < DIR_SIZE; ++i)
		{
			me.verify[i+1] = verify_ptr + sum;
			sum += verify_count[i];
		}

		//3: second scan: fill verify and shift
		unsigned int i = 0;
		for (pit = begin(needle(me)); pit != pit_end; ++pit)
		{
			unsigned short hash_plus_1;
			hash_plus_1 = _WuManber_Hash<TNeedle, Q>::hash(begin(*pit) + me.lmin-Q) + 1;

			//write into verify_tab
			*(me.verify[hash_plus_1]) = i;
			++i;			

			//update verify
			++me.verify[hash_plus_1];
		}
	}
//____________________________________________________________________________

	template <typename TFinder>
	static inline bool
	find(TFinder & finder, TPattern & me) 
	{
		typedef typename Haystack<TFinder>::Type THaystack;
		typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;

		THaystackIterator tit;
		THaystackIterator haystack_end = end(haystack(finder));
		THaystackIterator tit_end = haystack_end - Q + 1;
		unsigned short hash;

		if (empty(finder)) 
		{
//START
			_patternInit(me);
			_finderSetNonEmpty(finder);
			tit = hostIterator(finder) + me.lmin-Q;
		} 
		else 
		{
//RESUME
			tit = hostIterator(finder) + me.lmin-Q;
			goto VERIFY;
		}

//SEARCH
		while (tit < tit_end)
		{
			hash = _WuManber_Hash<TNeedle, Q>::hash(tit);

			if (me.shift[hash])
			{
//SHIFT
				tit += me.shift[hash];
			}
			else
			{
				me.to_verify_begin = me.verify[hash];
				me.to_verify_end = me.verify[hash+1];
//VERIFY
VERIFY:
				while (me.to_verify_begin != me.to_verify_end)
				{
					me.position = *me.to_verify_begin;
					++me.to_verify_begin;

					TKeyword & kw = needle(me)[me.position];

					TIterator pit = begin(kw, Standard());
					TIterator pit_end = end(kw, Standard());
					THaystackIterator tit2 = tit - (me.lmin-Q);

					if ((pit_end - pit) > (haystack_end - tit2)) continue; //rest of haystack too short

					while (true)
					{
						if (pit == pit_end)
						{
//FOUND
							setPosition(finder, tit - (me.lmin-Q) - begin(haystack(finder), Standard()));
							_setFinderLength(finder, length(kw));
							_setFinderEnd(finder, position(finder) + length(finder));
							return true;
						}
						if (*pit != *tit2) break;
						++pit;
						++tit2;
					}
				}

				++tit;
			}
		}
//END
		return false;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
//implementation of hash function

template <typename TNeedle>
struct _WuManber_Hash<TNeedle, 1>
{
	template <typename TIterator>
	inline static unsigned short
	hash(TIterator vals)
	{
		return ordValue(*vals);
	}
};
template <typename TNeedle>
struct _WuManber_Hash<TNeedle, 2>
{
	template <typename TIterator>
	inline static unsigned short
	hash(TIterator vals)
	{
		return ordValue(*vals)
			+ (ordValue(*(vals+1)) << _WuManber_Imp<TNeedle, 2>::SHIFT);
	}
};
template <typename TNeedle>
struct _WuManber_Hash<TNeedle, 3>
{
	template <typename TIterator>
	inline static unsigned short
	hash(TIterator vals)
	{
		return ordValue(*vals) 
			+ (ordValue(*(vals+1)) << _WuManber_Imp<TNeedle, 3>::SHIFT_MIDDLE)
			+ (ordValue(*(vals+2)) << _WuManber_Imp<TNeedle, 3>::SHIFT_LEFT);
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void _setHost_WuManber(Pattern<TNeedle, WuManber> & me, 
					   TNeedle2 const & needle_)
{
SEQAN_CHECKPOINT
SEQAN_ASSERT(!empty(needle_));

	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Value<TKeyword>::Type TValue;
	typedef typename Size<TKeyword>::Type TSize;

	//me.needle
	setValue(me.needle, needle_);

	//determine lmin
	me.lmin = length(needle(me)[0]);
	for (TNeedleIterator it = begin(needle(me)) + 1; it != end(needle(me)); ++it)
	{
		TSize len = length(*it);
		if (len < me.lmin)
		{
			me.lmin = len;
		}
	}

	if (me.lmin == 0) return;

	//compute q:
	unsigned int C = BitsPerValue<TValue>::VALUE;
	if (C > 12)
	{
		me.q = 1;
	}
	else 
	{
		//according to Wu & Manber: B = log_c(2mk) "is a good value"
		//i.e. C^B = 2mk
		//m = lmin
		//k = length(needle)
		//our heuristic: take B = 2 if C^2 >= mk, else B = 3

		if (C * C >= me.lmin * length(needle(me)))
		{
			me.q = 2;
		}
		else
		{
			me.q = 3;
		}
	}
	if (me.q > me.lmin) 
	{
		me.q = me.lmin;
	}

	//rest of preprocessing is done in _WuManber_Imp
	if (me.q == 2) _WuManber_Imp<TNeedle, 2>::initialize(me);
	else if (me.q == 3) _WuManber_Imp<TNeedle, 3>::initialize(me);
	else _WuManber_Imp<TNeedle, 1>::initialize(me);
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, WuManber> & me, 
			  TNeedle2 const & needle) 
{
	_setHost_WuManber(me, needle);
}

template <typename TNeedle, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, WuManber> & me, 
		TNeedle2 & needle)
{
	_setHost_WuManber(me, needle);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, WuManber> >::Type & 
host(Pattern<TNeedle, WuManber> & me)
{
SEQAN_CHECKPOINT
	return value(me.needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, WuManber> const>::Type & 
host(Pattern<TNeedle, WuManber> const & me)
{
SEQAN_CHECKPOINT
	return value(me.needle);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, WuManber> & me)
{
	return me.position;
}

//////////////////////////////////////////////////////////////////////////////

//called when search begins
template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, WuManber> & me) 
{
SEQAN_CHECKPOINT
	me.to_verify_begin = 0;
	me.to_verify_end = 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, 
				 Pattern<TNeedle, WuManber> & me) 
{
SEQAN_CHECKPOINT

	if (me.lmin == 0) return false;

	if (me.q == 2) return _WuManber_Imp<TNeedle, 2>::find(finder, me);
	else if (me.q == 3) return _WuManber_Imp<TNeedle, 3>::find(finder, me);
	else return _WuManber_Imp<TNeedle, 1>::find(finder, me);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
