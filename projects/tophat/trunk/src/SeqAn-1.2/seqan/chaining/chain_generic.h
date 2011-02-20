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
  $Id: chain_generic.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CHAIN_GENERIC_H
#define SEQAN_HEADER_CHAIN_GENERIC_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TWeight>
struct _Chain_Generic_Entry
{
	TPos me;		//position of fragment here (within Source)
	TPos pre;		//position of precursor or -1 for top (within Frags)
	TWeight weight;	//weight for best chain until here
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
struct _Chain_Generic_SortFragsPredFunctional
{
	TSource & source;
	_Chain_Generic_SortFragsPredFunctional(TSource & src)
		: source(src)
	{
	}
	_Chain_Generic_SortFragsPredFunctional(_Chain_Generic_SortFragsPredFunctional const & other)
		: source(other.source)
	{
	}
	inline _Chain_Generic_SortFragsPredFunctional &
	operator = (_Chain_Generic_SortFragsPredFunctional const & other)
	{
		source = other.source;
		return *this;
	}
	~_Chain_Generic_SortFragsPredFunctional()
	{
	}

	template <typename TFrags>
	inline bool 
	operator() (TFrags & left, TFrags & right)
	{
		return leftPosition(source[left.me], 0) < leftPosition(source[right.me], 0);
	}
};

//____________________________________________________________________________


template <typename TSource, typename TFrags, typename TScoring>
inline void
_chain_generic_initFrags(TSource & source,
						 TFrags & frags,
						 TScoring scoring)
{
	SEQAN_ASSERT(length(source))

	typedef typename Position<TSource>::Type TPos;
	typedef typename Value<TSource>::Type TFragment;
	typedef typename Value<TFrags>::Type TFrag;

	//create top fragment
	TFragment top(dimension(source[0]));
	makeBeginFragment(top);

	//create entry in frags for each item in source
	resize(frags, length(source));
	for (TPos i = beginPosition(source); i < endPosition(source); ++i)
	{
		TFrag & frag = frags[i];
		frag.me = i;
		frag.pre = ~0UL; //link it with the top fragment
		frag.weight = scoreChainGap(scoring, top, source[i]) + weight(source[i]);
	}

	std::sort(begin(frags, Standard()), end(frags, Standard()), _Chain_Generic_SortFragsPredFunctional<TSource>(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFrag>
inline bool
_chain_generic_chainable(TFrag & f1,
						 TFrag & f2)
{
	SEQAN_ASSERT(dimension(f1) == dimension(f2))

	unsigned int dim = dimension(f1); 
	while (dim > 0)
	{
		--dim;
		if (rightPosition(f1, dim) > leftPosition(f2, dim))
		{
			return false;
		}
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TFrags, typename TIterator, typename TScoring>
inline void
_chain_generic_findBest(TSource & source,
						TFrags & frags,
						TIterator & it_act,
						TScoring scoring)
{
	typedef typename Iterator<TFrags, Standard>::Type TFragsIterator;
	typedef typename Value<TScoring>::Type TWeight;
	typedef typename Reference<TIterator>::Type TFragRef;

	TFragRef act = *it_act;
	TWeight act_weight = weight(source[act.me]);
	TFragsIterator it_begin = begin(frags, Standard());
	for (TFragsIterator it = it_begin; it < it_act; ++it)
	{
		TFragRef frag = *it;
		if (_chain_generic_chainable(source[frag.me], source[act.me]))
		{
			TWeight score = frag.weight + scoreChainGap(scoring, source[frag.me], source[act.me]) + act_weight;
			if (score > act.weight)
			{//better predecessor found
				act.pre = it - it_begin;
				act.weight = score;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TFrags, typename TIterator, typename TDest, typename TFragment>
inline void
_chain_generic_Backtrace(TSource & source,
						 TFrags & frags,
						 TIterator & it_best,
						 TDest & dest,
						 TFragment & bottom)
{
	typedef typename Position<TSource>::Type TPos;

	clear(dest);
	//build chain in reverse order
	appendValue(dest, bottom); //chain will end with bottom fragment

	appendValue(dest, source[(*it_best).me]);
	for (TPos pos = (*it_best).pre; pos != ~0UL; pos = frags[pos].pre)
	{
		appendValue(dest, source[frags[pos].me]);
	}

	//chain will start with top fragment
	TFragment top(dimension(bottom));
	makeBeginFragment(top);
	appendValue(dest, top);

	//reverse chain
	std::reverse(begin(dest, Standard()), end(dest, Standard()));
}


//////////////////////////////////////////////////////////////////////////////

//spec for GenericChaining
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
globalChaining(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring, 
	  GenericChaining)
{
	typedef typename Value<TSource>::Type TFragment;
	typedef typename Weight<TFragment>::Type TWeight;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef _Chain_Generic_Entry<TSourcePosition, TWeight> TFrag;
	typedef String<TFrag> TFrags;
	typedef typename Iterator<TFrags, Standard>::Type TFragsIterator;

	//initialize fragments
	TFrags frags;
	_chain_generic_initFrags(source, frags, scoring);

	TFragsIterator it_begin = begin(frags, Standard());
	TFragsIterator it_end = end(frags, Standard());
	TFragsIterator it_best = it_begin;
	TWeight weight_best = InfimumValue<TWeight>::VALUE;


	//create bottom fragment
	unsigned int dim = dimension(source[0]);
	TFragment bottom(dim);
	makeEndFragment(bottom, source);

	//iterate all fragments
	for (TFragsIterator it = it_begin; it < it_end; ++it)
	{
		//find best predecessor for *it
		_chain_generic_findBest(source, frags, it, scoring);

		//determine heaviest fragment
		TWeight weight_it = (*it).weight + scoreChainGap(scoring, source[(*it).me], bottom);

		if (weight_it > weight_best)
		{
			it_best = it;
			weight_best = weight_it;
		}
	}

	//follow best fragment back to the beginning of the chain
	_chain_generic_Backtrace(source, frags, it_best, dest, bottom);

	return weight_best;
}



//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
