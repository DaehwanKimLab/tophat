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
  $Id: gaps_iterator_base.h 2431 2008-07-09 12:51:10Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GAPS_ITERATOR_BASE_H
#define SEQAN_HEADER_GAPS_ITERATOR_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Source
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<TGaps>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type Type;
};
template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef typename Source<TGaps>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Value
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type TSource;
	typedef typename Value<TSource>::Type TSourceValue;
	typedef typename GappedValueType<TSourceValue>::Type Type;

};
template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > const>::Type TSource;
	typedef typename Value<TSource>::Type TSourceValue;
	typedef typename GappedValueType<TSourceValue>::Type Type;

};

//////////////////////////////////////////////////////////////////////////////
// GetValue
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct GetValue<Iter<TGaps, GapsIterator<TSpec> > >:
	Value<Iter<TGaps, GapsIterator<TSpec> > > //no reference
{
};
template <typename TGaps, typename TSpec>
struct GetValue<Iter<TGaps, GapsIterator<TSpec> > const>:
	Value<Iter<TGaps, GapsIterator<TSpec> > const> //no reference
{
};

//////////////////////////////////////////////////////////////////////////////
// Reference
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};
template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator; //remove const
	typedef Proxy<IteratorProxy<TIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////
//specialization of basic_proxy.h/CompareType implementation

/*
template <typename TGaps, typename TSpec, typename T>
struct CompareType<Proxy<IteratorProxy<Iter<TGaps, GapsIterator<TSpec> > > >, T>
{
	typedef T Type;
};
*/

//???TODO: Symmetrie von CompareType herstellen


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const & 
operator ++(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	goNext(me);
	return me;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const 
operator ++(Iter<TGaps, GapsIterator<TSpec> > const & me, int)
{
SEQAN_CHECKPOINT
	Iter<TGaps, GapsIterator<TSpec> > ret = me;
	goNext(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const & 
operator --(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	goPrevious(me);
	return me;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > const 
operator --(Iter<TGaps, GapsIterator<TSpec> > const & me, int)
{
SEQAN_CHECKPOINT
	Iter<TGaps, GapsIterator<TSpec> > ret = me;
	goPrevious(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////
//todo???: weitere operatoren


//////////////////////////////////////////////////////////////////////////////

// insert a single gap at current iterator position
template <typename TGaps, typename TSpec>
inline void
insertGap(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	insertGaps(me, 1);
}

//////////////////////////////////////////////////////////////////////////////

// remove a single gap at current iterator position
template <typename TGaps, typename TSpec>
inline void
removeGap(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	removeGaps(me, 1);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
inline typename Reference< Iter<TGaps, GapsIterator<TSpec> > >::Type
value(Iter<TGaps, GapsIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Reference< Iter<TGaps, GapsIterator<TSpec> > >::Type TReference;
	return TReference(me);
}
template <typename TGaps, typename TSpec>
inline typename Reference< Iter<TGaps, GapsIterator<TSpec> > const>::Type
value(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Reference< Iter<TGaps, GapsIterator<TSpec> > const>::Type TReference;
	return TReference(me);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > & me,
			TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
	//else insert??? das waere neu.
}
template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > const & me,
			TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
}
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec, typename TValue>
inline void
moveValue(Iter<TGaps, GapsIterator<TSpec> > & me,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		moveValue(source(me), val);
	}
}
template <typename TGaps, typename TSpec, typename TValue>
inline void
moveValue(Iter<TGaps, GapsIterator<TSpec> > const & me,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	if (!isGap(me)) 
	{
		moveValue(source(me), val);
	}
}

//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________

template <typename TGaps, typename TSpec>
inline TGaps &
container(Iter<TGaps, GapsIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return *me.data_container;
}

template <typename TGaps, typename TSpec>
inline TGaps &
container(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return *me.data_container;
}

//____________________________________________________________________________

template <typename TGaps, typename TSpec>
inline typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TGaps, typename TSpec>
inline typename Source<Iter<TGaps, GapsIterator<TSpec> > const>::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//todo??? setSource? setContainer?
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec, typename TSize>
inline void
goFurther(Iter<TGaps, GapsIterator<TSpec> > & me,
		  TSize steps)
{
SEQAN_CHECKPOINT
	while (steps)
	{
		goNext(me);
		--steps;
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
