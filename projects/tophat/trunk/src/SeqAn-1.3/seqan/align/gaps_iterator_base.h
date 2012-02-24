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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

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
	typedef typename Source<TGaps>::Type TSource_;
	typedef typename Iterator<TSource_, Rooted>::Type Type;
};
template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef typename Source<TGaps>::Type TSource_;
	typedef typename Iterator<TSource_, Rooted>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Value
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type TSource_;
	typedef typename Value<TSource_>::Type TSourceValue_;
	typedef typename GappedValueType<TSourceValue_>::Type Type;

};
template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > const>::Type TSource_;
	typedef typename Value<TSource_>::Type TSourceValue_;
	typedef typename GappedValueType<TSourceValue_>::Type Type;

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
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};
template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator_; //remove const
	typedef Proxy<IteratorProxy<TIterator_> > Type;
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
