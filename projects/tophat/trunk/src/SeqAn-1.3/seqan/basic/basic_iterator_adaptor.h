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

#ifndef SEQAN_HEADER_BASIC_ITERATOR_ADAPTOR_H
#define SEQAN_HEADER_BASIC_ITERATOR_ADAPTOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tag

//An iterator that adapts a std iterator to a default seqan iterator
template <typename TIterator, typename TSpec = Default>
struct AdaptorIterator;

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Iterator_Default_Imp<T, Rooted>
{
	typedef typename Iterator<T, Standard>::Type TStandardIterator;
	typedef Iter<T, AdaptorIterator<TStandardIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// Adaptor Iterator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Adaptor Iterator:
..cat:Iterators
..general:Class.Iter
..summary:Adapts iterators to @Concept.Rooted Iterator@.
..signature:Iter<TContainer, AdaptorIterator<TIterator [, TSpec]> >
..param.TContainer:Type of the container that can be iterated by $TIterator$.
...remarks:Use @Metafunction.Container@ to get the container type for a given iterator.
..param.TIterator:Type of the iterator that is adapted to @Concept.Rooted Iterator@.
..remarks.text:Adaptor iterators can implicitly converted to $TIterator$.
..include:seqan/basic.h
*/

template <typename TContainer, typename TIterator, typename TSpec>
class Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
{
public:
	typename Pointer_<TContainer>::Type data_container;
	TIterator data_iterator;
//____________________________________________________________________________

/**
.Memfunc.AdaptorIterator#Iter:
..class:Spec.Adaptor Iterator
..summary:Constructor
..signature:Iter()
..signature:Iter(iter)
..signature:Iter(container [, iterator])
..param.iter:Another adaptor iterator object.
..param.container:The corresponding container object.
..param.iterator:A iterator of $container$. (optional)
...remarks.text:If this argument is omitted, the adaptor iterator is initialized to the @Function.begin.begin iterator@ of $container$.
*/

public:
	Iter():
		data_container(0)
	{
SEQAN_CHECKPOINT
		data_iterator = TIterator();
	}
/*//TODO: welches "begin" zur initialisierung von "data_iterator" aufrufen?
	Iter(typename Parameter_<TContainer>::Type container_):
		data_container(_toPointer(container_)),
		data_iterator(begin(container_))
	{
SEQAN_CHECKPOINT
	}
*/
	Iter(typename Parameter_<TContainer>::Type container_, TIterator it_):
		data_container(_toPointer(container_)),
		data_iterator(it_)
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_iterator(other_.data_iterator)
	{
SEQAN_CHECKPOINT
	}
/*
	template <typename TSource>
	Iter(TSource & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
	template <typename TSource>
	Iter(TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
	}
*/

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & 
	operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		data_iterator = other_.data_iterator;
		return *this;
	}
/*
	template <typename TSource>
	Iter const & 
	operator = (TSource & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
	template <typename TSource>
	Iter const & 
	operator = (TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
*/
//____________________________________________________________________________


//____________________________________________________________________________

	operator TIterator () const
	{
SEQAN_CHECKPOINT
		return data_iterator;
	}

//____________________________________________________________________________
};

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return _toParameter<TContainer>(me.data_container);
}
template <typename TContainer, typename TIterator, typename TSpec>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return _toParameter<TContainer>(me.data_container);
}
//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec>
inline void
setContainer(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,	typename Parameter_<TContainer>::Type container_)
{
SEQAN_CHECKPOINT
   typedef Iter<TContainer, AdaptorIterator<TIterator, TSpec> > TIter;
	if (me.data_container && me.data_iterator != TIterator())
	{
		typename Position<TIter>::Type pos = position(me);
		me.data_container = _toPointer(container_);
		setPosition(me, pos);
	}
	else
	{	
		me.data_container = _toPointer(container_);
	}
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec>
inline TIterator &
hostIterator(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}
template <typename TContainer, typename TIterator, typename TSpec>
inline TIterator const &
hostIterator(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}
//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////
/*
template <typename TContainer, typename TIterator, typename TSpec>
struct Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>:
	Reference<TIterator const>
{
};
*/

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// position
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return hostIterator(me) - begin(container(me), Standard());
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec, typename TContainer2>
inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
		 TContainer2 const &)
{
SEQAN_CHECKPOINT
	return hostIterator(me) - begin(container(me), Standard());
}

//////////////////////////////////////////////////////////////////////////////
// setPosition
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TPosition>
inline void 
setPosition(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
			TPosition pos_)
{
SEQAN_CHECKPOINT
	hostIterator(me) = begin(container(me), Standard()) + pos_;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return value(hostIterator(me));
}
template <typename TContainer, typename TIterator, typename TSpec>
inline typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type 
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return value(hostIterator(me));
}

/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(hostIterator(me), _value);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
			TValue const & _value)
{
SEQAN_CHECKPOINT
	assignValue(hostIterator(me), _value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(hostIterator(me), _value);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(hostIterator(me), _value);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator == (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) == hostIterator(right);
}

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator == (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 typename IterComplementConst<TIterator>::Type const & right)
{
SEQAN_CHECKPOINT
 	return hostIterator(left) == right;
}

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator == (typename IterComplementConst<TIterator>::Type const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return left == hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator != (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) != hostIterator(right);
}

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator != (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			 typename IterComplementConst<TIterator>::Type const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) != right;
}

template <typename TContainer, typename TIterator, typename TSpec>
inline bool 
operator != (typename IterComplementConst<TIterator>::Type const & left,
			 Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return left != hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline void
goNext(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	goNext(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline void
goPrevious(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	goPrevious(hostIterator(me));
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) + right);
}
// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) + right);
}

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (TIntegral left,
			Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(right), hostIterator(right) + left);
}
// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator + (int left,
			Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(right), hostIterator(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator += (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}
// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator += (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 int right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator - (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) - right);
}
// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >  
operator - (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) - right);
}

//____________________________________________________________________________

template <typename TContainer, typename TIterator, typename TSpec, typename TContainer2, typename TIterator2, typename TSpec2>
inline typename Difference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type  
operator - (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
			Iter<TContainer2, AdaptorIterator<TIterator2, TSpec2> > const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator -= (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}
// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator -= (Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
			 int right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// atEnd
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIterator, typename TSpec>
inline bool
atEnd(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
SEQAN_CHECKPOINT
	return atEnd(me, container(me));
}

template <typename TContainer, typename TIterator, typename TSpec>
inline bool
atEnd(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return atEnd(me, container(me));
}

//////////////////////////////////////////////////////////////////////////////
// assign (Conversion)
//////////////////////////////////////////////////////////////////////////////

template <typename TTargetContainer, typename TIterator, typename TSpec, typename TSource>
inline void
assign(Iter<TTargetContainer, AdaptorIterator<TIterator, TSpec> > & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.data_container = container(source);
	target.data_iterator = begin(container(source)) + position(source);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
