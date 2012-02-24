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

#ifndef SEQAN_HEADER_BASIC_ITERATOR_POSITION_H
#define SEQAN_HEADER_BASIC_ITERATOR_POSITION_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tag

struct PositionIterator;

//////////////////////////////////////////////////////////////////////////////
// Position Iterator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Position Iterator:
..cat:Iterators
..general:Class.Iter
..summary:Adapts @Metafunction.Position.position@ to @Concept.Rooted Iterator.iterator@.
..signature:Iter<TContainer, PositionIterator>
..param.TContainer:Type of the container.
...metafunction:Metafunction.Container
..remarks
...text:Position Iterators provide the concept @Concept.Rooted Iterator@.
..see:Metafunction.Position
..include:seqan/basic.h
*/

template <typename TContainer>
class Iter<TContainer, PositionIterator>
{
public:
	typedef typename Position<TContainer>::Type TPosition;

	typename Pointer_<TContainer>::Type data_container;
	TPosition data_position;
//____________________________________________________________________________

public:
/**
.Memfunc.PositionIterator#Iter:
..class:Spec.Position Iterator
..summary:Constructor
..signature:Iter()
..signature:Iter(iter)
..signature:Iter(container [, position])
..param.iter:Another position iterator object.
..param.container:The corresponding container object.
...metafunction:Metafunction.Container
..param.position:A position in $container$. (optional)
...metafunction:Metafunction.Position
...remarks.text:If this argument is omitted, the adaptor iterator is initialized to the @Function.beginPosition.begin position@ of $container$.
*/
	Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter(typename Parameter_<TContainer>::Type container_, TPosition position_ = 0):
		data_container(_toPointer(container_)),
		data_position(position_)
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_position(other_.data_position)
	{
SEQAN_CHECKPOINT
	}
	template <typename TContainer2, typename TSpec2>
	Iter(Iter<TContainer2, TSpec2> const & other_):
		data_container(container(other_)),
		data_position(position(other_))
	{
SEQAN_CHECKPOINT
	}
	~Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter const & 
	operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		data_position = other_.data_position;
		return *this;
	}

//____________________________________________________________________________
};

//____________________________________________________________________________

template <typename TContainer>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	return _toParameter<TContainer>(me.data_container);
}
template <typename TContainer>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, PositionIterator> const & me)
{
SEQAN_CHECKPOINT
	return _toParameter<TContainer>(me.data_container);
}
//____________________________________________________________________________

template <typename TContainer>
inline void
setContainer(Iter<TContainer, PositionIterator> & me, typename Parameter_<TContainer>::Type container_)
{
SEQAN_CHECKPOINT
	typedef Iter<TContainer, PositionIterator> TIter;
	typename Position<TIter>::Type pos = position(me);
	me.data_container = _toPointer(container_);
	setPosition(me, pos);
}

//____________________________________________________________________________

template <typename TContainer>
inline typename Position<TContainer>::Type &
position(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	return me.data_position;
}
template <typename TContainer>
inline typename Position<TContainer>::Type const &
position(Iter<TContainer, PositionIterator> const & me)
{
SEQAN_CHECKPOINT
	return me.data_position;
}
//____________________________________________________________________________

template <typename TContainer, typename TPosition>
inline void
setPosition(Iter<TContainer, PositionIterator> & me, TPosition position_)
{
SEQAN_CHECKPOINT
	me.data_position = position_;
}

//////////////////////////////////////////////////////////////////////////////
// value
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type 
value(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	return value(container(me), position(me));
}
template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type 
value(Iter<TContainer, PositionIterator> const & me)
{
SEQAN_CHECKPOINT
	return value(container(me), position(me));
}

/////////////////////////////////////////////////////////////////////////////
// assignValue
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> & me,
			TValue _value)
{
SEQAN_CHECKPOINT
	assignValue(container(me), position(me), _value);
}
template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> const & me,
			TValue _value)
{
SEQAN_CHECKPOINT
	assignValue(container(me), position(me), _value);
}

/////////////////////////////////////////////////////////////////////////////
// moveValue
//////////////////////////////////////////////////////////////////////////////

// manual forwards
template <typename T, typename TValue, typename TPos>
inline void moveValue(T & me, TPos pos, TValue const & _value);
template <typename T, typename TValue, typename TPos>
inline void moveValue(T const & me, TPos pos, TValue const & _value);

template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(container(me), position(me), _value);
}
template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> const & me,
		  TValue const & _value)
{
SEQAN_CHECKPOINT
	moveValue(container(me), position(me), _value);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator == (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) == position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator != (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) != position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator < / >
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator < (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) < position(right);
}

template <typename TContainer>
inline bool 
operator > (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) > position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator <= / >=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator <= (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) <= position(right);
}

template <typename TContainer>
inline bool 
operator >= (Iter<TContainer, PositionIterator> const & left,
			 Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) >= position(right);
}

//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goNext(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	setPosition(me, position(me) + 1);
}

//////////////////////////////////////////////////////////////////////////////
// goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goPrevious(Iter<TContainer, PositionIterator> & me)
{
SEQAN_CHECKPOINT
	setPosition(me, position(me) - 1);
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////
template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator + (Iter<TContainer, PositionIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) + right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator + (Iter<TContainer, PositionIterator> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) + right);
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator + (TIntegral left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(right), position(right) + left);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator + (int left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(right), position(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator> &
operator += (Iter<TContainer, PositionIterator> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) + right);
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator> &
operator += (Iter<TContainer, PositionIterator> & left,
			 int right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) + right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator>  
operator - (Iter<TContainer, PositionIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) - right);
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator>  
operator - (Iter<TContainer, PositionIterator> const & left,
			int right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, PositionIterator>(container(left), position(left) - right);
}

//____________________________________________________________________________

template <typename TContainer>
inline typename Difference<TContainer>::Type  
operator - (Iter<TContainer, PositionIterator> const & left,
			Iter<TContainer, PositionIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) - position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, PositionIterator> &
operator -= (Iter<TContainer, PositionIterator> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) - right);
	return left;
}
// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, PositionIterator> &
operator -= (Iter<TContainer, PositionIterator> & left,
			 int right)
{
SEQAN_CHECKPOINT
	setPosition(left, position(left) - right);
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// assign (Conversion)
//////////////////////////////////////////////////////////////////////////////

template <typename TTargetContainer, typename TSource>
inline void
assign(Iter<TTargetContainer, PositionIterator> & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.data_container = _toPointer(container(source));
	target.data_position = position(source);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
