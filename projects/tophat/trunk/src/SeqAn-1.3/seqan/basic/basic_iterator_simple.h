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

#ifndef SEQAN_HEADER_BASIC_ITERATOR_SIMPLE_H
#define SEQAN_HEADER_BASIC_ITERATOR_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Iter
//////////////////////////////////////////////////////////////////////////////

struct SimpleIterator;

/**
.Spec.SimpleIterator:
..cat:Iterators
..summary:A simple iterator.
..signature:Iter<TContainer, SimpleIterator>
..param.TContainer:Type of the container that can be iterated.
...metafunction:Metafunction.Container
..general:Class.Iter
..include:seqan/basic.h
*/
template <typename TContainer>
class Iter<TContainer, SimpleIterator>
{
public:
	typedef typename Value<TContainer>::Type TValue;
	TValue * data_ptr;

	Iter()
	{
	}
	Iter(Iter const & other_):
		data_ptr(other_.data_ptr)
	{
	}
	Iter(TValue * other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TContainer2>
	Iter(Iter<TContainer2, SimpleIterator> const & other_):
		data_ptr(other_.data_ptr)
	{
	}
	~Iter()
	{
	}
	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
		return *this;
	}
	Iter const &
	operator = (TValue * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}
	template <typename TContainer2>
	Iter const &
	operator = (Iter<TContainer2, SimpleIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
		return *this;
	}

	operator TValue * ()
	{
		return data_ptr;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Iterator_Default_Imp<T, Standard>
{
	typedef typename Value<T>::Type * Type;
//	typedef Iter<T, SimpleIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TContainer2>
inline typename Position<Iter<TContainer, SimpleIterator> const>::Type 
position(Iter<TContainer, SimpleIterator> const & me,
		 TContainer2 const & cont)
{
SEQAN_CHECKPOINT
	return me.data_ptr - begin(cont);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
