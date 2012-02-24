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

#ifndef SEQAN_HEADER_BASIC_PROXY_H
#define SEQAN_HEADER_BASIC_PROXY_H


namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Proxy
//////////////////////////////////////////////////////////////////////////////
/**
.Class.Proxy:
..cat:Basic
..summary:Emulates object of another class.
..signature:Proxy<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..remarks.text:Use @Metafunction.Value@ to get the emulated type.
An instance of $Proxy$ behaves like an object of its value type.
$Proxy$ can be used as reference type (see @Metafunction.Reference@).
..remarks.text:Note that functions that are both general and specialized for 
the value type should be specialized for $Proxy<TSpec>$ too, 
since otherwise the general version will be called.
..include:seqan/basic.h
*/

template <typename TSpec>
struct Proxy;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Proxy

template <typename TSpec>
struct Spec< Proxy<TSpec> >
{
	typedef TSpec Type;
};
template <typename TSpec>
struct Spec< Proxy<TSpec> const>
{
	typedef TSpec Type;
};

template <typename TIterator>
struct IteratorProxy;

//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Proxy

template <typename TIterator>
struct Value< Proxy<IteratorProxy<TIterator> > >:
	Value<TIterator>
{
};
template <typename TIterator>
struct Value< Proxy<IteratorProxy<TIterator> > const>
{
	typedef typename Value<TIterator>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Proxy

template <typename TIterator>
struct GetValue< Proxy<IteratorProxy<TIterator> > >:
	GetValue<TIterator>
{
};
template <typename TIterator>
struct GetValue< Proxy<IteratorProxy<TIterator> > const>
{
	typedef typename GetValue<TIterator const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Proxy

template <typename TIterator>
struct Reference< Proxy<IteratorProxy<TIterator> > >
{
	typedef Proxy<IteratorProxy<TIterator> > Type;
};
template <typename TIterator>
struct Reference< Proxy<IteratorProxy<TIterator> > const >
{
	typedef Proxy<IteratorProxy<TIterator> > const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.Proxy

template <typename TIterator>
struct Size< Proxy<IteratorProxy<TIterator> > >:
	Size<TIterator>
{
};
template <typename TIterator>
struct Size< Proxy<IteratorProxy<TIterator> > const>:
	Size<TIterator>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Difference.param.T.type:Class.Proxy

template <typename TIterator>
struct Difference< Proxy<IteratorProxy<TIterator> > >:
	Difference<TIterator>
{
};
template <typename TIterator>
struct Difference< Proxy<IteratorProxy<TIterator> > const>:
	Difference<TIterator>
{
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
typename GetValue<Proxy<TSpec> >::Type
getValue(Proxy<TSpec> & me)
{
	return getValue(iter(me));
}
template <typename TSpec>
typename GetValue<Proxy<TSpec> const>::Type
getValue(Proxy<TSpec> const & me)
{
	return getValue(iter(me));
}

//////////////////////////////////////////////////////////////////////////////
// Iterator Proxy
//////////////////////////////////////////////////////////////////////////////
/**
.Spec.Iterator Proxy:
..cat:Proxies
..general:Class.Proxy
..summary:Proxy that is implemented by an iterator.
..signature:Proxy<IteratorProxy<TIterator> >
..param.TIterator:Iterator type.
..remarks.text:The value type of an iterator proxy is the value type of the
iterator $TIterator$.
..include:seqan/basic.h
*/

template <typename TIterator>
struct IteratorProxy;

//____________________________________________________________________________

template <typename TIterator>
struct Proxy<IteratorProxy<TIterator> >
{
public:
	typedef typename Value<Proxy>::Type TValue;
	typedef typename GetValue<Proxy>::Type TAccessor;

	typedef typename RemoveConst_<TAccessor>::Type TAccessor_NotConst;

public:
	TIterator data_iterator;

public:
	Proxy(TIterator const _it):
		data_iterator(_it)
	{
SEQAN_CHECKPOINT
	}
	Proxy(Proxy const & _other):
		data_iterator(_other.data_iterator)
	{
SEQAN_CHECKPOINT
	}

	~Proxy()
	{
SEQAN_CHECKPOINT
	}

	Proxy const &
	operator = (Proxy const & _other) const
	{
SEQAN_CHECKPOINT
		assignValue(data_iterator, getValue(_other.data_iterator));
		return *this;
	}

	Proxy const &
	operator = (TValue const & _value) const
	{
SEQAN_CHECKPOINT
		assignValue(data_iterator, _value);
		return *this;
	}

	operator TAccessor_NotConst()
	{
SEQAN_CHECKPOINT
		return getValue(data_iterator);
	}

	operator TAccessor_NotConst() const
	{
SEQAN_CHECKPOINT
		return getValue(data_iterator);
	}

//____________________________________________________________________________

	//not documented
//____________________________________________________________________________
};

template <typename TIterator>
inline TIterator &
iter(Proxy<IteratorProxy<TIterator> > & me)
{
	return me.data_iterator;
}
template <typename TIterator>
inline TIterator const &
iter(Proxy<IteratorProxy<TIterator> > const & me)
{
	return me.data_iterator;
}


//////////////////////////////////////////////////////////////////////////////
// Comparison
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename T>
struct CompareType <Proxy<TSpec>, T>
{
	typedef typename Value<Proxy<TSpec> >::Type TValue;
	typedef typename RemoveConst_<TValue>::Type TValue_NoConst;
	typedef typename CompareType<TValue_NoConst, T>::Type Type;
};

//???TODO: Symmetrie von CompareType herstellen
//____________________________________________________________________________

template <typename TTarget, typename T, typename TSpec>
inline typename Convert<TTarget, Proxy<TSpec> >::Type
convertImpl(Convert<TTarget, T> const,
			Proxy<TSpec> & source)
{
	return convert<TTarget>(getValue(source));
}
template <typename TTarget, typename T, typename TSpec>
inline typename Convert<TTarget, Proxy<TSpec> const>::Type
convertImpl(Convert<TTarget, T> const,
			Proxy<TSpec> const & source)
{
	return convert<TTarget>(getValue(source));
}
//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename TSpec, typename TRight>
inline bool
operator == (Proxy<TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator == (TLeft const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator == (Proxy<TLeftSpec> const & left_, 
			 Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator == (Proxy<TSpec> const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) == convert<TAccessor>(right_);
}

//____________________________________________________________________________
// operator !=

template <typename TSpec, typename TRight>
inline bool
operator != (Proxy<TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator != (TLeft const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator != (Proxy<TLeftSpec> const & left_, 
			 Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator != (Proxy<TSpec> const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) != convert<TAccessor>(right_);
}


//____________________________________________________________________________
// operator <

template <typename TSpec, typename TRight>
inline bool
operator < (Proxy<TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator < (TLeft const & left_, 
			Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator < (Proxy<TLeftSpec> const & left_, 
			Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator < (Proxy<TSpec> const & left_, 
			Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) < convert<TAccessor>(right_);
}

//____________________________________________________________________________
// operator <=

template <typename TSpec, typename TRight>
inline bool
operator <= (Proxy<TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator <= (TLeft const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator <= (Proxy<TLeftSpec> const & left_, 
			 Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator <= (Proxy<TSpec> const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) <= convert<TAccessor>(right_);
}


//____________________________________________________________________________
// operator >

template <typename TSpec, typename TRight>
inline bool
operator > (Proxy<TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator > (TLeft const & left_, 
			Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator > (Proxy<TLeftSpec> const & left_, 
			Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator > (Proxy<TSpec> const & left_, 
			Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) > convert<TAccessor>(right_);
}


//____________________________________________________________________________
// operator >=

template <typename TSpec, typename TRight>
inline bool
operator >= (Proxy<TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator >= (TLeft const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator >= (Proxy<TLeftSpec> const & left_, 
			 Proxy<TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TLeftSpec> TLeft;
	typedef Proxy<TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator >= (Proxy<TSpec> const & left_, 
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
	return convert<TAccessor>(left_) >= convert<TAccessor>(right_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TSpec>
inline TStream &
operator >> (TStream & strm,
			 Proxy<TSpec> & proxy)
{
	typedef Proxy<TSpec> TProxy;
	typedef typename Value<TProxy>::Type TValue;
	TValue temp;
	strm >> temp;
	assignValue(iter(proxy), temp);
	return strm;
}
template <typename TStream, typename TSpec>
inline TStream &
operator >> (TStream & strm,
			 Proxy<TSpec> const& proxy)
{
	typedef Proxy<TSpec> TProxy;
	typedef typename Value<TProxy>::Type TValue;
	TValue temp;
	strm >> temp;
	assignValue(iter(proxy), temp);
	return strm;
}


template <typename TStream, typename TSpec>
inline TStream &
operator << (TStream & strm,
			 Proxy<TSpec> & proxy)
{
	return strm << getValue(proxy);
}
template <typename TStream, typename TSpec>
inline TStream &
operator << (TStream & strm,
			 Proxy<TSpec> const & proxy)
{
	return strm << getValue(proxy);
}


//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
