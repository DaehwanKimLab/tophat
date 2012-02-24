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

#ifndef SEQAN_HEADER_BASIC_DEFINITION_H
#define SEQAN_HEADER_BASIC_DEFINITION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Tag
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Length;

template <>
struct Length<void>
{
	enum { VALUE = 0 };
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.TagList:
..summary:A structure to represent a list of tags.
..signature:TagList<TTag1>
..signature:TagList<TTag1, TagList<TTag2> >
..signature:TagList<TTag1, TagList<TTag2, TagList<TTag3[...]> > >
..param.TTag1:The first tag of the list.
..param.TTag2:The second tag of the list.
..param.TTag3:The third tag of the list.
..include:seqan/basic.h
*/

template <typename TTag = void, typename TSubList = void>
struct TagList
{
	typedef TTag Type;
};

template <typename TTag>
struct Length< TagList<TTag, void> > {
	enum { VALUE = 1 };
};

template <typename TTag, typename TSubList>
struct Length< TagList<TTag, TSubList> > {
	enum { VALUE = Length<TSubList>::VALUE + 1 };
};

template <typename TTagList = void>
struct TagSelector
{
	int tagId;
	
	TagSelector()
	{
		tagId = 0;
	}

    inline bool operator==(TagSelector const & other) const
    {
        return other.tagId == tagId;
    }
};

/**
.Class.TagSelector:
..summary:A structure to select a tag from a @Tag.TagList@.
..signature:TagSelector<TTagList>
..param.TTagList:A tag list.
...type:Tag.TagList
.Memvar.TagSelector#tagId:
..class:Class.TagSelector
..type:nolink:int
..cat:Basic
..summary:Stores the index of a @Page.Glossary.Tag@ in the tag list.
..include:seqan/basic.h
*/

///

template <typename TTag, typename TSubList>
struct TagSelector< TagList<TTag, TSubList> >:
	TagSelector<TSubList>
{
	typedef TTag					Type;
	typedef TagSelector<TSubList>	Base;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Default:
..summary:Tag that specifies default behavior.
..tag.Default:Use default behavior. 
..include:seqan/basic.h
*/
struct Default_;
typedef Tag<Default_> const Default;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Move Switch:
..summary:Switch to force move.
..tag.Move:Move instead of assign. 
..remarks.text:The difference between move constructor and copy constructor
is that the source object is not copied but moved into the target object.
The source object can lose its content and will be empty after
this operation in this case.
A move constructor can sigificantly faster than a copy constructor.
..example.code:String source("hello");
String target(source, Move()); // source is moved to target
std::cout << source; //nothing printed since source lost content
std::cout << target; //"hello"
..see:Function.move
..include:seqan/basic.h
*/

struct Move_;
typedef Tag<Move_> const Move;

//////////////////////////////////////////////////////////////////////////////

//Pass to c'tor of iterator to move it to the end
struct GoEnd_;
typedef Tag<GoEnd_> const GoEnd;


//////////////////////////////////////////////////////////////////////////////

//construct without initializing
struct MinimalCtor_;
typedef Tag<MinimalCtor_> const MinimalCtor;

//construct with initializing
struct NonMinimalCtor_;
typedef Tag<NonMinimalCtor_> const NonMinimalCtor;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Nothing:
..summary:Tag that represents an absent parameter or an absent type.
..tag.Nothing:Omit parameter.
..include:seqan/basic.h
*/
///Empty Data Class.
struct Nothing {};



//////////////////////////////////////////////////////////////////////////////
// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct CopyConst_
{
	typedef TTo Type;
};
template <typename TFrom, typename TTo>
struct CopyConst_<TFrom const, TTo>
{
	typedef TTo const Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Internal.RemoveConst_:
..signature:RemoveConst_<T>
..returns:$t$ if $T$ is $t const$, otherwise $T$.
*/
template <typename T>
struct RemoveConst_
{
	typedef T Type;
};
template <typename T>
struct RemoveConst_<T const>:
	public RemoveConst_<T> {};

template <typename T>
struct RemoveConst_<T &>
{
	typedef typename RemoveConst_<T>::Type & Type;
};
template <typename T>
struct RemoveConst_<T *>
{
	typedef typename RemoveConst_<T>::Type * Type;
};
template <typename T, size_t I>
struct RemoveConst_<T const [I]>
{
	typedef T * Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Internal.MakeUnsigned_:
..signature:MakeUnsigned_<T>
..returns:$unsigned t$ if $T$ is not $unsigned t$, otherwise $T$.
*/
template <typename T>
struct MakeUnsigned_
{
	typedef
		typename If< IsSameType<T, char>::VALUE,         unsigned char,
		typename If< IsSameType<T, signed char>::VALUE,  unsigned char,
		typename If< IsSameType<T, signed short>::VALUE, unsigned short,
		typename If< IsSameType<T, signed int>::VALUE,   unsigned int,
		typename If< IsSameType<T, signed long>::VALUE,  unsigned long,
		typename If< IsSameType<T, __int64>::VALUE,      __uint64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeUnsigned_<T const> {
	typedef typename MakeUnsigned_<T>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Internal.MakeSigned_:
..signature:MakeSigned_<T>
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct MakeSigned_
{
	typedef
		typename If< IsSameType<T, char>::VALUE,           signed char,
		typename If< IsSameType<T, unsigned char>::VALUE,  signed char,
		typename If< IsSameType<T, unsigned short>::VALUE, signed short,
		typename If< IsSameType<T, unsigned int>::VALUE,   signed int,
		typename If< IsSameType<T, unsigned long>::VALUE,  signed long,
		typename If< IsSameType<T, __uint64>::VALUE,       __int64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeSigned_<T const> {
	typedef typename MakeSigned_<T>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Internal.ClassIdentifier_:
..signature:void * ClassIdentifier_<T>::getID()
..returns:A void * that identifies $T$.
...text:The returned values of two calls of $getID$ are equal if and only if
the used type $T$ was the same.
*/
template <typename T>
struct ClassIdentifier_
{
	static inline void *
	getID()
	{
SEQAN_CHECKPOINT
		static bool _id_dummy;
		return &_id_dummy;
	}
};

//////////////////////////////////////////////////////////////////////////////
/**
.Function.log2:
..cat:Miscellaneous
..summary:Computes logarithm of base 2 for integer types
..signature:unsigned int log2(i)
..param.i:An integer type.
..returns:The largest integer smaller or equal than
the logarithm of $i$.
..include:seqan/basic.h
*/

template <int BITS_MAX>
struct Log2Impl_
{
	template <typename T>
	static inline unsigned int
	log2(T val, unsigned int offset)
	{
		unsigned int val2 = val >> (BITS_MAX / 2);
		if (val2)
		{
			val = val2;
			offset += BITS_MAX / 2;
		}
		return Log2Impl_<BITS_MAX / 2>::log2(val, offset);
	}
};

template <>
struct Log2Impl_<1>
{
	template <typename T>
	static inline unsigned int
	log2(T /*val*/, unsigned int offset)
	{
		return offset;
	}
};


template <typename T>
inline unsigned int
log2(T val)
{
	enum
	{
//		BITS_PER_VALUE = BitsPerValue<T>::VALUE //TODO???
		BITS_PER_VALUE = sizeof(T) * 8
	};

	return Log2Impl_<BITS_PER_VALUE>::log2(val, 0);
}

template <typename TValue, typename TExponent>
inline TValue _intPow(TValue a, TExponent b)
{
SEQAN_CHECKPOINT
	TValue ret = 1;
	while (b != 0)
	{
		if (b & 1) ret *= a;
		a *= a;
		b >>= 1;
	}	
	return ret;
}

//////////////////////////////////////////////////////////////////////////////
// to avoid conflicts with non-standard macros and namespaces
// we define our own Min/Max functions

template<typename Tx_> inline
const Tx_& _min(const Tx_& _Left, const Tx_& Right_)
{	// return smaller of _Left and Right_
	if (_Left < Right_)
		return _Left;
	else
		return Right_;
}

template<typename Tx_, typename Ty_> inline
Tx_ _min(const Tx_& _Left, const Ty_& Right_)
{	// return smaller of _Left and Right_
    return (Right_ < _Left ? Right_ : _Left);
}

template<typename Ty_> inline
const Ty_& _max(const Ty_& _Left, const Ty_& Right_)
{	// return larger of _Left and Right_
	if (_Left < Right_)
		return Right_;
	else
		return _Left;
}

template<typename Tx_, typename Ty_> inline
Tx_ _max(const Tx_& _Left, const Ty_& Right_)
{	// return smaller of _Left and Right_
    return (Right_ < _Left ? _Left : Right_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline
T _abs(T const & x)
{
    if (x < static_cast<T>(0))
        return -x;
    else
        return x;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
inline bool 
_isSameType()
{
	return IsSameType<T1, T2>::VALUE;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

