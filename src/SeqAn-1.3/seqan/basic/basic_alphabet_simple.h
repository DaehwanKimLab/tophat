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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of SimpleType, the biological sequence types and conversion
// tables.
// ==========================================================================

#ifndef SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_H
#define SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Class that is used for various simple value types
//////////////////////////////////////////////////////////////////////////////

/**
.Class.SimpleType:
..cat:Basic
..summary:Implementation for "simple" types.
..signature:SimpleType<TValue, TSpec>
..param.TValue:Type that stores the values of an instance.
...remarks:TValue must be a simple type.
...metafunction:Metafunction.Value
..param.TSpec:Specialization tag.
...metafunction:Metafunction.Spec
..remarks:
...text:A "simple type" is a C++ type that can be constructed without constructor,
destructed without destructor and copied without copy constructor or assignment operator.
All basic types (like $char$, $int$ or $float$) are simple. Pointers, references and arrays of
simple types are simple.
POD types ("plain old data types"), that are - simplified spoken - C++-types that already existed in C,
are simple too. 
...text:Arrays of simple types can be copied very fast by memory manipulation routines, 
but the default implementation of functions like @Function.arrayCopyForward@ and @Function.arrayCopy@
are not optimized for simple types this way.
But for classes derived from $SimpleType$, optimized variants of array manipulation functions are applied. 
...text:Note that simple types need not to be derived or specialized from $SimpleType$, but
it could be convenient to do so.
..implements:Concept.Simple Type
..include:seqan/basic.h
*/
// TODO(holtgrew): This should actually be a class.
template <typename TValue, typename TSpec>
struct SimpleType
{
//____________________________________________________________________________

	TValue value;

//____________________________________________________________________________

	SimpleType() 
	{
SEQAN_CHECKPOINT
	}

//____________________________________________________________________________

	SimpleType(SimpleType const & other)
	{
SEQAN_CHECKPOINT
		assign(*this, other);
	}

	template <typename T> 
	SimpleType(T const & other) 
	{
SEQAN_CHECKPOINT
		assign(*this, other);
	}


//____________________________________________________________________________

	SimpleType & operator=(SimpleType const & other) 
	{ 
SEQAN_CHECKPOINT
		assign(*this, other);
		return *this;
	}

	template <typename T>
	SimpleType & operator=(T const & other) 
	{ 
SEQAN_CHECKPOINT
		assign(*this, other);
		return *this;
	}
//____________________________________________________________________________

	~SimpleType()
	{
SEQAN_CHECKPOINT
	}
//____________________________________________________________________________


    // Class.SimpleType specifies type conversion operators for all built-in
    // integer types since there is no way to extend the build-in types with
    // copy and assignment constructors in C++.
    //
    // This cannot be a template since it would conflict to the template
    // constructor.

    // TODO(holtgrew): These are candidates for breaking the style convention and simply write each function in one line.

	operator __int64() const
	{
SEQAN_CHECKPOINT
		__int64 c;
		assign(c, *this);
		return c;
	}

	operator __uint64() const
	{
SEQAN_CHECKPOINT
		__uint64 c;
		assign(c, *this);
		return c;
	}


	operator int() const
	{
SEQAN_CHECKPOINT
		int c;
		assign(c, *this);
		return c;
	}
	operator unsigned int() const
	{
SEQAN_CHECKPOINT
		unsigned int c;
		assign(c, *this);
		return c;
	}
	operator short() const
	{
SEQAN_CHECKPOINT
		short c;
		assign(c, *this);
		return c;
	}
	operator unsigned short() const
	{
SEQAN_CHECKPOINT
		unsigned short c;
		assign(c, *this);
		return c;
	}
	operator char() const
	{
SEQAN_CHECKPOINT
		char c;
		assign(c, *this);
		return c;
	}
	operator signed char() const
	{
SEQAN_CHECKPOINT
		signed char c;
		assign(c, *this);
		return c;
	}
	operator unsigned char() const
	{
SEQAN_CHECKPOINT
		unsigned char c;
		assign(c, *this);
		return c;
	}

//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSimple.param.T.type:Class.SimpleType

template <typename TValue, typename TSpec>
struct IsSimple<SimpleType<TValue, TSpec> >
{
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.SimpleType
template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> const >
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.SimpleType
template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> >
{
	typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> const >
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec>, Standard>
{
	typedef SimpleType<TValue, TSpec> * Type;
//	typedef Iter<SimpleType<TValue, TSpec>, SimpleIterator> * Type;
};

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec> const, Standard>
{
	typedef SimpleType<TValue, TSpec> const * Type;
//	typedef Iter<SimpleType<TValue, TSpec> const, SimpleIterator> * Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename TSourceValue, typename TSourceSpec>
inline typename RemoveConst_<TTarget>::Type
convertImpl(Convert<TTarget, T> const,
			SimpleType<TSourceValue, TSourceSpec> const & source_)
{
SEQAN_CHECKPOINT
	typename RemoveConst_<TTarget>::Type target_;
	assign(target_, source_);
	return target_;
}



//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator << (TStream & stream, 
			 SimpleType<TValue, TSpec> const & data)
{
SEQAN_CHECKPOINT
	stream << convert<char>(data);
	return stream;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator >> (TStream & stream, 
			 SimpleType<TValue, TSpec> & data)
{
SEQAN_CHECKPOINT
	char c;
	stream >> c;
	assign(data, c);
	return stream;
}

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:Class.SimpleType
///.Function.assign.param.source.type:Class.SimpleType


template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   SimpleType<TSourceValue, TSourceSpec> & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}
template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   SimpleType<TSourceValue, TSourceSpec> const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}

//____________________________________________________________________________

template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   TSource & source)
{
SEQAN_CHECKPOINT
	target.value = source;
}
template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target.value = source;
}

//____________________________________________________________________________
// Assign Proxy to SimpleType 
//??? Diese Funktionen wurden noetig wegen eines seltsamen VC++-Verhaltens

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   Proxy<TSourceSpec> & source)
{
SEQAN_CHECKPOINT
	target.value = getValue(source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   Proxy<TSourceSpec> const & source)
{
SEQAN_CHECKPOINT
	target.value = getValue(source);
}

//____________________________________________________________________________
//INTEGRAL TYPES
//note: it is not possible to write a single function here since "assign"
//must be specialized for the first argument at the first place

//__int64
template <typename TValue, typename TSpec>
inline void 
assign(__int64 & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(__int64 & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//__uint64
template <typename TValue, typename TSpec>
inline void 
assign(__uint64 & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(__uint64 & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//int
template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned int
template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//short
template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned short
template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//char
template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//signed char
template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//unsigned char
template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
	   SimpleType<TValue, TSpec> & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}
template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
	   SimpleType<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
	c_target = source.value;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// CompareType
//////////////////////////////////////////////////////////////////////////////

/**.Metafunction.CompareType:
..summary:Type to convert other types for comparisons.
..signature:CompareType<TLeft, TRight>::Type
..param.TLeft:Type of the left operand of a comparison.
..param.TRight:Type of the right operand of a comparison.
..return.Type:The Type in which the arguments are converted in order to compare them.
..remarks:Comparisons are for example operators like $==$ or $<$.
..remarks.text:Note that there is no rule that guarantees that $CompareType<T1, T2>::Type$
is the same as $CompareType<T2, T1>::Type$. It is also possible, that only one of these
two types is defined.
..remarks.text:This metafunction is used for the implementation of
comparisons that involve @Class.SimpleType@.
*/
//???TODO: muss geprueft werden, ob diese Metafunktion noch ausgeweitet oder aber versteckt wird.

template <typename TLeft, typename TRight>
struct CompareType;

template <typename T>
struct CompareType<T, T>
{
	typedef T Type;
};

//____________________________________________________________________________

template <typename TValue, typename TSpec, typename TRight>
struct CompareType<SimpleType<TValue, TSpec>, TRight>
{
	typedef TRight Type;
};

//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator == (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator == (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator == (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator == (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) == convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator == (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator == (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator !=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator != (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator != (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator != (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator != (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) != convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator != (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator != (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator <

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator < (SimpleType<TValue, TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator < (TLeft const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator < (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator < (SimpleType<TValue, TSpec> const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) < convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator < (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator < (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator <=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator <= (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator <= (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator <= (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator <= (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) <= convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator <= (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator <= (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}



//____________________________________________________________________________
// operator >

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator > (SimpleType<TValue, TSpec> const & left_, 
			TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator > (TLeft const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator > (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator > (SimpleType<TValue, TSpec> const & left_, 
			SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) > convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator > (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator > (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}


//____________________________________________________________________________
// operator >=

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator >= (SimpleType<TValue, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator >= (TLeft const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator >= (SimpleType<TLeftValue, TLeftSpec> const & left_, 
			 SimpleType<TRightValue, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
	typedef SimpleType<TRightValue, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator >= (SimpleType<TValue, TSpec> const & left_, 
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return convert<TValue>(left_) >= convert<TValue>(right_);
}


template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator >= (Proxy<TSpec> const & left_, 
			 SimpleType<TValue, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator >= (SimpleType<TValue, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}


//////////////////////////////////////////////////////////////////////////////

template<typename T_, typename TSpec> 
inline
bool lexLess(SimpleType<T_, TSpec> const &_Left, SimpleType<T_, TSpec> const &Right_)
{	// return lexicographical _Left < Right_
	typedef typename MakeUnsigned_<T_>::Type TUnsigned;
    return (TUnsigned)(_Left.value) < (TUnsigned)(Right_.value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator ++ (SimpleType<TValue, TSpec> & me)
{
	++me.value;
	return me;
}
template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator ++ (SimpleType<TValue, TSpec> & me,
			 int)
{
	SimpleType<TValue, TSpec> dummy = me;
	++me.value;
	return dummy;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator -- (SimpleType<TValue, TSpec> & me)
{
	--me.value;
	return me;
}
template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator -- (SimpleType<TValue, TSpec> & me,
			 int)
{
	SimpleType<TValue, TSpec> dummy = me;
	--me.value;
	return dummy;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Dna:
..cat:Alphabets
..summary:Alphabet for DNA.
..general:Class.SimpleType
..signature:Dna
..remarks:
...text:The @Metafunction.ValueSize@ of $Dna$ is 4. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$.
...text:Objects of type $Dna$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'T'}$ is converted to $'A'$.
...text:$Dna$ is typedef for $SimpleType<char,Dna_>$, while $Dna_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..see:Spec.Dna5
..include:seqan/basic.h
*/
struct Dna_ {};
typedef SimpleType<unsigned char,Dna_> Dna;

template <> struct ValueSize< Dna > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< Dna > { enum { VALUE = 2 }; };

//____________________________________________________________________________

/**
.Spec.Dna5:
..cat:Alphabets
..summary:Alphabet for DNA including 'N' character.
..general:Class.SimpleType
..signature:Dna5
..remarks:
...text:The @Metafunction.ValueSize@ of $Dna5$ is 5. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$. 
The 'N' character ("unkown nucleotide") is encoded by 4.
...text:Objects of type $Dna5$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'T'}$ is converted to $'N'$.
...text:$Dna5$ is typedef for $SimpleType<char,Dna5_>$, while $Dna5_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/
struct Dna5_ {};
typedef SimpleType<unsigned char, Dna5_> Dna5;

template <> struct ValueSize< Dna5 > { enum { VALUE = 5 }; };
template <> struct BitsPerValue< Dna5 > { enum { VALUE = 3 }; };

//____________________________________________________________________________

/**
.Spec.Rna:
..cat:Alphabets
..summary:Alphabet for RNA.
..general:Class.SimpleType
..signature:Rna
..remarks:
...text:The @Metafunction.ValueSize@ of $Rna$ is 4. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'U' = 3$.
...text:Objects of type $Rna$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'U'}$ is converted to $'A'$.
...text:$Rna$ is typedef for $SimpleType<char,Rna_>$, while $Rna_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..see:Spec.Rna5
..include:seqan/basic.h
*/
struct Rna_ {};
typedef SimpleType<unsigned char,Rna_> Rna;

template <> struct ValueSize< Rna > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< Rna > { enum { VALUE = 2 }; };

//____________________________________________________________________________

/**
.Spec.Rna5:
..cat:Alphabets
..summary:Alphabet for RNA including 'N' character.
..general:Class.SimpleType
..signature:Rna5
..remarks:
...text:The @Metafunction.ValueSize@ of $Rna5$ is 5. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'U' = 3$. 
The 'N' character ("unkown nucleotide") is encoded by 4.
...text:Objects of type $Rna5$ can be converted to various other types and vice versa. 
An object that has a value not in ${'A', 'C', 'G', 'U'}$ is converted to $'N'$.
...text:$Rna5$ is typedef for $SimpleType<char,Rna5_>$, while $Rna5_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/
struct Rna5_ {};
typedef SimpleType<unsigned char, Rna5_> Rna5;

template <> struct ValueSize< Rna5 > { enum { VALUE = 5 }; };
template <> struct BitsPerValue< Rna5 > { enum { VALUE = 3 }; };

//____________________________________________________________________________

/**
.Spec.Iupac:
..cat:Alphabets
..summary:Iupac code for DNA.
..general:Class.SimpleType
..signature:Iupac
..remarks:
...text:The @Metafunction.ValueSize@ of $Iupac$ is 16. 
The nucleotides are enumerated from 0 to 15 in this order: 
'U'=0, 'T', 'A', 'W', 'C', 'Y', 'M', 'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'=15. 
...text:Objects of type $Iupac$ can be converted to various other types and vice versa. 
Unkown values are converted to $'N'$.
...text:$Iupac$ is typedef for $SimpleType<char,Iupac_>$, while $Iupac_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/
struct Iupac_ {};
typedef SimpleType<unsigned char, Iupac_> Iupac;

template <> struct ValueSize< Iupac > { enum { VALUE = 16 }; };
template <> struct BitsPerValue< Iupac > { enum { VALUE = 4 }; };


//____________________________________________________________________________

/**
.Spec.AminoAcid:
..cat:Alphabets
..summary:Iupac code for amino acids.
..general:Class.SimpleType
..signature:AminoAcid
..remarks:
...text:The @Metafunction.ValueSize@ of $AminoAcid$ is 24. 
...text:The amino acids are enumerated from 0 to 15 in this order: 
...text:'A'=0, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=19.
...text:The remaining 4 symbols are:
...text: 'B'=20 (Aspartic Acid, Asparagine), 'Z'=21 (Glutamic Acid, Glutamine), 'X'=22 (unknown), '*'=23 (terminator)
...text:Objects of type $AminoAcid$ can be converted to $char$ and vice versa. 
Unkown values are converted to $'X'$.
...text:$AminoAcid$ is typedef for $SimpleType<char,AminoAcid_>$, while $AminoAcid_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/
struct AminoAcid_ {};
typedef SimpleType<unsigned char, AminoAcid_> AminoAcid;

template <> struct ValueSize< AminoAcid > { enum { VALUE = 24 }; };
template <> struct BitsPerValue< AminoAcid > { enum { VALUE = 5 }; };

//____________________________________________________________________________

/**
.Spec.Finite:
..cat:Alphabets
..summary:A finite alphabet of a fixed size.
..general:Class.SimpleType
..signature:SimpleType<TValue, Finite<SIZE> >
..param.TValue:The type that is use to store the values.
...default:$char$
..param.SIZE:The @Metafunction.ValueSize@ of the alphabet.
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/
template <unsigned SIZE>
struct Finite;

template <typename TValue, unsigned SIZE> 
struct ValueSize< SimpleType<TValue, Finite<SIZE> > > { enum { VALUE = SIZE }; };

template <typename TValue, unsigned SIZE> 
struct BitsPerValue< SimpleType<TValue, Finite<SIZE> > >: Log2<SIZE> {};

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign(Ascii & c_target, 
				   Dna const & source)
{
SEQAN_CHECKPOINT
	c_target = TranslateTableDna5ToAscii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, 
				   Dna5 const & source)
{
SEQAN_CHECKPOINT
	c_target = TranslateTableDna5ToAscii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii& target,
				   Rna const & source)
{
	SEQAN_CHECKPOINT
	target = TranslateTableRna5ToAscii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii& target,
				   Rna5 const & source)
{
	SEQAN_CHECKPOINT
	target = TranslateTableRna5ToAscii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, Iupac const & source)
{
SEQAN_CHECKPOINT
	c_target = TranslateTableIupacToAscii_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Ascii & c_target, AminoAcid const & source)
{
SEQAN_CHECKPOINT
	c_target = TranslateTableAAToAscii_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DNA (4 letters)

template <>
struct CompareType<Dna, __uint8> { typedef Dna Type; };
inline void assign(Dna & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableByteToDna_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Ascii> { typedef Dna Type; };
inline void assign(Dna & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToDna_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Unicode> { typedef Dna Type; };
inline void assign(Dna & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToDna_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Dna5> { typedef Dna Type; };
inline void assign(Dna & target, Dna5 const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value & 0x03;
}
//____________________________________________________________________________

template <>
struct CompareType<Dna, Iupac> { typedef Dna Type; };
inline void assign(Dna & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableIupacToDna_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//DNA (5 letters)

template <>
struct CompareType<Dna5, __uint8> { typedef Dna5 Type; };
inline void assign(Dna5 & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableByteToDna5_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Ascii> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToDna5_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Unicode> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToDna5_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Iupac> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableIupacToDna5_<>::VALUE[source.value];
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5, Dna> { typedef Dna Type; };
inline void assign(Dna5 & target, Dna const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value;
}

//////////////////////////////////////////////////////////////////////////////
//RNA (4 letters)

template <>
struct CompareType<Rna, __uint8> { typedef Rna Type; };
inline void assign(Rna & target, __uint8 c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableByteToRna_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna, Ascii> { typedef Rna Type; };
inline void assign(Rna & target, Ascii c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToRna_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna, Unicode> { typedef Rna Type; };
inline void assign(Rna & target, Unicode c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToRna_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna, Rna5> { typedef Rna Type; };
inline void assign(Rna & target, Rna5 const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value & 0x03;
}

//////////////////////////////////////////////////////////////////////////////
//RNA (5 letters)

template <>
struct CompareType<Rna5, __uint8> { typedef Rna5 Type; };
inline void assign(Rna5 & target, __uint8 c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableByteToRna5_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna5, Ascii> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Ascii c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToRna5_<>::VALUE[(unsigned char)c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna5, Unicode> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Unicode c_source)
{
	SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToRna5_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Rna5, Rna> { typedef Dna Type; };
inline void assign(Rna5 & target, Rna const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value;
}

//////////////////////////////////////////////////////////////////////////////
//IUPAC (4 bits)

template <>
struct CompareType<Iupac, __uint8> { typedef Iupac Type; };
inline void assign(Iupac & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableByteToIupac_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Iupac, Ascii> { typedef Iupac Type; };
inline void assign(Iupac & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToIupac_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<Iupac, Unicode> { typedef Iupac Type; };
inline void assign(Iupac & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToIupac_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

inline void assign(Iupac & target, Dna const & source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableDna5ToIupac_<>::VALUE[source.value];
}
//____________________________________________________________________________

inline void assign(Iupac & target, Dna5 const & source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableDna5ToIupac_<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//Amino Acid (5 bits)

template <>
struct CompareType<AminoAcid, __uint8> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableByteToAA_<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid, Ascii> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToAA_<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid, Unicode> { typedef AminoAcid Type; };
inline void assign(AminoAcid & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = TranslateTableAsciiToAA_<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
struct BaseAlphabet
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

// //DnaQ and Dna5Q

/**
.Spec.DnaQ:
..cat:Alphabets
..summary:Alphabet for DNA plus phred quality.
..general:Class.SimpleType
..signature:DnaQ
..remarks:
...text:The @Metafunction.ValueSize@ of $DnaQ$ is 4. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$.
...text:Objects of type $DnaQ$ can be converted to various other types and vice versa. 
...text:$DnaQ$ is typedef for $SimpleType<char,DnaQ_>$, while $DnaQ_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..see:Spec.Dna5Q
*/
struct DnaQ_ {};
typedef SimpleType <unsigned char, DnaQ_> DnaQ;

template <> struct ValueSize< DnaQ > { enum { VALUE = 4 }; };				// considering nucleotides
template <> struct InternalValueSize_< DnaQ > { enum { VALUE = 252 }; };	// considering nucleotides x Quality 0..62
template <> struct BitsPerValue< DnaQ > { enum { VALUE = 8 }; };

template <>
struct BaseAlphabet<DnaQ>
{
	typedef Dna Type;
};

//____________________________________________________________________________
 
/**
.Spec.Dna5Q:
..cat:Alphabets
..summary:Alphabet for DNA plus phred quality including 'N' character.
..general:Class.SimpleType
..signature:Dna5Q
..remarks:
...text:The @Metafunction.ValueSize@ of $Dna5Q$ is 5. 
The nucleotides are enumerated this way: $'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3$. 
The 'N' character ("unkown nucleotide") is encoded by 4.
...text:Objects of type $Dna5$ can be converted to various other types and vice versa. 
...text:$Dna5Q$ is typedef for $SimpleType<char,Dna5Q_>$, while $Dna5Q_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct Dna5Q_ {};
typedef SimpleType <unsigned char, Dna5Q_> Dna5Q;

static const unsigned char Dna5QValueN_ = 252;								// value representing N

template <> struct ValueSize< Dna5Q > { enum { VALUE = 5 }; };				// considering nucleotides + N
template <> struct InternalValueSize_< Dna5Q > { enum { VALUE = 253 }; };	// considering (nucleotides x Quality 0..62) + N
template <> struct BitsPerValue< Dna5Q > { enum { VALUE = 8 }; };

template <>
struct BaseAlphabet<Dna5Q>
{
	typedef Dna5 Type;
};

template <typename TValue>
struct QualityValueSize {
	enum { VALUE = ValueSize<TValue>::VALUE };
};
template <typename TValue>
struct QualityValueSize<TValue const>:
	public QualityValueSize<TValue> {};

template <> struct QualityValueSize< DnaQ >  { enum { VALUE = 63 }; };		// 64 - 1 (N)
template <> struct QualityValueSize< Dna5Q > { enum { VALUE = 63 }; };


// template <typename TValue, typename TValue2>
// struct CompareType<SimpleType<TValue,DnaQ_>, SimpleType<TValue2,Dna_> >
// {
// 	typedef SimpleType<TValue2,Dna_> Type;
// };
// 
// template <typename TValue, typename TValue2>
// struct CompareType<SimpleType<TValue,Dna_>, SimpleType<TValue2,DnaQ_> >
// {
// 	typedef SimpleType<TValue,Dna_> Type;
// };
// 
// 
// 
// template <typename TValue, typename TValue2>
// struct CompareType<SimpleType<TValue,Dna5Q_>, SimpleType<TValue2,Dna5_> >
// {
// 	typedef SimpleType<TValue2,Dna5_> Type;
// };
// 
// template <typename TValue, typename TValue2>
// struct CompareType<SimpleType<TValue,Dna5_>, SimpleType<TValue2,Dna5Q_> >
// {
// 	typedef SimpleType<TValue,Dna5_> Type;
// };


template <>
struct CompareType<Dna5Q, Dna5Q> { typedef Dna5 Type; };
template <>
struct CompareType<DnaQ, DnaQ> { typedef Dna Type; };

//////////////////////////////////////////////////////////////////////////////
//DNA (4 letters) with Qualities (0..60)


template <>
struct CompareType<DnaQ, Dna5Q> { typedef Dna Type; };
inline void assign(DnaQ & target, Dna5Q const & source)
{
    SEQAN_CHECKPOINT;

    // We perform the converstion from DNA5 to DNA5 with qualities by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
	// target.value = (source.value == Dna5QValueN_)? 0: source.value;

    static const unsigned table[] = {
          0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
         16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
         32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
         48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
         64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
         80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
         96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
        112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 0,   0,   0,   0
    };
    target.value = table[source.value];
}
//____________________________________________________________________________


template <>
struct CompareType<DnaQ, Dna> { typedef Dna Type; };
inline void assign(DnaQ & target, Dna const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value | (60 << 2);
}

//////////////////////////////////////////////////////////////////////////////
//DNA (5 letters)

template <>
struct CompareType<Dna5Q, DnaQ> { typedef Dna Type; };
inline void assign(Dna5Q & target, DnaQ const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}

//____________________________________________________________________________

template <>
struct CompareType<Dna5, Dna5Q> { typedef Dna5 Type; };
inline void assign(Dna5 & target, Dna5Q const & source)
{
    SEQAN_CHECKPOINT;

    // We perform the conversion from DNA5 to DNA5 with qualities by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
	// target.value = (source.value == Dna5QValueN_)? 4: source.value & 3;

    static const unsigned table[] = {
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 4 // <-- note the 4
    };
    target.value = table[source.value];
}


//____________________________________________________________________________

template <>
struct CompareType<Dna5Q, Dna5> { typedef Dna5 Type; };
inline void assign(Dna5Q & target, Dna5 const & source)
{
    SEQAN_CHECKPOINT;

    // We perform the conversion from DNA5 with qualities to DNA5 by a simple
    // table lookup.  The lookup below is equivalent to the following line:
    //
	// target.value = (source.value == 4)? Dna5QValueN_ : source.value | (60 << 2);

    static const unsigned table[] = {
        (60 << 2) + 0, (60 << 2) + 1, (60 << 2) + 2, (60 << 2) + 3, Dna5QValueN_
    };
    target.value = table[source.value];
}


//____________________________________________________________________________

template <>
struct CompareType<Dna, DnaQ> { typedef Dna Type; };
inline void assign(Dna & target, DnaQ const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value & 3;
}

//____________________________________________________________________________


template <>
struct CompareType<DnaQ, Iupac> { typedef Dna Type; };
inline void assign(DnaQ & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna) source);
}

//____________________________________________________________________________


template <>
struct CompareType<DnaQ, Dna5> { typedef Dna Type; };
inline void assign(DnaQ & target, Dna5 const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna) source);
}

template <>
struct CompareType<DnaQ, __uint8> { typedef Dna Type; };
inline void assign(DnaQ & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna) c_source);
}
//____________________________________________________________________________

template <>
struct CompareType<DnaQ, Ascii> { typedef Dna Type; };
inline void assign(DnaQ & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna) c_source);
}
//____________________________________________________________________________

template <>
struct CompareType<DnaQ, Unicode> { typedef Dna Type; };
inline void assign(DnaQ & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna) c_source);
}
//____________________________________________________________________________

inline void 
assign(DnaQ & target, DnaQ const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}

template <typename TSource>
inline void 
assign(DnaQ & target, TSource const & source)
{
SEQAN_CHECKPOINT
	target.value = (Dna)source;
}
//____________________________________________________________________________


template <>
struct CompareType<Dna5Q, Dna> { typedef Dna Type; };
inline void assign(Dna5Q & target, Dna const & source)
{
SEQAN_CHECKPOINT
	assign(target, (DnaQ) source);
}

//____________________________________________________________________________


template <>
struct CompareType<Dna, Dna5Q> { typedef Dna Type; };
inline void assign(Dna & target, Dna5Q const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)source);
}

//____________________________________________________________________________
template <>
struct CompareType<Dna5, DnaQ> { typedef Dna5 Type; };
inline void assign(Dna5 & target, DnaQ const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5Q)source);
}
//____________________________________________________________________________



template <>
struct CompareType<Dna5Q, __uint8> { typedef Dna5 Type; };
inline void assign(Dna5Q & target, __uint8 c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)c_source);
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5Q, Ascii> { typedef Dna5 Type; };
inline void assign(Dna5Q & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)c_source);
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5Q, Unicode> { typedef Dna5 Type; };
inline void assign(Dna5Q & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)c_source);
}
//____________________________________________________________________________

template <>
struct CompareType<Dna5Q, Iupac> { typedef Dna5 Type; };
inline void assign(Dna5Q & target, Iupac const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)source);
}

inline void 
assign(Dna5Q & target, Dna5Q const & source)
{
SEQAN_CHECKPOINT
	target.value = source.value;
}
template <typename TSource>
inline void 
assign(Dna5Q & target, TSource const & source)
{
SEQAN_CHECKPOINT
	assign(target, (Dna5)source);
}


//____________________________________________________________________________



//__int64

inline void 
assign(__int64 & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(__int64 & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//__uint64

inline void 
assign(__uint64 & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(__uint64 & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//int

inline void 
assign(int & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(int & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//unsigned int

inline void 
assign(unsigned int & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(unsigned int & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}


//short

inline void 
assign(short & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(short & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//unsigned short

inline void 
assign(unsigned short & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(unsigned short & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//char

inline void 
assign(char & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(char & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//signed char

inline void 
assign(signed char & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(signed char & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

//unsigned char

inline void 
assign(unsigned char & c_target, 
	   DnaQ & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}

inline void 
assign(unsigned char & c_target, 
	   DnaQ const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna(source);
}


//__int64

inline void 
assign(__int64 & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(__int64 & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//__uint64

inline void 
assign(__uint64 & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(__uint64 & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//int

inline void 
assign(int & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(int & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//unsigned int

inline void 
assign(unsigned int & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(unsigned int & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}


//short

inline void 
assign(short & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(short & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//unsigned short

inline void 
assign(unsigned short & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(unsigned short & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//char

inline void 
assign(char & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(char & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//signed char

inline void 
assign(signed char & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(signed char & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

//unsigned char

inline void 
assign(unsigned char & c_target, 
	   Dna5Q & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline void 
assign(unsigned char & c_target, 
	   Dna5Q const & source)
{
SEQAN_CHECKPOINT
	c_target = Dna5(source);
}

inline int getQualityValue(DnaQ const &c) 
{
	return c.value >> 2;
}

inline int getQualityValue(Dna5Q const &c) 
{
    // We use a lookup table to extract the qualities from DNA5Q.  The lookup
    // table based code is equivalent to the following line:
	// return (c.value == Dna5QValueN_)? 0: c.value >> 2;

    static const unsigned table[] = {
         0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,
         4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,
         8,  8,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12,
        12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16,
        17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21,
        21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 25, 25,
        25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29,
        29, 30, 30, 30, 30, 31, 31, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33,
        34, 34, 34, 34, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 37, 37, 38,
        38, 38, 38, 39, 39, 39, 39, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42,
        42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46,
        46, 47, 47, 47, 47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, 50, 50,
        51, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, 54, 54, 55,
        55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59,
        59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, 62, 62,
		0,  0,  0,  0};
    return table[c.value];
}

inline 
void convertQuality(Ascii & c, int q) 
{
	c = '!' + Ascii(q);
}


// TODO(holtgrew): What about different quality types? Guess scaling? Look at how other packages do this.
/**
.Function.assignQualityValue
..cat:Alphabets
..signature:assignQualityValue(c, q)
..summary:Assign quality to a character from an alphabet with integrated quality.
..param.c:Target character to assign quality to.
...type:Spec.DnaQ
..param.q:Quality to assign to the quality.
...type:nolink:int
...type:nolink:char
..remarks:If $q$ is a $char$ then $'!'$ is subtracted from $q$. This is useful for ASCII encoded PHRED scores.
 */
//set quality value
inline 
void assignQualityValue(DnaQ &c, int q)
{
	if (q < 0) q = 0;
    if (q >= QualityValueSize<DnaQ>::VALUE)
        q = QualityValueSize<DnaQ>::VALUE - 1;
	c.value = (c.value & 3) | (q << 2);
}

///.Function.assignQualityValue.param.c.type:Spec.Dna5Q
inline 
void assignQualityValue(Dna5Q &c, int q) 
{
	if (q < 0) q = 0;
    if (q >= QualityValueSize<Dna5Q>::VALUE)
        q = QualityValueSize<Dna5Q>::VALUE - 1;
	if (c.value != Dna5QValueN_)
		c.value = (c.value & 3) | (q << 2);
}

inline 
void assignQualityValue(DnaQ &c, Ascii q)
{
    int q1 = static_cast<int>(q - '!');
	if (q1 < 0) q1 = 0;
    if (q1 >= QualityValueSize<DnaQ>::VALUE)
        q1 = QualityValueSize<DnaQ>::VALUE - 1;
	assignQualityValue(c, q1);
}

inline 
void assignQualityValue(Dna5Q &c, Ascii q) 
{
    int q1 = static_cast<int>(q - '!');
	if (q1 < 0) q1 = 0;
	if (q1 >= QualityValueSize<Dna5Q>::VALUE)
        q1 = QualityValueSize<Dna5Q>::VALUE - 1;
	assignQualityValue(c, q1);
}

/**
.Function.assignQualities
..cat:Alphabets
..summary:Assign quality value between strings.
..signature:assignQualities(target, source)
..param.target:Target string
...type:nolink:@Class.String@ of any alphabet with qualities, e.g. @Spec.DnaQ@, @Spec.Dna5Q@
..param.source:Source string.
...type:nolink:@Class.String@ of $int$ or $char$.
..remarks:This funciton calls @Function.assignQualityValue@ for all entries of $target$ and $source$, look at the documentation of @Function.assignQualityValue@ on how the values of $source$ are interpreted.
..see:Function.assignQualityValue
 */
template <typename TDest, typename TSource>
void assignQualities(TDest &dst, TSource const &src)
{
	typedef typename Iterator<TDest>::Type TDestIter;
	typedef typename Iterator<TSource>::Type TSourceIter;

	TDestIter itDst = begin(dst, Standard());
	TDestIter itDstEnd = end(dst, Standard());
	TSourceIter itSrcEnd = end(src, Standard());
	
	for (TSourceIter itSrc = begin(src, Standard()); itDst != itDstEnd && itSrc != itSrcEnd; ++itDst, ++itSrc)
		assignQualityValue(*itDst, *itSrc);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
