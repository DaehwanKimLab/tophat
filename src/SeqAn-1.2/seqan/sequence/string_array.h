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
  $Id: string_array.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_STRING_ARRAY_H
#define SEQAN_HEADER_SEQUENCE_STRING_ARRAY_H


namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.Array String:
..cat:Strings
..general:Class.String
..summary:Fast but non-expandable string.
..signature:String<TValue, Array<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the capacity of the string.
...remarks:Note that the capacity of a stack string cannot be changed later.
*/
//////////////////////////////////////////////////////////////////////////////

template <unsigned int ISize>
struct Array;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int ISize>
class String<TValue, Array<ISize> >
{
public:
	mutable TValue data_begin[ISize];
	TValue * data_end;

//____________________________________________________________________________

public:
	String()
	{
SEQAN_CHECKPOINT
		data_end = data_begin;
	}

	template <typename TSource>
	String(TSource & source)
	{
SEQAN_CHECKPOINT
		data_end = data_begin;
		assign(*this, source);
	}
	template <typename TSource>
	String(TSource const & source)
	{
SEQAN_CHECKPOINT
		data_end = data_begin;
		assign(*this, source);
	}
	String(String const & source)
	{
SEQAN_CHECKPOINT
		data_end = data_begin;
		assign(*this, source);
	}

	template <typename TSource>
	String & operator =(TSource const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}
	String & operator =(String const & source)
	{
SEQAN_CHECKPOINT
		assign(*this, source);
		return *this;
	}

	~String()
	{
	}

//____________________________________________________________________________

	template <typename TPos>
	inline typename Reference<String>::Type
	operator [] (TPos pos)
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

	template <typename TPos>
	inline typename Reference<String const>::Type 
	operator [] (TPos pos) const
	{
SEQAN_CHECKPOINT
		return value(*this, pos);
	}

//____________________________________________________________________________


};
//////////////////////////////////////////////////////////////////////////////


template <typename TValue, unsigned int ISize>
inline typename Iterator<String<TValue, Array<ISize> >, Standard>::Type
begin(String<TValue, Array<ISize> > & me,
	Standard)
{
SEQAN_CHECKPOINT
	return me.data_begin;
}
template <typename TValue, unsigned int ISize>
inline typename Iterator<String<TValue, Array<ISize> > const, Standard>::Type
begin(String<TValue, Array<ISize> > const & me,
	Standard)
{
SEQAN_CHECKPOINT
	return me.data_begin;
}

//____________________________________________________________________________

template <typename TValue, unsigned int ISize>
inline typename Iterator<String<TValue, Array<ISize> >, Standard>::Type
end(String<TValue, Array<ISize> > & me,
	Standard)
{
SEQAN_CHECKPOINT
	return me.data_end;
}
template <typename TValue, unsigned int ISize>
inline typename Iterator<String<TValue, Array<ISize> > const, Standard>::Type
end(String<TValue, Array<ISize> > const & me,
	Standard)
{
SEQAN_CHECKPOINT
	return me.data_end;
}

//____________________________________________________________________________

template <typename TValue, unsigned int ISize>
inline size_t
capacity(String<TValue, Array<ISize> > &)
{
SEQAN_CHECKPOINT
	return ISize;
}

template <typename TValue, unsigned int ISize>
inline size_t
capacity(String<TValue, Array<ISize> > const &)
{
SEQAN_CHECKPOINT
	return ISize;
}
//____________________________________________________________________________

template <typename TValue, unsigned int ISize, typename TExpand>
inline size_t
reserve(
	String<TValue, Array<ISize> > & me, 
	size_t,
	Tag<TExpand> const)
{
SEQAN_CHECKPOINT
	return capacity(me);
}


//____________________________________________________________________________

/**
.Internal._setLength.param.object.type:Spec.Array String
*/
template <typename TValue, unsigned int ISize>
inline void 
_setLength(
	String<TValue, Array<ISize> > & me, 
	size_t new_length)
{
SEQAN_CHECKPOINT
	me.data_end = me.data_begin + new_length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int ISize>
struct DefaultOverflowImplicit<String<TValue, Array<ISize> > >
{
	typedef Limit Type;
};

template <typename TValue, unsigned int ISize>
struct DefaultOverflowImplicit<String<TValue, Array<ISize> > const >
{
	typedef Limit Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int ISize>
struct DefaultOverflowExplicit<String<TValue, Array<ISize> > >
{
	typedef Limit Type;
};

template <typename TValue, unsigned int ISize>
struct DefaultOverflowExplicit<String<TValue, Array<ISize> > const >
{
	typedef Limit Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int ISize>
struct IsContiguous< String<TValue, Array<ISize> > >
{
    typedef True Type;
	enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.LENGTH.param.T.type:Spec.Array String
template <typename TValue, unsigned int ISize>
struct LENGTH< String<TValue, Array<ISize> > >
{
	enum { VALUE = ISize };
};

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
