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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_ALPHABET_EXPANSION_H
#define SEQAN_HEADER_MODIFIER_ALPHABET_EXPANSION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Alphabet Expansion:
..summary:Modifier that adds a character to an alphabet.
..cat:Modifier
..signature:ModifiedAlphabet<TAlphabet, ModExpand<CHAR [,TSpec]> >
..param.TAlphabet:Original value type.
..param.CHAR:$char$ character that specifies, what value should added to the alphabet.
...remarks:$CHAR$ should not be a $char$ that already stands for a value in $TAlphabet$.
	For example, do not use $'A'$ or $'a'$ as $CHAR$ when expanding @Spec.Dna@.
...remarks:Some values of $CHAR$ have special meaning:
....table:$'-'$|A gap character. The value in the expanded alphabet that corresponds to $'-'$ will be returned by the @Function.gapValue@.
....table:$'\$'$|An end of string character.
..param.TSpec:Optional specialization tag.
...default:$Default$
...remarks:This modifier is intended to expand @Class.SimpleType@ classes.
*/

template <char CHAR, typename TSpec = Default>
struct ModExpand;


template <typename THost, char CHAR, typename TSpec>
class ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
{
public:
	typedef typename IntegralForValue<ModifiedAlphabet>::Type TData;
	TData data;

	ModifiedAlphabet() 
	{
	}
	ModifiedAlphabet(ModifiedAlphabet const & other)
		: data(other.data)
	{
	}
	template <typename TOther>
	ModifiedAlphabet(TOther const & other_data)
		: data(ordValue(convert<ModifiedAlphabet>(other_data)))
	{
	}
	~ModifiedAlphabet()
	{
	}
	ModifiedAlphabet const & 
	operator = (ModifiedAlphabet const & other)
	{
		data = other.data;
		return *this;
	}
	template <typename TOther>
	ModifiedAlphabet const & 
	operator = (TOther const & other_data)
	{
		data = ordValue(convert<ModifiedAlphabet>(other_data));
		return *this;
	}

/*	operator TData ()
	{
		return data;
	}

	operator THost()
	{
		return convert<THost>(data);
	}
*/
//____________________________________________________________________________

	//this cannot be a template since a template would be in conflict to
	//the template c'tor


	operator long() const
	{
SEQAN_CHECKPOINT
		return convert<long>(*this);
	}
	operator unsigned long() const
	{
SEQAN_CHECKPOINT
		return convert<unsigned long>(*this);
	}
	operator int() const
	{
SEQAN_CHECKPOINT
		return convert<int>(*this);
	}
	operator unsigned int() const
	{
SEQAN_CHECKPOINT
		return convert<unsigned int>(*this);
	}
	operator short() const
	{
SEQAN_CHECKPOINT
		return convert<short>(*this);
	}
	operator unsigned short() const
	{
SEQAN_CHECKPOINT
		return convert<unsigned short>(*this);
	}
	operator char() const
	{
SEQAN_CHECKPOINT
		return convert<char>(*this);
	}
	operator signed char() const
	{
SEQAN_CHECKPOINT
		return convert<signed char>(*this);
	}
	operator unsigned char() const
	{
SEQAN_CHECKPOINT
		return convert<unsigned char>(*this);
	}

};




//////////////////////////////////////////////////////////////////////////////
// sizes
//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
struct BitsPerValue<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >
{
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TValue;
	enum { VALUE = Log2< ValueSize<TValue>::VALUE >::VALUE };
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
struct ValueSize<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >
{
	enum { VALUE = ValueSize<THost>::VALUE + 1 };
};

//////////////////////////////////////////////////////////////////////////////
// conversions
//////////////////////////////////////////////////////////////////////////////


// some type => ModExpand
template <typename THost, char CHAR, typename TSpec, typename TSource>
inline void
_initializeAlphabetConversionTable(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > * buf,
								   TSource const &)
{
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;

	//assure that the conversion from TSource to THost is possible
//	_AlphabetConversionTable<THost, TSource>::initialize();
	
	//add the new character CHAR to the table
	buf[ordValue(convert<TSource>(CHAR))].data = ValueSize<THost>::VALUE;

	//copy the conversion table for converting TSouce => THost
	//maybe, if there is no CHAR in TSource, the entry for CHAR is overwritten now
	for (int i = ValueSize<TSource>::VALUE; i > 0; )
	{
		--i;
		buf[i].data = ordValue(convert<THost>(convert<TSource>(i)));
	}

}

template <int SIZE_OF_SOURCE>
struct _ConvertImpl_ModExpand
{
	//default implementation for large source types
	template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
	inline 
	static typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
	_convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
		TSource const & source_)
	{
		typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
		if (source_ == ValueSize<THost>::VALUE)
		{// the extra character
			return convert<TTarget>((int) ValueSize<THost>::VALUE);
		}
		return convert<TTarget>(convert<THost>(source_));
	}
};

//for 1 byte source: use translation table
template <>
struct _ConvertImpl_ModExpand<1>
{
	template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
	inline 
	static typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
	_convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
		TSource const & source_)
	{
		typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
		return _AlphabetConversionTable<TTarget, TSource>::table[ordValue(source_)];
	}
};

//generic source: dispatch for size of BytesPerValue
template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const convert_,
			TSource const & source_)
{
	return _ConvertImpl_ModExpand<BytesPerValue<TSource>::VALUE>::_convertImpl(convert_, source_);
}

//for SimpleType sources
template <typename THost, char CHAR, typename TSpec, typename T, typename TSourceValue, typename TSourceSpec>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, SimpleType<TSourceValue, TSourceSpec> >::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
			SimpleType<TSourceValue, TSourceSpec> const & source_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
	typedef SimpleType<TSourceValue, TSourceSpec> TSource;
	return _AlphabetConversionTable<TTarget, TSource>::table[ordValue(source_)];
}

//for Proxy sources
template <typename THost, char CHAR, typename TSpec, typename T, typename TSpec2>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, Proxy<TSpec2> >::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
			Proxy<TSpec2> const & source_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
	return convert<TTarget>(getValue(source_));
}


// ModExpand => some type

template <typename TTarget, typename THost, char CHAR, typename TSpec>
inline void
_initializeAlphabetConversionTable(TTarget * buf,
								   ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const &)
{
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TSource;

	//assure that the conversion from THost to TTarget is possible
	_AlphabetConversionTable<TTarget, THost>::initialize(); 
	
	//copy the conversion table for converting THost => TTarget
	for (int i = ValueSize<THost>::VALUE; i > 0; )
	{
		--i;
		buf[i] = convert<TTarget>(convert<THost>(i));
	}

	//add the new character CHAR to the table
	buf[ValueSize<THost>::VALUE] = convert<TTarget, char>(CHAR);
}

template <typename TTarget, typename T, typename THost, char CHAR, typename TSpec>
inline typename Convert<TTarget, ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >::Type
convertImpl(Convert<TTarget, T> const,
			ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & source_)
{
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TSource;
	return _AlphabetConversionTable<TTarget, TSource>::table[ordValue(source_)];
}


//////////////////////////////////////////////////////////////////////////////


/*
template 
<
	typename TTargetHost, char TARGET_CHAR, typename TTargetSpec, typename T, 
	typename TSourceHost, char SOURCE_CHAR, typename TSourceSpec
>
inline typename Convert<ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> > , ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > >::Type
convertImpl(Convert<ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> >, T> const,
			ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > const & source_)
{
	ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> > TTarget;
	ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > TSource;
	return convert<TTarget>(convert<TTargetHost>(convert<TSourceHost>(source_)));
}
*/

//no conversion 
template <typename THost, char CHAR, typename TSpec, typename T>
inline ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
			ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & source_)
{
	return source_;
}


//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
inline unsigned
ordValue(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & c) 
{
	return c.data;
}

//////////////////////////////////////////////////////////////////////////////
// comparisons
//////////////////////////////////////////////////////////////////////////////

template <typename TModExpand, typename THost, typename TRight, typename TCompareHostRight>
struct _CompareType_ModExpand_Impl
{
	typedef TCompareHostRight Type; //fallback
};
template <typename TModExpand, typename THost, typename TRight>
struct _CompareType_ModExpand_Impl<TModExpand, THost, TRight, THost>
{
	typedef TModExpand Type;
};
template <typename TModExpand, typename THost, typename TRight>
struct _CompareType_ModExpand_Impl<TModExpand, THost, TRight, TRight>
{
	typedef TRight Type;
};


template <typename THost, char CHAR, typename TSpec, typename TRight>
struct CompareType<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TRight>
{
	typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TModExpand;
	typedef typename CompareType<THost, TRight>::Type TCompareHostRight;
	typedef typename _CompareType_ModExpand_Impl<TModExpand, THost, TRight, TCompareHostRight>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace 

#endif
