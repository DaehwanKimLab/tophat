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

#ifndef SEQAN_HEADER_MODIFIER_ALPHABET_H
#define SEQAN_HEADER_MODIFIER_ALPHABET_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.ModifiedAlphabet:
..summary:Modifies value types.
..cat:Modifier
..signature:ModifiedAlphabet<TAlphabet, TSpec>
..param.TAlphabet:Original value type.
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
...remarks:There is no default specialization.
*/

template <typename THost, typename TSpec>
class ModifiedAlphabet;


//////////////////////////////////////////////////////////////////////////////
// sizes
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct BitsPerValue<ModifiedAlphabet<THost, TSpec> >:
	BitsPerValue<THost>
{};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct ValueSize<ModifiedAlphabet<THost, TSpec> >:
	ValueSize<THost>
{};

//////////////////////////////////////////////////////////////////////////////
// conversions
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename THost, typename TSpec>
inline typename Convert<TTarget, THost>::Type
convertImpl(Convert<TTarget, T> const convert_,
			ModifiedAlphabet<THost, TSpec> const & source_)
{
	return convertImpl( convert_, static_cast<THost const &>(source_));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline unsigned ordValue(ModifiedAlphabet<THost,TSpec> const &c) 
{
	return ordValue(static_cast<THost const &>(c));
}



//////////////////////////////////////////////////////////////////////////////
// comparisons
//////////////////////////////////////////////////////////////////////////////


template <typename THost, typename TSpec, typename TRight>
struct CompareType<ModifiedAlphabet<THost, TSpec>, TRight>
{
	typedef typename CompareType<THost, TRight>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// operator ==

template <typename THost, typename TSpec, typename TRight>
inline bool
operator == (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator == (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator == (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator == (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) == ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator == (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator == (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator == (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator == (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator !=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator != (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator != (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator != (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator != (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) != ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator != (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator != (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator != (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator != (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator <=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator <= (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator <= (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator <= (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator <= (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) <= ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator <= (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator <= (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator <= (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator <= (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator <

template <typename THost, typename TSpec, typename TRight>
inline bool
operator < (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator < (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator < (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator < (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) < ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator < (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator < (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator < (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator < (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator >=

template <typename THost, typename TSpec, typename TRight>
inline bool
operator >= (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator >= (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator >= (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator >= (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) >= ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator >= (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator >= (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator >= (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator >= (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
// operator >

template <typename THost, typename TSpec, typename TRight>
inline bool
operator > (ModifiedAlphabet<THost, TSpec> const & left_, 
			 TRight const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename THost, typename TSpec>
inline bool
operator > (TLeft const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftHost, typename TLeftSpec, typename TRightHost, typename TRightSpec>
inline bool
operator > (ModifiedAlphabet<TLeftHost, TLeftSpec> const & left_, 
			 ModifiedAlphabet<TRightHost, TRightSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<TLeftHost, TLeftSpec> TLeft;
	typedef ModifiedAlphabet<TRightHost, TRightSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename THost, typename TSpec>
inline bool
operator > (ModifiedAlphabet<THost, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	return ordValue(left_) > ordValue(right_);
}

//SimpleType

template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator > (SimpleType<TValue, TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef SimpleType<TValue, TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TValue, typename TSpec, typename THost, typename TSpec2>
inline bool
operator > (ModifiedAlphabet<THost, TSpec2> const & left_,
			 SimpleType<TValue, TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef SimpleType<TValue, TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

//Proxy

template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator > (Proxy<TSpec> const & left_, 
			 ModifiedAlphabet<THost, TSpec2> const & right_)
{
SEQAN_CHECKPOINT
	typedef Proxy<TSpec> TLeft;
	typedef ModifiedAlphabet<THost, TSpec> TRight;
	typedef typename CompareType<TRight, TLeft>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TSpec, typename THost, typename TSpec2>
inline bool
operator > (ModifiedAlphabet<THost, TSpec2> const & left_,
			 Proxy<TSpec> const & right_)
{
SEQAN_CHECKPOINT
	typedef ModifiedAlphabet<THost, TSpec> TLeft;
	typedef Proxy<TSpec> TRight;
	typedef typename CompareType<TLeft, TRight>::Type TCompareType;
	return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource>
inline void
_initializeAlphabetConversionTable(TTarget *,
								   TSource const &)
{
	//default: do nothing (because this array is not used)
	//define this function for each conversion table
}

template <typename TTarget, typename TSource>
struct _AlphabetConversionTable
{
	enum { SIZE = ValueSize<TSource>::VALUE };
private:
	static TTarget table_store[SIZE];
public:
	static TTarget * table;
	static TTarget * initialize()
	{
		static bool _is_initialized = false;
		if (! _is_initialized)
		{
			_initializeAlphabetConversionTable(table_store, TSource());
		}
		_is_initialized = true;
		return table_store;
	}
};

template <typename TTarget, typename TSource>
TTarget _AlphabetConversionTable<TTarget, TSource>::table_store[_AlphabetConversionTable<TTarget, TSource>::SIZE];

template <typename TTarget, typename TSource>
TTarget * _AlphabetConversionTable<TTarget, TSource>::table = _AlphabetConversionTable<TTarget, TSource>::initialize();


//////////////////////////////////////////////////////////////////////////////

}// namespace 

#endif
