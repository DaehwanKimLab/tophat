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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_ONEGAPPED_H
#define SEQAN_HEADER_SHAPE_ONEGAPPED_H

namespace SEQAN_NAMESPACE_MAIN
{


	//////////////////////////////////////////////////////////////////////////////
	// shape with one gap
	//////////////////////////////////////////////////////////////////////////////

	struct OneGappedShape;

	template <typename TValue>
	class Shape<TValue, OneGappedShape>
	{
	public:
//____________________________________________________________________________

		unsigned blockLen1;
		unsigned gapLen;
		unsigned blockLen2;
	
		typename Value<Shape>::Type	hValue;		// current hash value

		TValue						leftChar;
		typename Value<Shape>::Type	factor1;
		typename Value<Shape>::Type	factor2;

//____________________________________________________________________________

		Shape()	{}

		// c'tor for ungapped shapes
		Shape(unsigned _blockLen1, unsigned _gapLen, unsigned _blockLen2):
			blockLen1(_blockLen1),
			gapLen(_gapLen),
			blockLen2(_blockLen2)
		{
		SEQAN_CHECKPOINT
			typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;
			factor1 = _intPow((THValue)ValueSize<TValue>::VALUE, weight(*this) - 1);
			factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, blockLen2);
		}

		Shape(Shape const &other):
			blockLen1(other.blockLen1),
			gapLen(other.gapLen),
			blockLen2(other.blockLen2),
			hValue(other.hValue),
			leftChar(other.leftChar),
			factor1(other.factor1), 
			factor2(other.factor2) {}

		template <typename TSpec>
		Shape(Shape<TValue, TSpec> const &other)
		{
			*this = other;
		}	

		template <typename TSpec>
		Shape(GappedShape<TSpec> const &other)
		{
			*this = other;
		}

		template <typename TStringValue, typename TSpec>
		Shape(String<TStringValue, TSpec> const &bitmap)
		{
			*this = bitmap;
		}	

//____________________________________________________________________________


		template <typename TStringValue, typename TSpec>
		inline Shape &
		operator=(String<TStringValue, TSpec> const &bitmap)
		{
			stringToShape(*this, bitmap);
			return *this;
		}
	};

//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	inline typename Size< Shape<TValue, OneGappedShape> >::Type
	length(Shape<TValue, OneGappedShape> const &me)
	{
	SEQAN_CHECKPOINT
		return me.blockLen1 + me.gapLen + me.blockLen2;
	}

//____________________________________________________________________________

	template <typename TValue>
	inline typename Size< Shape<TValue, OneGappedShape> >::Type
	weight(Shape<TValue, OneGappedShape> const & me)
	{
	SEQAN_CHECKPOINT
		return me.blockLen1 + me.blockLen2;
	}

//____________________________________________________________________________

	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, OneGappedShape> >::Type
	hash(Shape<TValue, OneGappedShape> &me, TIter it)	
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, OneGappedShape> >::Type	TSize;

		me.hValue = ordValue(me.leftChar = *it);
		for(TSize i = 1; i < me.blockLen1; ++i) {
			goNext(it);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		goFurther(it, me.gapLen);
		for(TSize i = 0; i < me.blockLen2; ++i) {
			goNext(it);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		return me.hValue;
	}

	template <typename TValue, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, OneGappedShape> >::Type
	hash(Shape<TValue, OneGappedShape> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;
		TSize blockLen1 = me.blockLen1;
		TSize blockLen2 = me.blockLen2;

		if ((TSize)length(me) > charsLeft)
		{
			if (blockLen1 > charsLeft)
			{
				blockLen1 = charsLeft;
				blockLen2 = 0;
				if (blockLen1 == 0) return me.hValue = 0;
			} else
				if (blockLen1 + (TSize)me.gapLen > charsLeft)
					blockLen2 = 0;
				else
					blockLen2 = charsLeft - (blockLen1 + me.gapLen);
		}
		
		me.hValue = ordValue(me.leftChar = *it);
		for(TSize i = 1; i < blockLen1; ++i) {
			goNext(it);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		goFurther(it, me.gapLen);
		for(TSize i = 0; i < blockLen2; ++i) {
			goNext(it);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
		}
		// fill shape with zeros
		TSize w = (TSize)weight(me);
		for(TSize i = blockLen1 + blockLen2; i < w; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

	template <typename TValue, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, OneGappedShape> >::Type
	hashUpper(Shape<TValue, OneGappedShape> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;
		TSize blockLen1 = me.blockLen1;
		TSize blockLen2 = me.blockLen2;

		if ((TSize)length(me) > charsLeft)
		{
			if (blockLen1 > charsLeft)
			{
				blockLen1 = charsLeft;
				blockLen2 = 0;
			} else
				if (blockLen1 + (TSize)me.gapLen > charsLeft)
					blockLen2 = 0;
				else
					blockLen2 = charsLeft - (blockLen1 + me.gapLen);
		}

		me.hValue = 0;
		me.leftChar = *it;
		for(TSize i = 0; i < blockLen1; ++i) {
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			goNext(it);
		}
		goFurther(it, me.gapLen);
		for(TSize i = 0; i < blockLen2; ++i) {
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
			goNext(it);
		}
		++me.hValue;

		// fill shape with zeros
		TSize w = (TSize)weight(me);
		for(TSize i = blockLen1 + blockLen2; i < w; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, OneGappedShape> >::Type
	hashNext(Shape<TValue, OneGappedShape> &me, TIter &_it)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;
		TIter it(_it);

		// remove leftmost character
		me.hValue -= ordValue(me.leftChar) * me.factor1;
		me.leftChar = *it;

		// shift
		me.hValue *= ValueSize<TValue>::VALUE;

		// add emerging character
		goFurther(it, me.blockLen1 - 1);
		me.hValue += ordValue((TValue)*it) * me.factor2;

		// subtract vanishing character
		goFurther(it, me.gapLen);
		me.hValue -= ordValue((TValue)*it) * me.factor2;

		// add rightmost emerging character
		goFurther(it, me.blockLen2);
		me.hValue += ordValue((TValue)*it);

		return me.hValue;
	}


//____________________________________________________________________________

	template <typename TValue, typename TShapeString>
	inline bool
	stringToShape(
		Shape<TValue, OneGappedShape> &me, 
		TShapeString const &bitmap)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TShapeString const>::Type				TIter;
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;

		TIter it = begin(bitmap, Standard());
		TIter itEnd = end(bitmap, Standard());

		me.blockLen1 = 0;
		me.gapLen = 0;
		me.blockLen2 = 0;

		for(; it != itEnd && *it == '0' ; ++it) ;

		for(; it != itEnd && *it == '1' ; ++it)
			++me.blockLen1;

		for(; it != itEnd && *it == '0' ; ++it)
			++me.gapLen;

		for(; it != itEnd && *it == '1' ; ++it)
			++me.blockLen2;

		for(; it != itEnd && *it == '0' ; ++it) ;

		me.factor1 = _intPow((THValue)ValueSize<TValue>::VALUE, weight(me) - 1);
		me.factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, me.blockLen2);

		return it == itEnd && me.blockLen1 > 0;
	}

	template <typename TShapeString, typename TValue>
	inline void
	shapeToString(
		TShapeString &bitmap,
		Shape<TValue, OneGappedShape> const &me)
	{
	SEQAN_CHECKPOINT

		clear(bitmap);
		fill(bitmap, me.blockLen1, '1');
		fill(bitmap, me.blockLen1 + me.gapLen, '0');
		fill(bitmap, me.blockLen1 + me.gapLen + me.blockLen2, '1');
	}

//____________________________________________________________________________
	
	template <typename TValue>
	inline void
	reverse(Shape<TValue, OneGappedShape> &me)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, OneGappedShape> >::Type	THValue;

		unsigned temp = me.blockLen1;
		me.blockLen1 = me.blockLen2;
		me.blockLen2 = temp;
		me.factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, me.blockLen2);
	}
	
}	// namespace seqan

#endif
