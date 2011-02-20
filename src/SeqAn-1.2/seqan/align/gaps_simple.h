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
  $Id: gaps_array.h 1432 2007-12-19 15:11:24Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GAPS_SIMPLE_H
#define SEQAN_HEADER_GAPS_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//given a value type T, returns a value type that contains a "-" and 
//is also able to store all values of T

template <typename T>
struct GapAlphabet
{
	typedef char Type;
};
template <>
struct GapAlphabet<wchar_t>
{
	typedef wchar_t Type;
};


//////////////////////////////////////////////////////////////////////////////
// Tag

struct SequenceGaps;


//////////////////////////////////////////////////////////////////////////////
// Gaps - SequenceGaps Spec
//////////////////////////////////////////////////////////////////////////////
/**
.Spec.SequenceGaps:
..cat:Alignments
..general:Class.Gaps
..summary:Stores gapped sequences as sequences including blank signs.
..signature:Gaps<TSource, SequenceGaps>
..param.TSource:Type of the ungapped sequence.
...metafunction:Metafunction.Source
*/

template <typename TSource>
class Gaps<TSource, SequenceGaps>
{
public:
	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Position<Gaps>::Type TViewPosition;
	typedef typename Position<TSource>::Type TSourcePosition;

public:
	TString data; 
	TViewPosition begin_position;
	TSourcePosition source_begin_position;
	TSourcePosition source_end_position;

public:
	Gaps()
		: begin_position(0)
		, source_begin_position(0)
		, source_end_position(0)
	{
SEQAN_CHECKPOINT
	}
	Gaps(TSource & source_)
		: data(source_)
		, begin_position(0)
		, source_begin_position(0)
		, source_end_position(length(source_))
	{
SEQAN_CHECKPOINT
	}

	template <typename TSource2>
	Gaps(TSource2 const & source_)
		: data(source_)
		, begin_position(0)
		, source_begin_position(0)
		, source_end_position(length(source_))
	{
SEQAN_CHECKPOINT
	}

	Gaps(Gaps const & other_)
		: data(other_.data)
		, begin_position(other_.begin_position)
		, source_begin_position(other_.source_begin_position)
		, source_end_position(other_.source_end_position)
	{
SEQAN_CHECKPOINT
	}
	Gaps & operator = (Gaps const & other_)
	{
SEQAN_CHECKPOINT
		data = other_.data;
		begin_position = other_.begin_position;
		source_begin_position = other_.source_begin_position;
		source_end_position = other_.source_end_position;
		return *this;
	}

	~Gaps()
	{
SEQAN_CHECKPOINT
	}

	inline typename Reference<Gaps>::Type
	operator[](TViewPosition view_pos)
	{
SEQAN_CHECKPOINT
		return value(*this, view_pos);
	}
	inline typename Reference<Gaps const>::Type
	operator[](TViewPosition view_pos) const
	{
SEQAN_CHECKPOINT
		return value(*this, view_pos);
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
struct Value<Gaps<TSource, SequenceGaps> >
{
	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type Type;
};
template <typename TSource>
struct Value<Gaps<TSource, SequenceGaps> const>
{
	typedef typename Value<TSource const>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
struct GetSource<Gaps<TSource, SequenceGaps> >
{
	typedef TSource Type;
};

template <typename TSource>
struct GetSource<Gaps<TSource, SequenceGaps> const>
{
	typedef TSource Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename GetSource< Gaps<TSource, SequenceGaps> >::Type //returns temporary!
_source_SequenceGaps(Gaps<TSource, SequenceGaps> & me)
{
	typedef typename Iterator<TSource, Standard>::Type TSourceIterator;

	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Iterator<TString, Standard>::Type TStringIterator;

	TSource src;
	resize(src, sourceEndPosition(me) - sourceBeginPosition(me));

	TStringIterator it_from = begin(me.data);
	TStringIterator it_from_end = end(me.data);
	TSourceIterator it_to = begin(src);

	TGapValue gap_value = gapValue<TGapValue>();

	for (; it_from != it_from_end; goNext(it_from))
	{
		if (*it_from != gap_value)
		{
			*it_to = *it_from;
			goNext(it_to);
		}
	}

	return src;
}

template <typename TSource>
inline typename GetSource< Gaps<TSource, SequenceGaps> >::Type
source(Gaps<TSource, SequenceGaps> & me)
{
SEQAN_CHECKPOINT
	return _source_SequenceGaps(me);
}
template <typename TSource>
inline typename GetSource< Gaps<TSource, SequenceGaps> const >::Type
source(Gaps<TSource, SequenceGaps> const & me)
{
SEQAN_CHECKPOINT
	return _source_SequenceGaps(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename GetSource< Gaps<TSource, SequenceGaps> >::Type
sourceSegment(Gaps<TSource, SequenceGaps> & me)
{
SEQAN_CHECKPOINT
	return _source_SequenceGaps(me);
}
template <typename TSource>
inline typename GetSource< Gaps<TSource, SequenceGaps> const >::Type
sourceSegment(Gaps<TSource, SequenceGaps> const & me)
{
SEQAN_CHECKPOINT
	return _source_SequenceGaps(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSource2, typename TPosition1, typename TPosition2>
inline void
assignSource(Gaps<TSource, SequenceGaps> & me,
			 TSource2 const & source_,
			 TPosition1 source_begin_pos,
			 TPosition2 source_end_pos)
{
SEQAN_CHECKPOINT
	me.data = infix(source_, source_begin_pos, source_end_pos);
	me.begin_position = 0;
	me.source_begin_position = source_begin_pos;
	me.source_end_position = source_end_pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSource2, typename TPosition1, typename TPosition2>
inline void
moveSource(Gaps<TSource, SequenceGaps> & me,
		   TSource2 const & source_,
		   TPosition1 source_begin_pos,
		   TPosition2 source_end_pos)
{
SEQAN_CHECKPOINT
	assignSource(me, source_, source_begin_pos, source_end_pos);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline bool
dependentSource(Gaps<TSource, SequenceGaps> & )
{
SEQAN_CHECKPOINT
	return false;
}
template <typename TSource>
inline bool
dependentSource(Gaps<TSource, SequenceGaps> const & )
{
SEQAN_CHECKPOINT
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position< Gaps<TSource, SequenceGaps> >::Type
beginPosition(Gaps<TSource, SequenceGaps> & gaps)
{
SEQAN_CHECKPOINT
	return gaps.begin_position;
}
template <typename TSource>
inline typename Position< Gaps<TSource, SequenceGaps> const>::Type
beginPosition(Gaps<TSource, SequenceGaps> const & gaps)
{
SEQAN_CHECKPOINT
	return gaps.begin_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<Gaps<TSource, SequenceGaps> >::Type
endPosition(Gaps<TSource, SequenceGaps> & me)
{
SEQAN_CHECKPOINT
	return me.begin_position + length(me.data);
}
template <typename TSource>
inline typename Position<Gaps<TSource, SequenceGaps> >::Type
endPosition(Gaps<TSource, SequenceGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.begin_position + length(me.data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
sourceBeginPosition(Gaps<TSource, SequenceGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.source_begin_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
sourceEndPosition(Gaps<TSource, SequenceGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.source_end_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearGaps(Gaps<TSource, SequenceGaps> & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Iterator<TString, Standard>::Type TIterator;

	TIterator it = ::std::remove(begin(me.data, Standard()), end(me.data, Standard()), gapValue<TGapValue>());
	resize(me.data, it - begin(me.data, Standard()));

	me.begin_position = 0;	
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clear(Gaps<TSource, SequenceGaps> & me)
{
SEQAN_CHECKPOINT
	clear(me.data);
	me.begin_position = 0;
	me.source_begin_position = 0;
	me.source_end_position = 0;
}

//////////////////////////////////////////////////////////////////////////////

// transforms source- to view-position
template <typename TSource>
inline typename Position< Gaps<TSource, SequenceGaps> >::Type
toViewPosition(Gaps<TSource, SequenceGaps> const & me,
			   typename Position<TSource>::Type pos)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SequenceGaps> const TGaps;
	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Iterator<TString, Standard>::Type TIterator;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename Position<TSource>::Type TSourcePosition;

	TGapValue gap_value = gapValue<TGapValue>();

	TIterator it = begin(me.data);
	TViewPosition view_pos = me.begin_position;
	TSourcePosition source_pos = me.source_begin_position;
	while (!atEnd(it))
	{
		if (*it != gap_value)
		{
			if (source_pos == pos) break;
			++source_pos;
		}
		++ view_pos;
		++it;
	}
	
	return view_pos;
}

//____________________________________________________________________________
// transformates view- to source-position. (aufgerundet!)
template <typename TSource>
inline typename Position<TSource>::Type
toSourcePosition(Gaps<TSource, SequenceGaps> const & me,
				 typename Position< Gaps<TSource, SequenceGaps> >::Type pos)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SequenceGaps> const TGaps;
	typedef typename Value<TSource>::Type TValue;
	typedef typename GapAlphabet<TValue>::Type TGapValue;
	typedef String<TGapValue> const TString;
	typedef typename Iterator<TString, Standard>::Type TIterator;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename Position<TSource>::Type TSourcePosition;

	TGapValue gap_value = gapValue<TGapValue>();

	TIterator it = begin(me.data, Standard());
	TIterator it_end = end(me.data, Standard());
	TViewPosition view_pos = me.begin_position;
	TSourcePosition source_pos = me.source_begin_position;
	if (pos <= view_pos)
	{
		return source_pos;
	}
	while (it != it_end)
	{
		if (*it != gap_value) ++source_pos;
		++view_pos;
		if (view_pos == pos)
		{
			break;
		}
		++it;
	}

	return source_pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, SequenceGaps>, Tag<TTag> const>::Type
iter(Gaps<TSource, SequenceGaps> & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SequenceGaps> TGaps;
	typedef typename Iterator<TGaps, Tag<TTag> const>::Type TGapsIterator;
	return TGapsIterator(gaps, view_pos);
}
template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, SequenceGaps> const, Tag<TTag> const>::Type
iter(Gaps<TSource, SequenceGaps> const & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SequenceGaps> TGaps;
	typedef typename Iterator<TGaps const, Tag<TTag> >::Type TGapsIterator;
	return TGapsIterator(const_cast<TGaps &>(gaps), view_pos);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


template <typename TSource, typename TPosition>
inline void
setBeginPosition(Gaps<TSource, SequenceGaps> & me,
				 TPosition view_position)
{
SEQAN_CHECKPOINT
	me.begin_position = view_position;
}

//////////////////////////////////////////////////////////////////////////////

//DISABLED/.Function.setSourceBeginPosition.remarks.note: 
///Since @Spec.SequenceGaps@ stores the source implicitely,  
///setting the source begin or source end position has no effect for that specialization.
template <typename TSource, typename TPosition>
inline void
setSourceBeginPosition(Gaps<TSource, SequenceGaps> & ,
					   TPosition )
{
	//me.source_begin_position = source_position;
	//dummy implementation: do nothing
}

//////////////////////////////////////////////////////////////////////////////

//DISABLED/.Function.setSourceEndPosition.remarks.note: 
///Since @Spec.SequenceGaps@ stores the source implicitely,  
///setting the source begin or source end position has no effect for that specialization.
template <typename TSource, typename TPosition>
inline void
setSourceEndPosition(Gaps<TSource, SequenceGaps> & ,
					 TPosition )
{
	//me.source_end_position = source_position;
	//dummy implementation: do nothing
}

//////////////////////////////////////////////////////////////////////////////
// Gaps Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
class Iter<TGaps, GapsIterator<SequenceGaps> >
{
public:
	typedef typename Value<TGaps>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Position<TGaps>::Type TViewPosition;
	//typedef typename Position<TSource>::Type TSourcePosition;

	typedef typename Iterator<TString, Standard>::Type TStringIterator;

	TGaps * data_container;			//the gaps object
	mutable TStringIterator iter;	//iterator, valid if pos == 0
	mutable TViewPosition pos;		//steps before begin position

public:
	Iter() 
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		iter(other_.iter),
		pos(other_.pos)
	{
SEQAN_CHECKPOINT
	}
	Iter(TGaps & container_, TViewPosition view_position):
		data_container(& container_)
	{
SEQAN_CHECKPOINT
		if (view_position < container_.begin_position)
		{
			pos = container_.begin_position - view_position;
			iter = begin(container_.data, Standard());
		}
		else
		{
			pos = 0;
			iter = begin(container_.data, Standard()) + (view_position - container_.begin_position);
		}
	}

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		iter = other_.iter;
		pos = other_.pos;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<SequenceGaps> > >::Type
getValue(Iter<TGaps, GapsIterator<SequenceGaps> > & me)
{
SEQAN_CHECKPOINT

	typedef typename Value<TGaps>::Type TGapValue;
	if (me.pos)
	{
		return gapValue<TGapValue>();
	}
	return getValue(me.iter);
}
template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<SequenceGaps> > const>::Type
getValue(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT

	typedef typename Value<TGaps>::Type TGapValue;
	if (me.pos)
	{
		return gapValue<TGapValue>();
	}
	return getValue(me.iter);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename Position<TGaps>::Type
viewPosition(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT
	if (me.pos)
	{
		return me.data_container->begin_position - me.pos;
	}
	return me.iter - begin(me.data_container->data, Standard());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename Position<TGaps>::Type
sourcePosition(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT
	return toSourcePosition(*me.data_container, viewPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline bool 
isGap(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<TGaps>::Type TGapValue;
	return (me.pos || (getValue(me.iter) == gapValue<TGapValue>()));
}

//____________________________________________________________________________

template <typename TGaps>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
	if (me.pos)
	{
		return me.pos;
	}

	typedef typename Size<TGaps>::Type TSize;
	typedef typename Value<TGaps>::Type TGapValue;
	Iter<TGaps, GapsIterator<SequenceGaps> > it = me;
	TSize count = 0;
	for (; getValue(it) == gapValue<TGapValue>(); ++count)
	{
		++it;
	}
	return count;
}

//____________________________________________________________________________


template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<SequenceGaps> > & me)
{
	if (me.pos)
	{
		--me.pos;
	}
	else
	{
		goNext(me.iter);
	}
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
	if (me.pos)
	{
		--me.pos;
	}
	else
	{
		goNext(me.iter);
	}
}

//____________________________________________________________________________

template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<SequenceGaps> > & me)
{
	if (me.pos || (me.iter == begin(me.data_container->data)))
	{
		++me.pos;
	}
	else
	{
		goPrevious(me.iter);
	}
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
	if (me.pos || (me.iter == begin(me.data_container->data)))
	{
		++me.pos;
	}
	else
	{
		goPrevious(me.iter);
	}
}
//____________________________________________________________________________

template <typename TGaps>
inline bool 
atBegin(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT
	return (me.pos == 0);
}

//____________________________________________________________________________

template <typename TGaps>
inline bool 
atEnd(Iter<TGaps, GapsIterator<SequenceGaps> > const & me)
{
SEQAN_CHECKPOINT
	return (atEnd(me.iter));
}

//____________________________________________________________________________

template <typename TGaps, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<SequenceGaps> > const & me,
		   TCount size)
{
	typedef typename Size<TGaps>::Type TSize;

	if (me.pos)
	{
		me.pos += size;
		me.data_container->begin_position += size;
	}
	else
	{
		typedef typename Value<TGaps>::Type TGapValue;
		typedef String<TGapValue> TString;
		typedef typename Position<TString>::Type TPosition;

		TString blanks;
		fill(blanks, size, gapValue<TGapValue>());
		TPosition pos = me.iter - begin(me.data_container->data, Standard());
		replace(me.data_container->data, pos, pos, blanks);
	}
}

//____________________________________________________________________________

//delete up to size gaps 
	
template <typename TGaps, typename TCount>
inline void
removeGaps(Iter<TGaps, GapsIterator<SequenceGaps> > const & me,
		   TCount _size)
{
	typedef typename Value<TGaps>::Type TGapValue;
	typedef String<TGapValue> TString;
	typedef typename Position<TGaps>::Type TViewPosition;

	if (me.pos)
	{
		if (me.pos < (TViewPosition) _size)
		{
			_size = me.pos;
		}
		me.pos -= _size;
		me.data_container->begin_position -= _size;
	}
	else
	{

		if ((end(me.data_container->data, Standard()) - me.iter) < _size)
		{
			_size = end(me.data_container->data, Standard()) - me.iter;
		}
		TViewPosition pos = me.iter - begin(me.data_container->data);
		replace(me.data_container->data, pos, pos + _size, TString());
	}

}

//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SequenceGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SequenceGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SequenceGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SequenceGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SequenceGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SequenceGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SequenceGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SequenceGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SequenceGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
#endif //#ifndef SEQAN_HEADER_...
