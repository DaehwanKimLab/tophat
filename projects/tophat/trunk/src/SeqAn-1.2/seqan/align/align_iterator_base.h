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
  $Id: align_iterator_base.h 1045 2007-08-21 16:41:41Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_ITERATOR_BASE_H
#define SEQAN_HEADER_ALIGN_ITERATOR_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Align Iterator for Gaps alignment
//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
class Iter<TAlign, AlignColIterator<TSpec> >
{
public:
	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TRowIterator;
	typedef typename Position<TRow>::Type TRowPosition;
	typedef String<TRowIterator> TIterators;

	TAlign * data_host;
	TIterators data_iterators;

public:
	Iter()
	{
SEQAN_CHECKPOINT
	}
	Iter(TAlign & _align):
		data_host(& _align)
	{
SEQAN_CHECKPOINT
		typename Position<TRows>::Type _i = length(rows(_align));
		resize(data_iterators, _i, Exact());
	}
	Iter(TAlign & _align, TRowPosition _pos):
		data_host(& _align)
	{
SEQAN_CHECKPOINT
		typename Position<TRows>::Type _i = length(rows(_align));
		resize(data_iterators, _i, Exact());

		while (_i > 0)
		{
			--_i;
			data_iterators[_i] = iter(row(_align, _i), _pos);
		}
	}
	Iter(Iter const & _other):
		data_host(_other.data_host),
		data_iterators(_other.data_iterators)
	{
	}
	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const &
	operator = (Iter const & _other)
	{
SEQAN_CHECKPOINT
		data_host = _other.data_host;
		data_iterators = _other.data_iterators;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline TAlign &
host(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return *me.data_host;
}
template <typename TAlign, typename TSpec>
inline TAlign &
host(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return *me.data_host;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
setHost(Iter<TAlign, AlignColIterator<TSpec> > & me, TAlign & _host)
{
SEQAN_CHECKPOINT
	me.data_host = & _host;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline typename Cols<TAlign>::Type
container(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return cols(*me.data_host);
}
template <typename TAlign, typename TSpec>
inline typename Cols<TAlign>::Type
container(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return cols(*me.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
goNext(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TRowIterator;
	typedef String<TRowIterator> TIterators;
	typedef typename Iterator<TIterators, Standard>::Type TIteratorsIterator;

	TIteratorsIterator _it = begin(me.data_iterators);
	TIteratorsIterator _it_end = end(me.data_iterators);

	while (_it != _it_end)
	{
		goNext(*_it);
		++_it;
	}
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > & 
operator ++(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	goNext(me);
	return me;
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator ++(Iter<TAlign, AlignColIterator<TSpec> > & me, int)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	goNext(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
goPrevious(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TRowIterator;
	typedef String<TRowIterator> TIterators;
	typedef typename Iterator<TIterators, Standard>::Type TIteratorsIterator;

	TIteratorsIterator _it = begin(me.data_iterators);
	TIteratorsIterator _it_end = end(me.data_iterators);

	while (_it != _it_end)
	{
		goPrevious(*_it);
		++_it;
	}
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > & 
operator --(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	goPrevious(me);
	return me;
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator --(Iter<TAlign, AlignColIterator<TSpec> > & me, int)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	goPrevious(me);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
	return getValue(_left.data_iterators, 0) == getValue(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
			Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
	return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition>
inline typename Reference<TAlign>::Type
value(Iter<TAlign, AlignColIterator<TSpec> > & me,
	  TPosition pos_)
{
SEQAN_CHECKPOINT
	return value(me.data_iterators[pos_]);
}
template <typename TAlign, typename TSpec, typename TPosition>
inline typename Reference<TAlign>::Type
value(Iter<TAlign, AlignColIterator<TSpec> > const & me,
	  TPosition pos_)
{
SEQAN_CHECKPOINT
	return value(me.data_iterators[pos_]);
}
//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition>
inline typename GetValue<TAlign>::Type
getValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
		 TPosition pos_)
{
SEQAN_CHECKPOINT
	return getValue(me.data_iterators[pos_]);
}
template <typename TAlign, typename TSpec, typename TPosition>
inline typename GetValue<TAlign>::Type
getValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
		 TPosition pos_)
{
SEQAN_CHECKPOINT
	return getValue(me.data_iterators[pos_]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
			TPosition pos_,
			TValue & val)
{
SEQAN_CHECKPOINT
	return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
			TPosition pos_,
			TValue const & val)
{
SEQAN_CHECKPOINT
	return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
			TPosition pos_,
			TValue & val)
{
SEQAN_CHECKPOINT
	return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
			TPosition pos_,
			TValue const & val)
{
SEQAN_CHECKPOINT
	return assignValue(me.data_iterators[pos_], val);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
		  TPosition pos_,
		  TValue & val)
{
SEQAN_CHECKPOINT
	return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
		  TPosition pos_,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
		  TPosition pos_,
		  TValue & val)
{
SEQAN_CHECKPOINT
	return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
		  TPosition pos_,
		  TValue const & val)
{
SEQAN_CHECKPOINT
	return moveValue(me.data_iterators[pos_], val);
}

//////////////////////////////////////////////////////////////////////////////

//??? TODO
//disabled since GapsIterator has no operator - and +
/*
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > & 
operator +=(Iter<TAlign, AlignColIterator<TSpec> > & me,
			TSize size)
{
SEQAN_CHECKPOINT
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;
	typedef String<TRowIterator> TIterators;
	typedef typename Iterator<TIterators>::Type TIteratorsIterator;

	TIteratorsIterator _it = begin(me.data_iterators);
	TIteratorsIterator _it_end = end(me.data_iterators);

	while (_it != _it_end)
	{
		*_it += size;
		++_it;
	}
	return me;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator +(Iter<TAlign, AlignColIterator<TSpec> > & me,
		   TSize size)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	me += size;
	return me;
}
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator +(Iter<TAlign, AlignColIterator<TSpec> > const & me,
		   TSize size)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	me += size;
	return me;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > & 
operator -=(Iter<TAlign, AlignColIterator<TSpec> > & me,
			TSize size)
{
SEQAN_CHECKPOINT
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;
	typedef String<TRowIterator> TIterators;
	typedef typename Iterator<TIterators>::Type TIteratorsIterator;

	TIteratorsIterator _it = begin(me.data_iterators);
	TIteratorsIterator _it_end = end(me.data_iterators);

	while (_it != _it_end)
	{
		*_it -= size;
		++_it;
	}
	return me;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator -(Iter<TAlign, AlignColIterator<TSpec> > & me,
		   TSize size)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	me -= size;
	return me;
}
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > 
operator -(Iter<TAlign, AlignColIterator<TSpec> > const & me,
		   TSize size)
{
SEQAN_CHECKPOINT
	Iter<TAlign, AlignColIterator<TSpec> > ret = me;
	me -= size;
	return me;
}

//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline typename Difference<TAlign>::Type
operator -(Iter<TAlign, AlignColIterator<TSpec> > const & left,
		   Iter<TAlign, AlignColIterator<TSpec> > const & right)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(length(left.data_iterators))
	SEQAN_ASSERT(length(right.data_iterators))

	return (left.data_iterators[0] - right.data_iterators[0]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline typename Position<TAlign>::Type
position(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
	return position(me.data_iterators[0], row(host(me), 0));
}
template <typename TAlign, typename TSpec>
inline typename Position<TAlign>::Type
position(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return position(me.data_iterators[0], row(host(me), 0));
}
*/
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
