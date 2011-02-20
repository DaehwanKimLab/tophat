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
  $Id: gaps_array.h 3224 2009-01-14 12:25:06Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GAPS_ARRAY_H
#define SEQAN_HEADER_GAPS_ARRAY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tag

struct ArrayGaps;


//////////////////////////////////////////////////////////////////////////////
// Gaps - ArrayGaps Spec
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ArrayGaps:
..cat:Alignments
..general:Class.Gaps
..summary:Stores gaps sizes in an array.
..signature:Gaps<TSource, ArrayGaps >
..param.TSource:Type of the ungapped sequence.
...metafunction:Metafunction.Source
*/
template <typename TSource>
class Gaps<TSource, ArrayGaps>
{
public:
	typedef typename Size<Gaps>::Type TSize;
	typedef String<TSize> TArr;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef typename Position<Gaps>::Type TViewPosition;

public:
	String<TSize> data_arr; //a list of gap and non-gap region lengths
	TViewPosition data_end_position;

	Holder<TSource> data_source;
	TSourcePosition data_source_begin_position;
	TSourcePosition data_source_end_position;

public:
	Gaps():
		data_source_begin_position(0),
		data_source_end_position(0)
	{
SEQAN_CHECKPOINT
	}
	Gaps(TSize _size):
		data_source_begin_position(0),
		data_source_end_position(0)
	{
		_init_to_resize(*this, _size);
	}
	Gaps(TSource & source_):
		data_source(source_),
		data_source_begin_position(beginPosition(source_)),
		data_source_end_position(endPosition(source_))
	{
SEQAN_CHECKPOINT
		_init_to_resize(*this, length(source_));
	}

	template <typename TSource2>
	Gaps(TSource2 const & source_):
		data_source_begin_position(0),
		data_source_end_position(length(source_))
		//data_source_begin_position(beginPosition(source_)),
		//data_source_end_position(endPosition(source_))
	{
SEQAN_CHECKPOINT
		data_source = source_;
		_init_to_resize(*this, length(source_));
	}

	Gaps(Gaps const & other_):
		data_arr(other_.data_arr),
		data_end_position(other_.data_end_position),
//		data_source(value(other_.data_source)), //variant: setValue => Align benutzen gemeinsame Sources
		data_source(other_.data_source),		//variant: assignValue => Align kopieren Sources
		data_source_begin_position(other_.data_source_begin_position),
		data_source_end_position(other_.data_source_end_position)
	{
SEQAN_CHECKPOINT
	}
	Gaps & operator = (Gaps const & other_)
	{
SEQAN_CHECKPOINT
		data_arr = other_.data_arr;
		data_end_position = other_.data_end_position;
		setValue(data_source, source(other_));
		data_source_begin_position = other_.data_source_begin_position;
		data_source_end_position = other_.data_source_end_position; 
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
inline String<typename Size<Gaps<TSource, ArrayGaps> >::Type> &
_dataArr(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return me.data_arr;
}
template <typename TSource>
inline String<typename Size<Gaps<TSource, ArrayGaps> >::Type> const &
_dataArr(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_arr;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline Holder<TSource> &
_dataSource(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TSource>
inline Holder<TSource> const &
_dataSource(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<Gaps<TSource, ArrayGaps> >::Type
endPosition(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return me.data_end_position;
}
template <typename TSource>
inline typename Position<Gaps<TSource, ArrayGaps> >::Type
endPosition(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_end_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition>
inline void 
_setEndPosition(Gaps<TSource, ArrayGaps> & me, TPosition _pos)
{
SEQAN_CHECKPOINT
	me.data_end_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
sourceBeginPosition(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_source_begin_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setSourceBeginPosition(Gaps<TSource, ArrayGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.data_source_begin_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
sourceEndPosition(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_source_end_position;
}
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setSourceEndPosition(Gaps<TSource, ArrayGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.data_source_end_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSize>
inline void
_init_to_resize(Gaps<TSource, ArrayGaps> & me,
				TSize _size)
{
SEQAN_CHECKPOINT
	resize(_dataArr(me), 2);
	_dataArr(me)[0] = 0;
	_dataArr(me)[1] = _size;
	_setEndPosition(me, _size);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearGaps(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	_init_to_resize(me, sourceEndPosition(me) - sourceBeginPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clear(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	_init_to_resize(me, 0);
	_setSourceBeginPosition(me, 0);
	_setSourceEndPosition(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

// transforms source- to view-position
template <typename TSource>
inline typename Position< Gaps<TSource, ArrayGaps> >::Type
toViewPosition(Gaps<TSource, ArrayGaps> const & gaps,
			   typename Position<TSource>::Type pos)
{
SEQAN_CHECKPOINT

	SEQAN_ASSERT(pos >= sourceBeginPosition(gaps))
	pos -= sourceBeginPosition(gaps);

	typedef Gaps<TSource, ArrayGaps> const TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef String<TSize> const TArr;

	typename Iterator<TArr, Standard>::Type arr_begin = begin(_dataArr(gaps));
	typename Iterator<TArr, Standard>::Type arr_end = end(_dataArr(gaps));
	typename Position<TGaps>::Type view_pos = pos;
	typename Position<TSource>::Type source_pos = pos;
	
	while (true)
	{
		if (arr_begin == arr_end) return view_pos;
		TSize step = *arr_begin;
		view_pos += step;
	    
		++arr_begin;
	    
		step = *arr_begin;
		if (source_pos < step) return view_pos;
		source_pos -= step;
	    
		++arr_begin;
	}
}

//____________________________________________________________________________
// transformates view- to source-position. (aufgerundet!)
template <typename TSource>
inline typename Position<TSource>::Type
toSourcePosition(Gaps<TSource, ArrayGaps> const & gaps,
				 typename Position< Gaps<TSource, ArrayGaps> >::Type pos)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, ArrayGaps> const TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef String<TSize> const TArr;

	typename Iterator<TArr, Standard>::Type arr_begin = begin(_dataArr(gaps));
	typename Iterator<TArr, Standard>::Type arr_end = end(_dataArr(gaps));
	typename Position<TSource>::Type source_pos = sourceBeginPosition(gaps);
	typename Position<TGaps>::Type view_pos = pos;
	
	while (true)
	{
		if (arr_begin == arr_end) return source_pos;

		TSize step = *arr_begin;
		if (view_pos <= step) return source_pos;
		view_pos -= step;
	    
		++arr_begin;
		step = *arr_begin;
		if (view_pos <= step ) return source_pos + view_pos;
		source_pos += step;
		view_pos -= step;
	    
		++arr_begin;
	}
}

//////////////////////////////////////////////////////////////////////////////
// returns an iterator to view-Position 

template <typename TGaps, typename TPosition>
inline typename Iterator<TGaps, Standard>::Type
_iterator_gaps_array(TGaps & gaps, 
					 TPosition view_position)
{
	typedef typename Size<TGaps>::Type TSize;
	typedef typename TGaps::TArr const TArr;
	typedef typename Source<TGaps>::Type TSource;
	typedef typename Iterator<TArr, Standard>::Type TArrIterator;
	typedef typename Value<TArrIterator>::Type TArrIteratorType;

 	TArrIterator arr_begin = begin(_dataArr(gaps), Standard());
	TArrIterator arr_end = end(_dataArr(gaps), Standard());
	typename Position<TGaps>::Type block = 0;

	TArrIteratorType view_pos = view_position;

	if (emptySource(gaps))
	{
		for (; arr_begin != arr_end; ++arr_begin)
		{
			if (view_pos < *arr_begin) break;
			++block;
			view_pos -= *arr_begin;
		}
		return typename Iterator<TGaps, Standard>::Type(
			gaps, 
			block, 
			view_pos);
	}
	else
	{
		typename Position<TSource>::Type source_pos = sourceBeginPosition(gaps);

		while (true)
		{
			if ((arr_begin == arr_end) || (view_pos < *arr_begin))
			{
				break;
			}
			++block;
			view_pos -= *arr_begin;
			++arr_begin;

			if (view_pos < *arr_begin)
			{
				source_pos += view_pos;
				break;
			}
			++block;
			view_pos -= *arr_begin;
			source_pos += *arr_begin;
			++arr_begin;
		}

		return typename Iterator<TGaps, Standard>::Type(
			gaps,
            iter(source(gaps), source_pos),
			block, 
			view_pos);
	}

}

template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, ArrayGaps>, Tag<TTag> const>::Type
iter(Gaps<TSource, ArrayGaps> & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	return _iterator_gaps_array(gaps, view_pos);
}
template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, ArrayGaps> const, Tag<TTag> const>::Type
iter(Gaps<TSource, ArrayGaps> const & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	return _iterator_gaps_array<Gaps<TSource, ArrayGaps> const>(gaps, view_pos);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position< Gaps<TSource, ArrayGaps> >::Type
beginPosition(Gaps<TSource, ArrayGaps> & gaps)
{
SEQAN_CHECKPOINT
	return value(_dataArr(gaps), 0);
}
template <typename TSource>
inline typename Position< Gaps<TSource, ArrayGaps> const>::Type
beginPosition(Gaps<TSource, ArrayGaps> const & gaps)
{
SEQAN_CHECKPOINT
	return value(_dataArr(gaps), 0);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


template <typename TSource, typename TPosition>
inline void
setBeginPosition(Gaps<TSource, ArrayGaps> & me,
				 TPosition view_position)
{
SEQAN_CHECKPOINT
	TPosition old_pos = beginPosition(me);
	if (length(_dataArr(me)) == 0)
	{
		_init_to_resize(me, length(source(me)));
	}
	_dataArr(me)[0] = view_position;
	_setEndPosition(me, endPosition(me) + view_position - old_pos);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition>
inline void
setSourceBeginPosition(Gaps<TSource, ArrayGaps> & me,
					   TPosition source_position)
{
	SEQAN_ASSERT(length(_dataArr(me)))

	typedef Gaps<TSource, ArrayGaps> TGaps;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename TGaps::TArr TArr;

	TPosition old_source_begin_pos = sourceBeginPosition(me);
	if (old_source_begin_pos == source_position) return;
	else if (source_position < old_source_begin_pos)
	{
SEQAN_CHECKPOINT
		TViewPosition delta = old_source_begin_pos - source_position;
		TViewPosition view_begin_pos = beginPosition(me);
		if (view_begin_pos <= delta)
		{
			setBeginPosition(me, 0);
		}
		else
		{
			setBeginPosition(me, view_begin_pos - delta);
		}
		_setEndPosition(me, endPosition(me) + delta);
		_dataArr(me)[1] += delta;
		_setSourceBeginPosition(me, source_position);
	}
	else //(source_position > old_source_begin_pos)
	{
SEQAN_CHECKPOINT
		TViewPosition view_pos = toViewPosition(me, source_position);
		TViewPosition source_pos_left = source_position - sourceBeginPosition(me);
		TViewPosition gaps_count = 0;

		typename Iterator<TArr, Standard>::Type it_arr_begin = begin(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr_end = end(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr = it_arr_begin;

		while (it_arr != it_arr_end)
		{
			gaps_count += *it_arr;
			++it_arr;
			if (*it_arr > source_pos_left) 
			{
				*it_arr_begin = view_pos;
				*it_arr -= source_pos_left;
				if (it_arr - it_arr_begin > 1)
				{
					replace(_dataArr(me), 1, it_arr - it_arr_begin, "");
				}
				_setSourceBeginPosition(me, source_position);
				return;
			}
			gaps_count += *it_arr;
			source_pos_left -= *it_arr;
			++it_arr;
		}
		//alignment is empty
		clear(me);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition>
inline void
setSourceEndPosition(Gaps<TSource, ArrayGaps> & me,
					 TPosition source_position)
{
	typedef Gaps<TSource, ArrayGaps> TGaps;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename TGaps::TArr TArr;

	TArr arr = _dataArr(me);
	SEQAN_ASSERT(length(arr));

	TPosition old_end_begin_pos = sourceEndPosition(me);
	if (old_end_begin_pos == source_position) return;
	else if (source_position < old_end_begin_pos)
	{
SEQAN_CHECKPOINT
		typename Iterator<TArr, Standard>::Type it_arr_begin = begin(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr_end = end(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr = it_arr_begin;
		TViewPosition end_pos = 0;
		TViewPosition chars_to_scan = source_position - sourceBeginPosition(me);

		while (it_arr != it_arr_end)
		{
			end_pos += *it_arr;
			++it_arr;
			if (*it_arr >= chars_to_scan)
			{
				resize(_dataArr(me), it_arr - it_arr_begin + 1);
				*it_arr = chars_to_scan;
				_setEndPosition(me, end_pos + chars_to_scan);
				break;
			}
			end_pos += *it_arr;
			chars_to_scan -= *it_arr;
			++it_arr;
		}
	}
	else //(source_position > old_end_begin_pos)
	{
SEQAN_CHECKPOINT
		*(end(_dataArr(me)) - 1) += (source_position - old_end_begin_pos);
		_setEndPosition(me, endPosition(me) + source_position - old_end_begin_pos);
	}
	_setSourceEndPosition(me, source_position);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Size<TSource>::Type
sourceLength(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return sourceEndPosition(me) - sourceBeginPosition(me);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Gaps Iterator
//////////////////////////////////////////////////////////////////////////////


template <typename TGaps>
class Iter<TGaps, GapsIterator<ArrayGaps> >
{
public:
	typedef typename Size<TGaps>::Type TGapsSize;
	typedef String<TGapsSize> TArr;
	typedef typename Position<TArr>::Type TArrPosition;
	typedef typename Source<Iter>::Type TSourceIterator;

	TGaps * data_container;								//the gaps object
	mutable TSourceIterator data_source;			
	mutable TArrPosition data_block;					//block in array of container
	mutable TGapsSize data_sub;							//position block

public:
	Iter() 
	{
SEQAN_CHECKPOINT
	}
	Iter(Iter const & other_):
		data_container(other_.data_container),
		data_source(other_.data_source),
		data_block(other_.data_block),
		data_sub(other_.data_sub)
	{
SEQAN_CHECKPOINT
	}
	Iter(TGaps & container_,
			TSourceIterator source_iterator,
			TArrPosition block_ = 0,
			TGapsSize sub_ = 0):
		data_container(&container_),
		data_source(source_iterator),
		data_block(block_),
		data_sub(sub_)
	{
SEQAN_CHECKPOINT
	}

	Iter(TGaps & container_,
			TArrPosition block_,
			TGapsSize sub_ = 0):
		data_container(&container_),
		data_block(block_),
		data_sub(sub_)
	{
SEQAN_CHECKPOINT
	}

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		data_block = other_.data_block;
		data_sub = other_.data_sub;
		data_source = other_.data_source;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
getValue(Iter<TGaps, GapsIterator<ArrayGaps> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}
template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type
getValue(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}

//____________________________________________________________________________

template <typename TGaps>
inline bool 
isGap(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
SEQAN_CHECKPOINT
	return !(me.data_block & 1);
}

//____________________________________________________________________________

template <typename TGaps>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	if (isGap(me))
	{
		typedef typename Size<TGaps>::Type TGapsSize;
		typedef String<TGapsSize> const TArr;

		TArr & arr = _dataArr(container(me));
		typename Position<TArr>::Type pos = me.data_block;

		if (length(arr) == pos)
		{//counting trailing gaps
SEQAN_CHECKPOINT
			return 0;
		}
		else
		{
SEQAN_CHECKPOINT
			return arr[pos] - me.data_sub;
		}
	}
	else
	{
SEQAN_CHECKPOINT
		return 0;
	}
}

//____________________________________________________________________________

template <typename T>
inline void 
_goNext_arrayGapsIterator(T const & me)
{
	if (!isGap(me))
	{
SEQAN_CHECKPOINT
		goNext(me.data_source);
	}

	++me.data_sub;
	if ((me.data_block < length(_dataArr(container(me)))) && (me.data_sub >= _dataArr(container(me))[me.data_block]))
	{
SEQAN_CHECKPOINT
		++me.data_block;
		me.data_sub = 0;
	}
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<ArrayGaps> > & me)
{
	_goNext_arrayGapsIterator(me);
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	_goNext_arrayGapsIterator(me);
}

//____________________________________________________________________________

template <typename T>
inline void 
_goPrevious_arrayGapsIterator(T const & me)
{
	if (me.data_sub > 0)
	{
SEQAN_CHECKPOINT
		--me.data_sub;
	}
	else
	{
		if (me.data_block > 0)
		{
SEQAN_CHECKPOINT
			--me.data_block;
			me.data_sub = _dataArr(container(me))[me.data_block] - 1;
		}
	}
	if (!isGap(me))
	{
SEQAN_CHECKPOINT
		goPrevious(me.data_source);
	}
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<ArrayGaps> > & me)
{
	_goPrevious_arrayGapsIterator(me);
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	_goPrevious_arrayGapsIterator(me);
}
//____________________________________________________________________________

template <typename TGaps>
inline bool 
atBegin(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
SEQAN_CHECKPOINT
	return ((me.data_block == 1) && (me.data_sub == 0));
}

//____________________________________________________________________________

template <typename TGaps>
inline bool 
atEnd(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
SEQAN_CHECKPOINT
	return ((me.data_block == length(_dataArr(container(me)))) && (me.data_sub == 0));
}

//____________________________________________________________________________

template <typename TGaps, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & me,
		   TCount size)
{
	typedef typename Size<TGaps>::Type TGapsSize;
	typedef String<TGapsSize> TArr;

	TArr & arr = _dataArr(container(me));
	typename Position<TArr>::Type pos = me.data_block;
	typename Iterator<TArr, Rooted>::Type it = iter(arr, pos);

	if (isGap(me))
	{
		if (pos < length(arr))
		{//here is a gap: expand it
SEQAN_CHECKPOINT
			value(it) += size;
		}
		else
		{//gap at end: nothing to do
SEQAN_CHECKPOINT
			return;
		}
	}
	else
	{
		if (me.data_sub)
		{//here is no gap: insert one
SEQAN_CHECKPOINT
			_clearSpace(arr, 2, pos + 1, pos + 1, Generous());

			it = iter(arr, pos); //reload, since iterator could be invalid
			it[2] = it[0] - me.data_sub;
			it[0] = me.data_sub;
			it[1] = size;

			++me.data_block;
			me.data_sub = 0;
		}
		else
		{//insert gaps at begin of a char block
SEQAN_CHECKPOINT
			me.data_sub = arr[pos - 1];
			arr[pos - 1] += size; //note: pos > 0 because this is a char block
			--me.data_block;
		}
	}

	//adjust right position
	_setEndPosition(container(me), endPosition(container(me)) + size);
}

//____________________________________________________________________________

//delete up to size gaps 
	
template <typename TGaps, typename TCount>
inline void
removeGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & me,
		   TCount _size)
{
	typedef typename Size<TGaps>::Type TGapsSize;
	typedef String<TGapsSize> TArr;

	TGapsSize size = _size;

	if (isGap(me))
	{//here is a gap
		TArr & arr = _dataArr(container(me));
		if (me.data_block < length(arr))
		{//not trailing gaps
			typename Iterator<TArr, Standard>::Type it = iter(arr, me.data_block);

			TGapsSize rest_size = *it - me.data_sub;

			if (rest_size <= size)
			{//remove rest of gap block
				if (me.data_sub || !me.data_block)
				{//keep this gap block
SEQAN_CHECKPOINT
					*it = me.data_sub;
					++me.data_block;
					me.data_sub = 0;
				}
				else
				{//remove complete gap block, merge two char blocks
SEQAN_CHECKPOINT
					TGapsSize data_sub = *(it-1); //(note: this well defined, since "it" was not the first block)

					*(it-1) += *(it+1);
					replace(_dataArr(container(me)), me.data_block, me.data_block + 2, "");

					--me.data_block;
					me.data_sub = data_sub;
				}
				_setEndPosition(container(me), endPosition(container(me)) - rest_size);
			}
			else
			{//remove a part of this block
SEQAN_CHECKPOINT
				*it -= size;
				_setEndPosition(container(me), endPosition(container(me)) - size);
			}
		}
		//else: trailing gaps: infinite gaps here, so do nothing
	}
	//else: here is no gap - nothing to do
}

//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<ArrayGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block == _right.data_block) && (_left.data_sub == _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<ArrayGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block == _right.data_block) && (_left.data_sub == _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<ArrayGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block == _right.data_block) && (_left.data_sub == _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<ArrayGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block == _right.data_block) && (_left.data_sub == _right.data_sub);
}
//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<ArrayGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block != _right.data_block) || (_left.data_sub != _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<ArrayGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block != _right.data_block) || (_left.data_sub != _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<ArrayGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block != _right.data_block) || (_left.data_sub != _right.data_sub);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<ArrayGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<ArrayGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.data_block != _right.data_block) || (_left.data_sub != _right.data_sub);
}


//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
#endif //#ifndef SEQAN_HEADER_...
