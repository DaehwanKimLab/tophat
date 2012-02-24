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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

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
..include:seqan/align.h
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
	TViewPosition data_unclipped_end_position;

	Holder<TSource> data_source;
	TSourcePosition clipped_source_begin;
	TSourcePosition clipped_source_end;

public:
	Gaps():
		data_unclipped_end_position(0),
		clipped_source_begin(0),
		clipped_source_end(0)
	{
SEQAN_CHECKPOINT
	}
	Gaps(TSize _size):
		data_unclipped_end_position(0),
		clipped_source_begin(0),
		clipped_source_end(0)
	{
		_initToResize(*this, _size);
	}
	Gaps(TSource & source_):
		data_unclipped_end_position(0),
		data_source(source_),
		clipped_source_begin(beginPosition(source_)),
		clipped_source_end(endPosition(source_))
	{
SEQAN_CHECKPOINT
		_initToResize(*this, length(source_));
	}

	template <typename TSource2>
	Gaps(TSource2 const & source_):
		data_unclipped_end_position(0),
		clipped_source_begin(0),
		clipped_source_end(length(source_))
		//clipped_source_begin(beginPosition(source_)),
		//clipped_source_end(endPosition(source_))
	{
SEQAN_CHECKPOINT
		data_source = source_;
		_initToResize(*this, length(source_));
	}

	Gaps(Gaps const & other_):
		data_arr(other_.data_arr),
		data_end_position(other_.data_end_position),
		data_unclipped_end_position(other_.data_unclipped_end_position),
//		data_source(value(other_.data_source)), //variant: setValue => Align benutzen gemeinsame Sources
		data_source(other_.data_source),		//variant: assignValue => Align kopieren Sources
		clipped_source_begin(other_.clipped_source_begin),
		clipped_source_end(other_.clipped_source_end)
	{
SEQAN_CHECKPOINT
	}
	Gaps & operator = (Gaps const & other_)
	{
SEQAN_CHECKPOINT
		data_arr = other_.data_arr;
		data_end_position = other_.data_end_position;
		data_unclipped_end_position = other_.data_unclipped_end_position,
		setValue(data_source, source(other_));
		clipped_source_begin = other_.clipped_source_begin;
		clipped_source_end = other_.clipped_source_end; 
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

template <typename TSource>
inline typename Position<Gaps<TSource, ArrayGaps> >::Type
_getTrailingGaps(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_unclipped_end_position;
}

template <typename TSource, typename TSize>
inline void
_setTrailingGaps(Gaps<TSource, ArrayGaps> & me, TSize const & size)
{
	TSize zero = 0;
	if (size >= zero)
	{
		me.data_unclipped_end_position = size;
	}
	else
	{
		me.data_unclipped_end_position = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////
/**
 * Returns the length of the gaps sequence including the trailing gaps
 */
template <typename TSource>
inline typename Position<Gaps<TSource, ArrayGaps> >::Type
_unclippedLength(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return me.data_end_position + me.data_unclipped_end_position;
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

/**
.Function.clippedBeginPosition:
..summary:Begin position of the source segment.
..cat:Alignments
..signature:clippedBeginPosition(object)
..param.object:An object that has a source
...type:Class.Gaps
..returns:The position of the first item in @Function.source.source(object)@ that is used in object.
..see:Function.begin
..see:Function.source
..see:Function.clippedEndPosition
..include:seqan/align.h
*/
template <typename TSource>
inline typename Position<TSource>::Type
clippedBeginPosition(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.clipped_source_begin;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setClippedBeginPosition(Gaps<TSource, ArrayGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.clipped_source_begin = _pos;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.clippedEndPosition:
..summary:Position of the end of the source segment.
..cat:Alignments
..signature:clippedEndPosition(object)
..param.object:An object that has a source
...type:Class.Gaps
..returns:The position behind the last element of the source segment.
..see:Function.end
..see:Function.sourceEnd
..see:Function.clippedBeginPosition
..include:seqan/align.h
*/
template <typename TSource>
inline typename Position<TSource>::Type
clippedEndPosition(Gaps<TSource, ArrayGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.clipped_source_end;
}
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setClippedEndPosition(Gaps<TSource, ArrayGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.clipped_source_end = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSize>
inline void
_initToResize(Gaps<TSource, ArrayGaps> & me,
				TSize _size)
{
SEQAN_CHECKPOINT
	resize(_dataArr(me), 2);
	_dataArr(me)[0] = 0;
	_dataArr(me)[1] = _size;
	_setEndPosition(me, _size);
	_setTrailingGaps(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearGaps(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	_initToResize(me, clippedEndPosition(me) - clippedBeginPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clear(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	_initToResize(me, 0);
	_setClippedBeginPosition(me, 0);
	_setClippedEndPosition(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearClipping(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	_setEndPosition(me, length(source(me)));
	_setTrailingGaps(me, 0);
	_setClippedBeginPosition(me, 0);
	_setClippedEndPosition(me, length(source(me)));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.toViewPosition:
..summary:Transforms source to view position.
..cat:Alignments
..signature:toViewPosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the original sequence to get the view position for.
..returns:The position in the view/gaps position.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position and must be greater or equal the clipped begin position of $gap$.
..see:Function.toSourcePosition
..include:seqan/align.h
*/
template <typename TSource>
inline typename Position< Gaps<TSource, ArrayGaps> >::Type
toViewPosition(Gaps<TSource, ArrayGaps> const & gaps,
			   typename Position<TSource>::Type pos)
{
SEQAN_CHECKPOINT

	SEQAN_ASSERT_GEQ(pos, clippedBeginPosition(gaps));
	pos -= clippedBeginPosition(gaps);

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

/**
.Function.toSourcePosition:
..summary:Transforms view to source position, if the view position is a gap, the original position of the next non-gap entry is returned.
..cat:Alignments
..signature:toSourcePosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the view sequence (this includes gaps) to get the original position for.
..returns:The position in the source sequence.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position.
..see:Function.toViewPosition
..include:seqan/align.h
*/
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
	typename Position<TSource>::Type source_pos = clippedBeginPosition(gaps);
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
_iteratorGapsArray(TGaps & gaps, 
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
		typename Position<TSource>::Type source_pos = clippedBeginPosition(gaps);

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
	return _iteratorGapsArray(gaps, view_pos);
}
template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, ArrayGaps> const, Tag<TTag> const>::Type
iter(Gaps<TSource, ArrayGaps> const & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	return _iteratorGapsArray<Gaps<TSource, ArrayGaps> const>(gaps, view_pos);
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
		_initToResize(me, length(source(me)));
	}
	_dataArr(me)[0] = view_position;
	_setEndPosition(me, endPosition(me) + view_position - old_pos);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition>
inline void
setClippedBeginPosition(Gaps<TSource, ArrayGaps> & me,
					   TPosition source_position)
{
	SEQAN_ASSERT(length(_dataArr(me)))

	typedef Gaps<TSource, ArrayGaps> TGaps;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename TGaps::TArr TArr;

	TPosition old_clipped_begin_pos = clippedBeginPosition(me);
	if (old_clipped_begin_pos == source_position) return;
	else if (source_position < old_clipped_begin_pos)
	{
SEQAN_CHECKPOINT
		TViewPosition delta = old_clipped_begin_pos - source_position;
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
		_setClippedBeginPosition(me, source_position);
	}
	else //(source_position > old_clipped_begin_pos)
	{
SEQAN_CHECKPOINT
		TViewPosition view_pos = toViewPosition(me, source_position);
		TViewPosition source_pos_left = source_position - old_clipped_begin_pos;
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
				_setClippedBeginPosition(me, source_position);
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
setClippedEndPosition(Gaps<TSource, ArrayGaps> & me,
					 TPosition source_position)
{
	typedef Gaps<TSource, ArrayGaps> TGaps;
	typedef typename Position<TGaps>::Type TViewPosition;
	typedef typename TGaps::TArr TArr;

	TArr arr = _dataArr(me);
	SEQAN_ASSERT(length(arr));

	TPosition old_end_begin_pos = clippedEndPosition(me);
	if (old_end_begin_pos == source_position) return;
	else if (source_position < old_end_begin_pos)
	{
SEQAN_CHECKPOINT
		typename Iterator<TArr, Standard>::Type it_arr_begin = begin(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr_end = end(_dataArr(me));
		typename Iterator<TArr, Standard>::Type it_arr = it_arr_begin;
		TViewPosition end_pos = 0;
		TViewPosition chars_to_scan = source_position - clippedBeginPosition(me);

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
	_setClippedEndPosition(me, source_position);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Size<TSource>::Type
sourceLength(Gaps<TSource, ArrayGaps> & me)
{
SEQAN_CHECKPOINT
	return clippedEndPosition(me) - clippedBeginPosition(me);
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
		typename Position<TArr>::Type block_pos = me.data_sub;

		if (length(arr) == pos)
		{//counting trailing gaps
			if (block_pos <= _getTrailingGaps(container(me)))
			{//is in range of trailing gaps
SEQAN_CHECKPOINT
				return _getTrailingGaps(container(me)) - block_pos;
			}
			else
			{//no trailing gaps here
SEQAN_CHECKPOINT
				return 0;
			}
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

template <typename TGaps>
inline typename Size<TGaps>::Type
countCharacters(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	if (isGap(me))
	{// only count chars
SEQAN_CHECKPOINT
		return 0;
	}
	else
	{//count chars
SEQAN_CHECKPOINT
		typedef typename Size<TGaps>::Type TGapsSize;
		typedef String<TGapsSize> const TArr;

		TArr & arr = _dataArr(container(me));
		typename Position<TArr>::Type pos = me.data_block;
		typename Position<TArr>::Type block_pos = me.data_sub;

		return arr[pos] - block_pos;
	}
}

//____________________________________________________________________________

template <typename T>
inline void 
_goNextArrayGapsIterator(T const & me)
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
	_goNextArrayGapsIterator(me);
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	_goNextArrayGapsIterator(me);
}

//____________________________________________________________________________

template <typename T>
inline void 
_goPreviousArrayGapsIterator(T const & me)
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
	_goPreviousArrayGapsIterator(me);
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<ArrayGaps> > const & me)
{
	_goPreviousArrayGapsIterator(me);
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
	typename Position<TArr>::Type block_pos = me.data_sub;
	typename Iterator<TArr, Rooted>::Type it = iter(arr, pos);

	if (isGap(me))
	{
		if (pos < length(arr))
		{//here is a gap: expand it
SEQAN_CHECKPOINT
			value(it) += size;
		}
		else if (pos == length(arr))
		{//gap at end: add trailing gaps
			if (block_pos <= _getTrailingGaps(container(me)))
			{
SEQAN_CHECKPOINT
				_setTrailingGaps(container(me), _getTrailingGaps(container(me)) + size);
			}
			return; // do not adjust right position
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
		else if (me.data_block == length(arr))
		{
			if (me.data_sub <= _getTrailingGaps(container(me)))
			{//remove trailing gaps
				if (size > (_getTrailingGaps(container(me)) - me.data_sub))
				{//remove more than exists
SEQAN_CHECKPOINT
					_setTrailingGaps(container(me), me.data_sub);
				}
				else
				{//remove less or equal than exists
SEQAN_CHECKPOINT
					_setTrailingGaps(container(me), _getTrailingGaps(container(me)) - size);
				}
			}
		}//else: trailing gaps: infinite gaps here - nothing to do
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
