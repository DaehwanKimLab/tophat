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

#ifndef SEQAN_HEADER_GAPS_SUMLIST_H
#define SEQAN_HEADER_GAPS_SUMLIST_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tag

struct SumlistGaps;


//////////////////////////////////////////////////////////////////////////////
// Gaps - SumlistGaps Spec
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.SumlistGaps:
..cat:Alignments
..general:Class.Gaps
..summary:Stores gaps in a Sumlist.
..signature:Gaps<TSource, SumlistGaps>
..param.TSource:Type of the ungapped sequence.
...metafunction:Metafunction.Source
..include:seqan/align.h
*/

template <typename TSource>
class Gaps<TSource, SumlistGaps>
{
public:
	typedef typename Size<Gaps>::Type TSize;
	typedef SumList<2, TSize> TSumlist;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef typename Position<Gaps>::Type TViewPosition;

public:
	TSumlist data_sumlist; 

	Holder<TSource> data_source;
	TSourcePosition clipped_source_begin;
	TSourcePosition clipped_source_end;

public:
	Gaps():
		clipped_source_begin(0),
		clipped_source_end(0)
	{
SEQAN_CHECKPOINT
	}
	Gaps(TSize _size):
		clipped_source_begin(0),
		clipped_source_end(0)
	{
		_initToResize(*this, _size);
	}
	Gaps(TSource & source_):
		data_source(source_),
		clipped_source_begin(beginPosition(source_)),
		clipped_source_end(endPosition(source_))
	{
SEQAN_CHECKPOINT
		_initToResize(*this, length(source_));
	}

	template <typename TSource2>
	Gaps(TSource2 const & source_):
		clipped_source_begin(beginPosition(source_)),
		clipped_source_end(endPosition(source_))
	{
SEQAN_CHECKPOINT
		data_source = source_;
		_initToResize(*this, length(source_));
	}

	Gaps(Gaps const & other_):
		data_sumlist(other_.data_sumlist),
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
		data_sumlist = other_.data_sumlist;
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
inline Holder<TSource> &
_dataSource(Gaps<TSource, SumlistGaps> & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TSource>
inline Holder<TSource> const &
_dataSource(Gaps<TSource, SumlistGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position< Gaps<TSource, SumlistGaps> >::Type
beginPosition(Gaps<TSource, SumlistGaps> & gaps)
{
SEQAN_CHECKPOINT
	if (length(gaps.data_sumlist))
	{
		return getValue(begin(gaps.data_sumlist), 0) - getValue(begin(gaps.data_sumlist), 1);
	}
	return 0;
}
template <typename TSource>
inline typename Position< Gaps<TSource, SumlistGaps> const>::Type
beginPosition(Gaps<TSource, SumlistGaps> const & gaps)
{
SEQAN_CHECKPOINT
/*	
typedef Gaps<TSource, SumlistGaps> const TGaps;
typedef typename Size<TGaps>::Type TSize;
typedef SumList<2, TSize> const TSumlist;
typedef typename Iterator<TSumlist>::Type TSumlistIterator;

TSumlistIterator it = begin(gaps.data_sumlist);
*/
	if (length(gaps.data_sumlist))
	{
		return getValue(begin(gaps.data_sumlist), 0) - getValue(begin(gaps.data_sumlist), 1);
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<Gaps<TSource, SumlistGaps> >::Type
endPosition(Gaps<TSource, SumlistGaps> & me)
{
SEQAN_CHECKPOINT
	return getSum(me.data_sumlist, 0);
}
template <typename TSource>
inline typename Position<Gaps<TSource, SumlistGaps> >::Type
endPosition(Gaps<TSource, SumlistGaps> const & me)
{
SEQAN_CHECKPOINT
	return getSum(me.data_sumlist, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
clippedBeginPosition(Gaps<TSource, SumlistGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.clipped_source_begin;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setClippedBeginPosition(Gaps<TSource, SumlistGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.clipped_source_begin = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline typename Position<TSource>::Type
clippedEndPosition(Gaps<TSource, SumlistGaps> const & me)
{
SEQAN_CHECKPOINT
	return me.clipped_source_end;
}
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSourcePosition>
inline void
_setClippedEndPosition(Gaps<TSource, SumlistGaps> & me, TSourcePosition _pos)
{
SEQAN_CHECKPOINT
	me.clipped_source_end = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSize2>
inline void
_initToResize(Gaps<TSource, SumlistGaps> & me,
				TSize2 _size)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Size<TGaps>::Type TSize;

	clear(me.data_sumlist);
	TSize values[2];
	values[0] = _size;
	values[1] = _size;
	appendValues(me.data_sumlist, (TSize *) values);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearGaps(Gaps<TSource, SumlistGaps> & me)
{
SEQAN_CHECKPOINT
	_initToResize(me, clippedEndPosition(me) - clippedBeginPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clear(Gaps<TSource, SumlistGaps> & me)
{
SEQAN_CHECKPOINT
	_initToResize(me, 0);
	_setClippedBeginPosition(me, 0);
	_setClippedEndPosition(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
inline void
clearClipping(Gaps<TSource, SumlistGaps> & me)
{
SEQAN_CHECKPOINT
	_setClippedBeginPosition(me, 0);
	_setClippedEndPosition(me, length(source(me)));
}

//////////////////////////////////////////////////////////////////////////////

// transforms source- to view-position
template <typename TSource>
inline typename Position< Gaps<TSource, SumlistGaps> >::Type
toViewPosition(Gaps<TSource, SumlistGaps> const & gaps_,
			   typename Position<TSource>::Type pos)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef SumList<2, TSize> TSumlist;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;
	TGaps & gaps = const_cast<TGaps &>(gaps_);

	SEQAN_ASSERT(pos >= clippedBeginPosition(gaps))
	pos -= clippedBeginPosition(gaps);

	if (pos >= getSum(gaps.data_sumlist, 1))
	{
		return endPosition(gaps);
	}

	TSumlistIterator it(gaps.data_sumlist);
	searchSumList(it, pos, 1);
	return getSum(it, 0) + getValue(it, 0) - (getSum(it, 1) + getValue(it, 1)) + pos;
}

//____________________________________________________________________________
// transformates view- to source-position. (aufgerundet!)
template <typename TSource>
inline typename Position<TSource>::Type
toSourcePosition(Gaps<TSource, SumlistGaps> const & gaps,
				 typename Position< Gaps<TSource, SumlistGaps> >::Type pos)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> const TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef SumList<2, TSize> const TSumlist;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;
	TSumlistIterator it(gaps.data_sumlist);
	searchSumList(it, pos, 0);

	TSize dif = getValue(it, 0) - getValue(it, 1);

	if (pos - getSum(it, 0) > dif)
	{//view position within a non-gap
		return clippedBeginPosition(gaps) + getSum(it, 1) + (pos - getSum(it, 0) - dif);
	}
	else
	{//view position within a gap
		return clippedBeginPosition(gaps) + getSum(it, 1);
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
template <typename TSource, typename TTag>
inline typename Iterator<Gaps<TSource, SumlistGaps>, Tag<TTag> const>::Type
begin(Gaps<TSource, SumlistGaps> & gaps,
	  Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Iterator<TGaps, Tag<TTag> const>::Type TGapsIterator;
	return TGapsIterator(gaps);
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, SumlistGaps>, Tag<TTag> const>::Type
iter(Gaps<TSource, SumlistGaps> & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Iterator<TGaps, Tag<TTag> const>::Type TGapsIterator;
	return TGapsIterator(gaps, view_pos);
}
template <typename TSource, typename TPosition, typename TTag>
inline typename Iterator<Gaps<TSource, SumlistGaps> const, Tag<TTag> const>::Type
iter(Gaps<TSource, SumlistGaps> const & gaps,
	 TPosition view_pos,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Iterator<TGaps const, Tag<TTag> >::Type TGapsIterator;
	return TGapsIterator(const_cast<TGaps &>(gaps), view_pos);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


template <typename TSource, typename TPosition>
inline void
setBeginPosition(Gaps<TSource, SumlistGaps> & me,
				 TPosition view_position)
{
SEQAN_CHECKPOINT
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef SumList<2, TSize> TSumlist;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;
	TSumlistIterator it(me.data_sumlist);

	assignValue(it, 0, getValue(it, 1) + view_position);
}

//////////////////////////////////////////////////////////////////////////////
//Funktion gibt die Zahl, um die sich die view_positions verschieben zurueck.

template <typename TSource, typename TPosition>
inline typename Size< Gaps<TSource, SumlistGaps> >::Type
setClippedBeginPosition(Gaps<TSource, SumlistGaps> & me,
					   TPosition source_position)
{
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef SumList<2, TSize> TSumlist;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;
	typedef SumListValues<2, TSize> TValues;

	TSourcePosition old_begin_position = clippedBeginPosition(me);

	TSize ret_value = 0;

	if ((TSourcePosition) source_position < old_begin_position)
	{//expand to left
		TSumlistIterator it(me.data_sumlist);
		TSize new_block_size = getValue(it, 1) + (old_begin_position - source_position);
		//are there enough trailing blanks?
		if (getValue(it, 0) < new_block_size)
		{//not enough trailing blanks: move view positions!
			ret_value = new_block_size - getValue(it, 0); //the amount the view positions are moved
			assignValue(it, 0, new_block_size);
		}
		assignValue(it, 1, new_block_size);
	}
	else if ((TSourcePosition) source_position > old_begin_position)
	{//shrink left

		TSize pos_change = source_position - old_begin_position;

		//sum all values up to the new source_position
		TSumlistIterator it(me.data_sumlist);
		TValues vals;

		while (vals[1] + getValue(it, 1) < pos_change)
		{
			vals += getValues(it);
			removeValues(it);
		}
		//add the collected gaps to the first block
		assignValue(it, 0, getValue(it, 0) + vals[0] );
		assignValue(it, 1, getValue(it, 1) - (pos_change - vals[1]));
	}
	else
	{//nothing to do
		return 0;
	}

	_setClippedBeginPosition(me, source_position);

	return ret_value;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TPosition>
inline void
setClippedEndPosition(Gaps<TSource, SumlistGaps> & me,
					 TPosition source_position)
{
	typedef Gaps<TSource, SumlistGaps> TGaps;
	typedef typename Size<TGaps>::Type TSize;
	typedef SumList<2, TSize> TSumlist;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;
	typedef SumListValues<2, TSize> TValues;

	TSourcePosition old_end_position = clippedEndPosition(me);

	if (old_end_position < (TSourcePosition) source_position)
	{//expand end to right: just add the difference to the last block
		TSumlistIterator it(me.data_sumlist);
		goBeforeEnd(it);
		assignValue(it, 0, getValue(it, 0) + source_position - old_end_position);
		assignValue(it, 1, getValue(it, 1) + source_position - old_end_position);
	}
	else if (old_end_position > (TSourcePosition) source_position)
	{//shrink end to left
		TSumlistIterator it(me.data_sumlist);

		//ajust the new last block
		searchSumList(it, source_position - clippedBeginPosition(me), 1);
		TSize diff = getValue(it, 1) + getSum(it, 1) - source_position + clippedBeginPosition(me);
		assignValue(it, 0, getValue(it, 0) - diff);
		assignValue(it, 1, getValue(it, 1) - diff);

		//remove all follwing blocks
		goNext(it);
		while (!atEnd(it))
		{
			removeValues(it);
		}
	}
	else 
	{//nothing to do
		return;
	}

	_setClippedEndPosition(me, source_position);

}

//////////////////////////////////////////////////////////////////////////////


template <typename TGaps>
struct SumListOfGaps_
{
	typedef typename Size<TGaps>::Type TSize_;
	typedef SumList<2, TSize_> Type;
};

template <typename TGaps>
struct SumListOfGaps_<TGaps const>
{
	typedef typename Size<TGaps const>::Type TSize_;
	typedef SumList<2, TSize_> const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Gaps Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
class Iter<TGaps, GapsIterator<SumlistGaps> >
{
public:
	typedef typename Size<TGaps>::Type TSize;
	typedef typename SumListOfGaps_<TGaps>::Type TSumlist;
//	typedef SumList<2, TSize> TSumlist;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;

	TGaps * data_container;							//the gaps object
	mutable TSumlistIterator iter;			
	mutable TSize pos;

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
	Iter(TGaps & container_, TSize view_position):
		data_container(& container_),
		iter(container_.data_sumlist)
	{
SEQAN_CHECKPOINT
		searchSumList(iter, view_position, 0);
		pos = view_position - getSum(iter, 0);
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
struct Source< Iter<TGaps, GapsIterator<SumlistGaps> > >
{
	typedef typename Source<TGaps>::Type TGapsSource_;
	typedef typename Iterator<TGapsSource_>::Type Type;
};
template <typename TGaps>
struct Source< Iter<TGaps, GapsIterator<SumlistGaps> > const>
{
	typedef typename Source<TGaps>::Type TGapsSource_;
	typedef typename Iterator<TGapsSource_>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<SumlistGaps> > >::Type
getValue(Iter<TGaps, GapsIterator<SumlistGaps> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<SumlistGaps> > >::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(* me.data_container), sourcePosition(me));
}
template <typename TGaps>
inline typename GetValue< Iter<TGaps, GapsIterator<SumlistGaps> > const>::Type
getValue(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<SumlistGaps> > const>::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(* me.data_container), sourcePosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename Position<TGaps>::Type
viewPosition(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	return getSum(me.iter, 0) + me.pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename Position<TGaps>::Type
sourcePosition(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	if (isGap(me))
	{
		return getSum(me.iter, 1) + me.data_container->clipped_source_begin;
	}
	else
	{
		return getSum(me.iter, 1) + (me.pos - (getValue(me.iter, 0) - getValue(me.iter, 1))) + me.data_container->clipped_source_begin;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline typename Source<Iter<TGaps, GapsIterator<SumlistGaps> > >::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<SumlistGaps> > & me)
{
SEQAN_CHECKPOINT
	return begin(source(*me.data_container)) +  sourcePosition(me);
}
template <typename TGaps>
inline typename Source<Iter<TGaps, GapsIterator<SumlistGaps> > const>::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	return begin(source(*me.data_container)) + sourcePosition(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps>
inline bool 
isGap(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	return (atEnd(me.iter) || (me.pos < (getValue(me.iter, 0) - getValue(me.iter, 1))));
}

//____________________________________________________________________________

template <typename TGaps>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
	if (!atEnd(me) && isGap(me))
	{
		return getValue(me.iter, 0) - getValue(me.iter, 1) - me.pos;
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
_goNextSumlistGapsIterator(T & me)
{
	++me.pos;
	if (atEnd(me.iter)) return;

	if (getValue(me.iter, 0) > me.pos) return;

	me.pos = 0;
	goNext(me.iter);
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<SumlistGaps> > & me)
{
	_goNextSumlistGapsIterator(me);
}
template <typename TGaps>
inline void 
goNext(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
	_goNextSumlistGapsIterator(me);
}

//____________________________________________________________________________

template <typename T>
inline void 
_goPreviousSumlistGapsIterator(T & me)
{
	if (me.pos > 0)
	{
SEQAN_CHECKPOINT
		--me.pos;
	}
	else
	{
		searchSumList(me.iter, getSum(me.iter, 0) - 1, 0);
		me.pos = getValue(me.iter, 0) - 1;
	}
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<SumlistGaps> > & me)
{
	_goPreviousSumlistGapsIterator(me);
}
template <typename TGaps>
inline void 
goPrevious(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
	_goPreviousSumlistGapsIterator(me);
}
//____________________________________________________________________________

template <typename TGaps>
inline bool 
atBegin(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	return (atBegin(me.iter) && (me.pos == 0));
}

//____________________________________________________________________________

template <typename TGaps>
inline bool 
atEnd(Iter<TGaps, GapsIterator<SumlistGaps> > const & me)
{
SEQAN_CHECKPOINT
	return (atEnd(me.iter));
}

//____________________________________________________________________________

template <typename TGaps, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<SumlistGaps> > const & me,
		   TCount size)
{
	typedef typename Size<TGaps>::Type TSize;

	if (atEnd(me)) return;

	if (getValue(me.iter, 0) - getValue(me.iter, 1) >= me.pos)
	{//gap verlaengern
		assignValue(me.iter, 0, getValue(me.iter, 0) + size); 
	}
	else
	{
		//new left block: 
		TSize vals[2];
		vals[0] = me.pos;
		vals[1] = me.pos - (getValue(me.iter, 0) - getValue(me.iter, 1));

		//adjust the current block
		assignValue(me.iter, 0, getValue(me.iter, 0) - vals[0] + size);
		assignValue(me.iter, 1, getValue(me.iter, 1) - vals[1]);

		//insert the new block before the current block
		insertValues(me.iter, vals);

		//move iterator again to beginning of next block
		goNext(me.iter);
		me.pos = 0;
	}

}

//____________________________________________________________________________

//delete up to size gaps 
	
template <typename TGaps, typename TCount>
inline void
removeGaps(Iter<TGaps, GapsIterator<SumlistGaps> > const & me,
		   TCount _size)
{
	typedef typename Size<TGaps>::Type TSize;
	typedef typename SumListOfGaps_<TGaps>::Type TSumlist;
	typedef typename Iterator<TSumlist>::Type TSumlistIterator;

	TCount gaps_here = countGaps(me);
	if (gaps_here < _size) _size = gaps_here;
	if (_size)
	{
		TSize new_size = getValue(me.iter, 0) - _size;
		if (new_size > getValue(me.iter, 1))
		{//ok, there are still blanks at the beginning of this block. Just assign new block width
			assignValue(me.iter, 0, new_size);
			return;
		}
		//there are not enough blanks, i.e. the complete gap was deleted
		//is this the first block? if yes, nothing more to do
		if (me == begin(container(me))) return; 

		//this is not the first block, so it must be merged with the block before this block
		//find the preceding block
		TSumlistIterator it_before = me.iter;
		goPrevious(it_before); //this involves a search, alas!

		//adjust the pos-member of me:
		me.pos += getValue(it_before, 0);

		//add the current block size to it
		assignValue(it_before, 0, getValue(it_before, 0) + new_size);
		assignValue(it_before, 1, getValue(it_before, 1) + new_size);

		//remove the current block
		removeValues(me.iter);

		//adjust me (me.pos was already adjusted)
		me.iter = it_before;
	}
}

//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SumlistGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SumlistGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SumlistGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator == (Iter<TGaps1, GapsIterator<SumlistGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter == _right.iter) && (_left.pos == _right.pos);
}
//____________________________________________________________________________

template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SumlistGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SumlistGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SumlistGaps> > & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}
template <typename TGaps1, typename TGaps2>
inline bool
operator != (Iter<TGaps1, GapsIterator<SumlistGaps> > const & _left, 
			 Iter<TGaps2, GapsIterator<SumlistGaps> > const & _right)
{
SEQAN_CHECKPOINT
	return (_left.iter != _right.iter) || (_left.pos != _right.pos);
}


//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
#endif //#ifndef SEQAN_HEADER_...
