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

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

template <typename _TPos, typename _TGapAnchor, typename _TSpec = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type	TId;
	typedef _TPos										TPos;
	typedef _TGapAnchor									TGapAnchor;
	typedef _TSpec										TSpec;
	typedef String<TGapAnchor>							TGapAnchors;

	static const TId INVALID_ID;
	
	TId			id;
	TId			readId;
	TId			contigId;
	TId			pairMatchId;	// unique id. for multiple mate-pair matches (not matePairId)
	TPos		beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos		endPos;			// end position of ..., for reverse aligned reads holds end < begin
	TGapAnchors	gaps;

	AlignedReadStoreElement() : id(INVALID_ID), readId(INVALID_ID), contigId(INVALID_ID), pairMatchId(INVALID_ID), beginPos(0), endPos(0) {}
};

template <typename TScore, typename TSpec = void>
struct AlignQualityStoreElement
{
	TScore				pairScore;		// score of the mate-pair alignment (this read is part of)
	TScore				score;			// score of the single read alignment
	unsigned char		errors;			// absolute number of errors (Hamming or edit distance)
};


//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec>
const typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = SupremumValue<typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE;


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////


struct _SortContigId;
typedef Tag<_SortContigId> const SortContigId;


struct _SortId;
typedef Tag<_SortId> const SortId;

struct _SortBeginPos;
typedef Tag<_SortBeginPos> const SortBeginPos;

struct _SortEndPos;
typedef Tag<_SortEndPos> const SortEndPos;

struct _SortPairMatchId;
typedef Tag<_SortPairMatchId> const SortPairMatchId;

struct _SortReadId;
typedef Tag<_SortReadId> const SortReadId;


//////////////////////////////////////////////////////////////////////////////
// Sorting functors
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlignedRead, typename TTag>
struct _LessAlignedRead;

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return (a1.id) < (a2.id);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortContigId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.contigId < a2.contigId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortBeginPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _min(a1.beginPos, a1.endPos) < _min(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortEndPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _max(a1.beginPos, a1.endPos) < _max(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortPairMatchId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.pairMatchId < a2.pairMatchId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortReadId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.readId < a2.readId;
	}
};

//////////////////////////////////////////////////////////////////////////////
// Sorting function
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const &) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard()), 
		end(const_cast<TAlign&>(alignStore), Standard()), 
		less);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////


template <typename TGapAnchors = String<GapAnchor<unsigned> > >
struct AnchorGaps;

template <typename TSource, typename TGapAnchors>
class Gaps<TSource, AnchorGaps<TGapAnchors> >
{
public:
	typedef typename Value<TGapAnchors>::Type TGapAnchor;
	typedef typename Position<TGapAnchor>::Type TViewPosition;

public:
	Holder<TSource>		data_source;
	Holder<TGapAnchors>	data_gaps;
	int					data_viewCutBegin;	// how many alignment chars should be clipped at the beginning (can be negative too)
	int					data_viewCutEnd;	// how ...                                           end ...
	
	Gaps() {}

	Gaps(TSource &source): 
		data_source(source),
		data_viewCutBegin(0),
		data_viewCutEnd(0) {}

	Gaps(TSource &source, TGapAnchors &anchors):
		data_source(source),
		data_gaps(anchors),
		data_viewCutBegin(0),
		data_viewCutEnd(0) {}

	Gaps(TSource const &source, TGapAnchors &anchors):
		data_source(source),
		data_gaps(anchors),
		data_viewCutBegin(0),
		data_viewCutEnd(0) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline Holder<TSource> &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}
template <typename TSource, typename TGapAnchors>
inline Holder<TSource> const &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline TGapAnchors &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_gaps);
}

template <typename TSource, typename TGapAnchors>
inline TGapAnchors const &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_gaps);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAnchor, typename TSource, typename TGapAnchors, typename TIdx>
inline void
_getAnchor(TAnchor &anchor, Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TIdx idx)
{
SEQAN_CHECKPOINT
	if (idx > (TIdx)length(_dataAnchors(me)))
	{
		anchor.seqPos = length(source(me));
		if (empty(_dataAnchors(me)) && idx == 1)
			anchor.gapPos = anchor.seqPos;
		else
			anchor.gapPos = supremumValue(anchor.gapPos);
	}
	else if (idx > 0)
		anchor = _dataAnchors(me)[idx - 1];
	else
	{
		anchor.seqPos = 0;
		if (idx == 0)
			anchor.gapPos = 0;
		else
			anchor.gapPos = infimumValue(anchor.gapPos);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Size< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
_unclippedLength(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
	typedef typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type TSize;
	typedef typename Value<TGapAnchors>::Type TAnchor;
	TSize len = length(source(me));
	if (!empty(_dataAnchors(me)))
	{
		TAnchor const &last = back(_dataAnchors(me));
		len += last.gapPos - last.seqPos;
	}
	return len;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Size< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
length(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return _unclippedLength(me) - (me.data_viewCutBegin + me.data_viewCutEnd);
}

template <typename TSource, typename TGapAnchors>
inline typename Size< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
length(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
SEQAN_CHECKPOINT
	return _unclippedLength(me) - (me.data_viewCutBegin + me.data_viewCutEnd);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me);
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Standard)
{
SEQAN_CHECKPOINT
	return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Rooted)
{
SEQAN_CHECKPOINT
	return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me);
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Rooted)
{
SEQAN_CHECKPOINT
	return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Position< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
beginPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return me.data_viewCutBegin;
}

template <typename TSource, typename TGapAnchors>
inline typename Position< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
endPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
SEQAN_CHECKPOINT
	return _unclippedLength(me) - me.data_viewCutEnd;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors, typename TPosition>
inline void
setBeginPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, TPosition view_position)
{
SEQAN_CHECKPOINT
	me.data_viewCutBegin = view_position;
}

template <typename TSource, typename TGapAnchors, typename TPosition>
inline void
setEndPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, TPosition view_position)
{
SEQAN_CHECKPOINT
	me.data_viewCutEnd = _unclippedLength(me) - view_position;
}

//////////////////////////////////////////////////////////////////////////////

// to remove '... < 0 is always false' warning
template <typename T>
inline bool
_helperIsNegative(T, False)
{
	return false;
}
template <typename T>
inline bool
_helperIsNegative(T t, True)
{
	return t < 0;
}

template <typename TSource, typename TGapAnchors, typename TPosition>
inline TPosition
positionGapToSeq(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TPosition pos)
{
	typedef typename Position<typename Value<TGapAnchors>::Type >::Type TAnchorPos;	

	GapAnchor<int>	prevAnchor, nextAnchor;
	TPosition		seqPos;
	int				anchorIdx;

	if (_helperIsNegative(pos, typename TYPECMP<TPosition, typename _MakeSigned<TPosition>::Type>::Type()))
		anchorIdx = -1;
	else
	{
		TGapAnchors const & anchors = _dataAnchors(me);
		TAnchorPos seqLength = length(source(me));
		if (!empty(anchors))
		{
			anchorIdx = upperBoundGapAnchor(anchors, pos, SortGapPos()) - begin(anchors, Standard());
			if (anchorIdx < (int)length(anchors))
				if (anchors[anchorIdx].gapPos == (TAnchorPos)pos && anchors[anchorIdx].seqPos != seqLength)
					++anchorIdx;
		}
		else 
			anchorIdx = ((TAnchorPos)pos < seqLength)? 0: 1;
	}
	_getAnchor(prevAnchor, me, anchorIdx);
	_getAnchor(nextAnchor, me, anchorIdx + 1);

	if (nextAnchor.seqPos - prevAnchor.seqPos > (int)pos - prevAnchor.gapPos)
		seqPos = prevAnchor.seqPos + (pos - prevAnchor.gapPos);
	else
		seqPos = nextAnchor.seqPos;
	return seqPos;
}

template <typename TSource, typename TGapAnchors, typename TPosition>
inline TPosition
positionSeqToGap(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TPosition pos)
{
	typedef typename Position<typename Value<TGapAnchors>::Type >::Type TAnchorPos;	

	GapAnchor<int>	prevAnchor, nextAnchor;
	TPosition		gapPos;
	int				anchorIdx;

	if (_helperIsNegative(pos, typename TYPECMP<TPosition, typename _MakeSigned<TPosition>::Type>::Type()))
		anchorIdx = -1;
	else
	{
		TGapAnchors const & anchors = _dataAnchors(me);
		TAnchorPos seqLength = length(source(me));
		if (!empty(anchors))
		{
			anchorIdx = upperBoundGapAnchor(anchors, pos, SortSeqPos()) - begin(anchors, Standard());
			if (anchorIdx < (int)length(anchors))
				if (anchors[anchorIdx].seqPos == (TAnchorPos)pos)
					++anchorIdx;
		}
		else 
			anchorIdx = ((TAnchorPos)pos < seqLength)? 0: 1;
	}
	_getAnchor(prevAnchor, me, anchorIdx);
	_getAnchor(nextAnchor, me, anchorIdx + 1);

	if (nextAnchor.gapPos - prevAnchor.gapPos > (int)pos - prevAnchor.seqPos)
		gapPos = prevAnchor.gapPos + (pos - prevAnchor.seqPos);
	else
		gapPos = nextAnchor.gapPos;
	return gapPos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGaps_, typename TGapAnchors_> //Gaps<TSource_, AnchorGaps<TGapAnchors_> >
class Iter<TGaps_, GapsIterator<AnchorGaps<TGapAnchors_> > >
{
public:
	typedef TGaps_											TGaps;
	typedef typename Source<TGaps>::Type					TSource;
	typedef TGapAnchors_									TGapAnchors;

//	typedef typename Value<TGapAnchors>::Type				TGapAnchor;
	typedef GapAnchor<int>									TGapAnchor;
	typedef typename Size<TGapAnchor>::Type					TGapSize;
	typedef typename Iterator<TGapAnchors, Standard>::Type	TAnchorIter;

	TGaps *					data_container;							//the gaps object
	TGapSize				seqLength;
	mutable TGapAnchor		current;
	mutable TGapAnchor		prevAnchor;
	mutable TGapAnchor		nextAnchor;
	mutable TGapAnchor		viewBegin;
	mutable TGapAnchor		viewEnd;
	mutable int				anchorIdx;

public:
	Iter() 
	{
SEQAN_CHECKPOINT
		data_container = NULL;
		seqLength = 0;
	}
/*	Iter(Iter const & other_):
		data_container(other_.data_container),
		seqLength(other_.seqLength),
		current(other_.current),
		prevAnchor(other_.prevAnchor),
		nextAnchor(other_.nextAnchor),
		anchorIdx(other_.anchorIdx)
	{
SEQAN_CHECKPOINT
	}
*/	Iter(TGaps & container_):
		data_container(&container_),
		seqLength(length(source(container_)))
	{
SEQAN_CHECKPOINT
		_goTo_gapAnchorIterator(*this, data_container->data_viewCutBegin);
		viewBegin = current;
		viewEnd.gapPos = _unclippedLength(*data_container) - data_container->data_viewCutEnd;
		viewEnd.seqPos = positionGapToSeq(*data_container, viewEnd.gapPos);
	}
	Iter(TGaps & container_, TGapSize position):
		data_container(&container_),
		seqLength(length(source(container_)))
	{
SEQAN_CHECKPOINT
		_goTo_gapAnchorIterator(*this, position + data_container->data_viewCutBegin);
		viewBegin.gapPos = data_container->data_viewCutBegin;
		viewEnd.gapPos   = _unclippedLength(*data_container) - data_container->data_viewCutEnd;
		viewBegin.seqPos = positionGapToSeq(*data_container, viewBegin.gapPos);
		viewEnd.seqPos   = positionGapToSeq(*data_container, viewEnd.gapPos);
	}
	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & operator = (Iter const & other_)
	{
SEQAN_CHECKPOINT
		data_container = other_.data_container;
		seqLength = other_.seqLength;
		current = other_.current;
		prevAnchor = other_.prevAnchor;
		nextAnchor = other_.nextAnchor;
		anchorIdx = other_.anchorIdx;
		viewBegin = other_.viewBegin;
		viewEnd = other_.viewEnd;
		return *this;
	}
};

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline typename Source<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
source(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	return begin(source(*me.data_container), Rooted()) + me.current.seqPos;
}
template <typename TGaps, typename TGapAnchors>
inline typename Source<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
source(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	return begin(source(*me.data_container), Rooted()) + me.current.seqPos;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline typename GetValue< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
getValue(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}
template <typename TGaps, typename TGapAnchors>
inline typename GetValue< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
getValue(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else return getValue(source(me));
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline bool 
isGap(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	return me.current.seqPos == me.nextAnchor.seqPos;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
isClipped(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	return me.current.gapPos == me.nextAnchor.gapPos;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	typedef typename Size<TGaps>::Type TGapsSize;
	
SEQAN_CHECKPOINT
	if (isGap(me))
		return me.nextAnchor.gapPos - me.current.gapPos;
	else
		return 0;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline bool 
atBegin(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//	return me.current.seqPos == 0 && me.current.gapPos == 0;
	return me.current == me.viewBegin;
}
template <typename TGaps, typename TGapAnchors>
inline bool 
atBegin(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current.seqPos == 0 && me.current.gapPos == 0;
	return me.current == me.viewBegin;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//	return me.current == me.nextAnchor;
	return me.current == me.viewEnd;
}
template <typename TGaps, typename TGapAnchors>
inline bool 
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current == me.nextAnchor;
	return me.current == me.viewEnd;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline bool 
operator == (
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
	return left.current == right.current;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
operator != (
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
	return left.current != right.current;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
operator < (
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
	return left.current < right.current;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
operator > (
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
	return left.current > right.current;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me,
		   TCount size)
{
	TGapAnchors & anchors = _dataAnchors(*me.data_container);
	typedef typename Iterator<TGapAnchors, Standard>::Type TIter;
	typedef typename Value<TGapAnchors>::Type TGapAnchor;

	if (size <= 0) return;
	
	// insert a new anchor
	if (!isGap(me))
	{
		if (me.prevAnchor.gapPos == me.current.gapPos)
		{
			me.nextAnchor = me.prevAnchor;
			_getAnchor(me.prevAnchor, *me.data_container, --me.anchorIdx);
		}
		else
		{
			me.nextAnchor = me.current;
			insertValue(anchors, me.anchorIdx, me.nextAnchor, Generous());
		}
	}
	else
	{
		if (me.anchorIdx >= (int)length(anchors)) return;
		if (empty(anchors))
			appendValue(anchors, me.nextAnchor, Generous());
	}
	if (me.anchorIdx < (int)length(anchors))
	{
		if (me.anchorIdx >= 0)
			me.nextAnchor.gapPos += size;
		TIter it = begin(anchors, Standard());
		TIter itEnd = end(anchors, Standard());
		if (me.anchorIdx >= 0)
			it += me.anchorIdx;
		for (; it != itEnd; ++it)
			(*it).gapPos += size;
	}
	if (me.current.gapPos <= me.viewEnd.gapPos)
		me.viewEnd.gapPos += size;
/*
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > it2 = begin(*me.data_container) + me.current.gapPos;
	if (me.current != it2.current || me.prevAnchor != it2.prevAnchor || me.nextAnchor != it2.nextAnchor || me.anchorIdx != it2.anchorIdx)
		std::cout<<"*";
*/
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors, typename TCount>
inline void
removeGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me,
		   TCount size)
{
	TGapAnchors & anchors = _dataAnchors(*me.data_container);
	typedef typename Iterator<TGapAnchors, Standard>::Type TIter;
	typedef typename Value<TGapAnchors>::Type TGapAnchor;

	if (size <= 0 || !isGap(me)) return;
	
	if (me.current.gapPos + size > me.nextAnchor.gapPos)
		size = me.nextAnchor.gapPos - me.current.gapPos;
	
	if (me.prevAnchor.gapPos + me.current.seqPos == me.current.gapPos + me.prevAnchor.seqPos &&
		me.current.gapPos + size == me.nextAnchor.gapPos)
	{	
		// remove the gap
		if (me.anchorIdx < (int)length(anchors))
			erase(anchors, me.anchorIdx);
		_getAnchor(me.nextAnchor, *me.data_container, me.anchorIdx);
	}

	// shift anchors
	if (me.anchorIdx < (int)length(anchors))
	{
		me.nextAnchor.gapPos -= size;	
		if (me.anchorIdx >= 0)
			me.nextAnchor.gapPos += size;
		TIter it = begin(anchors, Standard());
		TIter itEnd = end(anchors, Standard());
		if (me.anchorIdx >= 0)
			it += me.anchorIdx;
		for (; it != itEnd; ++it)
			(*it).gapPos -= size;
	}
	if (me.current.gapPos <= me.viewEnd.gapPos)
		me.viewEnd.gapPos -= size;
/*
	Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > it2 = begin(*me.data_container) + me.current.gapPos;
	if (me.current != it2.current || me.prevAnchor != it2.prevAnchor || me.nextAnchor != it2.nextAnchor || me.anchorIdx != it2.anchorIdx)
		std::cout<<"*";
*/
}

//____________________________________________________________________________

template <typename T>
inline void 
_goNext_gapAnchorIterator(T & me)
{
	if (me.current.gapPos < me.nextAnchor.gapPos)
	{
		++me.current.gapPos;
		if (me.current.seqPos < me.nextAnchor.seqPos)
			++me.current.seqPos;
	}
	while (me.current.gapPos == me.nextAnchor.gapPos)
	{
		me.current = me.prevAnchor = me.nextAnchor;
		_getAnchor(me.nextAnchor, *me.data_container, ++me.anchorIdx + 1);
	}
}

template <typename T>
inline void 
_goPrevious_gapAnchorIterator(T & me)
{	
	while (me.current.gapPos == me.prevAnchor.gapPos)
	{
		me.current = me.nextAnchor = me.prevAnchor;
		_getAnchor(me.prevAnchor, *me.data_container, --me.anchorIdx);
	}
	--me.current.gapPos;
	if (me.nextAnchor.seqPos - me.prevAnchor.seqPos > me.current.gapPos - me.prevAnchor.gapPos)
		me.current.seqPos = me.prevAnchor.seqPos + (me.current.gapPos - me.prevAnchor.gapPos);
	else
		me.current.seqPos = me.nextAnchor.seqPos;
}


template <typename T, typename TPos>
inline void 
_goTo_gapAnchorIterator(T & me, TPos pos)
{
	typedef typename T::TGapAnchors	TGapAnchors;
	typedef typename T::TGapAnchor	TGapAnchor;
	typedef typename Position<typename Value<TGapAnchors>::Type >::Type TAnchorPos;

	if (_helperIsNegative(pos, typename TYPECMP<TPos, typename _MakeSigned<TPos>::Type>::Type()))
		me.anchorIdx = -1;
	else
	{
		TGapAnchors const & anchors = _dataAnchors(*me.data_container);
		if (!empty(anchors))
		{
			me.anchorIdx = upperBoundGapAnchor(anchors, pos, SortGapPos()) - begin(anchors, Standard());
			if (me.anchorIdx < (int)length(anchors))
				if (anchors[me.anchorIdx].gapPos == (TAnchorPos)pos && anchors[me.anchorIdx].seqPos != (TAnchorPos)me.seqLength)
					++me.anchorIdx;
		}
		else 
			me.anchorIdx = (pos < (TPos)me.seqLength)? 0: 1;
	}
	_getAnchor(me.prevAnchor, *me.data_container, me.anchorIdx);
	_getAnchor(me.nextAnchor, *me.data_container, me.anchorIdx + 1);

	me.current.gapPos = pos;
	if (me.nextAnchor.seqPos - me.prevAnchor.seqPos > (int)pos - me.prevAnchor.gapPos)
		me.current.seqPos = me.prevAnchor.seqPos + ((int)pos - me.prevAnchor.gapPos);
	else
		me.current.seqPos = me.nextAnchor.seqPos;
}

//____________________________________________________________________________

template <typename TGaps, typename TGapAnchors>
inline void
goNext(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	_goNext_gapAnchorIterator(me);
}

template <typename TGaps, typename TGapAnchors>
inline void
goPrevious(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	_goPrevious_gapAnchorIterator(me);
}

template <typename TGaps, typename TGapAnchors, typename TSize>
inline void
goFurther(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me, TSize steps)
{
	_goTo_gapAnchorIterator(me, me.current.gapPos + steps);
}

//____________________________________________________________________________

template <
	typename TCigar,
	typename TGaps1,
	typename TGaps2,
	typename TThresh>
inline void
getCigarString(
	TCigar &cigar,
	TGaps1 &gaps1,
	TGaps2 &gaps2,
	TThresh splicedGapThresh)
{
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	clear(cigar);
	char op, lastOp = ' ';
	unsigned numOps = 0;

//	std::cout << gaps1 << std::endl;
//	std::cout << gaps2 << std::endl;
	for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
	{
		if (isGap(it1))
		{
			if (isGap(it2))
				op = 'P';
			else if (isClipped(it2))
				op = '?';
			else
				op = 'I';
		} 
		else if (isClipped(it1))
		{
			op = '?';
		}
		else 
		{
			if (isGap(it2))
				op = 'D';
			else if (isClipped(it2))
				op = 'S';
			else
				op = 'M';
		}
		if (lastOp != op)
		{
			if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
				lastOp = 'N';
			if (numOps > 0)
			{
				std::stringstream num;
				num << numOps;
				append(cigar, num.str());
				appendValue(cigar, lastOp);
			}
			numOps = 0;
			lastOp = op;
		}
		++numOps;
	}
	if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
		lastOp = 'N';
	if (numOps > 0)
	{
		std::stringstream num;
		num << numOps;
		append(cigar, num.str());
		appendValue(cigar, lastOp);
	}
}

template <
	typename TCigar,
	typename TGaps1,
	typename TGaps2>
inline void
getCigarString(
	TCigar &cigar,
	TGaps1 &gaps1,
	TGaps2 &gaps2)
{
	return getCigarString(cigar, gaps1, gaps2, 20);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
