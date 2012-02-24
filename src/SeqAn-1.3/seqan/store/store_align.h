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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.AlignedReadStoreElement
..summary:Represents an alignment between read and contig.
..cat:Fragment Store
..signature:AlignedReadStoreElement<>
..signature:AlignedReadStoreElement<TPos[, TGapAnchor[, TSpec]]>
..param.TPos:Type to store gap-space positions.
..param.TGapAnchor:Type of a read gap anchor.
...type:Class.GapAnchor
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#alignedReadStore@ string.

.Typedef.AlignedReadStoreElement#TPos
..summary:Type of the $beginPos$ and $endPos$ members.
..class:Class.AlignedReadStoreElement
.Typedef.AlignedReadStoreElement#TGapAnchors
..summary:Type of the $gaps$ member.
..class:Class.AlignedReadStoreElement
.Typedef.AlignedReadStoreElement#TSpec
..summary:The specialization type.
..class:Class.AlignedReadStoreElement

.Memfunc.AlignedReadStoreElement#AlignedReadStoreElement
..summary:Constructor
..signature:AlignedReadStoreElement<> ()
..signature:AlignedReadStoreElement<TPos[, TGapAnchor[, TSpec]]> ()
..signature:AlignedReadStoreElement<TPos[, TGapAnchor[, TSpec]]> (id, readId, contigId, beginPos, endPos[, gaps])
..param.id:The alignment id refers to associated alignment information in @Memvar.FragmentStore#alignQualityStore@ or @Memvar.FragmentStore#alignedReadTagStore@.
..param.readId:Refers to the aligned read in the @Memvar.FragmentStore#readStore@.
..param.contigId:Refers to the contig in the @Memvar.FragmentStore#contigStore@ the read is aligned with.
..param.beginPos:Begin position of the alignment in gap-space.
..param.endPos:End position of the alignment in gap-space.
..param.gaps:Read gap anchors.
..remarks:The default constructor sets all ids to $INVALID_ID$ and $beginPos$ and $endPos$ to $0$.

..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#id
..summary:The alignment id refers to associated alignment information in @Memvar.FragmentStore#alignQualityStore@ or @Memvar.FragmentStore#alignedReadTagStore@.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#readId
..summary:Refers to the aligned read in the @Memvar.FragmentStore#readStore@.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#contigId
..summary:Refers to the contig in the @Memvar.FragmentStore#contigStore@ the read is aligned with.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#pairMatchId
..summary:Two read alignments having the same $pairMatchId$ form a valid pair match.
If $INVALID_ID$ the read is either not paired or could not be aligned as part of a pair match.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#beginPos
..summary:Begin position of the alignment in gap-space.
..type:Typedef.AlignedReadStoreElement#TPos
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#endPos
..summary:End position of the alignment in gap-space.
..type:Typedef.AlignedReadStoreElement#TPos
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#gaps
..summary:String of read gap anchors. Can be used to create a $Spec.AnchorGaps$ alignment row.
..type:Typedef.AlignedReadStoreElement#TGapAnchors
..class:Class.AlignedReadStoreElement
.Memvar.AlignedReadStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.AlignedReadStoreElement
..include:seqan/store.h
*/

template <typename TPos_, typename TGapAnchor_, typename TSpec_ = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type	TId;
	typedef TPos_										TPos;
	typedef TGapAnchor_									TGapAnchor;
	typedef TSpec_										TSpec;
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

	AlignedReadStoreElement(TId _id, TId _readId, TId _contigId, TPos _beginPos, TPos _endPos) :
		id(_id), 
		readId(_readId), 
		contigId(_contigId), 
		pairMatchId(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos) {}

	AlignedReadStoreElement(TId _id, TId _readId, TId _contigId, TPos _beginPos, TPos _endPos, TGapAnchors const &_gaps) :
		id(_id), 
		readId(_readId), 
		contigId(_contigId), 
		pairMatchId(INVALID_ID), 
		beginPos(_beginPos), 
		endPos(_endPos),
		gaps(_gaps) {}

    inline bool operator==(AlignedReadStoreElement const & other) const
    {
        return id == other.id &&
                readId == other.readId &&
                contigId == other.contigId &&
                pairMatchId == other.pairMatchId &&
                beginPos == other.beginPos &&
                endPos == other.endPos &&
                gaps == other.gaps;
    }
};


// TODO(holtgrew): I find this useful for debugging purposes. Keep it?
template <typename TStream, typename TPos, typename TGapAnchor, typename TSpec>
TStream & operator<<(TStream & stream, AlignedReadStoreElement<TPos, TGapAnchor, TSpec> const & e) {
    return stream << "AlignedReadStore(id=" << e.id << ", readId=" << e.readId << ", contigId=" << e.contigId << ", pairMatchId=" << e.pairMatchId << ", beginPos=" << e.beginPos << ", endPos=" << e.endPos << ", {gaps})";
}


//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec> 
const typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type 
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = MaxValue<typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

/**
.Class.AlignQualityStoreElement
..summary:Stores alignment qualities.
..cat:Fragment Store
..signature:AlignQualityStoreElement<TScore[, TSpec]>
..param.TScore:Type to store align and pair score values.
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#alignQualityStore@ string.

.Memfunc.AlignQualityStoreElement#AlignQualityStoreElement
..summary:Constructor
..signature:AlignQualityStoreElement<TScore[, TSpec]> ()
..remarks:Sets all members to $0$.

..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#pairScore
..summary:Combined score of both alignments of a pair match.
..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#score
..summary:Score of the alignment.
..class:Class.AlignQualityStoreElement
.Memvar.AlignQualityStoreElement#errors
..summary:Absolute number of errors in the alignment.
..type:nolink:unsigned char
..class:Class.AlignQualityStoreElement
..include:seqan/store.h
*/

template <typename TScore, typename TSpec = void>
struct AlignQualityStoreElement
{
	TScore				pairScore;		// score of the mate-pair alignment (this read is part of)
	TScore				score;			// score of the single read alignment
	unsigned char		errors;			// absolute number of errors (Hamming or edit distance)
	
	AlignQualityStoreElement():
		pairScore(0),
		score(0),
		errors(0) {}

    inline bool operator==(AlignQualityStoreElement const & other)
    {
        return pairScore == other.pairScore &&
                score == other.score &&
                errors == other.errors;
    }
};


//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.sortAlignedRead Tags
..summary:Tag to select a specific field to stably sort the @Memvar.FragmentStore#alignedReadStore@ by.
..cat:Fragment Store
..see:Function.sortAlignedReads
..see:Function.lowerBoundAlignedReads
..see:Function.upperBoundAlignedReads
..tag.SortContigId:
...summary:Sort alignedReads by $contigId$.
...signature:SortContigId
..tag.SortId:
...summary:Sort alignedReads by $id$.
...signature:SortId
..tag.SortBeginPos:
...summary:Sort alignedReads by $beginPos$.
...signature:SortBeginPos
..tag.SortEndPos:
...summary:Sort alignedReads by $endPos$.
...signature:SortEndPos
..tag.SortPairMatchId:
...summary:Sort alignedReads by $pairMatchId$.
...signature:SortPairMatchId
..tag.SortReadId:
...summary:Sort alignedReads by $readId$.
...signature:SortReadId
..include:seqan/store.h
*/

struct SortContigId_;
typedef Tag<SortContigId_> const SortContigId;

struct SortId_;
typedef Tag<SortId_> const SortId;

struct SortBeginPos_;
typedef Tag<SortBeginPos_> const SortBeginPos;

struct SortEndPos_;
typedef Tag<SortEndPos_> const SortEndPos;

struct SortPairMatchId_;
typedef Tag<SortPairMatchId_> const SortPairMatchId;

struct SortReadId_;
typedef Tag<SortReadId_> const SortReadId;


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

// TODO(holtgrew): Maybe add equalRangeAlignedReads?

/**
.Function.sortAlignedReads
..summary:Stably sort aligned reads.
..cat:Fragment Store
..signature:sortAlignedReads(alignStore, sortTag)
..signature:sortAlignedReads(alignStore, lessFunctor)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be sorted, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.sortTag:Selects the field to sort by.
...type:Tag.sortAlignedRead Tags
..param.lessFunctor:STL-less functor to compare two @Class.AlignedReadStoreElement.AlignedReadStoreElements@.
..remarks:This function calls $std::stable_sort$ to sort $alignStore$.
..include:seqan/store.h
..see:Function.lowerBoundAlignedReads
..see:Function.upperBoundAlignedReads

.Function.lowerBoundAlignedReads
..summary:Performs a binary lower bound search on the aligned reads.
..cat:Fragment Store
..signature:lowerBoundAlignedReads(alignStore, value, sortTag)
..signature:lowerBoundAlignedReads(itBegin, itEnd, value, sortTag)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be searched through, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.itBegin:An iterator to the first element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.itEnd:An iterator behind the last element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.value:The value to use for the comparison.
..param.sortTag:Selects the field for the comparison in the binary search.
...type:Tag.sortAlignedRead Tags
..remarks:This is equivalent to calling $std::lower_bound$ on $alignStore$ with according parameters.
..include:seqan/store.h
..see:Function.sortAlignedReads
..see:Function.upperBoundAlignedReads

.Function.upperBoundAlignedReads
..summary:Performs a binary upper bound search on the aligned reads.
..cat:Fragment Store
..signature:upperBoundAlignedReads(alignStore, value, sortTag)
..signature:upperBoundAlignedReads(itBegin, itEnd, value, sortTag)
..param.alignStore:A sequence of @Class.AlignedReadStoreElement@ to be searched through, e.g. @Memvar.FragmentStore#alignedReadStore@.
..param.itBegin:An iterator to the first element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.itEnd:An iterator behind the last element of the sequence of @Class.AlignedReadStoreElement@ to be searched through.
..param.value:The value to use for the comparison.
..param.sortTag:Selects the field for the comparison in the binary search.
...type:Tag.sortAlignedRead Tags
..remarks:This is equivalent to calling $std::upper_bound$ on $alignStore$ with according parameters.
..include:seqan/store.h
..see:Function.sortAlignedReads
..see:Function.lowerBoundAlignedReads
*/

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
	el.endPos = val;
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
	el.endPos = val;
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
	el.beginPos = val;
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
	el.beginPos = val;
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

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
lowerBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename T, typename TSpec, typename TSearchValue>
inline Iter<T, TSpec>
upperBoundAlignedReads(Iter<T, TSpec> const & alignedReadsItBegin,
                       Iter<T, TSpec> const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<Iter<T, TSpec> >::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<Iter<T, TSpec> >::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	el.endPos = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
lowerBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TReadAlignElement, typename TSearchValue>
inline TReadAlignElement *
upperBoundAlignedReads(TReadAlignElement * const & alignedReadsItBegin,
                       TReadAlignElement * const & alignedReadsItEnd,
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TReadAlignElement *>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		alignedReadsItBegin, 
		alignedReadsItEnd, 
		el,
		_LessAlignedRead<typename Value<TReadAlignElement *>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AnchorGaps:
..cat:Alignments
..general:Class.Gaps
..summary:Stores gaps anchors of the first characters behind gaps.
..signature:Gaps<TSource, AnchorGaps<TGapAnchors> >
..param.TSource:Type of the ungapped sequence.
...metafunction:Metafunction.Source
..param.TGapAnchors:Type of the sequence of gap anchors, e.g. a string of $Class.GapAnchor$.
..include:seqan/store.h

.Memfunc.Gaps#Gaps
..summary:Constructor
..signature:Gaps<TSource, AnchorGaps<TGapAnchors> > ()
..signature:Gaps<TSource, AnchorGaps<TGapAnchors> > (source[, anchors])
..signature:Gaps<TSource, AnchorGaps<TGapAnchors> > (anchors)
..param.source:The underlying ungapped sequence.
..param.anchors:The sequence of gap anchors, e.g. the $gaps$ members in $Class.ReadStoreElement$ or $Class.ContigStoreElement$.
*/

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

	Gaps(TGapAnchors &anchors):
		data_gaps(anchors),
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

	Gaps(TSource &source, TGapAnchors const &anchors):
		data_source(source),
		data_gaps(anchors),
		data_viewCutBegin(0),
		data_viewCutEnd(0) {}

	Gaps(TSource const &source, TGapAnchors const &anchors):
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

template <typename TSize, typename TSource, typename TGapAnchors>
inline void
_assignSourceLength(TSize & size, Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
SEQAN_CHECKPOINT
	if (IsSameType<TSource, Nothing>::VALUE)
		size = maxValue<TSize>();
	else
		size = length(value(me.data_source));
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
    SEQAN_CHECKPOINT;

	typedef typename Value<TGapAnchors>::Type TGapAnchor;
	typedef typename Position<TGapAnchor>::Type TPos;

	if (idx > (TIdx)length(_dataAnchors(me)))
	{
		_assignSourceLength(anchor.seqPos, me);
		if (empty(_dataAnchors(me)) && idx == 1)
			anchor.gapPos = anchor.seqPos;
		else
		{
			// for the sick case that an anchor seq position is beyond the sequence end
			if (!empty(_dataAnchors(me)))
			{
				// if there is no sequence but anchors -> assume infinite sequence
				if (anchor.seqPos == 0)
					anchor.seqPos = maxValue(anchor.gapPos);
				// if the sequence has a length > 0, but there is an anchor behind the end
				// -> elongate sequence
				else if ((__int64)anchor.seqPos < (__int64)back(_dataAnchors(me)).seqPos)
					anchor.seqPos = back(_dataAnchors(me)).seqPos;
			}
			anchor.gapPos = maxValue(anchor.gapPos);
		}
	}
	else if (idx > 0)
		anchor = _dataAnchors(me)[idx - 1];
	else
	{
		anchor.seqPos = 0;
		if (idx == 0)
			anchor.gapPos = 0;
		else
			anchor.gapPos = minValue(anchor.gapPos);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TGapAnchors>
inline typename Size< Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
_unclippedLength(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
	typedef typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type TSize;
	typedef typename Value<TGapAnchors>::Type TAnchor;
	TSize len;
	_assignSourceLength(len, me);
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

	GapAnchor<__int64> prevAnchor, nextAnchor;
	TPosition			seqPos;
	int					anchorIdx;

	if (_helperIsNegative(pos, typename IsSameType<TPosition, typename MakeSigned_<TPosition>::Type>::Type()))
		anchorIdx = -1;
	else
	{
		TGapAnchors const & anchors = _dataAnchors(me);
		TAnchorPos seqLength;
		_assignSourceLength(seqLength, me);
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

	GapAnchor<__int64>	prevAnchor, nextAnchor;
	TPosition			gapPos;
	int					anchorIdx;

	if (_helperIsNegative(pos, typename IsSameType<TPosition, typename MakeSigned_<TPosition>::Type>::Type()))
		anchorIdx = -1;
	else
	{
		TGapAnchors const & anchors = _dataAnchors(me);
		TAnchorPos seqLength;
		_assignSourceLength(seqLength, me);
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
		data_container(&container_)
	{
SEQAN_CHECKPOINT
		_assignSourceLength(seqLength, container_);
		_goToGapAnchorIterator(*this, data_container->data_viewCutBegin);
		viewBegin = current;
		viewEnd.gapPos = _unclippedLength(*data_container) - data_container->data_viewCutEnd;
		viewEnd.seqPos = positionGapToSeq(*data_container, viewEnd.gapPos);
	}
	Iter(TGaps & container_, TGapSize position):
		data_container(&container_)
	{
SEQAN_CHECKPOINT
		_assignSourceLength(seqLength, container_);
		_goToGapAnchorIterator(*this, position + data_container->data_viewCutBegin);
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
	else if (isUnknown(me)) return unknownValue<TValue>();
	else return getValue(source(me));
}
template <typename TGaps, typename TGapAnchors>
inline typename GetValue< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
getValue(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
SEQAN_CHECKPOINT
	typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type TValue;
	if (isGap(me)) return gapValue<TValue>();
	else if (isUnknown(me)) return unknownValue<TValue>();
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
isUnknown(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
	int len;
	_assignSourceLength(len, *me.data_container);
	return me.current.seqPos < 0 || me.current.seqPos >= len;
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
	return me.current <= me.viewBegin;
}
template <typename TGaps, typename TGapAnchors>
inline bool 
atBegin(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current.seqPos == 0 && me.current.gapPos == 0;
	return me.current <= me.viewBegin;
}

template <typename TGaps, typename TGapAnchors>
inline bool 
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//	return me.current == me.nextAnchor;
	return me.current >= me.viewEnd;
}
template <typename TGaps, typename TGapAnchors>
inline bool 
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//	return me.current == me.nextAnchor;
	return me.current >= me.viewEnd;
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
_goNextGapAnchorIterator(T & me)
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
_goPreviousGapAnchorIterator(T & me)
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
_goToGapAnchorIterator(T & me, TPos pos)
{
	typedef typename T::TGapAnchors	TGapAnchors;
	typedef typename T::TGapAnchor	TGapAnchor;
	typedef typename Position<typename Value<TGapAnchors>::Type >::Type TAnchorPos;

	if (_helperIsNegative(pos, typename IsSameType<TPos, typename MakeSigned_<TPos>::Type>::Type()))
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
	_goNextGapAnchorIterator(me);
}

template <typename TGaps, typename TGapAnchors>
inline void
goPrevious(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
	_goPreviousGapAnchorIterator(me);
}

template <typename TGaps, typename TGapAnchors, typename TSize>
inline void
goFurther(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me, TSize steps)
{
	_goToGapAnchorIterator(me, me.current.gapPos + steps);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
