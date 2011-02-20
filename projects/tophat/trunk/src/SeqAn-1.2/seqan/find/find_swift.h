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

#ifndef SEQAN_HEADER_FIND_SWIFT_H
#define SEQAN_HEADER_FIND_SWIFT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SWIFT to search a text for
// - semi-global alignments of one/multiple short sequences
// - local epsilon matches of one/multiple short sequences
//////////////////////////////////////////////////////////////////////////////

template < typename TObject, typename TSpec > 
class Index;

struct _SwiftLocal;
typedef Tag<_SwiftLocal> SwiftLocal;

struct _SwiftSemiGlobal;
typedef Tag<_SwiftSemiGlobal> SwiftSemiGlobal;

struct _SwiftSemiGlobalHamming;
typedef Tag<_SwiftSemiGlobalHamming> SwiftSemiGlobalHamming;


template <typename TSpec = SwiftSemiGlobal>
struct Swift {
	enum { SEMIGLOBAL = 1 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };		// q-gram must match exactly
	enum { HAMMING_ONLY = 0 };		// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 1 };	// params are determined only by seq.length
};

template <>
struct Swift<SwiftSemiGlobalHamming> {
	enum { SEMIGLOBAL = 1 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };		// q-gram must match exactly
	enum { HAMMING_ONLY = 1 };		// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 1 };	// params are determined only by seq.length
};

template <>
struct Swift<SwiftLocal> {
	enum { SEMIGLOBAL = 0 };	// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };		// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };	// allow 0 errors per q-gram
	enum { HAMMING_ONLY = 0 };	// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 0 };
};

struct SwiftParameters {
	int minThreshold;
	int minLog2Delta;
	int tabooLength;

	SwiftParameters():
		minThreshold(1),		// set minimal threshold to 1
		minLog2Delta(4),		// set minimal delta to 16
		tabooLength(1) {}		// minimal genomic distance between q-gram hits
};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSpec, typename _TSize, typename _TShortSize = unsigned short>
	struct _SwiftBucket 
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize					firstIncrement;
		TSize					lastIncrement;
		TShortSize				counter;		// q-gram hits
		TShortSize				threshold;		// at least threshold q-gram hits induce an approx match
#ifdef SEQAN_DEBUG_SWIFT
		TSize					_lastIncDiag;
#endif
	};

	template <typename _TSize, typename _TShortSize>
	struct _SwiftBucket<SwiftSemiGlobal, _TSize, _TShortSize> 
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize					lastIncrement;
		TShortSize				counter;		// q-gram hits
		TShortSize				threshold;		// at least threshold q-gram hits induce an approx match
#ifdef SEQAN_DEBUG_SWIFT
		int						_lastIncDiag;
#endif
	};

	template <typename TSpec, typename _TSize, typename _TShortSize = unsigned short>
	struct _SwiftBucketParams 
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize			firstBucket;	// first _SwiftBucket entry in pattern.buckets
		TSize			reuseMask;		// 2^ceil(log2(x)) reuse every x-th bucket)
		TShortSize		threshold;		// at least threshold q-gram hits induce an approx match
//		TShortSize		distanceCut;	// if lastIncrement is this far or farer away, threshold can't be reached
		TShortSize		delta;			// buckets begin at multiples of delta
		TShortSize		overlap;		// number of diagonals/columns a bucket shares with its neighbor
		TShortSize		tabooLength;	// minimal genomic distance between q-gram hits
		unsigned char	logDelta;		// log2(delta)
	};

//____________________________________________________________________________


	template <typename THstkPos>
	struct _SwiftHit 
	{
		THstkPos	hstkPos;			// parallelogram begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
		unsigned	bucketWidth;		// (non-diagonal) bucket width (bktHeight + delta + overlap (for diagonals))
	};

//____________________________________________________________________________


	template < typename THaystack, typename TSpec >
	class Finder< THaystack, Swift<TSpec> >
	{
	public:
		typedef typename Iterator<THaystack, Rooted>::Type			TIterator;
		typedef typename Position<THaystack>::Type					THstkPos;
		typedef _SwiftHit<__int64>									TSwiftHit;
		typedef String<TSwiftHit>									THitString;
		typedef typename Iterator<THitString, Standard>::Type		THitIterator;
		typedef typename SAValue<THaystack>::Type					TSAValue;
		typedef Repeat<TSAValue, unsigned>							TRepeat;
		typedef String<TRepeat>										TRepeatString;
		typedef typename Iterator<TRepeatString, Standard>::Type	TRepeatIterator;

		TIterator		data_iterator;
		TIterator		haystackEnd;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		startPos, curPos, endPos;
		THstkPos		dotPos, dotPos2;
		TRepeatString	data_repeats;
		TRepeatIterator	curRepeat, endRepeat;

		Finder():
			_needReinit(true) { }

		Finder(THaystack &haystack):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) { }

		template <typename TRepeatSize, typename TPeriodSize>
		Finder(THaystack &haystack, TRepeatSize minRepeatLen, TPeriodSize maxPeriod):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) 
		{
			findRepeats(data_repeats, haystack, minRepeatLen, maxPeriod);
		}

		Finder(TIterator &iter):
			data_iterator(iter),
			_needReinit(true) { }

		Finder(TIterator const &iter):
			data_iterator(iter),
			_needReinit(true) { }

		Finder(Finder const &orig):
			data_iterator(orig.data_iterator),
			haystackEnd(orig.haystackEnd),
			_needReinit(orig._needReinit),
			hits(orig.hits),
            startPos(orig.startPos),
            curPos(orig.curPos),
            endPos(orig.endPos),
            dotPos(orig.dotPos),
            dotPos2(orig.dotPos2),
			data_repeats(orig.data_repeats)
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
			curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
			endRepeat = end(data_repeats, Standard());
		};

		inline typename Reference<TIterator>::Type 
		operator* () { return value(hostIterator(*this)); }

		inline typename Reference<TIterator const>::Type 
		operator* () const { return value(hostIterator(*this)); }

		operator TIterator () const	{ return data_iterator;	}
        
        Finder & operator = (Finder const &orig) const 
        {
            data_iterator = orig.data_iterator;
            haystackEnd = orig.haystackEnd;
            _needReinit = orig._needReinit;
            hits = orig.hits;
            startPos = orig.startPos;
            curPos = orig.curPos;
            endPos = orig.endPos;
            dotPos = orig.dotPos;
            dotPos2 = orig.dotPos2;
            data_repeats = orig.data_repeats;
            curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
            endHit = end(hits, Standard());
            curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
            endRepeat = end(data_repeats, Standard());
            return *this;
        }
    };


//____________________________________________________________________________

	// forward
    template < typename TInput, typename TSpec >
    struct Pipe;

	template < typename TTuples, typename TPipeSpec, typename TSpec >
	class Finder< Pipe<TTuples, TPipeSpec>, Swift<TSpec> >
	{
	public:
		typedef Pipe<TTuples, TPipeSpec>						TInput;
		typedef typename Size<TInput>::Type						THstkPos;
		typedef _SwiftHit<__int64>								TSwiftHit;
		typedef String<TSwiftHit>								THitString;
		typedef typename Iterator<THitString, Standard>::Type	THitIterator;

		TInput			&in;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		curPos, dotPos, dotPos2;

		Finder(TInput &_in):
			in(_in),
			_needReinit(true) {}

		Finder(Finder const &orig):
			in(orig.in),
			hits(orig.hits),
			_needReinit(orig._needReinit) 
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
		};
	};


//____________________________________________________________________________

	
	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack, Swift<TSpec> > & me)
	{
		return hostIterator(hostIterator(me)) == hostIterator(me.haystackEnd);
	}

	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack, Swift<TSpec> > & me)
	{
		hostIterator(me) = me.haystackEnd;
	}


//____________________________________________________________________________


	template < typename TNeedle, typename TIndexSpec, typename TSpec >
	class Pattern< Index<TNeedle, TIndexSpec>, Swift<TSpec> >
	{
	public:
		typedef Index<TNeedle, TIndexSpec>								TIndex;
		typedef typename Size<TIndex>::Type								TSize;
		typedef unsigned												TShortSize;
		typedef typename Fibre<TIndex, Tag<_Fibre_SA> const >::Type		TSA;
		typedef typename Fibre<TIndex, Tag<_Fibre_Shape> const >::Type	TShape;
		typedef typename Iterator<TSA const, Standard>::Type			TIterator;
		
		typedef _SwiftBucket<TSpec, TSize, TShortSize>					TBucket;
		typedef String<TBucket>											TBucketString;
		typedef _SwiftBucketParams<TSpec, TSize, TShortSize>			TBucketParams;
		typedef String<TBucketParams>									TBucketParamsString;
		
		TShape					shape;
		TBucketString			buckets;
		TBucketParamsString		bucketParams;
		SwiftParameters			params;
		unsigned				curSeqNo;

		double					_currentErrorRate;
		int						_currentMinLengthForAll;

		Holder<TIndex>	data_host;

		Pattern() 
		{
			clear(*this);
		}
		Pattern(TIndex &_index): data_host(_index) 
		{
			clear(*this);
		}
		Pattern(TIndex const &_index): data_host(_index)
		{
			clear(*this);
		}
	};
	
//____________________________________________________________________________


template <typename TParams>
inline void _printSwiftParams(TParams &bucketParams)
{
	::std::cout << "  firstBucket: " << bucketParams.firstBucket << ::std::endl;
	::std::cout << "  reuseMask:   " << bucketParams.reuseMask << ::std::endl;
//	::std::cout << "  distanceCut: " << bucketParams.distanceCut << ::std::endl;
	::std::cout << "  delta:       " << bucketParams.delta << ::std::endl;
	::std::cout << "  threshold:   " << bucketParams.threshold << ::std::endl;
	::std::cout << "  overlap:     " << bucketParams.overlap << ::std::endl;
	::std::cout << "  logDelta:    " << (int)bucketParams.logDelta << ::std::endl << ::std::endl;
}

template < typename TNeedle, typename TIndexSpec, typename TSpec >
inline void _printSwiftBuckets(Pattern< Index<TNeedle, TIndexSpec>, Swift<TSpec> > &p)
{
	typedef Index<TNeedle, TIndexSpec> TIndex;
	typedef typename Pattern<TIndex, Swift<TSpec> >::TBucketParams TParams;

	unsigned j = 0;
	TParams *bucketParams = &_swiftBucketParams(p, 0);

	for(unsigned i=0; i<length(p.buckets) && i<10; ++i) 
	{
		if ((i & bucketParams->reuseMask) == 0)
		{
			::std::cout << ::std::endl << "ReadBucket #" << j << "    " << '"';
			::std::cout << indexText(host(p))[j] << '"' << ::std::endl;
			::std::cout << "  length:      " << sequenceLength(j, host(p)) << ::std::endl;
			bucketParams = &_swiftBucketParams(p, j++);
			_printSwiftParams(*bucketParams);
		}

		::std::cout << "    lastInc: " << (int)p.buckets[i].lastIncrement;
		::std::cout << "  \tCounter: " << p.buckets[i].counter << ::std::endl;
	}
}

template <typename TIndex, typename TSpec, typename TSize>
inline typename Pattern<TIndex, Swift<TSpec> >::TBucketParams &
_swiftBucketParams(Pattern<TIndex, Swift<TSpec> > & pattern, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return pattern.bucketParams[sequenceLength(seqNo, host(pattern))];
	else
		return pattern.bucketParams[seqNo];
}

template <typename TIndex, typename TSpec, typename TParams, typename TSize>
inline unsigned
_swiftBucketNo(Pattern<TIndex, Swift<TSpec> > const &, TParams &bucketParams, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return (bucketParams.reuseMask + 1) * seqNo;
	else
		return bucketParams.firstBucket;
}

template <typename TIndex, typename TSpec, typename TSeqNo>
inline int
_qgramLemma(Pattern<TIndex, Swift<TSpec> > const & pattern, TSeqNo seqNo, int errors)
{
	// q-gram lemma: How many conserved q-grams we see at least?
	// each error destroys at most <weight> many (gapped) q-grams
	return 
		sequenceLength(seqNo, host(pattern)) - length(indexShape(host(pattern))) + 1 
		- errors * weight(indexShape(host(pattern)));
}

template <typename TIndex, typename TSpec, typename TSeqNo, typename TThreshold>
inline void
setMinThreshold(Pattern<TIndex, Swift<TSpec> > & pattern, TSeqNo seqNo, TThreshold thresh) 
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;

	TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
	TBucketIterator it = begin(pattern.buckets, Standard()) + bucketParams.firstBucket;
	TBucketIterator itEnd = it + bucketParams.reuseMask;

	for (; it != itEnd; ++it)
		if ((*it).threshold < thresh)
			(*it).threshold = thresh;
}


template <typename TIndex, typename TFloat, typename _TSize, typename TSpec>
inline void _patternInit(Pattern<TIndex, Swift<TSpec> > &pattern, TFloat errorRate, _TSize minLengthForAll) 
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;
	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucket							TBucket;
	typedef typename TBucket::TSize								TBucketSize;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;

	double _newErrorRate = errorRate;
	if (pattern._currentErrorRate == _newErrorRate &&
		pattern._currentMinLengthForAll == minLengthForAll)	return;

	indexRequire(host(pattern), QGram_SADir());
	pattern.shape = indexShape(host(pattern));

	TIndex const &index = host(pattern);
	TSize seqCount = countSequences(index);
	TSize span = length(pattern.shape);
	TSize count = 0;
	TSize bucketsPerCol2Max = 0;
	TSize maxLength = 0;

	pattern._currentErrorRate = _newErrorRate;
	pattern._currentMinLengthForAll = minLengthForAll;
	
	if (Swift<TSpec>::PARAMS_BY_LENGTH) {
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) {
			TSize length = sequenceLength(seqNo, host(pattern));
			if (maxLength < length)
				maxLength = length;
		}
		resize(pattern.bucketParams, maxLength + 1);
	} else
		resize(pattern.bucketParams, seqCount);
	
	if (minLengthForAll != 0) 
	{
		// global matches
		TSize minLength = minLengthForAll;
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
		{
			// swift q-gram lemma
			TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
			// n..next length that could decrease threshold
			TSize n = (TSize) ceil((floor(errorRate * minLength) + 1) / errorRate);
			// minimal threshold is minimum errors of minLength and n
			int threshold = (TSize) min(
				(n + 1) - span * (floor(errorRate * n) + 1),
				(minLength + 1) - span * (floor(errorRate * minLength) + 1));

			if (threshold > pattern.params.minThreshold)
				bucketParams.threshold = threshold;
			else
				bucketParams.threshold = pattern.params.minThreshold;

			TSize errors = (TSize) floor((2 * bucketParams.threshold + span - 3) / (1 / errorRate - span));
			bucketParams.overlap = errors;
//			bucketParams.distanceCut = (bucketParams.threshold - 1) + span * (errors + span);
			bucketParams.logDelta = (TSize) ceil(log((double)errors) / log(2.0));
			if (bucketParams.logDelta < pattern.params.minLog2Delta) 
				bucketParams.logDelta = pattern.params.minLog2Delta;
			bucketParams.delta = 1 << bucketParams.logDelta;
			bucketParams.tabooLength = pattern.params.tabooLength;
			// TODO: classical swift for rectangular buckets
		}
	} else
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
		{
			// get pattern length and max. allowed errors
			TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
			TSize length = sequenceLength(seqNo, host(pattern));
			TSize errors = (TSize) floor(errorRate * length);
			TSize errorsWC = errors / (1 + Swift<TSpec>::QGRAM_ERRORS);

			// q-gram lemma: How many conserved q-grams we see at least?
			// (define a minimal threshold of 1)
			int threshold = length - span + 1 - errorsWC * weight(pattern.shape);
			if (threshold > pattern.params.minThreshold)
				bucketParams.threshold = threshold;
			else
				bucketParams.threshold = pattern.params.minThreshold;
			
			if (Swift<TSpec>::HAMMING_ONLY != 0)
				errors = 0;			

			// a bucket has distanceCut different positions of q-grams
			// if a q-gram is this far or farer away it can't belong to the
			// same bucket
//			bucketParams.distanceCut = length - (span - 1) + errors;

			TSize bucketsPerCol2;
			if (Swift<TSpec>::DIAGONAL == 1)
			{
				// Use overlapping parallelograms				
				bucketParams.overlap = errors;
				
				// delta must be a power of 2 greater then errors (define a minimal delta of 8)
				bucketParams.logDelta = (TSize) ceil(log((double)(errors + 1)) / log(2.0));
				if (bucketParams.logDelta < pattern.params.minLog2Delta) 
					bucketParams.logDelta = pattern.params.minLog2Delta;
				bucketParams.delta = 1 << bucketParams.logDelta;

				// the formula for bucketsPerCol is (worst-case):
				// (height-(q-1) - 1 - (delta+1-e))/delta + 3
				//    ^-- full paral. in the middle --^     ^-- 2 at the bottom, 1 at the top
				TSize bucketsPerCol = (length - span + 2 * bucketParams.delta + errors - 1) / bucketParams.delta;
				bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0));
			}
			else
			{
				// Use overlapping rectangles
				bucketParams.overlap = length - span + errors;

				// delta must be a power of 2 greater then seq.length + errors (define a minimal delta of 32)
				bucketParams.logDelta = (TSize) ceil(log((double)(length - span + 1 + errors)) / log(2.0));
				if (bucketParams.logDelta < pattern.params.minLog2Delta) 
					bucketParams.logDelta = pattern.params.minLog2Delta;
				bucketParams.delta = 1 << bucketParams.logDelta;

				bucketsPerCol2 = 2;
			}

//			SEQAN_ASSERT(distanceCut <= bucketsPerCol * (TSize) delta);

			bucketParams.firstBucket = count;
			bucketParams.reuseMask = bucketsPerCol2 - 1;
			bucketParams.tabooLength = pattern.params.tabooLength;
			
			if (Swift<TSpec>::PARAMS_BY_LENGTH) {
				++count;
				if (bucketsPerCol2Max < bucketsPerCol2)
					bucketsPerCol2Max = bucketsPerCol2;
			} else
				count += bucketsPerCol2;
			
/*			if (seqNo<3)
				_printSwiftParams(bucketParams);
*/		}

	if (Swift<TSpec>::PARAMS_BY_LENGTH) {
		count *= bucketsPerCol2Max;
		for(unsigned i = 0; i < length(pattern.bucketParams); ++i)
			pattern.bucketParams[i].reuseMask = bucketsPerCol2Max - 1;
	}
	resize(pattern.buckets, count);

	TBucketIterator	bkt = begin(pattern.buckets, Standard());
	TBucketIterator	bktEnd;
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
	{
		TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
		bktEnd = bkt + bucketParams.reuseMask + 1;
		for(; bkt != bktEnd; ++bkt) 
		{
			(*bkt).lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
			(*bkt).counter = 0;
			(*bkt).threshold = bucketParams.threshold;
		}
	}
}

template <
	typename THaystack,
	typename TIndex, 
	typename TSpec,
	typename THValue
>
inline bool _swiftMultiProcessQGram(
	Finder<THaystack, Swift<TSpec> > &finder, 
	Pattern<TIndex, Swift<TSpec> > &pattern,
	THValue hash)
{
	typedef Finder<THaystack, Swift<TSpec> >					TFinder;
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIter;
	typedef typename Value<TBucketString>::Type					TBucket;
	typedef typename TBucket::TShortSize						TShortSize;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TFinder::TSwiftHit							THit;
	
	TIndex const &index = host(pattern);
	
	TSAIter saBegin = begin(indexSA(index), Standard());
	TSAIter occ = saBegin + indexDir(index)[hash];
	TSAIter occEnd = saBegin + indexDir(index)[hash + 1];
	TBucketIter bktBegin = begin(pattern.buckets, Standard());
	Pair<unsigned> ndlPos;
	
	for(; occ != occEnd; ++occ) 
	{
		posLocalize(ndlPos, *occ, stringSetLimits(index));
		TBucketParams &bucketParams = _swiftBucketParams(pattern, getSeqNo(ndlPos));

		__int64 diag = finder.curPos;
		if (Swift<TSpec>::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);
		
		unsigned bktNo = (diag >> bucketParams.logDelta) & bucketParams.reuseMask;
		unsigned bktOfs = diag & (bucketParams.delta - 1);
		__int64  bktBeginHstk = diag & ~(__int64)(bucketParams.delta - 1);

		TBucketIter bkt = bktBegin + (_swiftBucketNo(pattern, bucketParams, getSeqNo(ndlPos)) + bktNo);		
		TShortSize hitCount;

		do 
		{
			if ((__int64)((*bkt).lastIncrement + bktOfs) < diag)
			{
				// last increment was before the beginning of the current bucket
				// (we must ensure that bucketIdx doesn't collide)
				hitCount = 1;
			}
			else
			{
				if ((*bkt).lastIncrement + bucketParams.tabooLength > finder.curPos) 
					goto checkOverlap;	// increment only once per sequence			
				hitCount = (*bkt).counter + 1;
			}

			(*bkt).lastIncrement = finder.curPos;
			(*bkt).counter = hitCount;
#ifdef SEQAN_DEBUG_SWIFT
			(*bkt)._lastIncDiag = diag;
#endif

			if (hitCount == (*bkt).threshold)
			{

				TSize height = 0;
				if (Swift<TSpec>::DIAGONAL == 1)
					height = sequenceLength(getSeqNo(ndlPos), host(pattern)) - 1;

#ifdef SEQAN_DEBUG_SWIFT
				// upper bucket no. of lastIncr. q-gram
				__int64 upperBktNo = (*bkt).lastIncrement >> bucketParams.logDelta;

				// we must decrement bucket no. until (no. mod reuse == bktNo)
				__int64 _bktBeginHstk = 
					 (upperBktNo - ((upperBktNo - bktNo) & bucketParams.reuseMask)) << bucketParams.logDelta;

				if ((*bkt)._lastIncDiag - _bktBeginHstk >= bucketParams.delta + bucketParams.overlap || (*bkt)._lastIncDiag < _bktBeginHstk) {
					::std::cerr << "qgram stored in wrong bucket (diag:" << (*bkt)._lastIncDiag << ", begin:" << _bktBeginHstk;
					::std::cerr << ", delta:" << bucketParams.delta << ", overlap:" << bucketParams.overlap << ")" << ::std::endl;
				}
#endif
//				if (bktBeginHstk >= 0) 
//				{
					THit hit = {
						bktBeginHstk,										// bucket begin in haystack
						getSeqNo(ndlPos),									// needle seq. number
						height + bucketParams.delta + bucketParams.overlap	// bucket width (non-diagonal)
					};
					appendValue(finder.hits, hit);
//				} else {
//					// match begins left of haystack begin
//					THit hit = {
//						0,													// bucket begin in haystack
//						getSeqNo(ndlPos),									// needle seq. number
//						height + bucketParams.delta + bucketParams.overlap	// bucket width (non-diagonal)
//						+ (diag & ~(__int64)(bucketParams.delta - 1))
//					};
//					appendValue(finder.hits, hit);
//				}
			}

		checkOverlap:
			if (bktOfs >= bucketParams.overlap) break;

			// repeat with the previous overlapping bucket
			bktBeginHstk -= bucketParams.delta;
			bktOfs += bucketParams.delta;
			if (bktNo) {
				--bktNo;
				--bkt;
			} else {
				bktNo = bucketParams.reuseMask;
				bkt += bktNo;
			}
		} while (true);
	}

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

template <
	typename TIndex, 
	typename TSpec
>
inline void _swiftMultiFlushBuckets(Pattern<TIndex, Swift<TSpec> > &pattern)
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename TPattern::TBucket							TBucket;
	typedef typename TBucket::TSize								TBucketSize;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;
	typedef typename TPattern::TBucketParams					TBucketParams;

	typedef typename Size<TIndex>::Type							TSize;

	TBucketIterator	bkt = begin(pattern.buckets, Standard());
	TBucketIterator	bktEnd;
	TSize seqCount = countSequences(host(pattern));

	for(TSize ndlSeqNo = 0; ndlSeqNo < seqCount; ++ndlSeqNo) 
	{
		TBucketParams &bucketParams = _swiftBucketParams(pattern, ndlSeqNo);
		bktEnd = bkt + (bucketParams.reuseMask + 1);
		for(unsigned bktNo = 0; bkt != bktEnd; ++bkt, ++bktNo)
		{
			(*bkt).lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
			(*bkt).counter = 0;
		}
	}
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
empty(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > & me) 
{
	return empty(me.bucketParams);
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline void 
clear(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > & me) 
{
	me._currentErrorRate = -1;
	me._currentMinLengthForAll = -1;
	clear(me.bucketParams);
	clear(me.buckets);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > &finder)
{
	typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(haystack(finder));

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return TPair((TPosition)hitBegin, (TPosition)hitEnd);
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > const &finder)
{
	typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(haystack(finder));

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return TPair((TPosition)hitBegin, (TPosition)hitEnd);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
range(Finder<THaystack, Swift<TSpec> > &finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(haystack(finder));

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return infix(haystack(finder), hitBegin, hitEnd);
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
range(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(text);

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return infix(text, hitBegin, hitEnd);
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline typename Value<TNeedle>::Type &
range(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern)
{
	return indexText(needle(pattern))[pattern.curSeqNo];
}

template <typename THaystack, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
_nextNonRepeatRange(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern,
	bool /*printDots*/)
{
	typedef typename Finder<THaystack, Swift<TSpec> >::TRepeat	TRepeat;
	typedef typename Value<TRepeat>::Type						TPos;

	if (finder.curRepeat == finder.endRepeat) return false;

	do 
	{
		finder.startPos = (*finder.curRepeat).endPosition;
		if (++finder.curRepeat == finder.endRepeat) 
		{
			finder.endPos = length(host(finder));
			if (finder.startPos + length(pattern.shape) > finder.endPos)
				return false;
			else
				break;
		} else
			finder.endPos = (*finder.curRepeat).beginPosition;
		// repeat until the shape fits in non-repeat range
	} while (finder.startPos + length(pattern.shape) > finder.endPos);

	finder.curPos = finder.startPos;
	hostIterator(finder) = begin(host(finder)) + finder.startPos;
	finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);

//	if (printDots)
//		::std::cerr << ::std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") ";

	return true;
}

template <typename THaystack, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
_firstNonRepeatRange(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern,
	bool printDots)
{
	typedef typename Finder<THaystack, Swift<TSpec> >::TRepeat	TRepeat;
	typedef typename Value<TRepeat>::Type						TPos;

	finder.curRepeat = begin(finder.data_repeats, Standard());
	finder.endRepeat = end(finder.data_repeats, Standard());

	if (finder.curRepeat == finder.endRepeat)
		finder.endPos = length(host(finder));
	else
		finder.endPos = (*finder.curRepeat).beginPosition;
		
	if (length(pattern.shape) > finder.endPos)
		return _nextNonRepeatRange(finder, pattern, printDots);

	finder.curPos = finder.startPos = 0;
	hostIterator(finder) = begin(host(finder));
	finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);
	
//	if (printDots)
//		::std::cerr << ::std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") ";
	
	return true;
}

template <typename THaystack, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
find(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern, 
	double errorRate,
	bool printDots)
{
	typedef Index<TNeedle, TIndexSpec>					TIndex;
	typedef typename Fibre<TIndex, QGram_Shape>::Type	TShape;
	typedef	typename Value<TShape>::Type				THashValue;

	if (empty(finder)) 
	{
		_patternInit(pattern, errorRate, 0);
		_finderSetNonEmpty(finder);
		finder.dotPos = 100000;
		finder.dotPos2 = 10 * finder.dotPos;

		if (!_firstNonRepeatRange(finder, pattern, printDots)) return false;
		if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}
	} else
		if (++finder.curHit < finder.endHit) 
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}

	// all previous matches reported -> search new ones
	clear(finder.hits);

	// are we at the end of the text?
	if (atEnd(finder) && finder.curRepeat == finder.endRepeat) 
	{
		finder.curHit = finder.endHit;
		return false;
	}

	do {

		if (printDots)
			while (finder.curPos >= finder.dotPos) 
			{
				finder.dotPos += 100000;
				if (finder.dotPos >= finder.dotPos2)
				{
					::std::cerr << (finder.dotPos2 / 1000000) << "M" << ::std::flush;
					finder.dotPos2 += 1000000;
				} else
					::std::cerr << "." << ::std::flush;
			}

		if (atEnd(++finder)) 
		{
			if (!_nextNonRepeatRange(finder, pattern, printDots)) 
			{
				_swiftMultiFlushBuckets(pattern);
				return false;
			}
			hash(pattern.shape, hostIterator(hostIterator(finder)));
		}
		else
		{
			++finder.curPos;
			hashNext(pattern.shape, hostIterator(hostIterator(finder)));
		}
		
		if (_swiftMultiProcessQGram(finder, pattern, value(pattern.shape)))
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}

	} while (true);
}

template <typename THashes, typename TPipeSpec, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
find(
	Finder<Pipe<THashes, TPipeSpec>, Swift<TSpec> > &finder,
	Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern, 
	double errorRate,
	bool printDots)
{
	if (empty(finder)) 
	{
		_patternInit(pattern, errorRate, 0);
		_finderSetNonEmpty(finder);
		finder.dotPos = 100000;
		finder.dotPos2 = 10 * finder.dotPos;

		beginRead(finder.in);
		if (eof(finder.in)) 
		{
			endRead(finder.in);
			return false;
		}
		finder.curPos = (*finder.in).i1;
		if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)))
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}
	} else
		if (++finder.curHit != finder.endHit) 
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}

	clear(finder.hits);
	if (eof(finder.in)) return false;

	do {
		++finder.in;
		if (eof(finder.in)) {
			endRead(finder.in);
#ifdef SEQAN_DEBUG_SWIFT
			_printSwiftBuckets(pattern);
#endif
			_swiftMultiFlushBuckets(pattern);
			return false;
		}
		finder.curPos = (*finder.in).i1;
		if (printDots)
			while (finder.curPos >= finder.dotPos) 
			{
				finder.dotPos += 100000;
				if (finder.dotPos >= finder.dotPos2)
				{
					::std::cerr << (finder.dotPos2 / 1000000) << "M" << ::std::flush;
					finder.dotPos2 += 1000000;
				} else
					::std::cerr << "." << ::std::flush;
			}
	} while (!_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)));

	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
