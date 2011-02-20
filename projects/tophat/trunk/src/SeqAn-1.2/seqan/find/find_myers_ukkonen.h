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
  $Id: find_myers_ukkonen.h 3532 2009-02-24 05:42:36Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_H

namespace SEQAN_NAMESPACE_MAIN 
{
 
//////////////////////////////////////////////////////////////////////////////
// MyersUkkonen
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Myers:
..cat:Pattern Matching
..general:Class.Pattern
..summary:Provides fast approximate searching of one string in another using Myer's fast bit-parallel algorithm with application of the Ukkonen-trick.
..signature:Pattern<TNeedle, Myers< [TSpec [, TFindBeginPatternSpec] ]> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TSpec:Specialization tag.
...default:$FindInfix$
...remarks:This could be $FindInfix$ for infix search or $FindPrefix$ for prefix search.
..param.TFindBeginPatternSpec:Specialization of @Class.Pattern@ used to find the begin of matches.
...default:@Metafunction.DefaultFindBeginPatternSpec@
...metafunctin:@Metafunction.FindBeginPatternSpec@
...remarks:This must be a finder for prefix search, e.g. @DPSearch<TScore, FindPrefix>@ or @Myers<FindPrefix>@.
Specify $void$ to suppress prefix searching.
..remarks.text:The needle-length must be smaller than the highest number that can be stored in an unsigned int.
*/

///.Class.Pattern.param.TSpec.type:Spec.MyersUkkonen

template <typename TSpec = FindInfix, typename TFindBeginPatternSpec = typename DefaultFindBeginPatternSpec<>::Type>
struct Myers;

//FindInfix and FindPrefix are defined int find_base.h
struct AlignTextBanded; // search query in a parallelogram

//deprecated shortcuts:
typedef Myers<FindInfix, void>  MyersUkkonen;		// search semi-global (query-global, text-local)
typedef Myers<FindPrefix, void> MyersUkkonenGlobal;	// search global (query-global, text-global)
typedef Myers<AlignTextBanded, void> MyersUkkonenBanded;		// search query in a parallelogram



//____________________________________________________________________________
// bit 0 of the HP bit-vector
// 0 for begin-gap-free haystack
// 1 for global alignments of haystack

template <typename T>
struct _MyersUkkonenHP0 {
	enum { VALUE = 0 };
};

template <>
struct _MyersUkkonenHP0<FindPrefix> {
	enum { VALUE = 1 };
};


//////////////////////////////////////////////////////////////////////////////
//overwrite _FindBegin to define host member if find begin is switched on

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFindBeginPatternSpec2>
struct _FindBegin< Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >, TFindBeginPatternSpec2>
{
private:
	typedef Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > TPattern;
	typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;

public:
	TFindBeginPattern data_findBeginPattern;
	Holder<TNeedle>	data_host;	//defines the 
	typedef True HasHost;
};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct _FindBegin< Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >, void>
{
	typedef False HasHost;
//need no findBegin if FindBeginPatternSpec is void
};

//////////////////////////////////////////////////////////////////////////////

struct _MyersLargePattern
{
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif
	String<TWord> VP;
	String<TWord> VN;
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row
	unsigned blockCount;		// the number of blocks
	unsigned lastBlock;			// the block containing the last active cell	
};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >: 
	public _FindBegin<Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > >
{
//____________________________________________________________________________
public:

#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned score;				// the current score
	unsigned k;					// the maximal number of differences allowed

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks

//	Holder<TNeedle>		data_host;	// by now, this holder is not needed, would waste mem.
									// but is is added again via _FindBegin base class if TFindBeginPatternSpec is not void
	_MyersLargePattern	*large;	// extra preprocessing info for large patterns

//____________________________________________________________________________

	Pattern(int _limit = -1):
		k(- _limit),
		large(NULL)
	{}

	Pattern(Pattern const & other)
		: needleSize(other.needleSize)
		, score(other.score)
		, k(other.k)
		, VP0(other.VP0)
		, VN0(other.VN0)
		, bitMasks(other.bitMasks)
		, large(NULL)
	{
		if (other.large)
		{
			large = new _MyersLargePattern;
			(*large) = *(other.large);
		}
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 & ndl, int _limit = -1):
		k(- _limit),
		large(NULL)
	{
		setHost(*this, ndl);
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int _limit = -1):
		k(- _limit),
		large(NULL)
	{
		setHost(*this, ndl);
	}


	~Pattern()
	{
		delete large;
	}


	Pattern &
	operator = (Pattern const & other)
	{
		needleSize = other.needleSize;
		score = other.score;
		k = other.k;
		VP0 = other.VP0;
		VN0 = other.VN0;
		bitMasks = other.bitMasks;
		large = NULL;
		if (other.large)
		{
			large = new _MyersLargePattern;
			(*large) = *(other.large);
		}
		return *this;
	}

//____________________________________________________________________________
};


template <typename TNeedle, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > {
//____________________________________________________________________________
public:
#ifdef SEQAN_SSE2_INT128
	typedef SSE2_int128 TWord;
#else
	typedef unsigned long TWord;
#endif

	typedef typename Iterator<TNeedle, Standard>::Type TIter;
	enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

	unsigned needleSize;
	unsigned score;				// the current score
	unsigned blockCount;		// the number of blocks
	unsigned k;					// the maximal number of differences allowed
	unsigned lastBlock;			// the block containing the last active cell

	TWord VP0;					// VP[0] (saves one dereferentiation)
	TWord VN0;					// VN[0]
	int scoreBit;

	Holder<TNeedle>		data_host;

//	String<int> mat;

	String<TWord> VP;
	String<TWord> VN;
	String<TWord> bitMasks;		// encoding the alphabet as bit-masks
	TWord scoreMask;			// the mask with a bit set at the position of the last active cell
	TWord finalScoreMask;		// a mask with a bit set on the position of the last row
	
	TIter ndlIter;				// iterate through the pattern
	
//____________________________________________________________________________

	Pattern(int _limit = -1):
		k(- _limit)
	{}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl, int _limit = -1):
		k(- _limit)
	{
		setHost(*this, ndl);
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > > 
{
	typedef TFindBeginPatternSpec Type;
};
template <typename TNeedle, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<FindPrefix, TFindBeginPatternSpec> > >
{// no find begin for FindPrefix
	typedef void Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
							  TNeedle2 & needle)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	me.needleSize = length(needle);
	unsigned blockCount = (me.needleSize + me.MACHINE_WORD_SIZE - 1) / me.MACHINE_WORD_SIZE;

	if (blockCount > 1) 
	{
		if (me.large == NULL)
			me.large = new _MyersLargePattern();

		me.large->blockCount = blockCount;
		me.large->finalScoreMask = (TWord)1 << ((me.needleSize + me.MACHINE_WORD_SIZE - 1) % me.MACHINE_WORD_SIZE);
	} 
	else 
	{
		delete me.large;
		me.large = NULL;
	}

	clear(me.bitMasks);
	fill(me.bitMasks, (ValueSize<TValue>::VALUE + 1) * blockCount, 0, Exact());

	// encoding the letters as bit-vectors
    for (unsigned j = 0; j < me.needleSize; j++)
		me.bitMasks[
			blockCount * ordValue((typename Value<TNeedle>::Type) value(needle, j))
			+ j / me.MACHINE_WORD_SIZE
		] |= (TWord)1 << (j % me.MACHINE_WORD_SIZE);
		//me.bitMasks[me.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/me.MACHINE_WORD_SIZE] = me.bitMasks[me.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));

	_findBeginInit(me);
}

template <typename TNeedle, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > & me, 
							  TNeedle2 & ndl)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> >::TWord TWord;
	me.needleSize = length(ndl);
	me.finalScoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);

	_findBeginInit(me);
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline void _patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, bool match)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;
	unsigned blockCount = (me.large == NULL)? 1: me.large->blockCount;

	// letters are encoded as bit-vectors
	for (unsigned j = 0; j < me.needleSize; j++)
	{
		TWord bit = (TWord)1 << (j % me.MACHINE_WORD_SIZE);
		bool allNull = true;
		int idx = j / me.MACHINE_WORD_SIZE;

		for (int i = 0; i < 4; ++i, idx += blockCount)
			allNull &= (me.bitMasks[idx] & bit) == (TWord)0;

		if (allNull)
		{	// all bits are 0 => this letter must be 'N'
			if (match)
			{
				for (; idx >= 0; idx -= blockCount)
					me.bitMasks[idx] |= bit;
			} else
			{
				for (; idx >= 0; idx -= blockCount)
					me.bitMasks[idx] &= ~bit;
			}
		}
	}
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline void _patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, bool match)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	unsigned blockCount = (me.large == NULL)? 1: me.large->blockCount;

	// letters are encoded as bit-vectors
	if (match)
	{
		for (unsigned j = 0; j < me.needleSize; j++)
			me.bitMasks[blockCount * 4 + j / me.MACHINE_WORD_SIZE] |= (TWord)1 << (j % me.MACHINE_WORD_SIZE);
	} else {
		for (unsigned j = 0; j < me.needleSize; j++)
			me.bitMasks[blockCount * 4 + j / me.MACHINE_WORD_SIZE] &= ~((TWord)1 << (j % me.MACHINE_WORD_SIZE));
	}
}



template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void _myersSetHost(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
				   TNeedle2 const & ndl,
				   True) 
{
	setValue(me.data_host, ndl);
}
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void _myersSetHost(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > &, 
				   TNeedle2 const &,
				   False) 
{
}


template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
			 TNeedle2 & ndl)
{
SEQAN_CHECKPOINT
	typedef Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(me, ndl, typename TPattern::HasHost());
	_patternFirstInit(me, ndl);
}
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void setHost(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
			 TNeedle2 const & ndl) 
{
	typedef Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > TPattern;
	_myersSetHost(me, ndl, typename TPattern::HasHost());
	_patternFirstInit(me, ndl);
}

//____________________________________________________________________________
/*
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TNeedle temp;
	resize(temp, me.needleSize, Exact());

	TValue v = TValue();
	for (unsigned i = 0; i < length(me.bitMasks); i += me.blockCount)
	{
		for (unsigned j = 0; j < me.needleSize; j++)
			if ((me.bitMasks[i + j / me.MACHINE_WORD_SIZE] & (TWord)1 << (j % me.MACHINE_WORD_SIZE)) != (TWord)0)
				temp[j] = v;
		++v;
	}
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >  const & me)
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TNeedle temp;
	resize(temp, me.needleSize, Exact());

	TValue v = TValue();
	for (unsigned i = 0; i < length(me.bitMasks); i += me.blockCount)
	{
		for (unsigned j = 0; j < me.needleSize; j++)
			if ((me.bitMasks[i + j / me.MACHINE_WORD_SIZE] & (TWord)1 << (j % me.MACHINE_WORD_SIZE)) != (TWord)0)
				temp[j] = v;
		++v;
	}
}
*/

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > >::Type & 
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > const>::Type & 
host(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >  const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int 
scoreLimit(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > const & me)
{
SEQAN_CHECKPOINT
	return - (int) me.k;
}

//____________________________________________________________________________

///.Function.scoreLimit.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	me.k = (- _limit);
}

//____________________________________________________________________________

///.Function.getScore.param.pattern.type:Spec.MyersUkkonen

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
int getScore(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me) 
{
	return -(int)me.score;
}
//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________


template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFinder>
void _patternInit(Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > &me, TFinder &)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;

	if (me.large == NULL)
	{
		me.score = me.needleSize;
		me.VP0 = ~(TWord)0;
		me.VN0 = 0;
	} 
	else 
	{
		_MyersLargePattern &large = *me.large;
		me.score = me.k + 1;
		large.scoreMask = (TWord)1 << (me.k % me.MACHINE_WORD_SIZE);
		large.lastBlock = me.k / me.MACHINE_WORD_SIZE; 
		if (large.lastBlock >= large.blockCount)
			large.lastBlock = large.blockCount - 1;

		clear(large.VP);
		fill(large.VP, large.blockCount, ~(TWord)0, Exact());

		clear(large.VN);
		fill(large.VN, large.blockCount, 0, Exact());
	}
}


template <typename TNeedle, typename TFinder, typename TFindBeginPatternSpec>
void _patternInit(Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > &me, TFinder &finder)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	me.ndlIter = begin(host(me), Standard());
	unsigned diagWidth = length(container(finder)) - me.needleSize;
	if (diagWidth >= me.needleSize)
		diagWidth = me.needleSize - 1;
	me.blockCount = diagWidth / me.MACHINE_WORD_SIZE + 1;

	clear(me.bitMasks);
	fill(me.bitMasks, ValueSize<TValue>::VALUE * me.blockCount, 0, Exact());

	if (me.blockCount == 1)
	{
		me.score = 0;
		me.scoreMask = 1;
		me.scoreBit = 0;
		me.VP0 = ~0;
		me.VN0 = 0;
	} 
	else 
	{
/*		me.score = me.k + 1;
		me.scoreMask = (TWord)1 << (me.k % me.MACHINE_WORD_SIZE);
		me.lastBlock = me.k / me.MACHINE_WORD_SIZE; 
		if (me.lastBlock >= me.blockCount)
			me.lastBlock = me.blockCount - 1;
*/
		clear(me.VP);
		fill(me.VP, me.blockCount, 0, Exact());

		clear(me.VN);
		fill(me.VN, me.blockCount, 0, Exact());
	}
}




//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////



template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder, 
									 Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;
	_MyersLargePattern &large = *me.large;

	while (position(finder) < haystack_length) 
	{
		carryD0 = carryHN = 0;
		carryHP = (int)_MyersUkkonenHP0<TSpec>::VALUE;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = large.lastBlock + (unsigned)(large.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == large.blockCount)
			limit--;

		shift = large.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) 
		{
			X = me.bitMasks[shift + currentBlock] | large.VN[currentBlock];
	
			temp = large.VP[currentBlock] + (X & large.VP[currentBlock]) + carryD0;
			if (carryD0 != (TWord)0)
				carryD0 = temp <= large.VP[currentBlock];
			else
				carryD0 = temp < large.VP[currentBlock];
			
			D0 = (temp ^ large.VP[currentBlock]) | X;
			HN = large.VP[currentBlock] & D0;
			HP = large.VN[currentBlock] | ~(large.VP[currentBlock] | D0);
			
			X = (HP << 1) | carryHP;
			carryHP = HP >> (me.MACHINE_WORD_SIZE - 1);
			
			large.VN[currentBlock] = X & D0;

			temp = (HN << 1) | carryHN;
			carryHN = HN >> (me.MACHINE_WORD_SIZE - 1);
								
		 	large.VP[currentBlock] = temp | ~(X | D0);

			//if the current block is the one containing the last active cell the new score is computed
			if (currentBlock == large.lastBlock) {
				if ((HP & large.scoreMask) != (TWord)0)
					me.score++;
				else if ((HN & large.scoreMask) != (TWord)0)
					me.score--;
			}
		}

		// updating the last active cell
		while (me.score > me.k) {
			if ((large.VP[large.lastBlock] & large.scoreMask) != (TWord)0)
				me.score--;
			else if ((large.VN[large.lastBlock] & large.scoreMask) != (TWord)0)
				me.score++;

			large.scoreMask >>= 1;
			if (large.scoreMask == (TWord)0) 
			{
				large.lastBlock--;
				if (TYPECMP<TSpec, FindPrefix>::VALUE && large.lastBlock == (unsigned)-1)
					break;
				large.scoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);
			}
		}

		if ((large.scoreMask == large.finalScoreMask) && (large.lastBlock == large.blockCount - 1))
		{
			_setFinderEnd(finder);
			if (TYPECMP<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
			return true;
		}
		else {
			large.scoreMask <<= 1;
			if (!large.scoreMask) {
				large.scoreMask = 1;
				large.lastBlock++;
			}
			
			if ((large.VP[large.lastBlock] & large.scoreMask) != (TWord)0)
				me.score++;
			else if ((large.VN[large.lastBlock] & large.scoreMask) != (TWord)0)
				me.score--;
		}

//		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}

//____________________________________________________________________________
// version for needles not longer than one machineword

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersSmallPatterns (TFinder & finder, 
									 Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me,
									 TSize haystack_length) 
{
SEQAN_CHECKPOINT

	typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;

	TWord X, D0, HN, HP;
	TWord lastBit = (TWord)1 << (me.needleSize - 1);

	// computing the blocks
	while (position(finder) < haystack_length) 
	{
		X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
		
		D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
		HN = me.VP0 & D0;
		HP = me.VN0 | ~(me.VP0 | D0);
		X = (HP << 1) | (TWord)(int)_MyersUkkonenHP0<TSpec>::VALUE;
		me.VN0 = X & D0;
		me.VP0 = (HN << 1) | ~(X | D0);

		if ((HP & lastBit) != (TWord)0)
			me.score++;
		else if ((HN & lastBit) != (TWord)0)
			me.score--;

		if (me.score <= me.k)
		{
			_setFinderEnd(finder);
			if (TYPECMP<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
			return true;
		}
/*
		if (TYPECMP<TSpec, FindPrefix>::VALUE)
		{//limit haystack length during prefix search

		}
*/		
		goNext(finder);
	}

	return false;
}




//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen as a banded alignment
// the band includes the main diagonal and the diagonals above
// the band width is (blockCount * MACHINE_WORD_SIZE)
//////////////////////////////////////////////////////////////////////////////



template <typename TFinder, typename TNeedle, typename TFindBeginPatternSpec, typename TSize>
inline bool 
_findMyersLargePatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > & me,
	TSize) 
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> >::TWord TWord;
	typedef typename Value<TNeedle>::Type TValue;

	TWord X, D0, HN, HP, temp;
	TWord carryD0, carryHP, carryHN;
	unsigned shift, limit, currentBlock;

	while (!atEnd(finder)) {
		// shift bitmasks and states
		if (!atEnd(me.ndlIter)) 
		{
			TWord carryVN = 0;
			TWord carryVP = 1;
			for(int j = me.blockCount - 1; j >= 0; --j) 
			{
				TWord newCarryVN = me.VN[j] & 1;
				TWord newCarryVP = me.VP[j] & 1;
				me.VN[j] = (me.VN[j] >> 1) | (carryVN << (me.MACHINE_WORD_SIZE - 1));
				me.VP[j] = (me.VP[j] >> 1) | (carryVP << (me.MACHINE_WORD_SIZE - 1));
				carryVN = newCarryVN;
				carryVP = newCarryVP;
			}
			for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i) 
			{
				TWord carry = 0;
				for(int j = me.blockCount - 1; j >= 0; --j)
				{
					unsigned pos = i * ValueSize<TValue>::VALUE + j;
					TWord newCarry = me.bitMasks[pos] & 1;
					me.bitMasks[pos] = (me.bitMasks[pos] >> 1) | (carry << (me.MACHINE_WORD_SIZE - 1));
					carry = newCarry;
				}
			}

			me.bitMasks[me.blockCount * (ordValue(*me.ndlIter) + 1) - 1]
				|= (TWord)1 << (me.MACHINE_WORD_SIZE - 1);

			goNext(me.ndlIter);
		}

		carryD0 = carryHN = 0;
		carryHP = _MyersUkkonenHP0<AlignTextBanded>::VALUE;

		// if the active cell is the last of it's block, one additional block has to be calculated
		limit = me.lastBlock + (me.scoreMask >> (me.MACHINE_WORD_SIZE - 1));

		if (limit == me.blockCount)
			limit--;

		shift = me.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

		// computing the necessary blocks, carries between blocks following one another are stored
		for (currentBlock = 0; currentBlock <= limit; currentBlock++) {
			X = me.bitMasks[shift + currentBlock] | me.VN[currentBlock];
	
			temp = me.VP[currentBlock] + (X & me.VP[currentBlock]) + carryD0;
			if (carryD0 != (TWord)0)
				carryD0 = temp <= me.VP[currentBlock];
			else
				carryD0 = temp < me.VP[currentBlock];
			
			D0 = (temp ^ me.VP[currentBlock]) | X;
			HN = me.VP[currentBlock] & D0;
			HP = me.VN[currentBlock] | ~(me.VP[currentBlock] | D0);
			
			X = (HP << 1) | carryHP;
			carryHP = HP >> (me.MACHINE_WORD_SIZE - 1);
			
			me.VN[currentBlock] = X & D0;

			temp = (HN << 1) | carryHN;
			carryHN = HN >> (me.MACHINE_WORD_SIZE - 1);
								
		 	me.VP[currentBlock] = temp | ~(X | D0);

			//if the current block is the one containing the last active cell the new score is computed
			if (currentBlock == me.lastBlock) {
				if ((HP & me.scoreMask) != (TWord)0)
					me.score++;
				else if ((HN & me.scoreMask) != (TWord)0)
					me.score--;
			}
		}

		// updating the last active cell
		while (me.score > me.k) {
			if ((me.VP[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score--;
			else if ((me.VN[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score++;

			me.scoreMask >>= 1;
			if (me.scoreMask == (TWord)0) {
				me.scoreMask = (TWord)1 << (me.MACHINE_WORD_SIZE - 1);
				me.lastBlock--;
			}
		}

		if ((me.lastBlock == me.blockCount-1) && (me.scoreMask == me.finalScoreMask))
			return true;
		else {
			me.scoreMask <<= 1;
			if (me.scoreMask == (TWord)0) {
				me.scoreMask = 1;
				me.lastBlock++;
			}
			
			if ((me.VP[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score++;
			else if ((me.VN[me.lastBlock] & me.scoreMask) != (TWord)0)
				me.score--;
		}

//		SEQAN_ASSERT (me.score >= 0);

		goNext(finder);
	}

	return false;
}
/*
template <typename TFinder, typename TNeedle>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Pattern<TNeedle, Myers<AlignTextBanded> >::TWord TWord;

	TWord X, D0, HN, HP;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	unsigned SHIFT=me.needleSize;
#endif

	if (!atEnd(me.ndlIter)) 
	{
		// Part 1: Go down the diagonal
		do
		{
			// shift bitmasks and states
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;
			me.bitMasks[ordValue(*me.ndlIter)] |= ((TWord)1 << (me.MACHINE_WORD_SIZE - 1));
			
			// Myers core
			X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
			D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
			HN = me.VP0 & D0;
			HP = me.VN0 | ~(me.VP0 | D0);
			X = (HP << 1);
			me.VN0 = X & D0;
			me.VP0 = (HN << 1) | ~(X | D0);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			--SHIFT;
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "   ";
			::std::cerr << "dD: ";
			for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
			{
				CharString vd = " 1 ";
				if ((D0 & ((TWord)1 << i)) != (TWord)0) vd = " 0 ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
#endif

			if ((D0 & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) == (TWord)0)
				++me.score;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			::std::cerr << me.score <<::std::endl;
#endif

			goNext(me.ndlIter);
//			if (atEnd(me.ndlIter))
			{
				if (me.score <= me.k)
					return true;
//				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}

	// Part 2: Go to the bottom-right of the parallelogram
	while (!atEnd(finder))
	{
		// Myers core
		X = me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | me.VN0;
		D0 = ((me.VP0 + (X & me.VP0)) ^ me.VP0) | X;
		HN = me.VP0 & D0;
		HP = me.VN0 | ~(me.VP0 | D0);
		X = (HP << 1) |1;
		me.VN0 = X & D0;
		me.VP0 = (HN << 1) | ~(X | D0);

#ifdef __SEQAN_DEBUG_MYERSBITVECTOR
		::std::cerr << "dH: ";
		for(int i=me.MACHINE_WORD_SIZE-1; i>=0 ;--i) 
		{
			CharString hd = " 0 ";
			if ((HP & ((TWord)1 << i)) != (TWord)0) hd = " 1 ";
			if ((HN & ((TWord)1 << i)) != (TWord)0) hd = "-1 ";
			::std::cerr << hd;
		}
		::std::cerr << "   ";
		::std::cerr << *finder;
#endif

		if ((HP & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) != (TWord)0)
			++me.score;
		else if ((HN & ((TWord)1 << (me.MACHINE_WORD_SIZE - 1))) != (TWord)0)
			--me.score;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
		::std::cerr << me.score <<::std::endl;
#endif

		if (me.score <= me.k)
			return true;
		goNext(finder);
	}
	return false;
}
*/

template <typename TWord, typename TAlignSpec, typename TFindBeginPatternSpec>
inline int
_myersCoreSmall(TWord &VP, TWord &VN, TWord const &bitmap, int scoreBit, Myers<TAlignSpec, TFindBeginPatternSpec> ) 
{
	// Myers core
	TWord X = bitmap | VN;
	TWord D0 = ((VP + (X & VP)) ^ VP) | X;
	TWord HN = VP & D0;
	TWord HP = VN | ~(VP | D0);
	X = HP << 1;
	VN = X & D0;
	VP = (HN << 1) | ~(X | D0);
	return (int)((HP >> scoreBit) & 1) - (int)((HN >> scoreBit) & 1);
}

template <typename TWord, typename TAlignSpec, typename TFindBeginPatternSpec>
inline int
_myersCoreSmallDiag(TWord &VP, TWord &VN, TWord const &bitmap, int scoreBit, Myers<TAlignSpec, TFindBeginPatternSpec> ) 
{
	// Myers core
	TWord X = bitmap | VN;
	TWord D0 = ((VP + (X & VP)) ^ VP) | X;
	TWord HN = VP & D0;
	TWord HP = VN | ~(VP | D0);
	X = HP << 1;
	VN = X & D0;
	VP = (HN << 1) | ~(X | D0);
	return (!D0 >> scoreBit) & 1;
}

template <typename TFinder, typename TNeedle, typename TFindBeginPatternSpec, typename TSize>
inline bool 
_findMyersSmallPatterns(
	TFinder & finder, 
	Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > & me,
	TSize)
{
SEQAN_CHECKPOINT
	typedef Pattern<TNeedle, Myers<AlignTextBanded, TFindBeginPatternSpec> > TPattern;
	typedef typename TPattern::TWord TWord;
	typedef typename TPattern::TIter TIter;

	TWord X, D0, HN, HP;
	//TWord maskMax = (TWord)1 << (length(container(finder)) - me.needleSize);
	int bitMax = length(container(finder)) - me.needleSize;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
	int SHIFT = me.needleSize - (me.ndlIter-begin(host(me), Standard()));
#endif

	TIter ndlEnd = end(host(me), Standard());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
	if (me.ndlIter != ndlEnd && me.scoreBit != bitMax)
#else
	if (me.ndlIter != ndlEnd)
#endif
	{
		// Part 1: left-upper triangle of parallelogram
		do
		{
			me.bitMasks[ordValue(*me.ndlIter)] |= me.scoreMask;
			me.score += _myersCoreSmallDiag(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			--SHIFT;
			::std::cerr << "1D:  ";
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "     ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
/*				char const *vd = "  1  ";
				if ((D0 & i) != (TWord)0) vd = "  0  ";
*/				char const *vd = "  0  ";
				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif

			me.scoreMask <<= 1;
			++me.scoreBit;
			goNext(me.ndlIter);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (me.scoreBit == bitMax || (me.ndlIter == ndlEnd))
			{
				if (me.score <= me.k)
					return true;
				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}


	if (me.ndlIter != ndlEnd)
	{
		// Part 2: go down the parallelogram
		do
		{
			me.bitMasks[ordValue(*me.ndlIter)] |= me.scoreMask;
			me.score += _myersCoreSmallDiag(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			--SHIFT;
			::std::cerr << "2D:  ";
			for(int i=0; i<SHIFT; ++i) 
				::std::cerr << "     ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
/*				char const *vd = "  0  ";
				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
*/				char const *vd = "  1  ";
				if ((D0 & i) != (TWord)0) vd = "  0  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif

			// shift bitmasks and states
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;

			goNext(me.ndlIter);
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (atEnd(me.ndlIter))
			{
				if (me.score <= me.k)
					return true;
				break;
			}
			goNext(finder);
		} while (true);
		
		goNext(finder);
	}

	if (!atEnd(finder))
	{
		// Part 3: go down the parallelogram
		do
		{
			me.score += _myersCoreSmall(
				me.VP0, me.VN0, me.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)],
				--me.scoreBit, Myers<AlignTextBanded, TFindBeginPatternSpec>());

			// shift bitmasks and states
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
			::std::cerr << "3H:  ";
			for(int i=me.scoreMask; i>0 ;i>>=1) 
			{
				char const *vd = "  0  ";
/*				if ((HP & i) != (TWord)0) vd = "  1  ";
				if ((HN & i) != (TWord)0) vd = " -1  ";
*/				if ((me.VP0 & i) != (TWord)0) vd = "  1  ";
				if ((me.VN0 & i) != (TWord)0) vd = " -1  ";
				::std::cerr << vd;
			}
			::std::cerr << "   ";
			::std::cerr << *finder;
			::std::cerr << me.score <<::std::endl;
#endif
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			resize(me.VP,1);
			resize(me.VN,1);
			me.VP[0] = me.VP0;
			me.VN[0] = me.VN0;
#endif
			me.VN0 >>= 1;
			me.VP0 = (me.VP0 >> 1) | ((TWord)1 << (me.MACHINE_WORD_SIZE - 1)); // ignore the field left of the leftmost diagonal
			for(unsigned i = 0; i < length(me.bitMasks); ++i)
				me.bitMasks[i] >>= 1;			

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
			return true;
#endif
			if (me.score <= me.k)
				return true;

			goNext(finder);
		} while (!atEnd(finder));
	}
    return false;
}


//////////////////////////////////////////////////////////////////////////////
// find
template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	int k = scoreLimit(me);

	TSize prefix_begin_position; //for prefix search: the position where the prefix begins

	if (empty(finder))
	{
		// in seqan k is treated as score, here we need it as penalty, that is why it is negated
		me.k = -k;

		_patternInit(me, finder);
		_finderSetNonEmpty(finder);

		prefix_begin_position = position(finder);

		//TODO: adapt myers-ukkonnen to dynamically change k
	}
	else
	{
		if (atEnd(finder)) return false;
		goNext(finder);

		prefix_begin_position = beginPosition(finder);
	}

	//limit search width for prefix search
	TSize haystack_length = length(container(hostIterator(finder)));
	if (TYPECMP<TSpec, FindPrefix>::VALUE)
	{
		TSize maxlen = prefix_begin_position + me.needleSize - k + 1;
		if (haystack_length > maxlen)
		{
			haystack_length = maxlen;
		}
	}

	// distinguish between the version for needles not longer than one machineword and the version for longer needles
	if (me.large == NULL) 
		return _findMyersSmallPatterns(finder, me, haystack_length);
	else
		return _findMyersLargePatterns(finder, me, haystack_length);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder, 
				  Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> > & me, 
				  int const k)
{
	setScoreLimit(me, k);
	return find(finder, me);
}

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
