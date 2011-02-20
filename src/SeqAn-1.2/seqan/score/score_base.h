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
  $Id: score_base.h 4504 2009-07-02 12:55:40Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_BASE_H
#define SEQAN_HEADER_SCORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct Simple;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.Score:
..cat:Miscellaneous
..summary:A scoring scheme.
..signature:Score<TValue, TSpec>
..param.TValue:The value type.
...default:int
..param.TSpec:The specializing type.
...default:@Spec.Simple Score@
*/
template <typename TValue = int, typename TSpec = Simple>
class Score;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value< Score<TValue, TSpec> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

//DEPRECATED
///*
//.Metafunction.ScoreValue:
//..summary:The ScoreValue type of a scoring scheme.
//..signature:Key<T>::Type
//..param.T:A type that holds a scoring scheme.
//...type:Class.Score
//..returns.param.Type:The value type of the scoring scheme.
//...default:$int$.
//*/
//
//template <typename T>
//struct ScoreValue
//{
//	typedef int Type;
//};
//template <typename TValue, typename TSpec>
//struct ScoreValue< Score<TValue, TSpec> >
//{
//	typedef TValue Type;
//};

//////////////////////////////////////////////////////////////////////////////

// ATTENTION: scoreGap...(TScore) is deprecated
// Better use the following functions:

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, TSpec> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return scoreGapOpen(me);
}

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, TSpec> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return scoreGapOpen(me);
}

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, TSpec> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return scoreGapExtend(me);
}

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, TSpec> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return scoreGapExtend(me);
}

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapHorizontal(
	Score<TValue, TSpec> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return scoreGap(me);
}

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapVertical(
	Score<TValue, TSpec> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &seq2)
{
	return scoreGap(me);
}

// ATTENTION: score(TScore, TVal1, TVal2) is deprecated
// Better use the following function:

template <typename TValue, typename TSpec, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, TSpec> const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	return score(me, seq1[pos1], seq2[pos2]);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
