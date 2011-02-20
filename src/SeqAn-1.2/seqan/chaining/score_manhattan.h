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
  $Id: score_manhattan.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_MANHATTAN_H
#define SEQAN_HEADER_SCORE_MANHATTAN_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Score spec for manhattan distance

/**.Spec.Score Manhattan
..summary:Scoring scheme for chaining that computes gap scores using manhattan distance.
..cat:Chaining
..general:Class.Score
..signature:Score<TValue, Manhattan>
..param.TValue:Type of the score values.
..remarks:The manhattan distance between two n-dimensional points is defined is the sum of the (absolute) differences of their coordinates. 
*/
template <typename TValue>
class Score<TValue, Manhattan>
{
public:
	TValue data_match;
	TValue data_mismatch;
	TValue data_gap;

public:
	Score( TValue _match = 0, TValue _misalign = 1 ):
		data_match( _match ),
		data_mismatch( _misalign ),
		data_gap( _misalign )
	{
	}

	Score( TValue score ):
		data_match( 0 ),
		data_mismatch( score ),
		data_gap( score )
	{
	}

	Score(Score const & other):
		data_match(other.data_match),
		data_mismatch(other.data_mismatch),
		data_gap(other.data_gap)
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		data_match = other.data_match;
		data_mismatch = other.data_mismatch;
		data_gap = other.data_gap;
		return *this;
	}



//____________________________________________________________________________
};


//____________________________________________________________________________

template <typename TValue>
inline TValue 
scoreMatch(Score<TValue, Manhattan> & me)
{
	return me.data_match;
}
template <typename TValue>
inline TValue
scoreMatch(Score<TValue, Manhattan> const & me)
{
	return me.data_match;
}

template <typename TValue>
inline TValue 
scoreMismatch(Score<TValue, Manhattan> & me)
{
	return me.data_mismatch;
}
template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, Manhattan> const & me)
{
	return me.data_mismatch;
}

template <typename TValue>
inline TValue 
scoreGap(Score<TValue, Manhattan> & me)
{
	return me.data_gap;
}
template <typename TValue>
inline TValue const &
scoreGap(Score<TValue, Manhattan> const & me)
{
	return me.data_gap;
}
//////////////////////////////////////////////////////////////////////////////

//Shortcut:


//template< typename TValue >
//typedef typename Score< TValue, Manhattan > ManhattanScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

//template <typename TValue, typename T>
//inline TValue
//score(Score<TValue, Manhattan> const & me,
//	  T const & left,
//	  T const & right)
//{
//	if (left == right) return scoreMatch(me);
//	else return scoreMismatch(me);
//}

//////////////////////////////////////////////////////////////////////////////
//compute score for chaining two fragments 
//return value this is only valid for f1 < f2, 
//that is f2 can be appended to f1

template <typename TValue, typename TFragment>
inline TValue
scoreChainGap(Score<TValue, Manhattan> const & me,
			  TFragment & f1,
			  TFragment & f2)
{
	SEQAN_ASSERT(dimension(f1) == dimension(f2))

	unsigned int dim = dimension(f1);
	TValue score = 0;
	TValue score_gap = scoreGap(me);
	for (unsigned int i = 0; i < dim; ++i)
	{
		score -= score_gap * (leftPosition(f2, i) - rightPosition(f1, i));
	}
	return score;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
