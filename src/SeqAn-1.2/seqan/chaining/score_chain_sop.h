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
  $Id: score_chain_sop.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_CHAIN_SOP_H
#define SEQAN_HEADER_SCORE_CHAIN_SOP_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Score spec for sum of pairs distance

/**.Spec.Score ChainSoP
..summary:Scoring scheme for chaining that uses a special method for scoring the gaps between two fragments.
..cat:Chaining
..general:Class.Score
..signature:Score<TValue, ChainSoP>
..param.TValue:Type of the score values.
..remarks:This scoring scheme is used to score gaps between to fragments in chaining.
*/

template <typename TValue>
class Score<TValue, ChainSoP>
{
private:
	TValue data_match;
	TValue data_mismatch;
	TValue data_gap;

public:
	Score( TValue _match = 0, TValue _mismatch = 3, TValue _gap = 2 ):
		data_match(_match),
		data_mismatch(_mismatch),
		data_gap(_gap)
	{
	}

	Score( TValue multiplier ):
		data_match( -1 * multiplier ),
		data_mismatch( -1 * multiplier ),
		data_gap( -1 * multiplier )
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

	friend inline TValue &
	scoreMatch(Score & me)
	{
		return me.data_match;
	}
	friend inline TValue const &
	scoreMatch(Score const & me)
	{
		return me.data_match;
	}

	friend inline TValue &
	scoreMismatch(Score & me)
	{
		return me.data_mismatch;
	}
	friend inline TValue const &
	scoreMismatch(Score const & me)
	{
		return me.data_mismatch;
	}

	friend inline TValue &
	scoreGap(Score & me)
	{
		return me.data_gap;
	}

	friend inline TValue const &
	scoreGap(Score const & me)
	{
		return me.data_gap;
	}
/*
	friend inline TValue &
	scoreGapOpen(Score<TValue, ChainSoP> & me)
	{
		return me.data_gap;
	}

	friend inline TValue const &
	scoreGapOpen(Score<TValue, ChainSoP> const & me)
	{
		return me.data_gap;
	}
*/
//____________________________________________________________________________
};


template <typename TValue>
inline TValue &
scoreMatch(Score<TValue, ChainSoP> & me)
{
	return me.data_match;
}
template <typename TValue>
inline TValue const &
scoreMatch(Score<TValue, ChainSoP> const & me)
{
	return me.data_match;
}

template <typename TValue>
inline TValue &
scoreMismatch(Score<TValue, ChainSoP> & me)
{
	return me.data_mismatch;
}
template <typename TValue>
inline TValue const &
scoreMismatch(Score<TValue, ChainSoP> const & me)
{
	return me.data_mismatch;
}

template <typename TValue>
inline TValue &
scoreGap(Score<TValue, ChainSoP> & me)
{
	return me.data_gap;
}

template <typename TValue>
inline TValue const &
scoreGap(Score<TValue, ChainSoP> const & me)
{
	return me.data_gap;
}
//////////////////////////////////////////////////////////////////////////////

//Shortcut:


//template< typename TValue >
//typedef typename Score< TValue, ChainSoP > ChainSoPScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

//template <typename TValue, typename T>
//inline TValue
//score( Score<TValue, ChainSoP> const & me,
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
scoreChainGap(Score<TValue, ChainSoP> const & me,
			  TFragment & f1,
			  TFragment & f2)
{
	SEQAN_ASSERT(dimension(f1) == dimension(f2))

	TValue score_gap = scoreGap(me);
	TValue score_miss = scoreMismatch(me);

	TValue sc = 0;
	unsigned int dim = dimension(f1);
	for (unsigned int d1 = 0; d1 < dim-1; ++d1)
	{
		TValue delta1 = (TValue) (leftPosition(f2, d1) - rightPosition(f1, d1));
		for (unsigned int d2 = d1+1; d2 < dim; ++d2)
		{
			TValue delta2 = (TValue) (leftPosition(f2, d2) - rightPosition(f1, d2));
			if (delta1 > delta2)
			{
				sc -= (score_miss*delta2 + score_gap*(delta1 - delta2));
			}
			else
			{
				sc -= (score_miss*delta1 + score_gap*(delta2 - delta1));
			}
		}
	}
	return sc;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
