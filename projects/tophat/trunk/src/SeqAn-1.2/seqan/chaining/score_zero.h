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
  $Id: score_zero.h 4615 2009-07-25 07:09:21Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_ZERO_H
#define SEQAN_HEADER_SCORE_ZERO_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
/**.Spec.Score Zero
..summary:Scoring scheme for chaining that set gap scores to 0
..cat:Chaining
..general:Class.Score
..signature:Score<TValue, Zero>
..param.TValue:Type of the score values.
*/

template <typename TValue>
class Score<TValue, Zero>
{
private:

public:
	Score()
	{
	}

	Score(Score const & )
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		if( this == & other )
			return *this;
		return *this;
	}

//____________________________________________________________________________
};

template <typename TValue>
inline TValue 
scoreMatch(Score<TValue, Zero> &)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreMatch(Score<TValue, Zero> const &)
{
	return 0;
}

template <typename TValue>
inline TValue 
scoreMismatch(Score<TValue, Zero> & me)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreMismatch(Score<TValue, Zero> const &)
{
	return 0;
}

template <typename TValue>
inline TValue 
scoreGap(Score<TValue, Zero> &)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreGap(Score<TValue, Zero> const &)
{
	return 0;
}
	
//////////////////////////////////////////////////////////////////////////////

//Shortcut:

//template< typename TValue >
//typedef typename Score<TValue, Zero> ZeroScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

template <typename TValue, typename T>
inline TValue
score(Score<TValue, Zero> const & me,
	  T const & left,
	  T const & right)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//compute score for chaining two fragments

template <typename TValue, typename TFragment>
inline TValue
scoreChainGap(Score<TValue, Zero> const &,
			  TFragment &,
			  TFragment &)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
