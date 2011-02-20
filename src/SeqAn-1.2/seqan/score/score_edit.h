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
  $Id:  $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_EDIT_H
#define SEQAN_HEADER_SCORE_EDIT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.EditDistance
..cat:Scoring
..summary:Edit distance scoring scheme.
..signature:Score<TValue, EditDistance>
..param.TValue:The value type.
...default:int
..general:Class.Score
..remarks:Semantically equivalent to a default contructed @Spec.Simple Score.Score<int, Simple>@.
..remarks:$EditDistance$ is a synonym for @Tag.LevenshteinDistance@.
*/

//EditDistance is defined in basic_tag.h

template <typename TValue>
class Score<TValue, EditDistance>
{
public:
	Score()
	{
	}

	Score(Score const &)
	{
	}
	~Score()
	{
	}

	Score & operator = (Score const &)
	{
		return *this;
	}

//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

/**
.Shortcut.EditDistanceScore:
..cat:scoring
..summary:Edit distance scoring scheme.
..signature:EditDistanceScore
..shortcutfor:Spec.EditDistance
...signature:Score<int, EditDistance>
..see:Spec.EditDistance
*/

typedef Score<int, EditDistance> EditDistanceScore;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline TValue
scoreMatch(Score<TValue, EditDistance> &)
{
	return 0;
}

template <typename TValue>
inline TValue
scoreMatch(Score<TValue, EditDistance> const &)
{
	return 0;
}

template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, EditDistance> &)
{
	return -1;
}

template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, EditDistance> const &)
{
	return -1;
}

template <typename TValue>
inline TValue
scoreGapExtend(Score<TValue, EditDistance> &)
{
	return -1;
}

template <typename TValue>
inline TValue
scoreGapExtend(Score<TValue, EditDistance> const &)
{
	return -1;
}

template <typename TValue>
inline TValue
scoreGapOpen(Score<TValue, EditDistance> &)
{
	return -1;
}

template <typename TValue>
inline TValue
scoreGapOpen(Score<TValue, EditDistance> const &)
{
	return -1;
}



//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
