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

#ifndef SEQAN_HEADER_SEQAN_CONSENSUS_SCORE_H
#define SEQAN_HEADER_SEQAN_CONSENSUS_SCORE_H


namespace SEQAN_NAMESPACE_MAIN
{


static const int SEQAN_CONSENSUS_UNITY = 1 << 20;

//////////////////////////////////////////////////////////////////////////////
// Consensus score tags
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

struct ConsensusScore_;
typedef Tag<ConsensusScore_> const ConsensusScore;

//////////////////////////////////////////////////////////////////////////////

struct FractionalScore_;
typedef Tag<FractionalScore_> const FractionalScore;

//////////////////////////////////////////////////////////////////////////////

template<typename TScore1, typename TScore2>
struct WeightedConsensusScore;




//////////////////////////////////////////////////////////////////////////////
// Scoring classes
//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// ConsensusScore
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, ConsensusScore>
{
public:
	String<TValue> consensus_set;		// Is the alphabet character part of the consensus set for each column

public:
	Score() {}

};

template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, ConsensusScore>& me,
			  TString const& profile)
{
	typedef typename Size<TString>::Type TSize;
	TSize alphSize = ValueSize<typename Value<TString>::Type>::VALUE;
	resize(me.consensus_set, alphSize * length(profile));

	typedef typename Iterator<TString, Standard>::Type TIter;
	typedef typename Iterator<String<TValue>, Standard>::Type TConsSetIter;
	TConsSetIter itConsSet = begin(me.consensus_set, Standard());
	TIter it = begin(profile, Standard());
	TIter itEnd = end(profile, Standard());
	TSize maxCount = 0;
	for(;it!=itEnd;++it) {
		maxCount = 0;
		for(TSize i = 0; i<alphSize; ++i)
			if ((*it).count[i] > maxCount) maxCount = (*it).count[i];
		for(TSize i = 0; i<alphSize; ++i, ++itConsSet)
			*itConsSet = ((*it).count[i] == maxCount)? 0 : (-SEQAN_CONSENSUS_UNITY);
	}
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, ConsensusScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return ((int) pos2 < 0) ? -SEQAN_CONSENSUS_UNITY : me.consensus_set[pos1 * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, ConsensusScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return ((int) pos2 < 0) ? -2 * SEQAN_CONSENSUS_UNITY : 2 * me.consensus_set[pos1 * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + (ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)];
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, ConsensusScore> const &,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, ConsensusScore> const &,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return -2 * SEQAN_CONSENSUS_UNITY;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ConsensusScore> const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &,
	  TSeq2 const &seq2)
{
	return me.consensus_set[pos1 * (ValueSize<typename Value<TSeq1>::Type>::VALUE) + seq2[pos2].count[0]];
}









//////////////////////////////////////////////////////////////////////////////
// FractionalScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, FractionalScore>
{
public:
	String<int> sum;		// Total number of profile characters in each column

public:
	Score() {}
};


template <typename TValue, typename TString>
inline void
assignProfile(Score<TValue, FractionalScore>& me,
			  TString const& profile)
{
	typedef typename Size<TString>::Type TSize;
	resize(me.sum, length(profile));
	typedef typename Iterator<TString, Standard>::Type TIter;
	typedef typename Iterator<String<int>, Standard>::Type TSumIter;
	TSumIter itSum = begin(me.sum, Standard());
	TIter it = begin(profile, Standard());
	TIter itEnd = end(profile, Standard());
	for(;it!=itEnd;++it, ++itSum) {
		*itSum = 0;
		for(TSize i = 0; i < (TSize) ValueSize<typename Value<TString>::Type>::VALUE; ++i) 
			*itSum += (*it).count[i];
	}
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, FractionalScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &)
{
	return (( (int) pos2 < 0) || (!me.sum[pos1])) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (( (int) seq1[pos1].count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[pos1]) * SEQAN_CONSENSUS_UNITY) / me.sum[pos1]);
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, FractionalScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &)
{
	return (( (int) pos2 < 0) || (!me.sum[pos1])) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (( (int) seq1[pos1].count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1] - me.sum[pos1]) * SEQAN_CONSENSUS_UNITY) / me.sum[pos1]);
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, FractionalScore> const &,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, FractionalScore> const &,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return -SEQAN_CONSENSUS_UNITY;
}

template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, FractionalScore> const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	return (!me.sum[pos1]) ? -SEQAN_CONSENSUS_UNITY : ((TValue) (((int) seq1[pos1].count[seq2[pos2].count[0]] - me.sum[pos1]) * SEQAN_CONSENSUS_UNITY) / me.sum[pos1]);
}








//////////////////////////////////////////////////////////////////////////////
// WeightedConsensusScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TScore1, typename TScore2>
class Score<TValue, WeightedConsensusScore<TScore1, TScore2> >
{
public:
	TScore1 sc1;
	TScore2 sc2;

public:
	Score() {}

};


template <typename TValue, typename TScore1, typename TScore2, typename TString>
inline void
assignProfile(Score<TValue, WeightedConsensusScore<TScore1, TScore2> >& me,
			  TString const& profile)
{
	assignProfile(me.sc1, profile);
	assignProfile(me.sc2, profile);
}



template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &seq2)
{
	return (scoreGapExtendHorizontal(me.sc1, pos1, pos2, seq1, seq2) + scoreGapExtendHorizontal(me.sc2, pos1, pos2, seq1, seq2)) / (TValue) 2;
}

template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &seq2)
{
	return (scoreGapOpenHorizontal(me.sc1, pos1, pos2, seq1, seq2) + scoreGapOpenHorizontal(me.sc2, pos1, pos2, seq1, seq2)) / (TValue) 2;
}


template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const & seq1,
	TSeq2 const & seq2)
{
	return (scoreGapExtendVertical(me.sc1, pos1, pos2, seq1, seq2) + scoreGapExtendVertical(me.sc2, pos1, pos2, seq1, seq2)) / (TValue) 2;
}

template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const & seq1,
	TSeq2 const & seq2)
{
	return (scoreGapOpenVertical(me.sc1, pos1, pos2, seq1, seq2) + scoreGapOpenVertical(me.sc2, pos1, pos2, seq1, seq2)) / (TValue) 2;
}


template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	return (score(me.sc1, pos1, pos2, seq1, seq2) + score(me.sc2, pos1, pos2, seq1, seq2)) / (TValue) 2;
}

}

#endif

