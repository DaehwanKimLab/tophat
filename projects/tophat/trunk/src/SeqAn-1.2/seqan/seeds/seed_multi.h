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
  $Id: seed_multi.h 3506 2009-02-20 08:14:53Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MultiSeed_H
#define SEQAN_HEADER_MultiSeed_H

namespace SEQAN_NAMESPACE_MAIN
{

struct _Seed_multi;
typedef Tag<_Seed_multi> const ChainedSeed;

/**
..Spec.ChainedSeed
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds. Additionaly diagonal segments
between start and end position2 are stored.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, ChainedSeed>
..param.TPosition:The type of number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.ChainedSeed#Seed:
..class:Spec.ChainedSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.length: Length of the seed.
*/

template<typename TPosition> 
class Seed<TPosition, ChainedSeed> {
	
public:
	std::list<Triple<TPosition, TPosition, TPosition> > seedSet;
	TPosition leftDiagonal;
	TPosition rightDiagonal;
 
	Seed(){
		SEQAN_CHECKPOINT
	}
	
	Seed(TPosition leftDim0,
		 TPosition leftDim1,
		 TPosition length)
	{
		SEQAN_CHECKPOINT
		seedSet.push_back(Triple<TPosition, TPosition, TPosition> (leftDim0, leftDim1, length));
		rightDiagonal = leftDiagonal = leftDim1 - leftDim0;
	}
	
	~Seed(){
	}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Standard Functions                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TPosition>
inline TPosition 
startDiagonal(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i2-seed.seedSet.front().i1;
}

template<typename TPosition>
inline TPosition 
endDiagonal(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.back().i2-seed.seedSet.back().i1;
}


template<typename TPosition>
inline TPosition 
leftDim0(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i1;
}


template<typename TPosition>
inline TPosition 
rightDim0(Seed<TPosition, ChainedSeed> const & seed)
{
	return seed.seedSet.back().i1+seed.seedSet.back().i3-1;
}


template<typename TPosition>
inline TPosition 
leftDim1(Seed<TPosition, ChainedSeed> const &seed)
{
	return seed.seedSet.front().i2;
}

template<typename TPosition>
inline TPosition 
rightDim1(Seed<TPosition, ChainedSeed> const & seed)
{
	return seed.seedSet.back().i2+seed.seedSet.back().i3-1;
}


template<typename TPosition>
inline TPosition 
length(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.seedSet.back().i1 + seed.seedSet.back().i3 - seed.seedSet.front().i1;	
}



template<typename TPosition>
inline void 
setLeftDim0(Seed<TPosition, ChainedSeed> &seed, 
			  TPosition start)
{
	SEQAN_CHECKPOINT
	TPosition lengthDiff = seed.seedSet.front().i1 - start;
	seed.seedSet.front().i3 += lengthDiff;
	seed.seedSet.front().i2 -= lengthDiff;
	seed.seedSet.front().i1 = start;
}


template<typename TPosition>
inline void 
setRightDim0(Seed<TPosition,ChainedSeed> & seed, 
			TPosition end)
{
	SEQAN_CHECKPOINT
	seed.seedSet.back().i3 = end - seed.seedSet.back().i1+1;
}


template<typename TPosition>
inline void 
setLeftDim1(Seed<TPosition, ChainedSeed> &seed, 
				 TPosition start)
{
	SEQAN_CHECKPOINT
	TPosition lengthDiff = seed.seedSet.front().i2 - start;
	seed.seedSet.front().i3 += lengthDiff;
	seed.seedSet.front().i1 -= lengthDiff;
	seed.seedSet.front().i2 = start;
}


template<typename TPosition>
inline void 
setRightDim1(Seed<TPosition,ChainedSeed> & seed, 
			   TPosition end)
{
	SEQAN_CHECKPOINT
	seed.seedSet.back().i3 = end - seed.seedSet.back().i2+1;
}

/*
.Function._getDiagSet:
..summary: Returns the set of matching diagonals.
..cat:Seed Handling
..signature:setRightDim1(seed)
..param.seed: The seed whose end position2 should be updated.
...type:Spec.ChainedSeed
..returns: A reference to the list of seed diagonals.
*/
template<typename TPosition>
inline const std::list<Triple<TPosition, TPosition, TPosition> >&
_getDiagSet(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet;
}

template<typename TPosition>
inline std::list<Triple<TPosition, TPosition, TPosition> >&
_getDiagSet(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet;
}

/**
.Function.appendDiag
..summary: Adds diagonal to the seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.diag: The diagonal to add.
...type:Class.Triple
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
*/
template<typename TPosition>
void
appendDiag(Seed<TPosition,ChainedSeed> & seed, 
		   Triple<TPosition, TPosition, TPosition> diag)
{
	SEQAN_CHECKPOINT
	seed.seedSet.push_back(diag);
}

template<typename TPosition>
Triple<TPosition, TPosition, TPosition>&
_getFirstDiag(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.front();
}

template<typename TPosition>
const Triple<TPosition, TPosition, TPosition>&
_getFirstDiag(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.front();
}

template<typename TPosition>
Triple<TPosition, TPosition, TPosition>&
_getLastDiag(Seed<TPosition,ChainedSeed> & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.back();
}

template<typename TPosition>
const Triple<TPosition, TPosition, TPosition>&
_getLastDiag(Seed<TPosition,ChainedSeed> const & seed){
	SEQAN_CHECKPOINT
	return seed.seedSet.back();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Merge Alogrithms                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, ChainedSeed> &firstSeed,
			   TPosition qPos,
			   TPosition dPos,
			   TPosition length,
			   Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator TIterator;
    TIterator begin1, end2, it;
	TPosition diag = dPos -qPos;
	//new seed would be longer?
	if (qPos+length-1 > rightDim0(firstSeed)){
		while ((_getLastDiag(firstSeed)).i1 > qPos)// || ((_getLastDiag(firstSeed)).i2 > dPos))
		{
			(firstSeed).seedSet.pop_back();
		}
		if ((rightDim0(firstSeed) < qPos) && (rightDim1(firstSeed) < dPos)) {
			appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
		} else {
			if (diag == endDiagonal(firstSeed)){
				setRightDim1(firstSeed,dPos+length-1);
			} else {
				TPosition tmp = diag - endDiagonal(firstSeed);
				if (tmp < 0){
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- dPos+1);
					if (tmp2 > 0){
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				} else {
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - qPos+1);
					if (tmp2 > 0){
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
				}
			}
		}
	}

	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, ChainedSeed> &firstSeed, 
			   Seed<TPosition, ChainedSeed> const &secondSeed, 
			   Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::const_iterator TIterator;
        TIterator begin1, end2, it;
	begin1 = _getDiagSet(secondSeed).begin();
	end2 = _getDiagSet(secondSeed).end();
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > _getFirstDiag(secondSeed).i1) //|| ((_getLastDiag(firstSeed).i2 > _getFirstDiag(secondSeed).i2)))
		{
			(firstSeed).seedSet.pop_back();
		}

		if ((rightDim0(firstSeed) < leftDim0(secondSeed)) && (rightDim1(firstSeed) < leftDim1(secondSeed))) {
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		} else {
		if (startDiagonal(secondSeed) == endDiagonal(firstSeed)){
			setRightDim1(firstSeed,(*begin1).i2+(*begin1).i3-1);
			++begin1;
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		}
		else {
			
			TPosition tmp = startDiagonal(secondSeed) - endDiagonal(firstSeed);
			if (tmp < 0){
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- leftDim1(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			} else {
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - leftDim0(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			}
		}	
		}
		
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}

template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void
_mergeTwoSeedsScore(Seed<TPosition, ChainedSeed> &firstSeed,
					TPosition3 &score1,
					TPosition qPos,
					TPosition dPos,
					TPosition length,
					TPosition3 score2,
					Score<TPosition2,Simple> const &scoreMatrix,
					TGapCost &,
				    Merge)
{
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator TIterator;
    TIterator begin1, end2, it;
	TPosition diag = dPos - qPos;
	score1 += score2;

	
	if (qPos+length-1 > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > qPos)// || (_getLastDiag(firstSeed).i2 > dPos))
		{
			score1 -= firstSeed.seedSet.back().i3*scoreMatch(scoreMatrix);
			TPosition x1 = firstSeed.seedSet.back().i1;
			TPosition x2 = firstSeed.seedSet.back().i2;
			firstSeed.seedSet.pop_back();
			score1 -=(abs(rightDim0(firstSeed)-x1)+abs(rightDim1(firstSeed)-x2))*scoreGap(scoreMatrix); 
		}
		score1 += abs(endDiagonal(firstSeed) - dPos + qPos)*scoreGap(scoreMatrix);
		score1 -=(max(abs(rightDim0(firstSeed)- qPos),abs(rightDim1(firstSeed)-dPos))+1)*scoreMatch(scoreMatrix);

		if ((rightDim0(firstSeed) < qPos) && (rightDim1(firstSeed) < dPos)) {
			appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
		} else {
			if (diag == endDiagonal(firstSeed))
			{
				setRightDim1(firstSeed,dPos+length-1);
			} 
			else 
			{
				TPosition tmp = diag - endDiagonal(firstSeed);
				if (tmp < 0)
				{
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- dPos+1);
					if (tmp2 > 0)
					{
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				} 
				else 
				{
					Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
					TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - qPos+1);
					if (tmp2 > 0)
					{
						_getLastDiag(firstSeed).i3 = tmp2;
						appendDiag(firstSeed,Triple<TPosition, TPosition, TPosition>(qPos, dPos, length));
					}
				}
			}
		}
	if (leftDiagonal(firstSeed) < diag)
		setLeftDiagonal(firstSeed, diag);
	if (rightDiagonal(firstSeed) > diag)
		setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void
_mergeTwoSeedsScore(Seed<TPosition, ChainedSeed> &firstSeed,
					TPosition3 &score1,
					Seed<TPosition, ChainedSeed> const &secondSeed,
					TPosition3 score2,
					Score<TPosition2,Simple> const &scoreMatrix,
					TGapCost &,
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	typedef typename std::list<Triple <TPosition, TPosition, TPosition> >::const_iterator TIterator;
        TIterator begin1, end2, it;
	begin1 = _getDiagSet(secondSeed).begin();
	end2 = _getDiagSet(secondSeed).end();
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		while (_getLastDiag(firstSeed).i1 > _getFirstDiag(secondSeed).i1)// || (_getLastDiag(firstSeed).i2 > _getFirstDiag(secondSeed).i2))
		{
			score1 -= firstSeed.seedSet.back().i3*scoreMatch(scoreMatrix);
			TPosition x1 = firstSeed.seedSet.back().i1;
			TPosition x2 = firstSeed.seedSet.back().i2;
			firstSeed.seedSet.pop_back();
			score1 -= (abs(rightDim0(firstSeed)-x1)+abs(rightDim1(firstSeed)-x2))*scoreGap(scoreMatrix); 
		}
		score1 += abs(endDiagonal(firstSeed) - startDiagonal(secondSeed))*scoreGap(scoreMatrix);
		score1 -= (max(abs(rightDim0(firstSeed)- leftDim0(secondSeed)),abs(rightDim1(firstSeed)-leftDim1(secondSeed)))+1)*scoreMatch(scoreMatrix);

		if ((rightDim0(firstSeed) < leftDim0(secondSeed)) && (rightDim1(firstSeed) < leftDim1(secondSeed))) {
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		} else {
		if (startDiagonal(secondSeed) == endDiagonal(firstSeed)){
			setRightDim1(firstSeed,(*begin1).i2+(*begin1).i3-1);
			++begin1;
			for (it = begin1; it != end2; it++){
				appendDiag(firstSeed,*it);
			}
		}
		else {
			
			TPosition tmp = startDiagonal(secondSeed) - endDiagonal(firstSeed);
			if (tmp < 0){
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2  = tmpSeed.i3 - ((tmpSeed.i2+tmpSeed.i3-1)- leftDim1(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			} else {
				Triple <TPosition, TPosition, TPosition> tmpSeed = _getLastDiag(firstSeed);
				TPosition tmp2 = tmpSeed.i3 -(tmpSeed.i1 +tmpSeed.i3-1 - leftDim0(secondSeed)+1);
				if (tmp2 > 0){
					_getLastDiag(firstSeed).i3 = tmp2;
					for (it = begin1; it != end2; it++){
						appendDiag(firstSeed,*it);
					}
				}
			}
		}	
		}
		
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}

template<typename TPosition, typename TText, typename TTPosition>
void 
extendSeed(Seed<TPosition,ChainedSeed> &seed, TPosition scoreDropOff, Score<TTPosition, Simple> const &scoreMatrix, String<TText> &query, String<TText> &database, TPosition direction, GappedXDrop){
	SEQAN_CHECKPOINT
	TPosition gapCost = scoreGap(scoreMatrix);
	//TPosition tmpScore = 0;
	TPosition infimum = infimumValue<TPosition>()+1-gapCost;
	
	//left extension
	if ((direction != 1)&&(leftDim0(seed)!=0)&&(leftDim1(seed)!=0)){
		TPosition upperBound = 0;
		TPosition lowerBound = 0;
		Segment<String<Dna>,PrefixSegment> dataSeg(database,leftDim1(seed));
		Segment<String<Dna>,PrefixSegment> querySeg(query,leftDim0(seed));
		TPosition xLength = length(querySeg);
		TPosition yLength = length(dataSeg);

		std::vector<TPosition> *antiDiag1 = new std::vector<TPosition>(1,0);		//smallest diagonal
		std::vector<TPosition> *antiDiag2 = new std::vector<TPosition>(2,infimum);
		std::vector<TPosition> *antiDiag3 = new std::vector<TPosition>(3,infimum);	//current diagonal
		std::vector<TPosition> *tmpDiag;

		//Matrix initialization
		if (gapCost >= (-1)*scoreDropOff){
			(*antiDiag2)[0] = gapCost;
			(*antiDiag2)[1] = gapCost;
		}
		if (2*gapCost >= (-1)*scoreDropOff){
			(*antiDiag3)[0] = 2*gapCost;
			(*antiDiag3)[2] = 2*gapCost;
		}
		
		TPosition b = 1; 
		TPosition u = 0;
		TPosition k = 1;
		TPosition tmp;
		TPosition tmpMax1 = 0;
		TPosition tmpMax2 = 0;
		
		//Extension as proposed by Zhang et al
		while(b<=u+1){
			++k;
			for (int i = b; i<= (u+1);++i){
				tmp = infimum;

				tmp = max((*antiDiag2)[i-1],(*antiDiag2)[i])+gapCost;
				tmp = max(tmp,(*antiDiag1)[i-1]+ score(scoreMatrix,xLength-i,yLength-(k-i),querySeg,dataSeg));
				tmpMax2 = max(tmpMax2,tmp);
				if (tmp < tmpMax1-scoreDropOff)
					(*antiDiag3)[i] = infimum;
				else
					(*antiDiag3)[i] = tmp;
			}
			while (((*antiDiag3)[b]  < tmpMax1-scoreDropOff) && (b<static_cast<TPosition>((*antiDiag3).size())-1)){
				
				++b;}
			++u;
			while (((*antiDiag3)[u]  < tmpMax1-scoreDropOff) && (u>0)){
				--u;}
			
			//borders for lower triangle of edit matrix
			b = max(b,k-yLength+1);
			u = min(u, xLength-1);
			
			if ((b < (k+1)/2)&&((k+1)/2-b>lowerBound))
				lowerBound = (k+1)/2-b;
			if ((u > k/2)&&(u-k/2>upperBound))
				upperBound = u-k/2;

			tmpDiag = antiDiag1;
			antiDiag1 = antiDiag2;
			antiDiag2 = antiDiag3;
			antiDiag3 = tmpDiag;

			int d = 0;
			while ((d<3) &&(static_cast<TPosition>((*antiDiag3).size())<=xLength)){
				(*antiDiag3).push_back(0);
				++d;
			}
			for (unsigned int eu = 0; eu < (*antiDiag3).size();++eu)
				(*antiDiag3)[eu] = infimum;

			if ((*antiDiag2)[0]+ gapCost >= tmpMax1-scoreDropOff)
				(*antiDiag3)[0] = (*antiDiag2)[0]+ gapCost;
			if ((*antiDiag2)[(*antiDiag2).size()-1]+ gapCost >=tmpMax1-scoreDropOff)
				(*antiDiag3)[(*antiDiag3).size()-1]=(*antiDiag2)[(*antiDiag2).size()-1]+ gapCost;
				tmpMax1 = tmpMax2;
		}
		
		//Calculate upper/lower bound for diagonals
		if (rightDiagonal(seed) > endDiagonal(seed)-upperBound)
			setRightDiagonal(seed, endDiagonal(seed)-upperBound);
	
		if (leftDiagonal(seed) < endDiagonal(seed)+lowerBound)
			setLeftDiagonal(seed, endDiagonal(seed)+lowerBound);
		


		// Find seed start
		TPosition tmpPos = 0;
		TPosition tmpMax = infimum;
		if ((k==xLength+yLength) &&((*antiDiag2)[xLength] >= tmpMax1-scoreDropOff)){
			tmpPos = xLength;
			tmpMax =0;
		} else{
			for (unsigned int eu = 0; eu < (*antiDiag1).size();++eu){
				if ((*antiDiag1)[eu] > tmpMax){
					tmpMax = (*antiDiag1)[eu];
					tmpPos = eu;		
				}
			}
			--k;
		}
		if(tmpMax != infimum){
			seed.seedSet.push_front(Triple<TPosition, TPosition, TPosition>(leftDim0(seed)-tmpPos,leftDim1(seed)-(k-tmpPos),1));
		}


		//free memory
		(*antiDiag1).clear();
		(*antiDiag2).clear();
		(*antiDiag3).clear();
		delete antiDiag1;
		delete antiDiag2;
		delete antiDiag3;
	}

	//right extension
	
	if ((direction != 0)&&(rightDim0(seed)< static_cast<TPosition>(length(query))-1)&&(rightDim1(seed)<static_cast<TPosition>(length(database))-1)){
		TPosition upperBound = 0;
		TPosition lowerBound = 0;
		Segment<String<Dna>,SuffixSegment> dataSeg(database,rightDim1(seed)+1);
		Segment<String<Dna>,SuffixSegment> querySeg(query,rightDim0(seed)+1);
		TPosition xLength = length(querySeg);
		TPosition yLength = length(dataSeg);

		std::vector<TPosition> *antiDiag1 = new std::vector<TPosition>(1,0);	//smallest diagonal
		std::vector<TPosition> *antiDiag2 = new std::vector<TPosition>(2,infimum);
		std::vector<TPosition> *antiDiag3 = new std::vector<TPosition>(3,infimum);	//current diagonal
		std::vector<TPosition> *tmpDiag;
		
		//Matrix initialization
		if (gapCost >= (-1)*scoreDropOff){
			(*antiDiag2)[0] = gapCost;
			(*antiDiag2)[1] = gapCost;
		}
		if (2*gapCost >= (-1)*scoreDropOff){
			(*antiDiag3)[0] = 2*gapCost;
			(*antiDiag3)[2] = 2*gapCost;
		}
		
		TPosition b = 1; 
		TPosition u = 0;
		TPosition k = 1;
		TPosition tmp;
		TPosition tmpMax1 = 0;
		TPosition tmpMax2 = 0;
	
		//Extension as proposed by Zhang
		while(b<=u+1){
			++k;
			for (int i = b; i<= (u+1);++i){
				tmp = infimum;
				tmp = max((*antiDiag2)[i-1],(*antiDiag2)[i])+gapCost;
				tmp = max(tmp,(*antiDiag1)[i-1]+ score(scoreMatrix,i-1,k-i-1,querySeg,dataSeg));
				tmpMax2 = max(tmpMax2,tmp);
				if (tmp < tmpMax1-scoreDropOff)
					(*antiDiag3)[i] = infimum;
				else
					(*antiDiag3)[i] = tmp;
			}
		
			while (((*antiDiag3)[b]  < tmpMax1-scoreDropOff) && (b< static_cast<TPosition>((*antiDiag3).size())-1)){
				++b;
			}
			++u;
			while (((*antiDiag3)[u]  < tmpMax1-scoreDropOff) && (u>0)){
				--u;}
			
			//borders for lower triangle of edit matrix
			b = max(b,k-yLength+1);
			u = min(u, xLength-1);
			

			if ((b < (k+1)/2)&&((k+1)/2-b>lowerBound)){
				lowerBound = (k+1)/2-b;
			}
		
			if ((u >= k/2)&&(u-k/2>=upperBound)){
				upperBound = u-k/2;
			}
			tmpDiag = antiDiag1;
			antiDiag1 = antiDiag2;
			antiDiag2 = antiDiag3;
			antiDiag3 = tmpDiag;

			int d = 0;

			while ((d<3) &&(static_cast<TPosition>((*antiDiag3).size())<=xLength)){
				(*antiDiag3).push_back(0);
				++d;
			}
			for (unsigned int eu = 0; eu < (*antiDiag3).size();++eu)
				(*antiDiag3)[eu] = infimum;

			if ((*antiDiag2)[0]+ gapCost >= tmpMax1-scoreDropOff)
				(*antiDiag3)[0] = (*antiDiag2)[0]+ gapCost;
			if ((*antiDiag2)[(*antiDiag2).size()-1]+ gapCost >=tmpMax1-scoreDropOff)
				(*antiDiag3)[(*antiDiag3).size()-1]=(*antiDiag2)[(*antiDiag2).size()-1]+ gapCost;
			tmpMax1 = tmpMax2;
		}
		
		//Calculate upper/lower bound for diagonals
		if (rightDiagonal(seed) > endDiagonal(seed)-upperBound)
			setRightDiagonal(seed, endDiagonal(seed)-upperBound);
	
		if (leftDiagonal(seed) < endDiagonal(seed)+lowerBound)
			setLeftDiagonal(seed, endDiagonal(seed)+lowerBound);

		//Find seed end
		TPosition tmpPos = 0;
		TPosition tmpMax = infimum;
		if ((k==xLength+yLength) && ((*antiDiag2)[xLength] >= tmpMax1-scoreDropOff)){
			tmpPos = xLength;
			tmpMax =0;
		} else{
			for (size_t eu = 0; eu < antiDiag1->size(); ++eu){
				if ((*antiDiag1)[eu] > tmpMax){
					tmpMax = (*antiDiag1)[eu];
					tmpPos = eu;
					
				}
			}
			--k;
		}

		if(tmpMax != infimum){
			seed.seedSet.push_back(Triple<TPosition, TPosition, TPosition>(rightDim0(seed)+tmpPos,rightDim1(seed)+k-tmpPos,1));
		}
		//free memory
		(*antiDiag1).clear();
		(*antiDiag2).clear();
		(*antiDiag3).clear();
		delete antiDiag1;
		delete antiDiag2;
		delete antiDiag3;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Alignment Construction Alogrithm										  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
.Function.getAlignment:
..summary: Constructs a alignment from a @Spec.ChainedSeed@.
..cat:Seed Handling
..signature:getAlignment(seed, align, query, database, scoreMatrix)
..param.seed: The alignment foundation.
...type:Spec.ChainedSeed
..param.align: An emtpy alignment object, that stores the constructed alignment.
...type:Class.Align
..param.query:The Query sequence.
...type:Class.String
..param.database:The database sequence.
...type:Class.String
..param.scoreMatrix:The scoring matrix.
...type:Spec.Simple Score
..returns: Score of the alignment.
*/
template<typename TPosition, typename TText, typename TPosition2>
int
getAlignment(Seed<TPosition,ChainedSeed> &seed,
			 Align<String<TText>, ArrayGaps> &aligned, 
			 String<TText> &query, 
			 String<TText> &database, 
			 Score<TPosition2, Simple> &scoreMatrix)
{
	SEQAN_CHECKPOINT
	int seedScore = 0;
	typename std::list<Triple< TPosition, TPosition, TPosition> > seedList = _getDiagSet(seed);
	typedef typename std::list<Triple< TPosition, TPosition, TPosition> >::iterator TIterator;

	resize(rows(aligned), 2);
	assignSource(row(aligned, 0), query, leftDim0(seed), rightDim0(seed)+1);
	assignSource(row(aligned, 1), database, leftDim1(seed), rightDim1(seed)+1);
	TIterator it1 = seedList.begin();
	TIterator it2 = ++seedList.begin();
	
	for (int i =0; i<(*it1).i3;++i){
		seedScore += score(scoreMatrix,(*it1).i1+i, (*it1).i2+i, query, database); 
	}
	
	if (seedList.size()>=2){
		TPosition gapLength;
		TPosition position1 =	it1->i1 - leftDim0(seed) + it1->i3;
		TPosition position2 =	it1->i2 - leftDim1(seed) + it1->i3;
		while (it2 != seedList.end()){
			if (it2->i1 == it1->i1 + it1->i3){ //query teile zusammen
				gapLength = it2->i2 - it1->i2 - it1->i3;
				insertGaps(row(aligned,0),position1, gapLength);
				seedScore += gapLength*scoreGap(scoreMatrix);
				position1 += gapLength;
			}
			else
				if ((*it2).i2 == it1->i2+it1->i3)
				{
					gapLength = it2->i1 - it1->i1 - it1->i3;
					insertGaps(row(aligned,1),position2,gapLength);
					position2 += gapLength;
					seedScore += gapLength*scoreGap(scoreMatrix);
				} 
				else 
				{
					Align<String<TText>, ArrayGaps> alignSeg;
					resize(rows(alignSeg), 2);
					assignSource(row(alignSeg, 0), query, (*it1).i1+(*it1).i3, (*it2).i1);
					assignSource(row(alignSeg, 1), database, (*it1).i2+(*it1).i3, (*it2).i2);
					seedScore += globalAlignment(alignSeg,scoreMatrix,NeedlemanWunsch());//needlemanWunsch(alignSeg,scoreMatrix);
			
					unsigned int j;
					bool gap;
					if (row(alignSeg,1).data_arr[0] == 0){
						j=1;
						gap=false;
					} else {
						j=0;
						gap=true;
					}
					while (j < length(row(alignSeg,1).data_arr)){
						if (gap){
							insertGaps(row(aligned,1),position2,row(alignSeg,1).data_arr[j]);
					
							gap = false;
							position2 += row(alignSeg,1).data_arr[j];
						} else {
							gap = true;
							position2 += row(alignSeg,1).data_arr[j];							
						}
						++j;
					}
					if (row(alignSeg,0).data_arr[0] == 0){
						j=1;
						gap=false;
					} else {
						j=0;
						gap=true;
					}
					 while (j < length(row(alignSeg,0).data_arr)){
						if (gap){
							insertGaps(row(aligned,0),position1,row(alignSeg,0).data_arr[j]);
							gap = false;
							position1+= row(alignSeg,0).data_arr[j];
						} else {
							gap = true;
							position1+= row(alignSeg,0).data_arr[j];							
						}
						++j;
					}
					TPosition tmp1, tmp2;
					tmp1 = position1;
					tmp2 = position2;
					if (tmp1 > tmp2){
						insertGaps(row(aligned,1),position2,tmp1-tmp2);
						position2 += tmp1-tmp2;
					} else 
						if (tmp2 > tmp1){
							insertGaps(row(aligned,0),position1,tmp2-tmp1);
							position1 += tmp2-tmp1;
						}
				}
				position1+= it2->i3;
				position2 += it2->i3;
				for (int i =0; i<it2->i3;++i)
					seedScore += score(scoreMatrix, it2->i1+i, it2->i2+i, query, database);

				++it1;
				++it2;
		}
	}
	return seedScore;
}


/**
.Function.scoreSeed:
..summary: Calculates the score of a seed. 
..cat:Seed Handling
..signature:scoreSeed(seed, query, database, scoreMatrix)
..param.seed: A seed.
...type:Spec.ChainedSeed
..param.query:The Query sequence.
...type:Class.String
..param.database:The database sequence.
...type:Class.String
..param.scoreMatrix:The scoring sheme.
...type:Spec.Simple Score
..returns: Score of the seed.
..remarks: Score has not the same value as the resulting alignment. Gaps between diagonals matches are scored as full length gaps.
*/
template<typename TPosition, typename TText, typename TScore>
TScore
scoreSeed(Seed<TPosition, ChainedSeed> &seed, String<TText> &query, String<TText> &database, Score<TScore, Simple> &matrix){
	SEQAN_CHECKPOINT
	typedef typename std::list<Triple< TPosition, TPosition, TPosition> >::iterator TIterator;
	int tmpScore =0;
	TIterator it = _getDiagSet(seed).begin();
	for (int i = 0; i < it->i3; ++i){
		tmpScore+=score(matrix, it->i1+i, it->i2+i, query, database);
	}

	if (_getDiagSet(seed).size()>=2){
		TIterator it_end = _getDiagSet(seed).end();
		for (TIterator it2 = ++_getDiagSet(seed).begin(); it2!= it_end; it2++){
			for (int i = 0; i < it2->i3; ++i){
				tmpScore+=score(matrix, it2->i1+i, it2->i2+i, query, database);
			}
			tmpScore += scoreGap(matrix)*(it2->i2-(it->i2+it->i3)+it2->i1-(it->i1+it->i3));
			++it;
		}
	}
	return tmpScore;
}


} //namespace Seqan

#endif //#ifndef SEQAN_HEADER_
