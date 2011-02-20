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
  $Id: seed_base.h 3739 2009-03-23 13:50:07Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/


#ifndef SEQAN_HEADER_SEED_H
#define SEQAN_HEADER_SEED_H

namespace SEQAN_NAMESPACE_MAIN
{


struct _Seed_simple;
typedef Tag<_Seed_simple> const SimpleSeed;



/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..see:Function.extendSeeds
..tag.MatchExtend:
	Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:
	Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:
	Gapped extension of a seed until score drops below a Value. Only @Spec.SimpleSeed@s.
*/


/**
.Tag.Seed Adding.tag.Merge:
	Merging of Seeds.
*/
struct _Chain_Merge;
typedef Tag<_Chain_Merge> const Merge;



struct _extendSeed_Match;
typedef Tag<_extendSeed_Match> const MatchExtend;


struct _extendSeed_UnGappedXDrop;
typedef Tag<_extendSeed_UnGappedXDrop> const UngappedXDrop;

struct _extendSeed_GappedXDrop;
typedef Tag<_extendSeed_GappedXDrop> const GappedXDrop;

//template<typename TPosition = int, typename TSpecSeed = SimpleSeed>class Seed;



/**
.Class.Seed:
..summary:Describes a seed.
..cat:Seed Handling
..signature:Seed<TPosition, TSpecSeed>
..param.TPosition:The type number that schuld be used. Must have negative numbers (e.g. int/long).
..param.TSpec:The seed type used.
*/

/**
.Spec.SimpleSeed:
..summary:Describes a seed with start and end position and diagonal upper and lower bounds.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, SimpleSeed>
..param.TPosition:The type number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.SimpleSeed#Seed:
..class:Spec.SimpleSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, qEndPos, dEndPos)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.qEndPos: End in query sequence.
..param.dEndPos: End in database sequence.
..param.length: Length of the seed.
*/
template<typename TPosition = int, typename TSpecSeed = SimpleSeed> 
class Seed{
public:
	TPosition leftDim0;
	TPosition leftDim1;
	TPosition rightDim0;
	TPosition rightDim1;
	TPosition leftDiagonal;
	TPosition rightDiagonal;

	Seed(){
		SEQAN_CHECKPOINT
	}

	Seed(TPosition leftDim0, TPosition leftDim1, TPosition length):leftDim0(leftDim0),leftDim1(leftDim1){
		SEQAN_CHECKPOINT
		rightDim0 = leftDim0 + length-1;
		rightDim1 = leftDim1 + length-1;
		rightDiagonal = leftDiagonal = leftDim1-leftDim0;
	}

	Seed(TPosition leftDim0, TPosition leftDim1, TPosition rightDim0, TPosition rightDim1):leftDim0(leftDim0),leftDim1(leftDim1),rightDim0(rightDim0), rightDim1(rightDim1){
		SEQAN_CHECKPOINT
		leftDiagonal = max(leftDim1 - leftDim0, rightDim1-rightDim0);
		rightDiagonal = min(leftDim1 - leftDim0, rightDim1-rightDim0);
	}


	~Seed(){
	}

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											Meta Functions		                                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///.Metafunction.Spec.param.T.type:Class.Seed

template <typename TPosition, typename TSpecSeed>
struct Spec<Seed<TPosition,TSpecSeed> >
{
	typedef TSpecSeed Type;
};

///.Metafunction.Value.param.T.type:Class.Seed
template <typename TPosition, typename TSpecSeed>
struct Value<Seed<TPosition,TSpecSeed> >
{
	typedef TPosition Type;
};



template< typename TBorder, typename TSpec >
struct Size< Seed< TBorder, TSpec > >
{
	typedef size_t Type;
};

template< typename TPosition, typename TSpec >
struct Key< Seed< TPosition, TSpec > >
{
	typedef TPosition Type;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Standard Functions                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.startDiagonal:
..summary: Returns the diagonal of the start point.
..cat:Seed Handling
..signature:startDiagonal(seed)
..param.seed: The seed whose start diagonal should be returned.
...type:Class.Seed
..returns: The diagonal of the start point.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
startDiagonal(Seed<TPosition, TSpecSeed> const &me)
{
	SEQAN_CHECKPOINT
	return me.leftDim1-me.leftDim0;
}

/**
.Function.endDiagonal:
..summary: Returns the diagonal of the end point.
..cat:Seed Handling
..signature:endDiagonal(seed)
..param.seed: The seed whose end diagonal should be returned.
...type:Class.Seed
..returns: The diagonal of the end point.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
endDiagonal(Seed<TPosition, TSpecSeed> const &me)
{
	SEQAN_CHECKPOINT
	return me.rightDim1-me.rightDim0;
}



/**.Function.leftPosition:
..summary:The begin position of segment in a seed.
..cat:Seed Handling
..signature:leftPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:Begin position of the $dim$-th segment in $seed$.
*/
template< typename TPosition, typename TSpecSeed, typename TSize> 
inline TPosition 
leftPosition(Seed<TPosition, TSpecSeed>  & me, 
			 TSize dim)
{
	SEQAN_CHECKPOINT
	return (dim)? me.leftDim1 : me.leftDim0;
}

/**.Function.setLeftPosition:
..summary:Sets begin position of segment in a seed.
..cat:Seed Handling
..signature:setLeftPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new begin position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
*/
template< typename TPosition, typename TSpecSeed, typename TSize, typename TPosition2> 
inline TPosition 
setLeftPosition(Seed<TPosition, TSpecSeed>  & me, 
				TSize dim,
				TPosition2 new_pos)
{
	SEQAN_CHECKPOINT
	if (dim) me.leftDim1 = new_pos;
	else me.leftDim0 = new_pos;
}

/**.Function.rightPosition:
..summary:The end position of segment in a seed.
..cat:Seed Handling
..signature:rightPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:End position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
*/

template< typename TPosition, typename TSpecSeed, typename TSize> 
inline TPosition 
rightPosition(Seed<TPosition, TSpecSeed>  & me, 
			  TSize dim)
{
	SEQAN_CHECKPOINT
	return (dim)? me.rightDim1 : me.rightDim0;
}

/**.Function.setRightPosition:
..summary:Sets end position of segment in a seed.
..cat:Seed Handling
..signature:setRightPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new end position of the $dim$-th segment in $seed$.
..see:Function.rightPosition
..see:Function.setLeftPosition
*/
template< typename TPosition, typename TSpecSeed, typename TSize, typename TPosition2> 
inline TPosition 
setRightPosition(Seed<TPosition, TSpecSeed>  & me, 
				TSize dim,
				TPosition2 new_pos)
{
	SEQAN_CHECKPOINT
	if (dim) me.rightDim1 = new_pos;
	else me.rightDim0 = new_pos;
}

/**.Function.dimension:
..summary:Dimension of a seed.
..cat:Seed Handling
..signature:dimension(seed)
..param.seed:A seed.
...type:Class.Seed
..returns:The number of segments in $seed$.
*/
template< typename TPosition, typename TSpec > inline
typename Size< Seed< TPosition, TSpec > >::Type
dimension( Seed< TPosition, TSpec > & me )
{
	return 2;
}

/**
.Function.leftDim0:
..summary: Returns the first position of the seed in the query.
..cat:Seed Handling
..signature:leftDim0(seed)
..param.seed: The seed whose query position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDim0(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDim0;
}

/**
.Function.rightDim0:
..summary: Returns the last position of the seed in the query.
..cat:Seed Handling
..signature:rightDim0(seed)
..param.seed: The seed whose last in the query position should be returned.
...type:Class.Seed
..returns: End of the seed.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDim0(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim0;
}

/**
.Function.leftDim1:
..summary: Returns the first position of the seed in the database.
..cat:Seed Handling
..signature:leftDim1(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDim1(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDim1;
}

/**
.Function.rightDim1:
..summary: Returns the last position of the seed in the database.
..cat:Seed Handling
..signature:rightDim1(seed)
..param.seed: The seed whose last in the database position should be returned.
...type:Class.Seed
..returns: End of the seed.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDim1(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim1;
}

/**
.Function.leftDiagonal:
..summary: Returns the most left diagonal of the seed (maximum diagonal value).
..cat:Seed Handling
..signature:leftDiagonal(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: The most left diagonal.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
leftDiagonal(Seed<TPosition, TSpecSeed> const &seed)
{
	SEQAN_CHECKPOINT
	return seed.leftDiagonal;
}

/**
.Function.rightDiagonal:
..summary: Returns the most right diagonal of the seed (minimum diagonal value).
..cat:Seed Handling
..signature:rightDiagonal(seed)
..param.seed: The seed whose database position should be returned.
...type:Class.Seed
..returns: The most right diagonal.
*/
template<typename TPosition, typename TSpecSeed>
inline TPosition 
rightDiagonal(Seed<TPosition,TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDiagonal;
}

///.Function.length.param.object.type:Class.Gaps
template<typename TPosition, typename TSpecSeed>
inline TPosition 
length(Seed<TPosition, TSpecSeed> const & seed)
{
	SEQAN_CHECKPOINT
	return seed.rightDim0-seed.leftDim0+1;
}


/**
.Function.setLeftDim0:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim0(seed, start)
..param.seed: The seed whose start position should be updated.
...type:Class.Seed
..param.start: The query position where the seed should start.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDim0(Seed<TPosition, TSpecSeed> &me, 
			  TPosition start)
{
	SEQAN_CHECKPOINT
	me.leftDim0 = start;
}

/**
.Function.setRightDim0:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim0(seed, end)
..param.seed: The seed whose end position should be updated.
...type:Class.Seed
..param.end: The query position where the seed should end.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDim0(Seed<TPosition,TSpecSeed> & me, 
			TPosition end)
{
	SEQAN_CHECKPOINT
	me.rightDim0 = end;
}

/**
.Function.setLeftDim1:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim1(seed, start)
..param.seed: The seed whose start position should be updated.
...type:Class.Seed
..param.start: The database position where the seed should start.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDim1(Seed<TPosition, TSpecSeed> &me, 
				 TPosition start)
{
	SEQAN_CHECKPOINT
	me.leftDim1 = start;
}

/**
.Function.setRightDim1:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim1(seed, end)
..param.seed: The seed whose end position should be updated.
...type:Class.Seed
..param.end: The database position where the seed should end.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDim1(Seed<TPosition,TSpecSeed> & me, 
			   TPosition end)
{
	SEQAN_CHECKPOINT
	me.rightDim1 = end;
}

/**
.Function.setLeftDiagonal:
..summary: Sets a new value for the most left diagonal.
..cat:Seed Handling
..signature:setLeftDiagonal(seed, diag)
..param.seed: The seed whose left diagonal value should be updated.
...type:Class.Seed
..param.diag: The new value for the most left diagonal.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setLeftDiagonal(Seed<TPosition, TSpecSeed> &me,
				TPosition diag)
{
	SEQAN_CHECKPOINT
	me.leftDiagonal = diag;
}

/**
.Function.setRightDiagonal:
..summary: Sets a new value for the most right diagonal.
..cat:Seed Handling
..signature:setRightDiagonal(seed, diag)
..param.seed: The seed whose right diagonal value should be updated.
...type:Class.Seed
..param.diag: The new value for the most right diagonal.
*/
template<typename TPosition, typename TSpecSeed>
inline void 
setRightDiagonal(Seed<TPosition,TSpecSeed> & seed, 
				 TPosition diag)
{
	SEQAN_CHECKPOINT
	seed.rightDiagonal = diag;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                      Merge Alogrithms                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename TPosition, typename TSpecSeed>
void 
_mergeTwoSeeds(Seed<TPosition, TSpecSeed> &firstSeed, 
			   Seed<TPosition, TSpecSeed> const &secondSeed, 
			   Merge)
{
	SEQAN_CHECKPOINT
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		setRightDim0(firstSeed,rightDim0(secondSeed));
		setRightDim1(firstSeed,rightDim1(secondSeed));
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}


template<typename TPosition, typename TSpecSeed>
void 
_mergeTwoSeeds(Seed<TPosition, TSpecSeed> &firstSeed, 
			   TPosition q,
			   TPosition d, 
			   TPosition l, 
			   Merge)
{
	SEQAN_CHECKPOINT
	setRightDim0(firstSeed,q+l-1);
	setRightDim1(firstSeed,d+l-1);
	if (leftDiagonal(firstSeed) < d-q)
		setLeftDiagonal(firstSeed, d-q);
	if (rightDiagonal(firstSeed) > d-q)
		setRightDiagonal(firstSeed, d-q);
}


template<typename TPosition>
void
_mergeTwoSeeds(Seed<TPosition, SimpleSeed> &firstSeed, 
			   TPosition qlPos, 
			   TPosition dlPos, 
			   TPosition qrPos, 
			   TPosition drPos, 
			   Merge)
{
	SEQAN_CHECKPOINT
	if (qrPos > rightDim0(firstSeed)){
		typename std::list<Triple <TPosition, TPosition, TPosition> >::iterator begin1, end2, it;
	
		setRightDim0(firstSeed,qrPos);
		setRightDim1(firstSeed,drPos);
	
		TPosition diag = dlPos -qlPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
		diag = drPos - qrPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
	}
}


template<typename TPosition, typename TSpecSeed, typename TPosition2, typename TPosition3, typename TGapCost>
void 
_mergeTwoSeedsScore(Seed<TPosition, TSpecSeed> &firstSeed, 
					TPosition3 &score1, 
					Seed<TPosition, TSpecSeed> const &secondSeed, 
					TPosition3 score2, 
					Score<TPosition2,Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	score1 += abs(endDiagonal(firstSeed)-startDiagonal(secondSeed))*scoreGap(scoreMatrix);
	score1 -= (max(abs(rightDim0(firstSeed)-leftDim0(secondSeed)),abs(rightDim1(firstSeed)-leftDim1(secondSeed)))+1)*scoreMatch(scoreMatrix);
	if (rightDim0(secondSeed) > rightDim0(firstSeed)){
		setRightDim0(firstSeed,rightDim0(secondSeed));
		setRightDim1(firstSeed,rightDim1(secondSeed));
		if (leftDiagonal(firstSeed) < leftDiagonal(secondSeed))
			setLeftDiagonal(firstSeed, leftDiagonal(secondSeed));
		if (rightDiagonal(firstSeed) > rightDiagonal(secondSeed))
			setRightDiagonal(firstSeed, rightDiagonal(secondSeed));
	}
}


template<typename TPosition, typename TPosition2, typename TSpecSeed, typename TPosition3, typename TGapCost>
void 
_mergeTwoSeedsScore(Seed<TPosition, TSpecSeed> &firstSeed, 
					TPosition3 &score1, 
					TPosition q, 
					TPosition d, 
					TPosition l, 
					TPosition3 score2, 
					Score<TPosition2, Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	score1 += score2;
	score1 += abs((TPosition2)(endDiagonal(firstSeed) - d + q))*scoreGap(scoreMatrix);
	score1 -= (::std::max<TPosition>(abs((TPosition2)(rightDim0(firstSeed)- q)),abs((TPosition2)(rightDim1(firstSeed)-d)))+1)*scoreMatch(scoreMatrix);
	setRightDim0(firstSeed,q+l-1);
	setRightDim1(firstSeed,d+l-1);
	if (leftDiagonal(firstSeed) < d-q)
		setLeftDiagonal(firstSeed, d-q);
	if (rightDiagonal(firstSeed) > d-q)
		setRightDiagonal(firstSeed, d-q);
}


template<typename TPosition, typename TPosition2, typename TPosition3, typename TGapCost>
void           
_mergeTwoSeedsScore(Seed<TPosition, SimpleSeed> &firstSeed, 
					TPosition3 &score1, 
					TPosition qlPos, 
					TPosition dlPos, 
					TPosition qrPos, 
					TPosition drPos, 
					TPosition3 score2, 
					Score<TPosition2,Simple> const &scoreMatrix, 
					TGapCost &, 
					Merge)
{
	SEQAN_CHECKPOINT
	if (qrPos > rightDim0(firstSeed)){
		
		score1 += score2;
		score1 += abs(endDiagonal(firstSeed) - dlPos + qlPos)*scoreGap(scoreMatrix);
		score1 -= (max(abs(rightDim0(firstSeed)- qlPos),abs(rightDim1(firstSeed)-dlPos))+1)*scoreMatch(scoreMatrix);

		setRightDim0(firstSeed,qrPos);
		setRightDim1(firstSeed,drPos);
		TPosition diag = dlPos -qlPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
		diag = drPos - qrPos;
		if (leftDiagonal(firstSeed) < diag)
			setLeftDiagonal(firstSeed, diag);
		if (rightDiagonal(firstSeed) > diag)
			setRightDiagonal(firstSeed, diag);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Extension Algorithms                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, tag)
..signature:extendSeed(seed, scoreDropOff, scoreMatrix, query, database, direction, tag)
..param.seed: The seed to extend.
...type:Class.Seed
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
*/


template<typename TPosition, typename TSpecSeed, typename TText>
void 
extendSeed(Seed<TPosition, TSpecSeed> &seed, 
		   String<TText> const &query, 
		   String<TText> const &database, 
		   TPosition direction, 
		   MatchExtend)
{
	SEQAN_CHECKPOINT
	//left extension
	if (direction != 1){
		TPosition queryPos =leftDim0(seed) ;
		TPosition dataPos = leftDim1(seed);
		while ((queryPos-1>=0) && (dataPos-1>=0) && (query[queryPos-1] == database[dataPos-1])){
			--queryPos;
			--dataPos;
		}
		setLeftDim0(seed,queryPos);
		setLeftDim1(seed,dataPos);

	}

	//right extension
	if (direction != 0){
		int queryLength = length(query);
		int databaseLength = length(database);
		TPosition queryPos =rightDim0(seed) ;
		TPosition dataPos = rightDim1(seed);
		while ((queryPos+1 < queryLength) && (dataPos+1 < databaseLength) && (query[queryPos+1] == database[dataPos+1])){
			++queryPos;
			++dataPos;
		}

		setRightDim0(seed,queryPos);
		setRightDim1(seed,dataPos);
	}
}


template<typename TPosition, typename TSpecSeed, typename TText, typename TScore>
void 
extendSeed(Seed<TPosition,TSpecSeed> &seed, 
		   TScore scoreDropOff, 
		   Score<TScore, Simple> const &scoreMatrix,
		   String<TText> const &query,
		   String<TText> const &database,
		   TPosition direction, 
		   UngappedXDrop)
{
	SEQAN_CHECKPOINT
	scoreDropOff *=-1;
	TPosition tmpScore = 0;

	//left extension
	if (direction != 1){
		TPosition xPos = leftDim0(seed)-1;
		TPosition yPos = leftDim1(seed)-1;
		TPosition last = 0;

		while ((tmpScore > scoreDropOff) && (xPos >= 0) && (yPos>=0)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmpScore += score(scoreMatrix, xPos, yPos, query, database);
				if (tmpScore > 0)
					tmpScore = 0;
			} else{
				tmpScore += score(scoreMatrix, xPos, yPos, query, database);
				++last;
			}
		--xPos;
		--yPos;
		}

		setLeftDim0(seed,xPos+last+1);
		setLeftDim1(seed,yPos+last+1);
	}

	//right extension
	if (direction != 0){
		TPosition xLength = length(query);
		TPosition yLength = length(database);
		TPosition xPos = rightDim0(seed)+1;
		TPosition yPos = rightDim1(seed)+1;

		TPosition last = 0;
		tmpScore= 0;
		while ((tmpScore > scoreDropOff) && (xPos < xLength) && (yPos < yLength)){
			if (query[xPos] == database[yPos]){
				last = 0;
				tmpScore += score(scoreMatrix, xPos, yPos, query, database);
				if (tmpScore > 0)
					tmpScore = 0;
			}else{
				tmpScore += score(scoreMatrix, xPos, yPos, query, database);
				++last;
			}
			++xPos;
			++yPos;
		}
		setRightDim0(seed,xPos-last-1);
		setRightDim1(seed,yPos-last-1);
	}
}


template<typename TPosition, typename TText, typename TScore>
void 
extendSeed(Seed<TPosition,SimpleSeed> &seed, 
		   TScore scoreDropOff, 
		   Score<TScore, Simple> const &scoreMatrix, 
		   String<TText> &query, 
		   String<TText> &database, 
		   TPosition direction, 
		   GappedXDrop)
{
	SEQAN_CHECKPOINT
	TPosition gapCost = scoreGap(scoreMatrix);
	//TPosition tmpScore = 0;
	TPosition infimum = infimumValue<TPosition>()+1-gapCost;
	
	//left extension
	if ((direction != 1)&&(leftDim0(seed)!=0)&&(leftDim1(seed)!=0)){
		TPosition upperBound = 0;
		TPosition lowerBound = 0;
		Segment<String<TText>,PrefixSegment> dataSeg(database,leftDim1(seed));
		Segment<String<TText>,PrefixSegment> querySeg(query,leftDim0(seed));
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

				tmp = ::std::max<TPosition>((*antiDiag2)[i-1],(*antiDiag2)[i])+gapCost;
				tmp = ::std::max<TPosition>(tmp,(*antiDiag1)[i-1]+ score(scoreMatrix, xLength-i, yLength-(k-i), querySeg, dataSeg));
				tmpMax2 = ::std::max<TPosition>(tmpMax2,tmp);
				if (tmp < tmpMax1-scoreDropOff)
					(*antiDiag3)[i] = infimum;
				else
					(*antiDiag3)[i] = tmp;
			}
			while (((*antiDiag3)[b]  < tmpMax1-scoreDropOff) && (b < static_cast<TPosition>((*antiDiag3).size())-1))
			{
				++b;
			}
			++u;
			while (((*antiDiag3)[u]  < tmpMax1-scoreDropOff) && (u>0)){
				--u;}
			
			//borders for lower triangle of edit matrix
			b = ::std::max<TPosition>(b,k-yLength+1);
			u = ::std::min<TPosition>(u, xLength-1);
			
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
			setLeftDim0(seed,leftDim0(seed)-tmpPos);
			setLeftDim1(seed,leftDim1(seed)-(k-tmpPos));
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
	
	if ((direction != 0)&&(rightDim0(seed)<static_cast<TPosition>(length(query))-1)&&(rightDim1(seed)<static_cast<TPosition>(length(database))-1)){
		TPosition upperBound = 0;
		TPosition lowerBound = 0;
		Segment<String<TText>,SuffixSegment> dataSeg(database,rightDim1(seed)+1);
		Segment<String<TText>,SuffixSegment> querySeg(query,rightDim0(seed)+1);
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
				tmp = ::std::max<TPosition>((*antiDiag2)[i-1],(*antiDiag2)[i])+gapCost;
				tmp = ::std::max<TPosition>(tmp,(*antiDiag1)[i-1]+ score(scoreMatrix,i-1,k-i-1,querySeg,dataSeg));
				tmpMax2 = ::std::max<TPosition>(tmpMax2,tmp);
				if (tmp < tmpMax1-scoreDropOff)
					(*antiDiag3)[i] = infimum;
				else
					(*antiDiag3)[i] = tmp;
			}
		
			while (((*antiDiag3)[b]  < tmpMax1-scoreDropOff) && (b<static_cast<TPosition>((*antiDiag3).size())-1)){
				++b;
			}
			++u;
			while (((*antiDiag3)[u]  < tmpMax1-scoreDropOff) && (u>0)){
				--u;}
			
			//borders for lower triangle of edit matrix
			b = ::std::max<TPosition>(b,k-yLength+1);
			u = ::std::min<TPosition>(u, xLength-1);
			

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
			for (size_t eu = 0; eu < (*antiDiag3).size();++eu)
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
			for (unsigned int eu = 0; eu < (*antiDiag1).size();++eu){
				if ((*antiDiag1)[eu] > tmpMax){
					tmpMax = (*antiDiag1)[eu];
					tmpPos = eu;
					
				}
			}
			--k;
		}

		if(tmpMax != infimum){
			setRightDim0(seed,rightDim0(seed)+tmpPos);
			setRightDim1(seed,rightDim1(seed)+k-tmpPos);
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





} //end of Seqan namespace

#endif //#ifndef SEQAN_HEADER_
