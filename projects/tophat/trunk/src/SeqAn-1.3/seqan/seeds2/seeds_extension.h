// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
//
// Based on the code by Carsten Kemena <carsten.kemena@crg.es>, debugged
// by Birte Kehr <birte.kehr@fu-berlin.de>.
// ==========================================================================
// Seed extension algorithms.
//
// The approach for gapped X-drop extension is based on the algorithm in
// Figure 2 from (Zhang et al., 2000).
//
//  Zhang Z, Schwartz S, Wagner L, Miller W.  A greedy algorithm for aligning
//  DNA sequences.  Journal of computational biologyA a journal of
//  computational molecular cell biology.  2000;7(1-2):203-14.
//  doi:10.1089/10665270050081478
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_EXTENSION_H_
#define SEQAN_SEEDS_SEEDS_EXTENSION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..see:Function.extendSeeds
..see:Function.extendSeedScore
..see:Function.extendSeedsScore
..tag.MatchExtend:Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:Gapped extension of a seed until score drops below a Value. Only @Spec.SimpleSeed@s.
..include:seqan/seeds.h
*/
struct MatchExtend_;
typedef Tag<MatchExtend_> const MatchExtend;

struct UngappedXDrop_;
typedef Tag<UngappedXDrop_> const UnGappedXDrop;

struct GappedXDrop_;
typedef Tag<GappedXDrop_> const GappedXDrop;

enum ExtensionDirection
{
    EXTEND_LEFT,
    EXTEND_RIGHT,
    EXTEND_BOTH,
	EXTEND_NONE
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, MatchExtend)
..signature:extendSeed(seed, query, database, direction, scoreDropOff, scoreMatrix, {UngappedXDrop, GappedXDrop})
..param.seed: The seed to extend.
...type:Class.Seed
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
..include:seqan/seeds.h
*/

// We need one specialization for each combination of the extension
// variants and seeds.  It is not worth to extract the common parts
// for simple and chained seeds.

template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Simple Seeds, we can simply update the
    // begin and end values in each dimension.
	SEQAN_CHECKPOINT;

    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
		}
		setBeginDim0(seed, posDim0);
		setBeginDim1(seed, posDim1);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
		TPosition posDim0 = getEndDim0(seed) ;
		TPosition posDim1 = getEndDim1(seed);
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
		}
		setEndDim0(seed, posDim0);
		setEndDim1(seed, posDim1);
	}
}


template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Chained Seeds, we extend the first and
    // the last Seed Diagonal.
	SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seed), 0u);

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TSeedDiagonal & diag = front(seed);
		TPosition posDim0 = diag.beginDim0;
		TPosition posDim1 = diag.beginDim1;
        TSize diagonalLength = diag.length;
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
            ++diagonalLength;
		}
        diag.beginDim0 = posDim0;
        diag.beginDim1 = posDim1;
        diag.length = diagonalLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posDim0 = diag.beginDim0 + diag.length;
		TPosition posDim1 = diag.beginDim1 + diag.length;
        TSize diagonalLength = diag.length;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
            ++diagonalLength;
		}
        diag.length = diagonalLength;
	}
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<Simple, TConfig> & seed,
           TQuery const & query,
           TDatabase const & database,
           ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
           UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Simple Seeds, we can simply
    // update the begin and end values in each dimension.
	SEQAN_CHECKPOINT;

    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posDim0 >= 1 && posDim1 >= 1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0 - 1] == database[posDim1 - 1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
			--posDim0;
			--posDim1;
		}
		setBeginDim0(seed, posDim0 + mismatchingSuffixLength);
		setBeginDim1(seed, posDim1 + mismatchingSuffixLength);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
		TPosition posDim0 = getEndDim0(seed) ;
		TPosition posDim1 = getEndDim1(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0] == database[posDim1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
            ++posDim0;
            ++posDim1;
		}
		setEndDim0(seed, posDim0 - mismatchingSuffixLength);
		setEndDim1(seed, posDim1 - mismatchingSuffixLength);
    }

    // TODO(holtgrew): Update score?!
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
		   UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Chained Seeds, we extend the
    // first and the last Seed Diagonal.
	SEQAN_CHECKPOINT;

    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
        TSeedDiagonal & diag = front(seed);
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
        TSize diagonalLength = diag.length;
		while (posDim0 >= 1 && posDim1 >= 1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0 - 1] == database[posDim1 - 1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
			--posDim0;
			--posDim1;
            ++diagonalLength;
		}
        diag.beginDim0 = posDim0 + mismatchingSuffixLength;
        diag.beginDim1 = posDim1 + mismatchingSuffixLength;
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posDim0 = diag.beginDim0 + diag.length;
		TPosition posDim1 = diag.beginDim1 + diag.length;
        TSize diagonalLength = diag.length;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0] == database[posDim1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
            ++posDim0;
            ++posDim1;
            ++diagonalLength;
		}
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

    // TODO(holtgrew): Update score?!
}

template<typename TAntiDiag, typename TDropOff, typename TScoreValue>
inline void
_initAntiDiags(TAntiDiag & ,
               TAntiDiag & antiDiag2,
               TAntiDiag & antiDiag3,
               TDropOff dropOff,
               TScoreValue gapCost,
               TScoreValue undefined) {
SEQAN_CHECKPOINT
	// antiDiagonals will be swaped in while loop BEFORE computation of antiDiag3 entries
	//  -> no initialization of antiDiag1 necessary

    resize(antiDiag2, 1);
    antiDiag2[0] = 0;

    resize(antiDiag3, 2);
    if (-gapCost > dropOff) {
        antiDiag3[0] = undefined;
        antiDiag3[1] = undefined;
    } else {
        antiDiag3[0] = gapCost;
        antiDiag3[1] = gapCost;
    }
}

template<typename TAntiDiag>
inline void
_swapAntiDiags(TAntiDiag & antiDiag1,
               TAntiDiag & antiDiag2,
			   TAntiDiag & antiDiag3) {
SEQAN_CHECKPOINT
    TAntiDiag temp;
    move(temp, antiDiag1);
    move(antiDiag1, antiDiag2);
    move(antiDiag2, antiDiag3);
    move(antiDiag3, temp);
}

template<typename TAntiDiag, typename TSize, typename TScoreValue>
inline TSize
_initAntiDiag3(TAntiDiag & antiDiag3,
			   TSize offset,
			   TSize maxCol,
			   TSize antiDiagNo,
			   TScoreValue minScore,
               TScoreValue gapCost,
               TScoreValue undefined) {
SEQAN_CHECKPOINT
	resize(antiDiag3, maxCol + 1 - offset);

    antiDiag3[0] = undefined;
	antiDiag3[maxCol - offset] = undefined;

	if ((int)antiDiagNo * gapCost > minScore) {
		if (offset == 0) {
			// init first column
			antiDiag3[0] = antiDiagNo * gapCost;
		}
		if (antiDiagNo - maxCol == 0) {
			// init first row
			antiDiag3[maxCol - offset] = antiDiagNo * gapCost;
		}
	}
	return offset;
}

template<typename TDiagonal, typename TSize>
inline void
_calcExtendedLowerDiag(TDiagonal & lowerDiag,
					   TSize minCol,
					   TSize antiDiagNo) {
SEQAN_CHECKPOINT
	TSize minRow = antiDiagNo - minCol;
	if ((TDiagonal)minCol - (TDiagonal)minRow < lowerDiag) {
		lowerDiag = (TDiagonal)minCol - (TDiagonal)minRow;
    }
}

template<typename TDiagonal, typename TSize>
inline void
_calcExtendedUpperDiag(TDiagonal & upperDiag,
					   TSize maxCol,
					   TSize antiDiagNo) {
SEQAN_CHECKPOINT
	TSize maxRow = antiDiagNo + 1 - maxCol;
	if ((TDiagonal)maxCol - 1 - (TDiagonal)maxRow > upperDiag) {
		upperDiag = maxCol - 1 - maxRow;
    }
}

template<typename TSeed, typename TSize, typename TDiagonal>
inline void
_updateExtendedSeed(TSeed & seed,
					ExtensionDirection direction,
					TSize cols,
					TSize rows,
					TDiagonal lowerDiag,
					TDiagonal upperDiag) {
SEQAN_CHECKPOINT
	if (direction == EXTEND_LEFT) {
		// Set lower and upper diagonals.
		TDiagonal startDiag = getStartDiagonal(seed);
		if (getLowerDiagonal(seed) > startDiag + lowerDiag)
			setLowerDiagonal(seed, startDiag + lowerDiag);
		if (getUpperDiagonal(seed) < startDiag + upperDiag)
			setUpperDiagonal(seed, startDiag + upperDiag);

		// Set new start position of seed.
		setBeginDim0(seed, getBeginDim0(seed) - cols);
		setBeginDim1(seed, getBeginDim1(seed) - rows);
	} else {  // direction == EXTEND_RIGHT
		// Set new lower and upper diagonals.
		TDiagonal endDiag = getEndDiagonal(seed);
		if (getUpperDiagonal(seed) < endDiag - lowerDiag)
			setUpperDiagonal(seed, endDiag - lowerDiag);
		if (getLowerDiagonal(seed) > endDiag - upperDiag)
			setLowerDiagonal(seed, endDiag - upperDiag);

		// Set new end position of seed.
		setEndDim0(seed, getEndDim0(seed) + cols);
		setEndDim1(seed, getEndDim1(seed) + rows);
	}
    SEQAN_ASSERT_GEQ(getUpperDiagonal(seed), getLowerDiagonal(seed));
    SEQAN_ASSERT_GEQ(getUpperDiagonal(seed), getStartDiagonal(seed));
    SEQAN_ASSERT_GEQ(getUpperDiagonal(seed), getEndDiagonal(seed));
    SEQAN_ASSERT_GEQ(getStartDiagonal(seed), getLowerDiagonal(seed));
    SEQAN_ASSERT_GEQ(getEndDiagonal(seed), getLowerDiagonal(seed));
}

template<typename TConfig, typename TQuerySegment, typename TDatabaseSegment, typename TScoreValue>
TScoreValue
_extendSeedGappedXDropOneDirection(
        Seed<Simple, TConfig> & seed,
        TQuerySegment const & querySeg,
        TDatabaseSegment const & databaseSeg,
        ExtensionDirection direction,
        Score<TScoreValue, Simple> const & scoringScheme,
		TScoreValue scoreDropOff) {
SEQAN_CHECKPOINT
    typedef typename Size<TQuerySegment>::Type TSize;
	typedef typename Seed<Simple,TConfig>::TDiagonal TDiagonal;

    TSize cols = length(querySeg)+1;
    TSize rows = length(databaseSeg)+1;
	if (rows == 1 || cols == 1) return 0;

    TScoreValue gapCost = scoreGap(scoringScheme);
    TScoreValue undefined = minValue<TScoreValue>() - gapCost;

    // DP matrix is calculated by anti-diagonals
    String<TScoreValue> antiDiag1;	//smallest anti-diagonal
	String<TScoreValue> antiDiag2;
	String<TScoreValue> antiDiag3;	//current anti-diagonal

	// Indices on anti-diagonals include gap column/gap row:
	//   - decrease indices by 1 for position in query/database segment
	//   - first calculated entry is on anti-diagonal n� 2

    TSize minCol = 1;
    TSize maxCol = 2;

	TSize offset1 = 0; // number of leading columns that need not be calculated in antiDiag1
	TSize offset2 = 0; //                                                       in antiDiag2
	TSize offset3 = 0; //                                                       in antiDiag3

    _initAntiDiags(antiDiag1, antiDiag2, antiDiag3, scoreDropOff, gapCost, undefined);
    TSize antiDiagNo = 1; // the currently calculated anti-diagonal

	TScoreValue best = 0; // maximal score value in the DP matrix (for drop-off calculation)

	TDiagonal lowerDiag = 0;
	TDiagonal upperDiag = 0;

    while (minCol < maxCol) {
        ++antiDiagNo;
		_swapAntiDiags(antiDiag1, antiDiag2, antiDiag3);
		offset1 = offset2;
		offset2 = offset3;
		offset3 = minCol-1;
		_initAntiDiag3(antiDiag3, offset3, maxCol, antiDiagNo, best - scoreDropOff, gapCost, undefined);

		TScoreValue antiDiagBest = antiDiagNo * gapCost;
        for (TSize col = minCol; col < maxCol; ++col) {
			// indices on anti-diagonals
			TSize i3 = col - offset3;
			TSize i2 = col - offset2;
			TSize i1 = col - offset1;

			// indices in query and database segments
			TSize queryPos, dbPos;
			if (direction == EXTEND_RIGHT) {
				queryPos = col - 1;
				dbPos = antiDiagNo - col - 1;
			} else { // direction == EXTEND_LEFT
				queryPos = cols - 1 - col;
				dbPos = rows - 1 + col - antiDiagNo;
			}

			// Calculate matrix entry (-> antiDiag3[col])
			TScoreValue tmp = _max(antiDiag2[i2-1], antiDiag2[i2]) + gapCost;
			tmp = _max(tmp, antiDiag1[i1-1] + score(scoringScheme, queryPos, dbPos, querySeg, databaseSeg));
			if (tmp < best - scoreDropOff) {
				antiDiag3[i3] = undefined;
			} else {
				antiDiag3[i3] = tmp;
				antiDiagBest = _max(antiDiagBest, tmp);
			}
        }
		best = _max(best, antiDiagBest);

		// Calculate new minCol and minCol
		while (minCol - offset3 < length(antiDiag3) && antiDiag3[minCol - offset3] == undefined &&
			   minCol - offset2 - 1 < length(antiDiag2) && antiDiag2[minCol - offset2 - 1] == undefined) {
			++minCol;
		}

		// Calculate new maxCol
		while (maxCol - offset3 > 0 && (antiDiag3[maxCol - offset3 - 1] == undefined) &&
			                           (antiDiag2[maxCol - offset2 - 1] == undefined)) {
			--maxCol;
		}
		++maxCol;

        // Calculate new lowerDiag and upperDiag of extended seed
		_calcExtendedLowerDiag(lowerDiag, minCol, antiDiagNo);
		_calcExtendedUpperDiag(upperDiag, maxCol - 1, antiDiagNo);

		// end of databaseSeg reached?
		minCol = _max((int)minCol, (int)antiDiagNo + 2 - (int)rows);
		// end of querySeg reached?
		maxCol = _min(maxCol, cols);
    }

	// find positions of longest extension

	// reached ends of both segments
	TSize longestExtensionCol = length(antiDiag3) + offset3 - 2;
	TSize longestExtensionRow = antiDiagNo - longestExtensionCol;
	TScoreValue longestExtensionScore = antiDiag3[longestExtensionCol - offset3];

	if (longestExtensionScore == undefined) {
		if (antiDiag2[length(antiDiag2)-2] != undefined) {
			// reached end of query segment
			longestExtensionCol = length(antiDiag2) + offset2 - 2;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
		}
		else if (antiDiag2[length(antiDiag2)-3] != undefined) {
			// reached end of database segment
			longestExtensionCol = length(antiDiag2) + offset2 - 3;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
		}
	}

	if (longestExtensionScore == undefined) {
		// general case
		for (TSize i = 0; i < length(antiDiag1); ++i) {
			if (antiDiag1[i] > longestExtensionScore) {
				longestExtensionScore = antiDiag1[i];
				longestExtensionCol = i + offset1;
				longestExtensionRow = antiDiagNo - 2 - longestExtensionCol;
			}
		}
	}

	// update seed
    if (longestExtensionScore != undefined) {
		_updateExtendedSeed(seed, direction, longestExtensionCol, longestExtensionRow, lowerDiag, upperDiag);
    }
	return longestExtensionScore;
}

template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
           Score<TScoreValue, Simple> const & scoringScheme,
           TScoreValue scoreDropOff,
		   GappedXDrop const &)
{
    // For gapped X-drop extension of Simple Seeds, we can simply
    // update the begin and end values in each dimension as well as the diagonals.
	SEQAN_CHECKPOINT;

    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;

    // The algorithm only works for linear gap scores < 0, mismatch scores < 0
    // and match scores > 0.
    SEQAN_ASSERT_GT(scoreMatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreMismatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapOpen(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapExtend(scoringScheme), 0);
    SEQAN_ASSERT_EQ(scoreGapExtend(scoringScheme), scoreGapOpen(scoringScheme));

	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        // Do not extend to the left if we are already at the beginning of an
        // infix or the sequence itself.

        typedef typename Prefix<TQuery const>::Type TQueryPrefix;
        typedef typename Prefix<TDatabase const>::Type TDatabasePrefix;

        TQueryPrefix queryPrefix = prefix(query, getBeginDim0(seed));
        TDatabasePrefix databasePrefix = prefix(database, getBeginDim1(seed));
        _extendSeedGappedXDropOneDirection(seed, queryPrefix, databasePrefix, EXTEND_LEFT, scoringScheme, scoreDropOff);
    }

	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        // Do not extend to the right if we are already at the beginning of an
        // infix or the sequence itself.

        typedef typename Suffix<TQuery const>::Type TQuerySuffix;
        typedef typename Suffix<TDatabase const>::Type TDatabaseSuffix;
        
        TQuerySuffix querySuffix = suffix(query, getEndDim0(seed));
        TDatabaseSuffix databaseSuffix = suffix(database, getEndDim1(seed));
		// std::cout << "database = " << database << std::endl;
		// std::cout << "database Suffix = " << databaseSuffix << std::endl;
		// std::cout << "query = " << query << std::endl;
		// std::cout << "query Suffix = " << querySuffix << std::endl;
        _extendSeedGappedXDropOneDirection(seed, querySuffix, databaseSuffix, EXTEND_RIGHT, scoringScheme, scoreDropOff);
    }

    // TODO(holtgrew): Update seed's score?!
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & /*seed*/,
		   TQuery const & /*query*/,
		   TDatabase const & /*database*/,
		   ExtensionDirection /*direction*/,
           Score<TScoreValue, Simple> const & /*scoringScheme*/,
           TScoreValue /*scoreDropOff*/,
		   GappedXDrop const &)
{
    // For ungapped X-drop extension of Chained Seeds, we have to append
    // diagonals to the front and end of the list of seed diagonals and modify
    // the first and last one of the current set of seed diagonals.
	SEQAN_CHECKPOINT;

    SEQAN_ASSERT_FAIL("Write me! Look into the function where this assertion fails for instructions on how to do this.");
    // TODO(holtgrew): Implement gapped X-drop extension with Chained seeds. As follows:
    //
    // Create a simple seed, copy over from chained seed.  Then,
    // performed gapped x-drop extension on the simple seed.  Perform
    // banded alignment on the left and right extended parts.  Use the
    // internal functions for this instead of the user-level functions
    // to initialize, fill the matrix, and compute the traceback
    // object.  Construct the correct SeedDiagonal objects from the
    // traceback objects and add them to the list of diagonals for the
    // diagonal seed.
    //
    // An alternative implementation with storing the banded extension
    // matrix would be too much work and it is questionable if this
    // was faster.  The banded seed alignment code from Tobias Rausch
    // is very optimized.

    // TODO(holtgrew): Update seed's score?!
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_EXTENSION_H_
