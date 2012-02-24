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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BANDED_SMITH_WATERMAN_CLUMP_H
#define SEQAN_HEADER_GRAPH_ALIGN_BANDED_SMITH_WATERMAN_CLUMP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template<typename TString>
inline void
_initAlign(Graph<Alignment<StringSet<TString, Dependent<> > > >& align,
           StringSet<TString, Dependent<> > const& str) {
SEQAN_CHECKPOINT
	assignStringSet(align, str);
}

template<typename TString, typename TPos>
inline void
_finishAlign(Graph<Alignment<StringSet<TString, Dependent<> > > >&,
             TPos,
             TPos,
             TPos,
             TPos) {
SEQAN_CHECKPOINT
    // Nothing to be done?
    //stringSet(g)[0] = infix(stringSet(g)[0], begin1, end1);
    //stringSet(g)[1] = infix(stringSet(g)[1], begin2, end2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TAlign, typename TTrace, typename TVal, typename TIndexPair, typename TDiagonal, typename TForbidden>
inline void
_alignBandedSmithWatermanTrace(TStringSet const& str,
                       TAlign& align,
                       TTrace const& trace,
                       TVal const initialDir,
                       TIndexPair const& indexPair,
                       TDiagonal diagL,
                       TDiagonal diagU,
                       TForbidden& forbidden) {
SEQAN_CHECKPOINT
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Id<TStringSet>::Type TId;
    typedef typename Size<TTrace>::Type TSize;
    typedef typename Value<TTrace>::Type TTraceValue;

    // Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // Initialization
    _initAlign(align, str);
    TString const& str1 = str[0];
    TString const& str2 = str[1];
    TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
    TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
    TSize len1 = length(str1);
    TSize len2 = length(str2);
    TSize lo_row = (diagU <= 0) ? -diagU : 0;
    TSize diagonalWidth = (TSize) (diagU - diagL + 1);

    // Start the trace from the cell with the max value
    TSize row = indexPair[0];
    TSize col = indexPair[1];
    TSize endRow = indexPair[0] + lo_row;
    TSize endCol = indexPair[1] + diagL + endRow;
	TSize actualRow = row + lo_row;
    TSize actualCol = col + diagL + actualRow;
    if ((actualCol == 0) || (actualRow == 0)) return;
	if (actualCol < len1) _alignTracePrint(align, str, id1, actualCol, id2, actualRow, len1 - actualCol, Horizontal);
	if (actualRow < len2) _alignTracePrint(align, str, id1, actualCol, id2, actualRow, len2 - actualRow, Vertical);
	_setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
	
    TTraceValue traceValue = initialDir;
    TTraceValue nextTraceValue = traceValue;
    if (traceValue == Diagonal) nextTraceValue = trace[(--row) * diagonalWidth + col];
    else if (traceValue == Vertical) nextTraceValue = trace[(--row) * diagonalWidth + (++col)];
    else if (traceValue == Horizontal) nextTraceValue = trace[row * diagonalWidth + (--col)];
    actualRow = row + lo_row;
    actualCol = col + diagL + actualRow;
	if (nextTraceValue == Diagonal) _setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
    TSize segLen = 1;
	
    while (nextTraceValue != Stop) {
        if (traceValue == nextTraceValue) {
            ++segLen;
        } else {
            _alignTracePrint(align, str, id1, actualCol, id2, actualRow, segLen, traceValue);
            segLen = 1;
        }
        traceValue = nextTraceValue;
        if (traceValue == Diagonal) nextTraceValue = trace[(--row) * diagonalWidth + col];
        else if (traceValue == Vertical) nextTraceValue = trace[(--row) * diagonalWidth + (++col)];
        else if (traceValue == Horizontal) nextTraceValue = trace[row * diagonalWidth + (--col)];
        actualRow = row + lo_row;
        actualCol = col + diagL + actualRow;
		if (nextTraceValue == Diagonal) _setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
    }
    if (segLen) _alignTracePrint(align, str, id1, actualCol, id2, actualRow, segLen, traceValue);
    
	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol, Horizontal);
	if (actualRow != 0) _alignTracePrint(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow, Vertical);

    _finishAlign(align, actualCol, endCol, actualRow, endRow);
}

////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TIndexPair, typename TDiagonal, typename TForbidden>
inline typename Value<TScore>::Type
_alignBandedSmithWaterman(TTrace& trace,
                 TStringSet const& str,
                 TScore const& sc,
                 typename Value<TTrace>::Type& initialDir,
                 TIndexPair& indexPair,
                 TDiagonal diagL,
                 TDiagonal diagU,
                 TForbidden forbidden) {
SEQAN_CHECKPOINT
    typedef typename Value<TTrace>::Type TTraceValue;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Size<TTrace>::Type TSize;

    // TraceBack values for Smith Waterman
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // Initialization
    TString const& str1 = str[0];
    TString const& str2 = str[1];
    TSize len1 = length(str1) + 1;
    TSize len2 = length(str2) + 1;
    TSize diagonalWidth = (TSize) (diagU - diagL + 1);
    TSize hi_diag = diagonalWidth;
    TSize lo_diag = 0;
    if (diagL > 0) lo_diag = 0;
    else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1 - diagL);
    TSize lo_row = (diagU <= 0) ? -diagU : 0;
    TSize hi_row = len2;
    if (len1 - diagL < hi_row) hi_row = len1 - diagL;
    TSize height = hi_row - lo_row;

    typedef String<TScoreValue> TColumn;
    TColumn mat;
    resize(mat, diagonalWidth, 0);
    resize(trace, height * diagonalWidth, Stop);
    
    // Record the max score
    TScoreValue score_max = 0;
    indexPair[0] = 0; indexPair[1] = 0;
    initialDir = Stop;

    // Classical DP
    //TScoreValue max_val = 0;
    TSize actualCol, actualRow;
    TScoreValue verti_val, hori_val;

    // Initialize first row
    typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
    typedef typename Iterator<TColumn, Standard>::Type TMatIter;
    TTraceIter traceIt;
    TMatIter matIt;
    
    if (lo_diag > 0) --lo_diag;
    if (lo_row >= len1 - diagU) --hi_diag;

    for (TSize row = 1; row < height; ++row) {
        actualRow = row + lo_row;
        if (lo_diag > 0) --lo_diag;
        if (actualRow >= len1 - diagU) --hi_diag;
        traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
        matIt = begin(mat, Standard()) + lo_diag;
        hori_val = 0;
        
        for (TSize col = lo_diag; col < hi_diag; ++col, ++matIt, ++traceIt) {
            actualCol = col + diagL + actualRow;
            if (actualCol != 0) {
                if (_isClumping(forbidden, col+1, row+1, diagonalWidth)) {
                    *matIt = 0;
                } else {
                    // Get the new maximum for mat
                    *matIt += score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                    *traceIt = Diagonal;
                }

                // Get the new maximum for vertical
                if (col < diagonalWidth - 1) {
                    verti_val = *(matIt+1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                } else {
                    verti_val = -2;
                }
                if (verti_val > *matIt) {
                    *matIt = verti_val;
                    *traceIt = Vertical;
                }

                // Get the new maximum for horizontal
                if (col > 0) {
                    hori_val = hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                } else {
                    hori_val = -2;
                }
                if (hori_val > *matIt) {
                    *matIt = hori_val;
                    *traceIt = Horizontal;
                }

                // Check if new maximum is greater than 0
                if (0 > *matIt) {
                    *matIt = 0;
                    *traceIt = Stop;
                    
                }

                // Record the new best score
                if (*matIt > score_max) {
                    indexPair[0] = row; indexPair[1] = col;
                    score_max = *matIt;
                    initialDir = *traceIt;
                }
            } else {
                // Init first col
                *matIt = 0;
                *traceIt = Stop;
            }
            hori_val = *matIt;
        }
    }
 //   // Debug code
 //   std::cerr << std::endl;
	//for(TSize i= 0; i<height;++i) {
	//	for(TSize j= 0; j<diagonalWidth;++j) {
	//		std::cerr << (TSize) getValue(trace, i*diagonalWidth + j) << ',';
	//	}
	//	std::cerr << std::endl;
	//}
	//std::cerr << "Max score: " << indexPair[0] << ',' << indexPair[1] << ':' << score_max << " (" << (TSize) initialDir << ")" << std::endl;
    return score_max;
}

////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TForbidden, typename TScore, typename TDiagonal>
inline typename Value<TScore>::Type
_localAlignment(TStringSet& str,
                TAlign& align,	
                TForbidden& forbidden,
                TScore const& sc,
                TDiagonal diag1,
                TDiagonal diag2,
                BandedSmithWatermanClump) {
SEQAN_CHECKPOINT
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Size<TStringSet>::Type TSize;
    typedef unsigned char TDir;

    // Maximum value
    TScoreValue maxScore;
    TSize indexPair[2];

    // Create the trace
    TDir initialDir;
    String<TDir> trace;
    maxScore = _alignBandedSmithWaterman(trace, str, sc, initialDir, indexPair, diag1, diag2, forbidden);

    // Follow the trace and create the alignment
    _alignBandedSmithWatermanTrace(str, align, trace, initialDir, indexPair, diag1, diag2, forbidden);

    return maxScore;
}

template<typename TString, typename TAlignments, typename TScores, typename TScoreValue, typename TSpec2, typename TDiagonal>
inline void
_localAlignment(StringSet<TString, Dependent<> > const& str,
				TAlignments& alignments,
				TScores& scores,
				Score<TScoreValue, TSpec2> const& sc,
				TScoreValue minScore,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedSmithWatermanClump) {
SEQAN_CHECKPOINT
	typedef typename Value<TAlignments>::Type TAlign;
	typedef typename Size<TString>::Type TSize;
  
	// For clumpping remember the used positions
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
    TSize diagonalWidth = (TSize) (diag2 - diag1 + 1);
    
    TSize lo_row = (diag2 < 0) ? -diag2 : 0;
    TSize hi_row = (len1 - diag1 < len2) ? len1 - diag1 : len2;
    TSize height = hi_row - lo_row;

	String<bool> forbidden;
	resize(forbidden, (height+1) * diagonalWidth, false);//(len1+1) * (len2+1), false);

	// Stop looking for local alignments, if their score is too low
	while (true) {
		// Create the local alignment
        TAlign align;
		TScoreValue local_score = _localAlignment(str, align, forbidden, sc, diag1, diag2, BandedSmithWatermanClump());

        //// Debug code
        //for (int i = 0; i < height; ++i) {
        //    for (int j = 0; j < diagonalWidth; ++j) {
        //        std::cerr << forbidden[i * diagonalWidth + j] << ", ";
        //    }
        //    std::cerr << std::endl;
        //}

		if (local_score >= minScore) {
            appendValue(alignments, align);
            appendValue(scores, local_score);
		} else break;
	}
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

