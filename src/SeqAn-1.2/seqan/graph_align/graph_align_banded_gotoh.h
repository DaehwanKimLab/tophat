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
  $Id: graph_align_banded_gotoh.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BANDED_GOTOH_H
#define SEQAN_HEADER_GRAPH_ALIGN_BANDED_GOTOH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Banded Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_align_banded_gotoh_trace(TAlign& align,
						  TStringSet const& str,
						  TTrace const& trace,
						  TValPair const& overallMaxValue,
						  TIndexPair const& overallMaxIndex,
						  TDiagonal const diagL,
						  TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TTrace>::Type TSize;
	typedef unsigned char TTraceValue;

	// Gotoh back-trace values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = 0;
	if (diagU <= 0) lo_row = -1 * diagU;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cout << count << ',';
	//		++count;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// Start the trace from the cell with the max value
	TSize row = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[0] : overallMaxIndex[2];
	TSize col = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[1] : overallMaxIndex[3];

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;
	if (actualCol + 1 < len1) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len1 - (actualCol + 1)),  Horizontal);
	if (actualRow + 1 < len2) _align_trace_print(align, str, id1, actualCol, id2, actualRow, (len2 - (actualRow + 1)),  Vertical);

	// Walk until we hit a border
	TTraceValue tv = (trace[row * diagonalWidth + col] & 3);
	TTraceValue oldTraceValue = tv;
	TSize seqLen = 0;
	bool hitBorder = false;
	do {
		actualRow = row + lo_row;
		actualCol = col + diagL + actualRow;

		// Direction changed, so make aligned segments
		if (oldTraceValue != tv) {
			_align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
			seqLen = 0;
		}

		// Check if we hit a border
		if ((actualRow == 0) || (actualCol == 0)) hitBorder = true;
		else {
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
			
			// Last value was diagonal
			if (tv == Diagonal) {
				oldTraceValue = Diagonal;
				tv = (trace[row * diagonalWidth + col] & 3);
				if (tv == Diagonal) {--row; ++seqLen;}
			} else if (tv == Horizontal) { // Last value was horizontal
				oldTraceValue = Horizontal;
				if ((col > 0) && ((trace[row * diagonalWidth + col] >> 2) & 1)) tv = Diagonal;
				--col; ++seqLen;
			} else { // Vertical
				oldTraceValue = Vertical;
				if ((col < diagonalWidth - 1) && ((trace[row * diagonalWidth + col] >> 3) & 1)) tv = Diagonal;
				--row; ++col; ++seqLen;
			}	
		}
	} while(!hitBorder);
	
	// Align left overs
	if (seqLen) _align_trace_print(align, str, id1, actualCol, id2, actualRow, seqLen, tv);

	// Handle the remaining sequence
	if (actualCol != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}


//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_banded_gotoh(TTrace& trace,
					TStringSet const& str,
					TScore const & sc,
					TValPair& overallMaxValue,
					TIndexPair& overallMaxIndex,
					TDiagonal diagL,
					TDiagonal diagU,
					TAlignConfig const)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TTrace>::Type TSize;

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	
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
	TColumn vertical;
	resize(mat, diagonalWidth);
	resize(vertical, diagonalWidth);
	resize(trace, height * diagonalWidth);
	overallMaxValue[0] = InfimumValue<TScoreValue>::VALUE;
	overallMaxValue[1] = InfimumValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = diagonalWidth;
	overallMaxIndex[1] = height;
	overallMaxIndex[2] = diagonalWidth;
	overallMaxIndex[3] = height;
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cerr << count << ',';
	//		++count;
	//	}
	//	std::cerr << std::endl;
	//}
	//std::cerr << std::endl;

	// Classical DP with affine gap costs
	typedef typename Iterator<TColumn, Standard>::Type TColIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceValue tvMat = 0;
	TSize actualRow = 0;
	TSize actualCol = 0;
	TScoreValue a = 0;
	TScoreValue b = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TColIter vertIt = begin(vertical, Standard()) + lo_diag;
		TColIter matIt = begin(mat, Standard()) + lo_diag;
		hori_val = InfimumValue<TScoreValue>::VALUE;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++vertIt, ++matIt, ++traceIt) {
			actualCol = col + diagL + actualRow;
			//std::cerr << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for vertical
				*traceIt = 0;
				if (col < diagonalWidth - 1) {
					a = *(matIt + 1) + scoreGapOpenVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
					b = (*(vertIt + 1) != InfimumValue<TScoreValue>::VALUE) ? *(vertIt + 1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE;
					if (a > b) {*vertIt = a; *traceIt = 1;}
					else *vertIt = b;
				} else *vertIt = InfimumValue<TScoreValue>::VALUE;

				// Get the new maximum for horizontal
				*traceIt <<= 1;
				if (col > 0) {
					a = *(matIt -1 ) + scoreGapOpenHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
					b = (hori_val != InfimumValue<TScoreValue>::VALUE) ? hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE;
					if (a > b) {hori_val = a; *traceIt |= 1;}
					else hori_val = b;
				} else hori_val = InfimumValue<TScoreValue>::VALUE;

				*traceIt <<= 2;
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), actualCol-1, actualRow-1, str1, str2);
				tvMat = Diagonal;
				if (*vertIt > *matIt) {
					*matIt = *vertIt;
					tvMat = Vertical;
				}
				if (hori_val > *matIt) {
					*matIt = hori_val;
					tvMat = Horizontal;
				}
				*traceIt |= tvMat;
			} else {
				// Usual initialization for first row and column
				if (actualRow == 0) {
					if (actualCol != 0) {
						_initFirstRow(TAlignConfig(), *matIt, scoreGapOpenHorizontal(sc, 0, -1, str1, str2) + (actualCol - 1) * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2));
						*vertIt = *matIt + scoreGapOpenVertical(sc, ((int) actualCol - 1), 0, str1, str2) - scoreGapExtendVertical(sc, ((int) actualCol - 1), 0, str1, str2);
						hori_val = InfimumValue<TScoreValue>::VALUE;
					} else {
						*matIt = 0;
						*vertIt = InfimumValue<TScoreValue>::VALUE;
						hori_val = InfimumValue<TScoreValue>::VALUE;
					}
				} else {
					_initFirstColumn(TAlignConfig(), *matIt, scoreGapOpenVertical(sc, -1, 0, str1, str2) + (actualRow - 1) * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2));
					hori_val = *matIt + scoreGapOpenHorizontal(sc, 0, ((int) actualRow - 1), str1, str2) - scoreGapExtendHorizontal(sc, 0, ((int) actualRow - 1), str1, str2);
					*vertIt = InfimumValue<TScoreValue>::VALUE;
				}
			}

			// Store the maximum
			if (actualCol == len1 - 1) _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			if (actualRow == len2 - 1) _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			//std::cerr << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
			//std::cerr << row << ',' << col << ':' << value(horizontal, row * diagonalWidth + col) << std::endl;
			//std::cerr << row << ',' << col << ':' << value(vertical, row * diagonalWidth + col) << std::endl;
		}
	}
	return (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxValue[0] : overallMaxValue[1];
}

////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 TDiagonal diag1,
				 TDiagonal diag2,
				 BandedGotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];

	// Create the trace
	String<unsigned char> trace;
	TScoreValue maxScore = _align_banded_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
	
	// Follow the trace and create the graph
	_align_banded_gotoh_trace(align, str, trace, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig, typename TDiagonal>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 TDiagonal diag1,
				 TDiagonal diag2,
				 BandedGotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];

	// Calculate the score
	String<unsigned char> trace;
	return _align_banded_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, (int) diag1, (int) diag2, TAlignConfig());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
