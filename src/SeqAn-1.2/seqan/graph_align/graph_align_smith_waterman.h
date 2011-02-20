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
  $Id: graph_align_smith_waterman.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize>
inline void
_setForbiddenCell(String<bool, TSpec>& forbidden,
				  TSize len1,
				  TSize len2,
				  TSize numRows)
{
	forbidden[(len1 - 1)*numRows + (len2 - 1)] = true;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void
_setForbiddenCell(Nothing&,
				  TSize,
				  TSize,
				  TSize)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize>
inline bool
_isClumping(String<bool, TSpec> const& forbidden,
			TSize row,
			TSize col,
			TSize len2)
{
	return forbidden[(col-1) * len2 + (row-1)];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline bool
_isClumping(Nothing&,
			TSize,
			TSize,
			TSize)
{
	return false;
}

////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TVal, typename TIndexPair, typename TForbidden>
inline void
_align_smith_waterman_trace(TAlign& align,
							TStringSet const& str,
							TTrace const& trace,
							TVal const initialDir,
							TIndexPair const& indexPair,
							TForbidden& forbidden)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// TraceBack values for Gotoh
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = indexPair[1];
	TSize len2 = indexPair[0];
	if ((indexPair[0] == 0) || (indexPair[1] == 0)) return;
	TSize numCols = length(str[0]);
	TSize numRowsOrig = length(str[1]);
	if (len1 < numCols) _align_trace_print(align, str, id1, len1, id2, len2, numCols - len1, Horizontal);
	if (len2 < numRowsOrig) _align_trace_print(align, str, id1, len1, id2, len2, numRowsOrig - len2, Vertical);
	TSize numRows = (numRowsOrig >> 1) + (numRowsOrig & 1);
	
	

	// Initialize everything	
	TTraceValue nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
	TTraceValue tv = Diagonal;
	if (initialDir == Diagonal) tv = (nextTraceValue & 3);
	else if (initialDir == Horizontal) {
		if ((nextTraceValue >> 2) & 1) _align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1, Horizontal);
		else tv = Horizontal;
	} else if (initialDir == Vertical) {
		if ((nextTraceValue >> 3) & 1) _align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1, Vertical);
		else tv = Vertical;
	}
	TSize segLen = 0;
	TTraceValue tvOld = tv;

	// Now follow the trace
	do {
		nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
		if ((nextTraceValue & 3) == Stop) break;
		_setForbiddenCell(forbidden, len1, len2, numRowsOrig);	
		if (tv == Diagonal) tv = (nextTraceValue & 3);
		else if (tv == Horizontal) {
			if ((nextTraceValue >> 2) & 1) tv = Diagonal; 
			else tv =  Horizontal;
		} else if (tv == Vertical) {
			if ((nextTraceValue >> 3) & 1) tv =  Diagonal; 
			else tv =  Vertical;
		}
		if (tv == Diagonal) {
			if (tv != tvOld) {
				if (tvOld == Vertical) --len2;
				else --len1;
				_align_trace_print(align, str, id1, len1, id2, len2, ++segLen, tvOld);
				tvOld = tv; segLen = 0;
			} else {
				++segLen;
				--len1; --len2;
			}
		} else if (tv == Horizontal) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if ((nextTraceValue >> 2) & 1) {
					_align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1, Horizontal);
					tv = Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len1;
				}
			} else {
				++segLen;
				--len1;
			}
		} else if (tv == Vertical) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				if ((nextTraceValue >> 3) & 1) {
					_align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1, Vertical);
					tv = Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len2;
				}
			} else {
				++segLen;
				--len2;
			}
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen) _align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, Horizontal);
	if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, Vertical);
}


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TIndexPair, typename TForbidden>
inline typename Value<TScore>::Type
_align_smith_waterman(TTrace& trace,
					  TStringSet const& str,
					  TScore const & sc,
					  typename Value<TTrace>::Type& initialDir,
					  TIndexPair& indexPair,
					  TForbidden& forbidden)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	
	// TraceBack values for Smith Waterman
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

	// The DP Matrix for diagonal walks
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	typedef typename Iterator<TColumn, Standard>::Type TMatIter;

	// Initialization
	typedef typename Value<TStringSet>::Type TString;
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	fill(trace, len1 * ((len2 >> 1) + (len2 & 1)), 0);
	TTraceValue tvMat= 0;

	// Record the max score
	TScoreValue score_max = 0;
	indexPair[0] = 0; indexPair[1] = 0;
	initialDir = Stop;
	
	// Classical DP
	TScoreValue max_val = 0;
	TScoreValue a = 0;
	TScoreValue b = 0;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	TMatIter matIt = begin(mat, Standard() );
	TMatIter horiIt = begin(horizontal, Standard() );
	*matIt = 0;
	for(TSize row = 1; row <= len2; ++row) {
		*(++matIt) = 0;
		*(++horiIt) = scoreGapOpenHorizontal(sc, 0, row-1, str1, str2) - scoreGapExtendHorizontal(sc, 0, row-1, str1, str2);
	}
	for(TSize col = 1; col <= len1; ++col) {
		matIt = begin(mat, Standard() );
		horiIt = begin(horizontal, Standard() );
		TScoreValue diagValMat = *matIt;
		*matIt = 0;
		vert = scoreGapOpenVertical(sc, col-1, 0, str1, str2) - scoreGapExtendVertical(sc, col-1, 0, str1, str2);
		TSize row = 1;
		while(row <= len2) {
			if (_isClumping(forbidden, row, col, len2)) {
				*it <<= 3;
				*it |= Stop;
				max_val = 0;
				vert = 0;
				*(++horiIt) = 0;
				++matIt;
			} else {
				// Get the new maximum for vertical
				a = *matIt + scoreGapOpenVertical(sc, col-1, row-1, str1, str2);
				b = vert + scoreGapExtendVertical(sc, col-1, row-1, str1, str2);
				if (a > b) { vert = a; *it |= 1;} 
				else vert = b;
	
				// Get the new maximum for horizontal
				*it <<= 1;
				a = *(++matIt) + scoreGapOpenHorizontal(sc, col-1, row-1, str1, str2);
				b = *(++horiIt) + scoreGapExtendHorizontal(sc, col-1, row-1, str1, str2);
				if (a > b) {*horiIt = a; *it |= 1; } 
				else *horiIt =  b;
	
				// Get the new maximum for mat
				*it <<= 2;
				max_val = diagValMat + score(const_cast<TScore&>(sc), col-1, row-1, str1, str2);
				tvMat =  Diagonal;
				if (vert > max_val) {
					max_val = vert;
					tvMat =  Vertical;
				}
				if (*horiIt > max_val) {
					max_val = *horiIt;
					tvMat =  Horizontal;
				}
				if (0 >= max_val) {
					max_val = 0;
					tvMat =  Stop;
				}
				*it |= tvMat;
			}

			// Assign the new diagonal values
			diagValMat = *matIt;
			*matIt = max_val;

			// Record the new best score
			if (max_val > score_max) {
				indexPair[0] = row; indexPair[1] = col;
				score_max = max_val;
				initialDir = tvMat;
			}

			if (row & 1) *it <<= 1; else ++it;
			++row;
		}
		if (!(row & 1)) {*it <<= 3; ++it; }
	}

	//// Debug code
	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "Max score: " << best_row << ',' << best_col << ':' << score_max << " (" << (TSize) initialDir << ")" << std::endl;

	return score_max;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet const& str,
				TScore const& sc,
				SmithWaterman)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize indexPair[2];
	Nothing forbidden;

	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;

	// Create the trace
	maxScore = _align_smith_waterman(trace, str, sc, initialDir, indexPair, forbidden);	
	// Follow the trace and create the graph
	_align_smith_waterman_trace(align, str, trace, initialDir, indexPair, forbidden);
	
	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TStringSet const& str,
				TScore const& sc,
				SmithWaterman)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	TSize indexPair[2];
	Nothing forbidden;
	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;
	return _align_smith_waterman(trace, str, sc, initialDir, indexPair, forbidden);	
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
