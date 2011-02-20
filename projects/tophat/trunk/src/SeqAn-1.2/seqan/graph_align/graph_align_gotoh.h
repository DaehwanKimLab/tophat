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
  $Id: graph_align_gotoh.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H
#define SEQAN_HEADER_GRAPH_ALIGN_GOTOH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair, typename TVal>
inline void
_align_gotoh_trace(TAlign& align,		 
				   TStringSet const& str,
				   TTrace const& trace,
				   TIndexPair const& overallMaxIndex,
				   TVal const initialDir)
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = overallMaxIndex[0];
	TSize len2 = overallMaxIndex[1];
	TSize numCols = length(str[0]);
	TSize numRows = length(str[1]);
	if (len1 < numCols) _align_trace_print(align, str, id1, len1, id2, len2, numCols - len1,  Horizontal);
	else if (len2 < numRows) _align_trace_print(align, str, id1, len1, id2, len2, numRows - len2,  Vertical);
	numRows = (numRows >> 1) + (numRows & 1);

	if ((len1 != 0) && (len2 !=0)) {

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
			} else if(tv == Horizontal) {
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					if ((nextTraceValue >> 2) & 1) {
						_align_trace_print(align, str, id1, --len1, id2, len2, (TSize) 1,  Horizontal);
						tv =  Diagonal; segLen = 0;
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
						_align_trace_print(align, str, id1, len1, id2, --len2, (TSize) 1,  Vertical);
						tv =  Diagonal; segLen = 0;
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
	}

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1,  Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2,  Vertical);
}



//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_gotoh(TTrace& trace,	     
			 TStringSet const& str,
			 TScore const & sc,
			 TValPair& overallMaxValue,
			 TIndexPair& overallMaxIndex,
			 typename Value<TTrace>::Type& initialDir,
			 TAlignConfig const)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;


	// The DP Matrix for diagonal walks
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	fill(trace, len1 * ((len2 >> 1) + (len2 & 1)), 0);
	TTraceValue tvMat = 0;
	
	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	overallMaxValue[0] = InfimumValue<TScoreValue>::VALUE;
	overallMaxValue[1] = InfimumValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = len1;
	overallMaxIndex[1] = len2;
	
	typedef typename Iterator<TColumn, Standard>::Type TColIter;
	TColIter matIt = begin(mat, Standard());
	*matIt = 0;
	TColIter horiIt = begin(horizontal, Standard());
	
	TScoreValue a = 0;
	TScoreValue b = 0;
	TScoreValue max_val = 0;
	for(TSize row = 1; row <= len2; ++row) {
		_initFirstColumn(TAlignConfig(), *(++matIt), scoreGapOpenVertical(sc, -1, 0, str1, str2) + (row - 1) * scoreGapExtendVertical(sc, -1, row - 1, str1, str2));
		*(++horiIt) = *matIt + scoreGapOpenHorizontal(sc, 0, row-1, str1, str2) - scoreGapExtendHorizontal(sc, 0, row-1, str1, str2);
	}
	_lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, 0);
	if (overallMaxIndex[0] == 0) initialDir = Vertical;
	for(TSize col = 1; col <= len1; ++col) {
		matIt = begin(mat, Standard());
		horiIt = begin(horizontal, Standard());
		TScoreValue diagValMat = *matIt;
		_initFirstRow(TAlignConfig(), *matIt, scoreGapOpenHorizontal(sc, 0, -1, str1, str2) + (col - 1) * scoreGapExtendHorizontal(sc, col-1, -1, str1, str2));
		vert = *matIt + scoreGapOpenVertical(sc, col-1, 0, str1, str2) - scoreGapExtendVertical(sc, col-1, 0, str1, str2);
		TSize row = 1; 
		while(row <= len2) {
			// Get the new maximum for vertical
			a = *matIt + scoreGapOpenVertical(sc, col-1, row-1, str1, str2);
			b = vert + scoreGapExtendVertical(sc, col-1, row-1, str1, str2);
			if (a > b) {vert = a; *it |= 1;}
			else vert = b;

			// Get the new maximum for left
			*it <<= 1;
			a = *(++matIt) + scoreGapOpenHorizontal(sc, col-1, row-1, str1, str2);
			b = *(++horiIt) + scoreGapExtendHorizontal(sc, col-1, row-1, str1, str2);
			if (a > b) {*horiIt = a; *it |= 1;}
			else *horiIt = b;
			
			// Get the new maximum for mat
			*it <<= 2;
			max_val = diagValMat + score(const_cast<TScore&>(sc), col-1, row-1, str1, str2);
			tvMat = Diagonal;
			if (vert > max_val) {
				max_val = vert;
				tvMat = Vertical;
			}
			if (*horiIt > max_val) {
				max_val = *horiIt;
				tvMat = Horizontal;
			}
			*it |= tvMat;

			// Assign the new diagonal values
			diagValMat = *matIt;
			*matIt = max_val;

			if (row & 1) *it <<= 1; else ++it;
			++row;
		}
		if (!(row & 1)) {*it <<= 3; ++it; }
		_lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, col);
		// If we got a new index, store direction
		if (overallMaxIndex[0] == col) initialDir = tvMat;
	}
	_lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, mat);
	
	// If we got a new index, store direction
	if ((overallMaxIndex[1] != len2)  && (overallMaxValue[1] > overallMaxValue[0])) 
		initialDir = (horizontal[overallMaxIndex[1]] == mat[overallMaxIndex[1]]) ? Horizontal : Diagonal;

	// If we end up in the bottom right corner, get direction
	if ((overallMaxIndex[0] == len1) && (overallMaxIndex[1] == len2)) {
		initialDir =  Diagonal;
		if (horizontal[len2] == mat[len2]) initialDir =  Horizontal;
		else if (vert == mat[len2]) initialDir =  Vertical;
	}
	return _maxOfAlignment<TScoreValue>(TAlignConfig(), overallMaxValue, overallMaxIndex, len1, len2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];

	// Create the trace
	TScoreValue maxScore = _align_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, initialDir, TAlignConfig());
	// Follow the trace and create the graph
	_align_gotoh_trace(align, str, trace, overallMaxIndex, initialDir);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	return _align_gotoh(trace, str, sc, overallMaxValue, overallMaxIndex, initialDir, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
