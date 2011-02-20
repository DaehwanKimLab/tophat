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
  $Id: graph_align_needleman_wunsch.h 2075 2008-05-21 11:13:28Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H
#define SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
// Gap extension score is taken as the constant gap score!!!
//////////////////////////////////////////////////////////////////////////////


template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair>
void
_align_needleman_wunsch_trace(TAlign& align,
							  TStringSet const& str,
							  TTrace const& trace,
							  TIndexPair const& overallMaxIndex)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Value<TTrace>::Type TTraceValue;


	// TraceBack values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	
	// Initialization
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = overallMaxIndex[0];
	TSize len2 = overallMaxIndex[1];
	TSize numCols = length(str[0]);
	TSize numRows = length(str[1]);
	if (len1 < numCols) _align_trace_print(align, str, id1, len1, id2, len2, numCols - len1, Horizontal);
	else if (len2 < numRows) _align_trace_print(align, str, id1, len1, id2, len2, numRows - len2, Vertical);
	
	if ((len1 != 0) && (len2 !=0)) {
	
		// Initialize everything
		TTraceValue tv = trace[(len1-1)*numRows + (len2-1)];
		TTraceValue tvOld = tv;  // We need to know when the direction of the trace changes

		TSize segLen = 1;
		if (tv == Diagonal) {
			--len1; --len2;
		}
		else if (tv == Horizontal) --len1;
		else if (tv == Vertical) --len2;

		// Now follow the trace
		if ((len1 != 0) && (len2 !=0)) {
			do {
				tv = value(trace, (len1-1)*numRows + (len2-1));
				if (tv == Diagonal) {
					if (tv != tvOld) {
						_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len1; --len2;
				} else if (tv == Horizontal) {
					//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
					if (tv != tvOld) {
						_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len1;
				} else if (tv == Vertical) {
					//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
					if (tv != tvOld) {
						_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len2;
				}
			} while ((len1 != 0) && (len2 !=0));
		}
		// Process left-overs
		_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
	}

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, Vertical);
}

	

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_align_needleman_wunsch(TTrace & trace,
						TStringSet const & str,
						TScore const & _sc,
						TValPair & overallMaxValue,
						TIndexPair & overallMaxIndex,
						TAlignConfig const) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TTrace>::Type TTraceValue;

	// TraceBack values
	TTraceValue Diagonal = 0;
	TTraceValue Horizontal = 1;
	TTraceValue Vertical = 2;

	// One DP Matrix column
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	TColumn column;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	overallMaxValue[0] = InfimumValue<TScoreValue>::VALUE;
	overallMaxValue[1] = InfimumValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = len1;
	overallMaxIndex[1] = len2;

	resize(column, len2 + 1);
	resize(trace, len1*len2);
	typedef typename Iterator<TColumn, Standard>::Type TColIterator;
	TColIterator coit = begin(column, Standard());
	TColIterator col_end = end(column, Standard());
	*coit = 0;
	for(TSize row = 1; row <= len2; ++row) _initFirstColumn(TAlignConfig(), *(++coit), (TScoreValue) (row) * scoreGapExtendVertical(_sc, -1, row - 1, str1, str2));
	_lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, column[len2], 0);
	//for(TSize i = 0; i <= len2; ++i) std::cout << value(column, i) << std::endl;
	
	// Classical DP
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	TScoreValue diagVal = 0;
	TScoreValue max_diag = 0;
	TScoreValue max_verti = 0;
	TScoreValue max_hori = 0;
	TSize col2 = 0;
	for(TSize col = 0; col < len1; ++col) 
	{
		coit = begin(column, Standard());
		diagVal = *coit;
		_initFirstRow(TAlignConfig(), *coit, (TScoreValue) (col+1) * scoreGapExtendHorizontal(_sc, col, -1, str1, str2));
		max_verti = *coit;
		col2 = 0;

		for(;++coit != col_end; ++it)
		{
			// Get max for vertical, horizontal and diagonal
			max_verti += scoreGapExtendVertical(_sc, col, col2, str1, str2);
			max_hori = *coit + scoreGapExtendHorizontal(_sc, col, col2, str1, str2);
			max_diag = diagVal + score(_sc, col, col2++, str1, str2); //compute the maximum in vertiVal 
			
			diagVal = *coit;
			// Choose the max
			if (max_diag >= _max(max_verti, max_hori)) {
				*it = Diagonal;
				max_verti = *coit = max_diag;
			} else if (max_hori >= max_verti) {
				*it = Horizontal;
				max_verti = *coit = max_hori;
			} else {
				*it = Vertical;
				*coit = max_verti;
			}
		}
		_lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, max_verti, col+1);
	}
	_lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, column);

	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return _maxOfAlignment<TScoreValue>(TAlignConfig(), overallMaxValue, overallMaxIndex, len1, len2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 NeedlemanWunsch)
				 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;

	// Create the trace
	String<TraceBack> trace;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	TScoreValue	maxScore = _align_needleman_wunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	

	// Follow the trace and create the graph
	_align_needleman_wunsch_trace(align, str, trace, overallMaxIndex);	

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TAlignConfig>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 TAlignConfig const,
				 NeedlemanWunsch)
{
	SEQAN_CHECKPOINT

	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	
	String<TraceBack> trace;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	return _align_needleman_wunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
