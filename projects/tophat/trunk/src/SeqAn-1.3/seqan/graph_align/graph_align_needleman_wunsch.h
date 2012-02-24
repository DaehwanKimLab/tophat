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

#ifndef SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H
#define SEQAN_HEADER_GRAPH_NEEDLEMAN_WUNSCH_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
// Gap extension score is taken as the constant gap score!!!
//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair>
inline void
_alignNeedlemanWunschMatrix(TAlign&,
							  TStringSet const&,
							  TTrace const&,
							  TIndexPair const&)
{
	// empty by default
}

template <typename TAlign, typename TStringSet, typename TTrace, typename TIndexPair>
void
_alignNeedlemanWunschTrace(TAlign& align,
							  TStringSet const& str,
							  TTrace const& trace,
							  TIndexPair const& overallMaxIndex)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Value<TTrace>::Type TTraceValue;
	

	_alignNeedlemanWunschMatrix(align, str, trace, overallMaxIndex);
	
	// TraceBack values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	static const TTraceValue tvmap[] = { Diagonal, Diagonal, Horizontal, Diagonal, Vertical, Diagonal, Horizontal, Diagonal };
	
	// Initialization
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = overallMaxIndex[0];
	TSize len2 = overallMaxIndex[1];
	TSize numCols = length(str[0]);
	TSize numRows = length(str[1]);
	if (len1 < numCols) _alignTracePrint(align, str, id1, len1, id2, len2, numCols - len1, Horizontal);
	else if (len2 < numRows) _alignTracePrint(align, str, id1, len1, id2, len2, numRows - len2, Vertical);
	
	if ((len1 != 0) && (len2 !=0)) {
	
		// Initialize everything
		TTraceValue tv = tvmap[(int)trace[(len1-1)*numRows + (len2-1)]];
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
				tv = tvmap[(int)trace[(len1-1)*numRows + (len2-1)]];
				if (tv == Diagonal) {
					if (tv != tvOld) {
						_alignTracePrint(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len1; --len2;
				} else if (tv == Horizontal) {
					//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
					if (tv != tvOld) {
						_alignTracePrint(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len1;
				} else if (tv == Vertical) {
					//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
					if (tv != tvOld) {
						_alignTracePrint(align, str, id1, len1, id2, len2, segLen, tvOld);
						tvOld = tv; segLen = 1;
					} else ++segLen;
					--len2;
				}
			} while ((len1 != 0) && (len2 !=0));
		}
		// Process left-overs
		_alignTracePrint(align, str, id1, len1, id2, len2, segLen, tvOld);
	}

	// Handle the remaining sequence
	if (len1 != 0) _alignTracePrint(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, Horizontal);
	else if (len2 != 0) _alignTracePrint(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, Vertical);
}

	

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TValPair, typename TIndexPair, typename TAlignConfig>
inline typename Value<TScore>::Type
_alignNeedlemanWunsch(TTrace & trace,
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
	const int Diagonal = 0;
	const int Horizontal = 1;
	const int Vertical = 2;

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
	overallMaxValue[0] = MinValue<TScoreValue>::VALUE;
	overallMaxValue[1] = MinValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = len1;
	overallMaxIndex[1] = len2;

	resize(column, len2 + 1);
//	resize(trace, len1*len2);
	resize(trace, len1*len2, 0);
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
	TScoreValue max_all;
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
			
			// Choose the max values
			max_all = _max(_max(max_verti, max_hori), max_diag);
			if (max_verti == max_all) *it |= 1 << Vertical;
			if (max_diag == max_all)  *it |= 1 << Diagonal;
			if (max_hori == max_all)  *it |= 1 << Horizontal;
			*coit = max_verti = max_all;
/*			if (max_diag >= _max(max_verti, max_hori)) {
				*it = Diagonal;
				max_verti = *coit = max_diag;
			} else if (max_hori >= max_verti) {
				*it = Horizontal;
				max_verti = *coit = max_hori;
			} else {
				*it = Vertical;
				*coit = max_verti;
			}
*/		}
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
	String<char> trace;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	TScoreValue	maxScore = _alignNeedlemanWunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	

	// Follow the trace and create the graph
	_alignNeedlemanWunschTrace(align, str, trace, overallMaxIndex);	

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
	
	String<char> trace;
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[2];
	return _alignNeedlemanWunsch(trace, str, sc, overallMaxValue, overallMaxIndex, TAlignConfig());	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
