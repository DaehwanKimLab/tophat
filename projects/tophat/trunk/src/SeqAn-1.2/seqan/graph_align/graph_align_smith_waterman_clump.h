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
  $Id: graph_align_smith_waterman_clump.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H
#define SEQAN_HEADER_GRAPH_SMITH_WATERMAN_CLUMP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Smith Waterman Alignment with Clumping, affine gap cost
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TForbidden, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
				TStringSet& str,
				TForbidden& forbidden,
				TScore const& sc,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize indexPair[2];
	
	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;

	// Create the trace
	maxScore = _align_smith_waterman(trace, str, sc, initialDir, indexPair, forbidden);	

	//// Debug code
	//for(TSize i= 0; i<length(str[1]);++i) {
	//	for(TSize j= 0; j<length(str[0]);++j) {
	//		std::cout << (TSize) getValue(forbidden, j*length(str[1]) + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	
	// Follow the trace and create the alignment
	_align_smith_waterman_trace(align, str, trace, initialDir, indexPair, forbidden);
	
	return maxScore;
}



//////////////////////////////////////////////////////////////////////////////
template<typename TString, typename TMatches, typename TScores, typename TScore, typename TSize1>
inline void
_localAlignment(StringSet<TString, Dependent<> > const& str,
				TMatches& matches,
				TScores& scores,
				TScore const& sc,
				TSize1 numAlignments,
				SmithWatermanClump)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TMatches>::Type TSize;
  
	// For clumpping remember the used positions
	TSize len0 = length(str[0]);
	TSize len1 = length(str[1]);
	String<bool> forbidden;
	fill(forbidden, len0 * len1, false);

	// Stop looking for local alignments, if there score is too low
	TScoreValue local_score = 0;
	TScoreValue last_score = 0;
	for(TSize count = 0; count < (TSize) numAlignments; ++count) {
		// Create the local alignment
		TSize from = length(matches);
		local_score = _localAlignment(matches, str, forbidden, sc, SmithWatermanClump());
		TSize to = length(matches);
		if (2 * local_score < last_score) {
			resize(matches, from, Generous());
			break;
		} 
		last_score = local_score;
		resize(scores, to);
		for(TSize k = from; k<to; ++k) scores[k] = local_score;
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
