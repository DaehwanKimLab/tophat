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
	maxScore = _alignSmithWaterman(trace, str, sc, initialDir, indexPair, forbidden);	

	//// Debug code
	//for(TSize i= 0; i<length(str[1]);++i) {
	//	for(TSize j= 0; j<length(str[0]);++j) {
	//		std::cout << (TSize) getValue(forbidden, j*length(str[1]) + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	
	// Follow the trace and create the alignment
	_alignSmithWatermanTrace(align, str, trace, initialDir, indexPair, forbidden);
	
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
	resize(forbidden, len0 * len1, false);

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
