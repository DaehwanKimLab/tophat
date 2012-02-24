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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment graph generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Segment Match Generation:
..summary:A tag that specifies how to generate segment matches.
..include:seqan/graph_msa.h
*/


/**
.Tag.Segment Match Generation.value.GlobalPairwiseLibrary:
	Segment matches from pairwise global alignments.
..include:seqan/graph_msa.h
*/

struct GlobalPairwiseLibrary_;
typedef Tag<GlobalPairwiseLibrary_> const GlobalPairwiseLibrary;


/**
.Tag.Segment Match Generation.value.LocalPairwiseLibrary:
	Segment matches from pairwise local alignments.
..include:seqan/graph_msa.h
*/

struct LocalPairwiseLibrary_;
typedef Tag<LocalPairwiseLibrary_> const LocalPairwiseLibrary;

/**
.Tag.Segment Match Generation.value.KmerLibrary:
	Segment matches from pairwise kmer alignments.
..include:seqan/graph_msa.h
*/

struct KmerLibrary_;
typedef Tag<KmerLibrary_> const KmerLibrary;


/**
.Tag.Segment Match Generation.value.LcsLibrary:
	Segment matches from pairwise longest common subsequence comparisons.
..include:seqan/graph_msa.h
*/

struct LcsLibrary_;
typedef Tag<LcsLibrary_> const LcsLibrary;





//////////////////////////////////////////////////////////////////////////////
// Pair selection to calculate alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// Dummy function selecting all pairs
template<typename TString, typename TSpec, typename TSize2, typename TSpec2>
inline void 
selectPairs(StringSet<TString, TSpec> const& str,
			String<TSize2, TSpec2>& pList)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, TSpec> TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Iterator<String<TSize2, TSpec2>, Standard>::Type TPairIter;

	TSize nseq = length(str);
	resize(pList, nseq * (nseq - 1));
	TPairIter itPair = begin(pList, Standard());
	for(TSize i=0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			*itPair = i; ++itPair;
			*itPair = j; ++itPair;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Alignment statistics
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TStringSet, typename TPos, typename TSize1>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1> const& matches,
					   TStringSet& str,
					   TPos const from,
					   TPos const to,
					   TSize1& matchLength,	// Number of identical characters
					   TSize1& overlapLength,	// Number of character in overlapping segments (with mismatches and gaps)
					   TSize1& alignLength)	// Length of the alignment
{
	SEQAN_CHECKPOINT
	typedef String<TFragment, TSpec1> TFragmentMatches;
	typedef typename Size<TFragmentMatches>::Type TSize;
	typedef typename Id<TFragmentMatches>::Type TId;
	typedef typename Iterator<TFragmentMatches, Standard>::Type TFragIter;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet; 
	matchLength = 0;
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);

	
	TSize minId1 = len1 + len2;
	TSize minId2 = len1 + len2;
	TSize maxId1 = 0;
	TSize maxId2 = 0;
	TSize matchMismatch_length = 0;

	TFragIter itFrag = begin(matches, Standard());
	TFragIter itFragEnd = itFrag;
	itFrag += from;
	itFragEnd += to;
	TId id1 = sequenceId(*itFrag, 0);
	TId id2 = sequenceId(*itFrag, 1);
	TSize fragLen = 0;
	TSize beginI = 0;
	TSize beginJ = 0;
	for(;itFrag != itFragEnd; ++itFrag) {
		fragLen = fragmentLength(*itFrag, id1);
		beginI = fragmentBegin(*itFrag, id1);
		beginJ = fragmentBegin(*itFrag, id2);
		if (beginI < minId1) minId1 = beginI;
		if (beginJ < minId2) minId2 = beginJ;
		if (beginI + fragLen > maxId1) maxId1 = beginI + fragLen;
		if (beginJ + fragLen > maxId2) maxId2 = beginJ + fragLen;
		typedef typename Infix<TString>::Type TInfix;
		typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
		TInfix inf1 = label(*itFrag, str, id1);
		TInfix inf2 = label(*itFrag, str, id2);
		TInfixIter sIt1 = begin(inf1, Standard());
		TInfixIter sIt2 = begin(inf2, Standard());
		TInfixIter sIt1End = end(inf1, Standard());
		matchMismatch_length += fragLen;
		for(;sIt1 != sIt1End; ++sIt1, ++sIt2) 
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) ++matchLength;
	}
	alignLength = matchMismatch_length + (len1 - matchMismatch_length) + (len2 - matchMismatch_length);
	overlapLength = alignLength -  minId1 - minId2 - (len1 + len2 - maxId1 - maxId2);
}

//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TSegmentMatches, typename TScores>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 String<TSize2, TSpec2> const& pList,
					 TSegmentMatches& matches,
					 TScores& scores,
					 LcsLibrary)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<TSpec> > TStringSet;
	typedef String<TSize2, TSpec2> TPairList;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Value<TSegmentMatches>::Type TFragment;
	typedef typename Value<TScores>::Type TScoreValue;
	typedef typename Iterator<TPairList, Standard>::Type TPairIter;

	// Pairwise longest common subsequence
	TPairIter itPair = begin(pList, Standard());
	TPairIter itPairEnd = end(pList, Standard());
	for(;itPair != itPairEnd; ++itPair) {
		TStringSet pairSet;
		TId id1 = positionToId(str, *itPair); ++itPair;
		TId id2 = positionToId(str, *itPair);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id2);

		// Lcs between first and second string
		TSize from = length(matches);
		globalAlignment(matches, pairSet, Lcs());
		TSize to = length(matches);

		// Record the scores
		resize(scores, to);
		typedef typename Iterator<TSegmentMatches, Standard>::Type TMatchIter;
		typedef typename Iterator<TScores, Standard>::Type TScoreIter;
		TScoreIter itScore = begin(scores, Standard());
		TScoreIter itScoreEnd = end(scores, Standard());
		TMatchIter itMatch = begin(matches, Standard());
		itScore+=from;
		itMatch+=from;
		for(;itScore != itScoreEnd; ++itScore, ++itMatch) *itScore = (*itMatch).len;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores, typename TAlphabet, typename TSize>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 TSegmentMatches& matches,
					 TScores& scores,
					 TSize ktup,
					 TAlphabet,
					 KmerLibrary)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<TSpec> > TStringSet;
	typedef typename Value<TScores>::Type TScoreValue;
	typedef typename Value<TSegmentMatches>::Type TFragment;
	typedef typename Id<TStringSet>::Type TId;
	typedef String<TSize> TTupelString;
	typedef String<TTupelString> TTupelStringSet;
	
	// Initialization
	TSize nseq = length(str);
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, nseq);
	for(TSize k=0;k<nseq;++k) {
		_getTupelString(str[k], tupSet[k], ktup, TAlphabet());
	}

	// Build one q-gram Index for all sequences
	typedef std::pair<TSize, TSize> TPosSeqPair;
	typedef std::set<TPosSeqPair> TQGramOcc;
	String<TQGramOcc> qIndex;
	TSize qIndexSize = 1;
	for(TSize i=0; i<(TSize) ktup;++i) qIndexSize *= alphabet_size;
	resize(qIndex, qIndexSize);
	for(TSize k=0;k<nseq;++k) {
		for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) {
			qIndex[ tupSet[k][i] ].insert(std::make_pair(i, k));
		}
	}
	for(TSize q=0;q< (TSize) qIndexSize;++q) {
		typename TQGramOcc::const_iterator pos = qIndex[q].begin();
		typename TQGramOcc::const_iterator posEnd = qIndex[q].end();
		while (pos != posEnd) {
			typename TQGramOcc::const_iterator pos2 = pos;
			++pos2;
			while (pos2 != posEnd) {
				if (pos->second != pos2->second) {
					appendValue(matches, TFragment(positionToId(str, pos->second), pos->first, positionToId(str, pos2->second), pos2->first, ktup));
					appendValue(scores, (TScoreValue) ktup);
				}
				++pos2;
			}
			++pos;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores, typename TSize>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 TSegmentMatches& matches,
					 TScores& scores,
					 TSize ktup,
					 KmerLibrary)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, matches, scores, ktup,  typename Value<TString>::Type(), KmerLibrary());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 TSegmentMatches& matches,
					 TScores& scores,
					 KmerLibrary)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, matches, scores, 3, KmerLibrary());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScores>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 String<TSize2, TSpec2> const& pList,
					 TScore const& score_type,
					 TSegmentMatches& matches,
					 TScores& scores,
					 LocalPairwiseLibrary)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<TSpec> > TStringSet;
	typedef String<TSize2, TSpec2> TPairList;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Iterator<TPairList, Standard>::Type TPairIter;

	// Pairwise alignments
	TPairIter itPair = begin(pList, Standard());
	TPairIter itPairEnd = end(pList, Standard());
	for(;itPair != itPairEnd; ++itPair) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = positionToId(str, *itPair); ++itPair;
		TId id2 = positionToId(str, *itPair);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id2);

		multiLocalAlignment(pairSet, matches, scores, score_type, 4, SmithWatermanClump() );
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize>
inline void
_resizeWithRespectToDistance(String<TValue, TSpec>& dist, 
							  TSize nseq)
{
	SEQAN_CHECKPOINT;
	resize(dist, nseq * nseq, 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void
_resizeWithRespectToDistance(Graph<Undirected<TCargo, TSpec> >& dist, TSize nseq)
{
	SEQAN_CHECKPOINT
	clear(dist);
	reserve(_getVertexString(dist), nseq);
	for(TSize i=0;i<nseq; ++i) addVertex(dist);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void
_resizeWithRespectToDistance(Nothing&, TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TValue,  typename TSpec, typename TSize>
inline void 
_setDistanceValue(String<TFragment, TSpec1>& matches,
				   StringSet<TString, TSpec2>& pairSet,			
				   String<TValue, TSpec>& dist,
				   TSize i,
				   TSize j,
				   TSize nseq,
				   TSize from)
{
	SEQAN_CHECKPOINT
	typedef typename Position<String<TFragment, TSpec1> >::Type TPos;
	
	// Determine a sequence weight
	TValue matchLen = 0;
	TValue overlapLen = 0;
	TValue alignLen = 0;
	getAlignmentStatistics(matches, pairSet, (TPos) from, (TPos) length(matches),  matchLen, overlapLen, alignLen);
			
	// Calculate sequence similarity
	TValue x = SEQAN_DISTANCE_UNITY - static_cast<TValue>(static_cast<double>(matchLen) / static_cast<double>(alignLen) * static_cast<double>(SEQAN_DISTANCE_UNITY));
	if (i < j)
	  dist[i * nseq + j] = x;
	else
	  dist[j * nseq + i] = x;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TCargo,  typename TSpec, typename TSize>
inline void 
_setDistanceValue(String<TFragment, TSpec1>& matches,
				   StringSet<TString, TSpec2>& pairSet,
				   Graph<Undirected<TCargo, TSpec> >& dist,
				   TSize i,
				   TSize j,
				   TSize,
				   TSize from)
{
	SEQAN_CHECKPOINT
		
	// Determine a sequence weight
	TCargo matchLen = 0;
	TCargo overlapLen = 0;
	TCargo alignLen = 0;
	getAlignmentStatistics(matches, pairSet, (TSize) from, (TSize) length(matches),  matchLen, overlapLen, alignLen);
			
	// Calculate sequence similarity
	TCargo normalizedSimilarity = SEQAN_DISTANCE_UNITY - (TCargo) (((double) matchLen / (double) overlapLen) * ((double) overlapLen / (double) alignLen) * (double) SEQAN_DISTANCE_UNITY);

	addEdge(dist, i, j, normalizedSimilarity);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec, typename TString, typename TSpec2, typename TSize>
inline void 
_setDistanceValue(String<TFragment, TSpec>&,
				   StringSet<TString, TSpec2>&,
				   Nothing&,
				   TSize,
				   TSize,
				   TSize,
				   TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScoreValues, typename TDistance, typename TAlignConfig>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 String<TSize2, TSpec2> const& pList,
					 TScore const& score_type,
					 TSegmentMatches& matches,
					 TScoreValues& scores,
					 TDistance& dist,
					 TAlignConfig const& ac,
					 GlobalPairwiseLibrary)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<TSpec> > TStringSet;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TScoreValues>::Type TScoreValue;
	typedef typename Iterator<String<TSize2, TSpec2>, Standard>::Type TPairIter;

	// Initialization
	TSize nseq = length(str);
	_resizeWithRespectToDistance(dist, nseq);
	
	// Pairwise alignments
	TPairIter itPair = begin(pList, Standard());
	TPairIter itPairEnd = end(pList, Standard());
	for(;itPair != itPairEnd; ++itPair) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = positionToId(str, *itPair); ++itPair;
		TId id2 = positionToId(str, *itPair);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id1);
		assignValueById(pairSet, const_cast<TStringSet&>(str), id2);
				
		// Alignment
		TSize from = length(matches);
		TScoreValue myScore = globalAlignment(matches, pairSet, score_type, ac, Gotoh() );
		TSize to = length(matches);

		// Record the scores
		resize(scores, to);
		typedef typename Iterator<TScoreValues, Standard>::Type TScoreIter;
		TScoreIter itScore = begin(scores, Standard());
		TScoreIter itScoreEnd = end(scores, Standard());
		itScore+=from;
		for(;itScore != itScoreEnd; ++itScore) *itScore = myScore;
			
		// Get the alignment statistics
		_setDistanceValue(matches, pairSet, dist, (TSize) *(itPair-1), (TSize) *itPair, (TSize) nseq, (TSize)from);
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScoreValues, typename TDistance>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 String<TSize2, TSpec2> const& pList,
					 TScore const& score_type,
					 TSegmentMatches& matches,
					 TScoreValues& scores,
					 TDistance& dist,					 
					 GlobalPairwiseLibrary)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, pList, score_type, matches, scores, dist, AlignConfig<>(), GlobalPairwiseLibrary() );
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
