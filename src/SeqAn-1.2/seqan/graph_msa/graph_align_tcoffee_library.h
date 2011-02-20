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
  $Id: graph_align_tcoffee_library.h 1901 2008-04-28 13:07:56Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

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
*/


/**
.Tag.Segment Match Generation.value.GlobalPairwise_Library:
	Segment matches from pairwise global alignments.
*/

struct GlobalPairwise_Library_;
typedef Tag<GlobalPairwise_Library_> const GlobalPairwise_Library;


/**
.Tag.Segment Match Generation.value.LocalPairwise_Library:
	Segment matches from pairwise local alignments.
*/

struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;

/**
.Tag.Segment Match Generation.value.Kmer_Library:
	Segment matches from pairwise kmer alignments.
*/

struct Kmer_Library_;
typedef Tag<Kmer_Library_> const Kmer_Library;


/**
.Tag.Segment Match Generation.value.Lcs_Library:
	Segment matches from pairwise longest common subsequence comparisons.
*/

struct Lcs_Library_;
typedef Tag<Lcs_Library_> const Lcs_Library;





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
					 Lcs_Library)
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
					 Kmer_Library)
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
					 Kmer_Library)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, matches, scores, ktup,  typename Value<TString>::Type(), Kmer_Library());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSegmentMatches, typename TScores>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 TSegmentMatches& matches,
					 TScores& scores,
					 Kmer_Library)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, matches, scores, 3, Kmer_Library());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TSize2, typename TSpec2, typename TScore, typename TSegmentMatches, typename TScores>
inline void 
appendSegmentMatches(StringSet<TString, Dependent<TSpec> > const& str,
					 String<TSize2, TSpec2> const& pList,
					 TScore const& score_type,
					 TSegmentMatches& matches,
					 TScores& scores,
					 LocalPairwise_Library)
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
__resizeWithRespectToDistance(String<TValue, TSpec>& dist, 
							  TSize nseq)
{
	SEQAN_CHECKPOINT
	resize(dist, nseq * nseq);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void
__resizeWithRespectToDistance(Graph<Undirected<TCargo, TSpec> >& dist, TSize nseq)
{
	SEQAN_CHECKPOINT
	clear(dist);
	reserve(_getVertexString(dist), nseq);
	for(TSize i=0;i<nseq; ++i) addVertex(dist);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void
__resizeWithRespectToDistance(Nothing&, TSize)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TValue,  typename TSpec, typename TSize>
inline void 
__setDistanceValue(String<TFragment, TSpec1>& matches,
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
	if (i < j) dist[i*nseq+j] = SEQAN_DISTANCE_UNITY - (TValue) (((double) matchLen / (double) overlapLen) * ((double) overlapLen / (double) alignLen) * (double) SEQAN_DISTANCE_UNITY);
	else dist[j*nseq+i] = SEQAN_DISTANCE_UNITY - (TValue) (((double) matchLen / (double) overlapLen) * ((double) overlapLen / (double) alignLen) * (double) SEQAN_DISTANCE_UNITY);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TString, typename TSpec2, typename TCargo,  typename TSpec, typename TSize>
inline void 
__setDistanceValue(String<TFragment, TSpec1>& matches,
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
__setDistanceValue(String<TFragment, TSpec>&,
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
					 GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<TSpec> > TStringSet;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TScoreValues>::Type TScoreValue;
	typedef typename Iterator<String<TSize2, TSpec2>, Standard>::Type TPairIter;

	// Initialization
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);
	
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
		__setDistanceValue(matches, pairSet, dist, (TSize) *(itPair-1), (TSize) *itPair, (TSize) nseq, (TSize)from);
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
					 GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	appendSegmentMatches(str, pList, score_type, matches, scores, dist, AlignConfig<>(), GlobalPairwise_Library() );
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
