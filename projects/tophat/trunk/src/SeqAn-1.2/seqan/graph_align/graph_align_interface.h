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
  $Id: graph_align_interface.h 1809 2008-03-31 12:57:59Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_INTERFACE_H
#define SEQAN_HEADER_GRAPH_ALIGN_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Function.globalAlignment:
..summary:Computes the best global alignment of the two sequences.
..cat:Alignments
..signature:globalAlignment(align, score [, align_config], tag)
..signature:globalAlignment(result, strings, score [, align_config], tag)
..signature:globalAlignment(result, strings, score [, align_config], diagLow, diagHigh, BandedGotoh() )
..param.align:An alignment data structure containing two sequences.
...type:Spec.Alignment Graph
...type:Class.Align
..param.result:A data structure that gets the result of the alignment procedure, 
 e.g., a file stream, or std::cout for a textual alignment, or a FragmentString for storing all the matches.
..param.strings:A string set with that contains two sequences.
...type:Class.StringSet
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.align_config:Alignment configuration options. (optional)
...type:Class.AlignConfig
...remarks:The class AlignConfig has four boolean parameters, i.e., TTop, TLeft, TRight, and TBottom.
If TTop is true the first row of the DP Matrix is initialized with 0's. If TLeft is true the first
column is initialized with 0's. If TRight is true, the maximum is search in the last column. If TBottom
is true, the maximum is searched in the last row. All options can be combined in all possible ways.
....text:This feature is not yet supported for all alignment algorithms (e.g. Hirschberg).
..param.tag:A tag indicating the alignment algorithm to use
...type:Tag.Global Alignment Algorithms
..returns:The maximum score of an global alignment between two sequences given in $align$ or $strings$.
...param.align:An optimal global alignment.
....remarks:If there was an alignment stored in $align$ before $globalAlignment$ was called, it will be replaced.
...param.result:An optimal global alignment.
*/
template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return globalAlignment(file,str,sc, AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return globalAlignment(str,sc,AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return globalAlignment(g,stringSet(g),sc, AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc, TAlignConfig(), diag1, diag2, BandedNeedlemanWunsch());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TDiagonal>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	return globalAlignment(g,sc, AlignConfig<>(), diag1, diag2, BandedNeedlemanWunsch());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc, TAlignConfig(), diag1, diag2, BandedNeedlemanWunsch());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TDiagonal>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	return globalAlignment(file,str,sc, AlignConfig<>(), diag1, diag2, BandedNeedlemanWunsch());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedNeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc, TAlignConfig(), diag1, diag2, BandedNeedlemanWunsch());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedGotoh)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc, TAlignConfig(), diag1, diag2, BandedGotoh());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TDiagonal>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedGotoh)
{
	SEQAN_CHECKPOINT
	return globalAlignment(g,sc, AlignConfig<>(), diag1, diag2, BandedGotoh());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedGotoh)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc, TAlignConfig(), diag1, diag2, BandedGotoh());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TDiagonal>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedGotoh)
{
	SEQAN_CHECKPOINT
	return globalAlignment(file,str,sc, AlignConfig<>(), diag1, diag2, BandedGotoh());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TDiagonal>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedGotoh)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc, TAlignConfig(), diag1, diag2, BandedGotoh());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline int
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Lcs)
{
	return globalAlignment(g, stringSet(g), Lcs());
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.localAlignment:
..summary:Computes the best local alignment of two sequences.
..cat:Alignments
..signature:
localAlignment(strSet, score, tag)
localAlignment(graph, score, tag)
localAlignment(file, strSet, score, tag)
..param.strSet:A string set with 2 sequences.
...type:Class.StringSet
...remarks: If an alignment graph is used that graph must contain a string set with two sequences
..param.graph:The alignment graph having 2 sequences.
...type:Spec.Alignment Graph
..param.file:A file stream or std::cout to write a textual alignment
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.tag:A tag indicating the alignment algorithm to use
...remarks:SmithWaterman
..returns:The maximum score of the best local alignment.
*/
template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
localAlignment(TAlign& file,
			   TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	return _localAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
localAlignment(TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	return _localAlignment(str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
inline TScoreValue
localAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
			   Score<TScoreValue, TSpec2> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _localAlignment(g,stringSet(g),sc,TTag());
}


//////////////////////////////////////////////////////////////////////////////


/**
.Function.multiLocalAlignment:
..summary:Computes multiple local alignments of two sequences.
..cat:Alignments
..signature:
multiLocalAlignment(strSet, matches, scores, scType, numAlign, tag)
..param.strSet:A string set of 2 sequences.
...type:Class.StringSet
..param.matches:The set of all local alignment matches.
..param.scores:The respective local alignment score that match was coming from.
..param.scType:A scoring object.
...type:Class.Score
..param.numAlign:The desired number of local alignments to compute.
..param.tag:A tag indicating the alignment algorithm to use
...remarks:SmithWatermanClump
..returns:void
*/

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TMatches, typename TScores, typename TScoreValue, typename TSpec2, typename TSize, typename TTag>
inline void
multiLocalAlignment(StringSet<TString, Dependent<> > const& str,
					TMatches& matches,
					TScores& scores,
					Score<TScoreValue, TSpec2> const& sc,
					TSize numAlignments,
					TTag)
{
	// Make a multiple local alignment and get all matches
	_localAlignment(str,matches,scores,sc,numAlignments,TTag());
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
