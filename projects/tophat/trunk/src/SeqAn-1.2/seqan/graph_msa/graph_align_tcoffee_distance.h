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
  $Id: graph_align_tcoffee_distance.h 1764 2008-03-07 10:28:01Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Distance matrix calculation
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Distance Calculation:
..summary:A tag to specify how to calculate distance matrices.
*/

/**
.Tag.Distance Calculation.value.LibraryDistance:
	Using the library itself and heaviest common subsequence to determine a distance matrix
*/
struct LibraryDistance_;
typedef Tag<LibraryDistance_> const LibraryDistance;


/**
.Tag.Distance Calculation.value.KmerDistance:
	Using a simple kmer count to determine a distance matrix
*/
struct KmerDistance_;
typedef Tag<KmerDistance_> const KmerDistance;






//////////////////////////////////////////////////////////////////////////////
// LibraryDistance
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  LibraryDistance)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TMatrix>::Type TValue;

	// Initialization
	clear(distanceMatrix);
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	resize(distanceMatrix, nseq * nseq);

	// All pairwise alignments
	typedef String<String<TVertexDescriptor> > TSegmentString;
	TValue maxScore = 0;
	for(TSize i=0; i<nseq; ++i) {
		TSegmentString seq1;
		TSize len1 = length(str[i]);
		_buildLeafString(g, i, seq1);
		for(TSize j=i+1; j<nseq; ++j) {
			// Align the 2 strings
			TSegmentString seq2;
			TSize len2 = length(str[j]);
			_buildLeafString(g, j, seq2);
			TSegmentString alignSeq;
			TValue score = heaviestCommonSubsequence(g,seq1,seq2,alignSeq);
			
			// Normalize by distance
			if (len1 > len2) score /= len1;
			else score /= len2;
			if (score > maxScore) maxScore = score;
			
			// Remember the value
			distanceMatrix[i*nseq+j] = score;
		}
	}

	// Normalize values
	for(TSize i=0; i<nseq; ++i) 
		for(TSize j=i+1; j<nseq; ++j) 
			distanceMatrix[i*nseq+j] = SEQAN_DISTANCE_UNITY - ((distanceMatrix[i*nseq+j] * SEQAN_DISTANCE_UNITY) / maxScore );
}


//////////////////////////////////////////////////////////////////////////////
// KmerDistance
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize, typename TAlphabet>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  TAlphabet,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, ktup, TAlphabet());
	
	// Similarity to distance conversion
	TMatrixIterator matIt = begin(distanceMatrix, Standard());
	TMatrixIterator endMatIt = end(distanceMatrix, Standard());
	for(;matIt != endMatIt;++matIt) 
		*matIt = SEQAN_DISTANCE_UNITY - (*matIt);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix, typename TSize>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  TSize ktup,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, ktup, typename Value<typename Value<TStringSet>::Type>::Type(), KmerDistance() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix,
				  KmerDistance)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, 3, KmerDistance() );
}


//////////////////////////////////////////////////////////////////////////////


/**
.Function.getDistanceMatrix:
..summary:Computes a pairwise distance matrix from an alignment graph.
..cat:Graph
..signature:
getDistanceMatrix(graph, mat [, tag])
getDistanceMatrix(graph, mat [, ktup] [, alphabet], KmerDistance)
..param.graph:An alignment graph containing the sequences and possible alignment edges.
...type:Spec.Alignment Graph
..param.mat:Out-parameter:Pairwise distance matrix.
...type:Class.String
..param.ktup:Length of k-mers.
...remarks:For KmerDistance the length of the k-mers.
..param.alphabet:Alphabet
...remarks:For KmerDistance the alphabet to use for k-mer counting (e.g., compressed alphabets).
..param.tag:Distance tag
...type:Tag.Distance Calculation
...remarks:Possible values are LibraryDistance or KmerDistance.
...default:KmerDistance
..returns:void
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void 
getDistanceMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				  TMatrix& distanceMatrix)
{
	SEQAN_CHECKPOINT
	getDistanceMatrix(g, distanceMatrix, KmerDistance() );
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
