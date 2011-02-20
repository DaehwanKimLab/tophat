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
  $Id: graph_align_tcoffee_kmer.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Simple k-mer counter
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TTupelString, typename TKTup, typename TAlphabet>
inline void
_getTupelString(TString const& str, 
				TTupelString& tupelString,
				TKTup const ktup, 
				TAlphabet) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<typename Value<TTupelString>::Type>::Type TWord;
	
	// Alphabet size
	TWord alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Assign a unique number to each k-tupel
	String<TWord> prod;  // Scaling according to position in k-tupel
	resize(prod,ktup);
	for (TWord i=0; i< (TWord) ktup;++i) {
		prod[ktup-i-1] = 1;
		for(TWord j=0;j<i;++j) prod[ktup-i-1] *= alphabet_size;
	}

	TWord len = length(str);
	clear(tupelString);
	resize(tupelString, len-(ktup - 1)); 
	TWord tupelIndex = 0;
	TWord endTupel = 0;
	tupelString[tupelIndex] = 0;
	for(;endTupel< (TWord) ktup;++endTupel) {
		tupelString[tupelIndex] += (TWord) (ordValue((TAlphabet) str[endTupel])) * prod[endTupel];
	}
	++tupelIndex;
	for(;endTupel<len;++endTupel) {
		tupelString[tupelIndex] = tupelString[tupelIndex - 1];
		tupelString[tupelIndex] -= (TWord) (ordValue((TAlphabet) str[endTupel - ktup])) * prod[0];
		tupelString[tupelIndex] *= alphabet_size;
		tupelString[tupelIndex] += (TWord) (ordValue((TAlphabet) str[endTupel]));
		++tupelIndex;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix, typename TSize, typename TAlphabet>
inline void
getKmerSimilarityMatrix(StringSet<TString, TSpec> const& strSet, 
						THitMatrix& mat, 
						TSize ktup, 
						TAlphabet) 
{
	SEQAN_CHECKPOINT
	typedef TSize TWord;
	typedef String<TWord> TTupelString;
	typedef String<TTupelString> TTupelStringSet;
	typedef typename Value<THitMatrix>::Type TValue;

	// Number of sequences
	TSize nseq = length(strSet);
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;
	TWord qIndexSize = (TWord) std::pow((double)alphabet_size, (double)ktup);

	// Initialization
	// Matrix for common k-tupels between sequence i and j
	resize(mat, nseq*nseq);

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, length(strSet));
	for(TSize k=0;k<(TSize) length(strSet);++k) _getTupelString(strSet[k], tupSet[k], ktup, TAlphabet());

	// Build for each sequence the q-gram Index and count common hits
	String<TWord> qIndex;
	String<TWord> compareIndex;
	for(TSize k=0;k<nseq;++k) {
		clear(qIndex);
		fill(qIndex, qIndexSize, (TWord) 0, Exact());
		for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) ++qIndex[ tupSet[k][i] ];
		TWord value;
	    for (TSize k2=k; k2<nseq; ++k2) {
			clear(compareIndex);
			fill(compareIndex, qIndexSize, (TWord) 0, Exact());
			value = 0;
			for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
				//std::cout << tupSet[k2][i] << "," << compareIndex[ tupSet[k2][i] ] << "," << qIndex[ tupSet[k2][i] ]<< std::endl;
				if (compareIndex[ tupSet[k2][i] ] < qIndex[ tupSet[k2][i] ]) ++value;
				++compareIndex[ tupSet[k2][i] ];
			}
			mat[k*nseq+k2] = value;
		}
	}

	// Scale counts
	for(TWord row = 0; row < (TWord) nseq; ++row) {
		for(TWord col = row+1; col < (TWord) nseq; ++col) {
			// Fractional common kmer count
			// = Number of common q-grams / Number of possible common q-grams
			TValue minVal = (mat[col*nseq+col] < mat[row*nseq+row]) ? mat[col*nseq+col] : mat[row*nseq+row];
			mat[row*nseq+col] = (mat[row*nseq+col] * SEQAN_DISTANCE_UNITY) / minVal;
			//std::cout << mat[row*nseq+col] << ",";
		}
		//std::cout << std::endl;
	}
}




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
