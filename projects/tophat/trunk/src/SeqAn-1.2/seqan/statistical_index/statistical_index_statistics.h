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
  $Id: statistics.h 2009-03-10 utro@math.unipa.it srombo@deis.unical.it $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STATISTICAL_INDEX_INTERFACE_H
#define SEQAN_HEADER_STATISTICAL_INDEX_INTERFACE_H

namespace seqan
{

template <typename TAlgorithm, typename TFloat, typename TAlphabet, typename TStringSet>
void _numOccurrences(TFloat &nW, String<TAlphabet>& haystack, TStringSet& needle);

/*
.Function._zscore:
..summary:Auxiliary function to compute the z-score index for a set of patterns w.r.t. a set of text strings and a MarkovModel
..signature:template <TFloat,TStringSet,TAlphabet,TSpec,tTAlgorithm>_zscore(W,X,M)
..param.TFloat:The type of the exploited arrays.
..param.TStringSet:A set of strings.
..param.TAlphabet:The type of the alphabet.
..param.TAlgorithm:The algorithm to exploit to compute the number of occurrences of patterns in the text strings.
..param.W:The set of patterns.
...type:Class.StringSet
..param.X:The text strings.
...type:Class.StringSet
..param.M:The @MarkovModel@ object.
...type:Class.MarkovModel
..returns:The z-score for W w.r.t. X and M.
..remarks:If the alphabet is Dna, then the suitable correction factors are computed.
*/

template <typename TAlgorithm, typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat _zscore( TStringSet& W,  TStringSet& X, MarkovModel<TAlphabet, TFloat, TSpec> & M)
{


 TFloat z_score=0;
 TFloat nW=0;
 //compute occurrances
 for(unsigned int i=0; i< length(X); i++)
 {
	 String<TAlphabet> temp = getValueById(X, i);
	_numOccurrences<TAlgorithm>(nW, temp, W);
 }

  //compute expectation
  TFloat E = expectation(W, X, M);

  //compute variance
  TFloat V = _computeVariance(W, X, M, E);

  //compute z-score
  z_score=(nW-E)/sqrt(V);

  return z_score;

}

/*
.Function._computevariance:
..summary:Auxiliary function to compute the variance for a set of patterns w.r.t. a set of text strings and a MarkovModel
..signature:template <TFloat,TStringSet,TAlphabet,TSpec>_computevariance(W,X,M)
..param.TFloat:The type of the exploited arrays.
..param.TStringSet:A set of strings.
..param.TAlphabet:The type of the alphabet.
..param.W:The set of patterns.
...type:nolink:TStringSet
..param.X:The text strings.
...type:nolink:TStringSet
..param.M:The @MarkovModel@ object.
...type:Class.MarkovModel
..returns:The variance for W w.r.t. X and M.
..remarks:If the alphabet is Dna, then the suitable correction factors are computed.
*/

template <typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat _computeVariance( TStringSet& W,  TStringSet& X, MarkovModel<TAlphabet, TFloat, TSpec> & M, TFloat &E)
{
	//V=B+2C-E^2
	TFloat V = E;

	//C=D+A

	//compute A and D

	TFloat A = 0;
	TFloat D = 0;
	TFloat tmpA, eQPPPe, eQPPQPPe;
	unsigned int sizeW=length(W);
	unsigned int n;

	String <TFloat> pStar;
	fill(pStar, sizeW, 0);

	Shape<TAlphabet, SimpleShape> orderShape;
	resize(orderShape, M.order);
	unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
	unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) M.order);

	for(unsigned int j=0; j<sizeW; j++){
		String<TAlphabet> string =getValueById(W, j);

		int row = hash(orderShape,begin(string));
		TFloat p = 1;
		for(unsigned int i=1; i<length(string)-M.order+1; i++)
		{
			int column=hash(orderShape,begin(string)+i);
			p*=value(M.transition,_getStringIndexFromMatrix(row,column,column_size));
			row = column;
		}
		value(pStar, j) = p;
	}



	for(unsigned int z=0; z<length(X); z++){

	  for(unsigned int i=0; i<length(X); i++){

	 	n = length(getValueById(X, i));

		 for(unsigned int j=0; j<sizeW; j++){

			String<TAlphabet> Wj =getValueById(W, j);

			TFloat q = (TFloat)(n-(2*length(Wj))+2);

			for(unsigned int k=0; k<sizeW; k++){

				tmpA=value(pStar,j)*value(pStar,k);

				unsigned int jfirst, jlast, kfirst;

				jfirst = hash(orderShape,begin(Wj));

				jlast = hash(orderShape,end(Wj)-M.order);

				kfirst = hash(orderShape,begin(getValueById(W, k)));

				eQPPPe = value(M._qppp, _getStringIndexFromMatrix(jlast,kfirst,column_size));

				eQPPQPPe = value(M._qppqpp, _getStringIndexFromMatrix(jlast,kfirst,column_size));

				tmpA  *= value(M.stationaryDistribution, jfirst) * ((q*(q+1)/2)* value(M.stationaryDistribution, kfirst) - (q-1)*eQPPPe - eQPPQPPe);

				A += tmpA;
			}
		 }
	  }

	  // Compute D
	  D+= _overlapExpectation(W,M,length(getValueById(X, z)));
	}




	//Compute Variance
	V += (2*A) + (2*D) -  std::pow((double) E, (int) 2);

	//return V;
	return V;
}

/*
.Function._overlapExpectation:
..summary:Auxiliary function necessary when correction factors have to be computed
..signature:template <TFloat,TStringSet,TAlphabet,TSpec>_overlapExpectation(W,X,M)
..param.TFloat:The type of the exploited arrays.
..param.TStringSet:A set of strings.
..param.TAlphabet:The type of the alphabet.
..param.W:The set of patterns.
...type:Class.StringSet
..param.X:The text strings.
...type:Class.StringSet
..param.M:The @MarkovModel@ object.
...type:Class.MarkovModel
..returns:A value of overlapping for the expectation.
*/

template <typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat _overlapExpectation(TStringSet& W, MarkovModel<TAlphabet, TFloat, TSpec> &M, unsigned int n)
{
	TFloat E_overlap = 0;
	unsigned int sizeW = length(W);
	for(unsigned int i=0; i<sizeW; i++)
	{
		String<TAlphabet> patt1 = getValueById(W, i);
		unsigned int size1 = length(patt1);
		for(unsigned int j=0; j<sizeW; j++)
		{
			String<TAlphabet> patt2 = getValueById(W, j);
			unsigned int k=1;
			unsigned int size2 = length(patt2);
			if(size1>size2)
			{
				k = size1 - size2 + 1;
			}
			for(; k<size1; k++)
			{
				if(isEqual(infix(patt1,begin(patt1)+k,end(patt1)),infix(patt2,begin(patt2),begin(patt2)+k-1)))
				{
					String<TAlphabet> temp = infix(patt1, begin(patt1),begin(patt1)+k-1);
					append(temp,infix(patt2,begin(patt2),end(patt2)));
					E_overlap += (n - size1 + 1)*M.emittedProbability(temp);
				}
			}
		}
	}
	return E_overlap;
}

/*
.Function._numOccurrences:
..summary:Auxiliary function to compute the number of occurrences of a set of patterns in a set of text strings
..signature:template <tTAlgorithm,TFloat,TAlphabet,TStringSet>_numOccurrences(W,haystack,needle)
..param.TAlgorithm:The algorithm to exploit to compute the number of occurrences of patterns in the text strings.
..param.TFloat:The type of the exploited arrays.
..param.TAlphabet:The type of the alphabet.
..param.TStringSet:A set of strings.
..param.W:The set of patterns.
...type:Class.StringSet
..param.haystack:The text strings.
...type:Metafunction.Haystack
..param.needle:The sequence that is searched in the @Metafunction.Haystack@.
*/

template <typename TAlgorithm, typename TFloat, typename TAlphabet, typename TStringSet>
void _numOccurrences(TFloat &nW, String<TAlphabet>& haystack, TStringSet& needle)
{
	Finder<String<TAlphabet> > finder(haystack);
	Pattern<TStringSet, AhoCorasick> pattern(needle);
	while (find(finder, pattern))
	{
		nW++;
	}
}

/*
.Function._addReveseComplements:
..summary:Computes the reverse complements of a set of strings in input.
..signature:<TStringSet> void _addReveseComplements(needle)
..param.needle:The sequence to be computed the reverse complement.
*/

template <typename TStringSet>
void _addReveseComplements(TStringSet& needle)
{
	unsigned int num= length(needle);

	for(unsigned int i=0; i< num; i++){
  	     DnaStringReverseComplement mycom(getValueById(needle, i));
		 appendValue(needle, mycom);
	}
}

/*
.Function._computeExpectation:
..summary:Auxiliary function to compute the expectation for a set of patterns w.r.t. a text string  and a MarkovModel
..signature:template <TAlphabet,TFloat,TSpec,TStringSet>_computeExpectation(mm,W,n)
..param.TAlphabet:The type of the alphabet.
..param.TFloat:The type of the exploited arrays.
..param.TStringSet:A set of strings.
..param.mm:The @MarkovModel@ object.
...type:Class.MarkovModel
..param.W:The set of patterns.
...type:Class.StringSet
..param.n:The length of the string.
...type:nolink:unsigned int

..returns:The expectation for W w.r.t. a string and M.
*/

template <typename TAlphabet, typename TFloat, typename TSpec, typename TStringSet>
TFloat _computeExpectation(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					 TStringSet & W, unsigned int n)
{
	TFloat E=0;
	for (unsigned int i=0; i<length(W); i++){
		String<TAlphabet> temp = getValueById(W, i);
		E += (n - length(temp) + 1)*mm.emittedProbability(temp);
	}
	return E;
}



///////////////////////////////////////////////////////////////////////
// Extern functions to be provided by SeqAn
///////////////////////////////////////////////////////////////////////

typedef Dna TDnaAlphabet;
typedef String<TDnaAlphabet> TDnaSequence;

/**
.Function.zscore:
..summary:Computes the z-score index for a set of patterns w.r.t. a set of text strings and a MarkovModel
..signature:zscore<TAlgorithm>(W,X,M)
..param.TAlgorithm:The algorithm to exploit to compute the number of occurrences of patterns in the text strings (see @Spec.AhoCorasick@ etc.).
..param.W:The set of patterns.
...type:Class.StringSet
..param.X:The set of text strings.
...type:Class.StringSet
..param.M:The MarkovModel object.
...type:Class.MarkovModel
..returns:The z-score for W w.r.t. X and M.
..remarks:If the alphabet is Dna, then the suitable correction factors are computed.
*/

template <typename TAlgorithm, typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat zscore(TStringSet W,  TStringSet &X, MarkovModel<TAlphabet, TFloat, TSpec> &M)
{
   return _zscore<TAlgorithm, TFloat, TStringSet, TAlphabet, TSpec>(W,X,M);
}

template <typename TAlgorithm, typename TFloat, typename TSpec>
TFloat zscore(StringSet<TDnaSequence> W,  StringSet<TDnaSequence> &X, MarkovModel<Dna, TFloat, TSpec> &M)
{
   //add-reverse complements
   _addReveseComplements(W);

   TFloat z_score=0;
   TFloat nW=0;
   //compute occurrences
   for(unsigned int i=0; i < length(X); i++)
   {
		 String<Dna> temp = getValueById(X, i);
		_numOccurrences<TAlgorithm>(nW, temp, W);
	}

	//compute expectation
	TFloat E = expectation(W, X, M);

	//compute variance
	TFloat V = _computeVariance(W, X, M, E);

	//compute correction factor
	TFloat correction = 0;

	unsigned int n;
	unsigned int sizeW= length(W);

	for(unsigned int j=0; j<length(X); j++){

	 	n = length(getValueById(X, j));

		for(unsigned int i=0; i<sizeW; i++)
		{
			String<Dna> patt = getValueById(W, i);
			DnaStringReverseComplement revpatt(patt);
			String<Dna> revc= revpatt;
			if (isEqual(patt,revc))
			{
				correction += (n-length(patt)+1)*M.emittedProbability(revc);
			}
		}
	}

	V+= correction;

	//compute z-score
	z_score=(nW-E)/sqrt(V);

	return z_score;
}

/**
.Function.variance:
..summary:Computes the variance for a set of patterns w.r.t. a set of text strings and a MarkovModel
..signature:variance(W,X,M)
..param.W:The set of patterns.
...type:Class.StringSet
..param.X:The set of text strings.
...type:Class.StringSet
..param.M:The MarkovModel object.
...type:Class.MarkovModel
..returns:The variance for W w.r.t. X and M.
..remarks:If the alphabet is Dna, then the suitable correction factors are computed.
*/

template <typename TFloat,  typename TStringSet, typename TAlphabet, typename TSpec>
TFloat variance(TStringSet W, TStringSet& X, MarkovModel<TAlphabet, TFloat, TSpec> & M)
{
   TFloat E = expectation(W, X, M);

   return _computeVariance(W,X,M,E);
}

template <typename TFloat, typename TSpec>
TFloat variance(StringSet<TDnaSequence> W, StringSet<TDnaSequence> &X, MarkovModel<Dna, TFloat, TSpec> & M)
{

   //add-reverse complements
	_addReveseComplements(W);

	TFloat E = expectation(W, X, M);

	TFloat var =  _computeVariance(W,X,M,E);

	//compute correction factor
	TFloat correction = 0;

	unsigned int n;
	unsigned int sizeW= length(W);


   for(unsigned int j=0; j<length(X); j++){

	 	n = length(getValueById(X, j));

		for(unsigned int i=0; i<sizeW; i++)
		{
			String<Dna> patt = getValueById(W, i);
			DnaStringReverseComplement revpatt(patt);
			String<Dna> revc= revpatt;
			if (isEqual(patt,revc))
			{
				correction += (n-length(patt)+1)*M.emittedProbability(revc);
			}
		}
	}
	var+=correction;

  return var;
}

/**
.Function.expectation:
..summary:Computes the expectation for a set of patterns w.r.t. a set of text strings and a MarkovModel
..signature:expectation(W,X,M)
..param.W:The set of patterns.
...type:Class.StringSet
..param.X:The set of text strings.
...type:Class.StringSet
..param.M:The MarkovModel object.
...type:Class.MarkovModel
..returns:The expectation for W w.r.t. X and M.
*/

template <typename TAlphabet, typename TFloat, typename TSpec, typename TStringSet>
TFloat expectation(TStringSet & W, TStringSet & X, MarkovModel<TAlphabet, TFloat, TSpec> & M)
{
	unsigned int n;
	TFloat E = 0;

	for(unsigned int i=0; i<length(X); i++){
	 	n = length(getValueById(X, i));
        E += _computeExpectation(M, W, n);
	}

    return E;
}


} //namespace seqan

#endif //#ifndef SEQAN_HEADER_...
