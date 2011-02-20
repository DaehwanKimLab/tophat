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
  $Id:
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MOTIF_BASE_H
#define SEQAN_HEADER_FIND_MOTIF_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

///////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Class.MotifFinder:
..summary:Holds the algorithm parameter values and the motif instance(s) found by the appropriate
          motif discovery algorithm.
..cat:Motif Search
..signature:MotifFinder<TValue, TSpec>
..param.TValue:The type of sequences to be analyzed.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:The motif finding algorithm to search with.
...type:Spec.Projection
...type:Spec.EPatternBranching
...type:Spec.PMS1
...type:Spec.PMSP
*/

template <typename TValue, typename TSpec>
class MotifFinder;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.MotifFinder

template<typename TValue, typename TSpec>
struct Value< MotifFinder<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< MotifFinder<TValue, TSpec> const>
{
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.findMotif:
..summary:Represents the main function which is used to start the search for noticeable motif patterns.
..cat:Motif Search
..signature:findMotif(finder,dataset,seq_model)
..param.finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.StringSet
..param.seq_model:The seq_model object.
...type:Tag.OOPS
...type:Tag.OMOPS
...type:Tag.ZOOPS
...type:Tag.TCM
...remarks:The sequence models rely on different assumptions about the distribution of motif occurrences
           across the sample sequences. 
..remarks:The PROJECTION algorithm is able to run in @Tag.OOPS@, @Tag.ZOOPS@ and @Tag.TCM@ mode.
..remarks:The ePatternBranching algorithm is able to run in @Tag.OOPS@ and @Tag.OMOPS@ mode.
..remarks:The PMS1 and PMSP algorithm is able to run in  @Tag.OOPS@,  @Tag.OMOPS@,  @Tag.ZOOPS@ and  
          @Tag.TCM@ mode.
*/

/**
.Function.factorial:
..summary:Calculates the factorial value of any integer number.
..cat:Motif Search
..signature:factorial(value)
..param.value:The value object.
...remarks:$value$ must be a positive integer.
..remarks:The factorial of a non-negative integer $value$ is 
          the product of all positive integers less than or equal to $value$.  
*/

template<typename TType>
TType factorial(TType n)
{
	TType result = 0;

	if(n==0)
	{
		result = 1;
	}
	else
	{
		result = n*factorial(n-1);
	}
   
	return result;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.binomialCoefficient:
..summary:Calculates the binomial coefficient C(n,k).
..cat:Motif Search
..signature:binomialCoefficient(n,k)
..param.n:The first parameter object.
...remarks:$n$ must be a positive integer.
..param.k:The second parameter object.
...remarks:$k$ must be a positive integer.
..remarks:The binomial coefficient of $n$ and $k$ is equal to zero 
          if $k$ is greater than $n$.   
*/

template<typename TType>
TType binomialCoefficient(TType n, TType k)
{
	//SEQAN_ASSERT(!(n<0) & !(k<0));
	TType result = 1;
	for(TType i=(n-k+1); i<=n; ++i)
	{
		result*=i;
	}
	result = result/factorial(k);
	
	return result;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hammingDistance:
..summary:Determines the Hamming distance between two sequences.
..cat:Motif Search
..signature:hammingDistance<TType>(begin1,end1,begin2)
..param.TType:Distance type.
..param.begin1:An iterator pointing to the beginning of the first sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end1:An iterator pointing to the end of the first sequence which is either
            a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.begin2:An iterator pointing to the beginning of the second sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
*/

template<typename TType, typename TStringIterator>
TType hammingDistance(TStringIterator start1, TStringIterator end1, TStringIterator start2)
{
	TType num_of_mismatches = 0;
	while(start1!=end1)
	{
		if(*start1!=*start2)
		{
			++num_of_mismatches;
		}
		++start1;
		++start2;
	}

	return num_of_mismatches;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.inverseHash:
..summary:Determines the corresponding sequence pattern given the hash value.
..cat:Motif Search
..signature:inverseHash<TValue>(hash_value,alphabet_size,seq_size)
..param.hash_value:The hash_value object.
..param.alphabet_size:The alphabet_size object.
...remarks:$alphabet_size$ is four for nucleotide sequences and twenty for amino acid sequences.
..param.seq_size:The seq_size object representing the size of the corresponding sequence.
*/

template<typename TValue, typename TType>
String<TValue>
inverseHash(TType const & hash_value, 
			typename Size<TValue>::Type const & alp_size, 
			typename Size< String<TValue> >::Type const & seq_size)
{
	typedef String<TValue> TString;
	TString seq;
	resize(seq, seq_size);

	TType hash_val = hash_value;
	typedef typename Position<TString>::Type TPos;
	for(TPos i=0; i<seq_size; ++i)
	{
		int letter = hash_val%alp_size;
		seq[i] = (TValue)letter;
		hash_val = (hash_val-letter)/alp_size;
	}

	std::reverse(begin(seq), end(seq));
	return seq;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.displayResult:
..summary:Displays all found motif candidates. In the case of the Projection Motif Finder
          the function displays the consensus pattern of the found motif candidate.
..cat:Motif Search
..signature:displayResult(motif_finder)
..param.motif_finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
*/

template<typename TValue, typename TAlgorithm>
void
displayResult(MotifFinder<TValue, TAlgorithm> & finder)
{
	typedef String<TValue> TString;
	typedef String<TString> TStrings;

	if(length(finder.set_of_motifs)!=0)
	{
		unsigned int counter = 0;
		typename Iterator<TStrings>::Type iter = begin(finder.set_of_motifs);
		for(; !atEnd(iter, finder.set_of_motifs); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "NO MOTIF HAS BEEN FOUND!!!\n";
	}
}

/////////////////////////////////////////////////////////////////////////
/**.Metafunction.Motif:
*/

template <typename T>
struct Motif;

template <typename TValue, typename TSpec>
struct Motif< MotifFinder<TValue,TSpec> >
{
	typedef String<TValue> Type;
};

/////////////////////////////////////////////////////////////////////////
/**.Function.getMotif:
*/

template <typename TValue, typename TSpec, typename TPosition>
inline typename Motif<MotifFinder<TValue, TSpec> >::Type &
getMotif(MotifFinder<TValue, TSpec> & me,
		 TPosition pos)
{
	return me.set_of_motifs[pos];
}
template <typename TValue, typename TSpec>
inline typename Motif<MotifFinder<TValue, TSpec> >::Type &
getMotif(MotifFinder<TValue, TSpec> & me)
{
	return me.set_of_motifs[0];
}

/////////////////////////////////////////////////////////////////////////
/**.Function.motifCount:
*/

template <typename TValue, typename TSpec>
inline size_t
motifCount(MotifFinder<TValue, TSpec> const & me)
{
	return length(me.set_of_motifs);
}


/////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
