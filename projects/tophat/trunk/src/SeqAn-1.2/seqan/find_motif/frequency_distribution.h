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

#ifndef SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H
#define SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.FrequencyDistribution:
..summary:Holds a collection of objects of a specific type, where each object represents
          the frequency (absolute or relative probability) of a particular residue which is a member
		  of a fixed sequence alphabet.
..cat:Motif Search
..signature:FrequencyDistribution<TValue[, TSpec]>
..param.TValue:The type of sequence which is considered.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:The type of probability distribution. 
...metafunction:Metafunction.Spec
...default: $double$
...remarks: It is preferable to use $double$.
..remarks:The number of objects in @Class.FrequencyDistribution@ equals the size of the sequence alphabet.
*/

template <typename TValue, typename TSpec = double>
class FrequencyDistribution
{
//____________________________________________________________________________________________

	enum { SIZE = ValueSize<TValue>::VALUE }; 

//____________________________________________________________________________________________

public:
	String<TSpec> frequency_list;

	// constructor & destructor
	FrequencyDistribution()
	{
		resize(frequency_list, (unsigned int) SIZE);
		std::fill(begin(frequency_list), end(frequency_list), (TSpec)0);
	}
	FrequencyDistribution(TValue const & letter_)
	{
		resize(frequency_list, (unsigned int) SIZE);
		convertResidueToFrequencyDist(*this, letter_);
	}
	FrequencyDistribution(FrequencyDistribution const & other_)
	{
		frequency_list = other_.frequency_list; 
	}
	~FrequencyDistribution()
	{
	}

	// overloading operators
	FrequencyDistribution & 
	operator = (FrequencyDistribution const & other_)
	{
		if(this!=&other_)
		{
			clear(frequency_list);
			frequency_list = other_.frequency_list; 
		}
		return *this;
	}

	FrequencyDistribution &
	operator += (FrequencyDistribution const & other_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]+=other_.frequency_list[i];
		}

		return *this;
	}



	FrequencyDistribution &
	operator -= (FrequencyDistribution const & other_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]-=other_.frequency_list[i];
		}

		return *this;
	}

	template<typename TType>
	FrequencyDistribution &
	operator *= (TType value_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]*= (TSpec)value_;
		}

		return *this;
	}

	template<typename TPos>
	inline TSpec &
	operator [] (TPos index_)
	{
		return frequency_list[index_];
	}

	template<typename TPosition>
	inline TSpec const & 
	operator [] (TPosition const index_) const
	{
		return frequency_list[index_];
	}


//____________________________________________________________________________________________

};

template <typename TValue, typename TSpec>
FrequencyDistribution<TValue, TSpec> 
operator + (FrequencyDistribution<TValue, TSpec> const & lhs_, FrequencyDistribution<TValue, TSpec> const & rhs_)
{
	FrequencyDistribution<TValue, TSpec> ret(lhs_);
	ret+=rhs_;

	return ret;
}


template <typename TValue, typename TSpec>
FrequencyDistribution<TValue, TSpec> 
operator - (FrequencyDistribution<TValue, TSpec> const & lhs_, FrequencyDistribution<TValue, TSpec> const & rhs_)
{
	FrequencyDistribution<TValue, TSpec> ret(lhs_);
	ret-=rhs_;

	return ret;
}

template <typename TValue, typename TSpec, typename TType>
FrequencyDistribution<TValue, TSpec> 
operator * (FrequencyDistribution<TValue, TSpec> const & fd_, TType value_)
{
	FrequencyDistribution<TValue, TSpec> ret(fd_);
	ret*=value_;

	return ret;
}

template <typename TValue, typename TSpec>
inline std::ostream & 
operator << (std::ostream & ostr, FrequencyDistribution<TValue, TSpec> & fd_) 
{ 
	for(unsigned int i=0; i<FrequencyDistribution<TValue, TSpec>::SIZE; ++i)
	{	
		ostr.width(15);
		ostr << std::left << fd_.frequency_list[i];
	}
	
	return ostr;  
} 


//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.FrequencyDistribution

template <typename TValue, typename TSpec>
struct Spec< FrequencyDistribution<TValue, TSpec> >
{
	typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TSpec Type;
};

//the following is a workaround for an error in GCC 4.1.2
template <typename TValue>
struct Spec< FrequencyDistribution<TValue, double> >
{
	typedef double Type;
};
template <typename TValue>
struct Spec< FrequencyDistribution<TValue, double> const>
{
	typedef double Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.FrequencyDistribution

template <typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator< FrequencyDistribution<TValue, TSpec>, TIteratorSpec >
{
	typedef String<TSpec> TString;
	typedef typename Iterator<TString, TIteratorSpec>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Position.param.T.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> >
{
	typedef String<TSpec> TString;
	typedef typename Position<TString>::Type Type;
};
template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> const>
{
	typedef String<TSpec> TString;
	typedef typename Position<TString const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
struct Size< FrequencyDistribution<TValue, TSpec> >
{
	typedef String<TSpec> TString;
	typedef typename Size<TString>::Type Type;
};
template<typename TValue, typename TSpec>
struct Size<FrequencyDistribution<TValue, TSpec> const>
{
	typedef String<TSpec> TString;
	typedef typename Size<TString const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.FrequencyDistribution
/*
.Metafunction.Value:
..summary:Returns the sequence type of a @Class.FrequencyDistribution@ type
	     (TValue for FrequencyDistribution<TValue, TSpec>).
*/

template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.absFreqOfLettersInSeq:
..summary:Counts the number of times each residue of a fixed sequence alphabet occurs in a given sequence.
..cat:Motif Search
..signature:absFreqOfLettersInSeq(frequencies,begin,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which will hold the calculated frequencies.
...type:Class.FrequencyDistribution
..param.begin:An iterator pointing to the beginning of a given sequence which is either
              a string of @Spec.Dna@ or a string of @Spec.AminoAcid@. 
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end:An iterator pointing to the end of a given sequence which is either
            a string of @Spec.Dna@ or a string of @Spec.AminoAcid@.  
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
*/

template<typename TValue, typename TSpec, typename TSeqIter> 
void 
absFreqOfLettersInSeq(FrequencyDistribution<TValue, TSpec> & fd,
					  TSeqIter seq_start,
					  TSeqIter seq_end) 
{	
	while(seq_start!=seq_end)
	{
		++fd[(int)*seq_start];
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.absFreqOfLettersInSetOfSeqs:
..summary:Counts the number of times each residue of a fixed sequence alphabet occurs in a given set of sequences.
..cat:Motif Search
..signature:absFreqOfLettersInSetOfSeqs(frequencies,begin,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.begin:An iterator pointing to the first sequence of a given set of sequences which is considered. 
...type:Concept.Iterator
..param.end:An iterator pointing to the last sequence of a given set of sequences which is considered. 
...type:Concept.Iterator
..remarks.text:This function is similar to @Function.absFreqOfLettersInSeq@ except that the function is performed
               on a set of sequences.
*/

template<typename TValue, typename TSpec, typename TIter>
void
absFreqOfLettersInSetOfSeqs(FrequencyDistribution<TValue, TSpec> & fd,
							TIter seq_start,
							TIter seq_end)
{
	while(seq_start!=seq_end)
	{
		absFreqOfLettersInSeq(fd, begin(*seq_start), end(*seq_start));
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addValue:
..summary:Adds a value of a specific type to each element of a given @Class.FrequencyDistribution@ object.
..cat:Motif Search
..signature:addValue(frequencies,value)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.value:The value object which is added to each element of a @Class.FrequencyDistribution@ object.
...remarks:The $value$ object should be identical in type to the elements of the @Class.FrequencyDistribution@ object.
*/

template<typename TValue, typename TSpec, typename TType>
void 
addValue(FrequencyDistribution<TValue, TSpec> & fd, TType const & val)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i]+= (TSpec)val;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.backgroundFrequency:
..summary:Determines the background letter frequencies in a given dataset
..cat:Motif Search
..signature:backgroundFrequency(frequencies,begin,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.begin:An iterator pointing to the first sequence of a given dataset (set of sequences) which is considered. 
...type:Concept.Iterator
..param.end:An iterator pointing to the last sequence of a given dataset (set of sequences) which is considered. 
...type:Concept.Iterator
*/

template<typename TValue, typename TSpec,typename TDatasetIter> 
void 
backgroundFrequency(FrequencyDistribution<TValue, TSpec> & fd,
					TDatasetIter dataset_start,
					TDatasetIter dataset_end)
{
	absFreqOfLettersInSetOfSeqs(fd, dataset_start, dataset_end);

	// check for zero entries
	if(std::find(begin(fd), end(fd), (TSpec)0)!= end(fd))
	{
		// add pseudocounts
		double epsilon = 0.1;
		seqan::Pseudocount<TValue, CMode> p(epsilon);
		addValue(fd, p.pseudocount);
	}
	normalize(fd);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.begin.param.object.type:Class.FrequencyDistribution

template <typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
begin(FrequencyDistribution<TValue, TSpec> & me)
{
	return begin(me.frequency_list);
}
template <typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
begin(FrequencyDistribution<TValue, TSpec> const & me)
{
	return begin(me.frequency_list);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
void 
clear(FrequencyDistribution<TValue, TSpec> & fd) 
{
	fd = FrequencyDistribution<TValue, TSpec>();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertResidueToFrequencyDist:
..summary:Coverts a residue to a frequency distribution (profile).
..cat:Motif Search
..signature:convertResidueToFrequencyDist(frequencies,residue)
..param.frequencies:The @Class.FrequencyDistribution@ object representing the profile for a specific residue.
...type:Class.FrequencyDistribution
..param.residue:The residue object which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:This function is used to convert a sequence pattern into a profile.
..see:Function.convertPatternToProfile
*/

template<typename TValue, typename TSpec>
void 
convertResidueToFrequencyDist(FrequencyDistribution<TValue, TSpec> & fd, TValue const & residue)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	TSpec probability = 
		(TSpec)(0.5/(ValueSize<TValue>::VALUE-1));

	for(TPos i=0; i<length(fd); ++i)
	{
		if(i==residue)
		{
			fd[i] = 0.5;
		}
		else
		{
			fd[i] = probability;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.end.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
end(FrequencyDistribution<TValue, TSpec> & me)
{
	return begin(me)+length(me);
}
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> const >::Type
end(FrequencyDistribution<TValue, TSpec> const & me)
{
	return begin(me)+length(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> & me)
{
	return length(me.frequency_list);
}

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> const & me)
{
	return length(me.frequency_list);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.logarithmize:
..summary:Logarithmizes each element of a given @Class.FrequencyDistribution@ object.
..cat:Motif Search
..signature:logarithmize(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
void 
logarithmize(FrequencyDistribution<TValue, TSpec> & fd) 
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = (TSpec)log(fd[i]);
	}
}

//////////////////////////////////////////////////////////////////////////////

/* s. normalize() (profile.h)
.Function.normalize:
..summary:Determines the normalized frequencies.
..cat:Motif Search
..signature:normalize(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
void 
normalize(FrequencyDistribution<TValue, TSpec> & fd)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	TSpec amount = sum(fd);
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = (TSpec)(fd[i]/amount);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.posOfMax:
..summary:Determines the residue position in a given @Class.FrequencyDistribution@ object with the maximum frequency.
..cat:Motif Search
..signature:posOfMax(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
typename Position< FrequencyDistribution<TValue, TSpec> >::Type
posOfMax(FrequencyDistribution<TValue, TSpec> & me)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution; 
	typedef typename Position<TFrequencyDistribution>::Type TPos;

	TPos position = 0;
	TSpec max_value = (TSpec)0;
	for(TPos i=0; i<length(me); ++i)
	{
		if(me[i]>max_value)
		{
			position = i;
			max_value = me[i];
		}
	}

	return position;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.sum:
..summary:Determines the sum of all frequencies in a given @Class.FrequencyDistribution@ object.
..cat:Motif Search
..signature:sum(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
TSpec 
sum(FrequencyDistribution<TValue, TSpec> & me)
{
	TSpec amount = 
		std::accumulate(begin(me), end(me), (TSpec)0);

	return amount;
}
								
//////////////////////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
