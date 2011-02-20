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

#ifndef SEQAN_HEADER_PSEUDOCOUNT_MODE_P_H
#define SEQAN_HEADER_PSEUDOCOUNT_MODE_P_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// PMode
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.PMode:
..summary: Represents the P computation scheme for handling "zero" probabilities.
..general:Class.Pseudocount
..cat:Motif Search
..signature:Pseudocount<TValue, PMode>
..param.TValue:The type of sequence which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The P mode computation scheme distributes the pseudocounts among the various residue
          according to their background probabilities.
*/

///.Class.Pseudocount.param.TSpec.type:Spec.PMode

struct _PMode;
typedef Tag<_PMode> PMode;


template<typename TValue>
class Pseudocount<TValue, PMode>
{
	enum { SIZE = ValueSize<TValue>::VALUE };

//_________________________________________________________________________________________________

public:
	double pseudocounts[SIZE];
	double epsilon;

//_________________________________________________________________________________________________

	template<typename TFrequencyDistribution>
	Pseudocount(double epsilon_, TFrequencyDistribution & background_):
		epsilon(epsilon_)
	{
		_computePseudocount(background_); 
	}
	Pseudocount(Pseudocount const & other_):
		pseudocounts(other_.pseudocounts),
		epsilon(other_.epsilon)
	{
	}
	~Pseudocount()
	{
	}

	Pseudocount const &
	operator = (Pseudocount const & other_)
	{
		this->pseudocount = other_.pseudocounts;
		this->epsilon = other_.epsilon;

		return *this;
	}

//_________________________________________________________________________________________________

private:

	template<typename TFrequencyDistribution>
	void
	_computePseudocount(TFrequencyDistribution & background_frequency) 
	{
		typedef typename Position< TFrequencyDistribution >::Type TPos;
		
		// calculating pseudocounts for each residue i.
		// pseudocount = epsilon*fi
		for(TPos i=0; i<length(background_frequency); ++i)
		{
			pseudocounts[i] = (double)(epsilon*background_frequency[i]);
		}
	}
//_________________________________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

// Function.normalize (s. profile.h)

template<typename TProfile, typename TValue>
void 
normalize(TProfile & profile, Pseudocount<TValue, PMode> & mode)
{
	typedef typename Value<TProfile>::Type TFreqDist;
	typedef typename Spec<TFreqDist>::Type TFrequencyType;

	typename Size<TProfile>::Type profile_size = length(profile);
	for(typename Position<TProfile>::Type i=0; 
		i<profile_size; 
		++i)
	{
		if(std::find(begin(profile[i]), end(profile[i]), 0)!=end(profile[i]))
		{
			// N:=row sum
			TFrequencyType N = sum(profile[i]);

			// add pseudocounts
			for(typename Position<TFreqDist>::Type j=0; 
				j<length(profile[i]); 
				++j)
			{
				profile[i][j] = 
					((TFrequencyType)(profile[i][j]+mode.pseudocounts[j]))/
					((TFrequencyType)(N+mode.epsilon));
			}
		}
		else
		{
			normalize(profile[i]);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
