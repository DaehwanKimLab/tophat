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

#ifndef SEQAN_HEADER_PSEUDOCOUNT_MODE_C_H
#define SEQAN_HEADER_PSEUDOCOUNT_MODE_C_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// CMode
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.CMode:
..summary: Represents the C ("constant") computation scheme for handling "zero" probabilities.
..general:Class.Pseudocount
..cat:Motif Search
..signature:Pseudocount<TValue, CMode>
..param.TValue:The type of sequence which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The pseudocount is identical for each residue (pseudocount = epsilon/alphabet_size).
*/

///.Class.Pseudocount.param.TSpec.type:Spec.CMode

struct _CMode;
typedef Tag<_CMode> CMode;


template<typename TValue>
class Pseudocount<TValue, CMode>
{

//_____________________________________________________________________________________________

public:
	double pseudocount;
	double epsilon;

//_____________________________________________________________________________________________

	Pseudocount():
		pseudocount(0),
		epsilon(0)
	{
	}
	Pseudocount(double epsilon_):
		pseudocount(0),
		epsilon(epsilon_)
	{
		_computePseudocount();
	}
	Pseudocount(Pseudocount const & other_):
		pseudocount(other_.pseudocount),
		epsilon(other_.epsilon)
	{
	}
	~Pseudocount()
	{
	}

	Pseudocount const &
	operator = (Pseudocount const & other_)
	{
		pseudocount = other_.pseudocount;
		epsilon = other_.epsilon;

		return *this;
	}

//_____________________________________________________________________________________________

private:
	inline void
	_computePseudocount() 
	{
		// alphabet_size = ValueSize<TValue>::VALUE
		pseudocount = 
			(double)(epsilon/ValueSize<TValue>::VALUE);
	}

//_____________________________________________________________________________________________
	
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

// Function.normalize (s. profile.h)

template<typename TProfile, typename TValue>
void 
normalize(TProfile & profile, 
		  Pseudocount<TValue, CMode> const & mode)
{
	typedef typename Value<TProfile>::Type TFreqDist;
	typedef typename Spec<TFreqDist>::Type TFrequencyType;

	typename Size<TProfile>::Type profile_size = length(profile);
	for(typename Position<TProfile>::Type i=0; 
		i<profile_size; 
		++i)
	{
		typename Iterator<TFreqDist>::Type fd_begin = begin(profile[i]);
		typename Iterator<TFreqDist>::Type fd_end = end(profile[i]);
		if(std::find(fd_begin, fd_end, 0)!= fd_end)
		{
			// N:=row sum
			TFrequencyType N = sum(profile[i]);

			// add pseudocounts
			for(typename Position<TFreqDist>::Type j=0; 
				j<length(profile[i]); 
				++j)
			{
				profile[i][j] = 
					((TFrequencyType)(profile[i][j]+mode.pseudocount))/
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
