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

#ifndef SEQAN_HEADER_PROFILE_H
#define SEQAN_HEADER_PROFILE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertPatternToProfile:
..summary:Converts a pattern into a profile which consists of a set of frequency distributions.
..cat:Motif Search
..signature:convertPatternToProfile(profile,begin,end)
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.begin:An iterator pointing to the beginning of a given sequence pattern which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@.
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end:An iterator pointing to the end of a given sequence pattern which is either
            a @Shortcut.DnaString@ or a @Shortcut.Peptide@.
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..remarks:The number of @Class.FrequencyDistribution@ objects which together form the profile
          equals the length of the given sequence.
..remarks:e.g.:$profile[0]$ represents the frequency distribution for the first residue of
          the given sequence.
..see:Function.convertResidueToFrequencyDist
*/

template<typename TProfile, typename TIterStr>
void 
convertPatternToProfile(TProfile & profile,
						TIterStr str_begin,
						TIterStr str_end)
{
	typedef typename Position<TProfile>::Type TPos;
	unsigned int str_size = str_end-str_begin;
	resize(profile, str_size);
	TPos pos = 0;
	while(str_begin!=str_end)
	{
		convertResidueToFrequencyDist(profile[pos], *str_begin);
		++str_begin;
		++pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertSetOfPatternsToProfile:
..summary:Converts a set of sequence patterns into a profile.
..cat:Motif Search
..signature:convertSetOfPatternsToProfile(profile,l_mers,pseudocount_mode)
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.l_mers:The set of sequence patterns.
...type:Class.StringSet
..param.pseudocount_mode:The @Class.Pseudocount@ object for determining the pseudocount method.
...type:Class.Pseudocount
..remarks:This function is used, for example, in the refinement step of the PROJECTION algorithm to convert
          the collection of l-mers inside the corresponding buckets into a profile. 
*/

template<typename TProfile, typename TStrings, typename TPseudocountMode>
void
convertSetOfPatternsToProfile(TProfile & profile,
					 TStrings & l_mers, 
					 TPseudocountMode & pseudocount)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TProfile>::Type TFrequencyDistribution;

	typename Size<TString>::Type l = length(l_mers[0]);
	resize(profile, l);

	typename Iterator<TStrings>::Type l_mers_iter, l_mers_end;
	typename Iterator<TString>::Type l_mer_iter, l_mer_end;
	l_mers_iter = begin(l_mers);
	l_mers_end = end(l_mers);
	while(l_mers_iter!=l_mers_end)
	{
		l_mer_iter = begin(*l_mers_iter);
		l_mer_end = end(*l_mers_iter);
		while(l_mer_iter!=l_mer_end)
		{
			++profile[(int)(l_mer_iter-begin(*l_mers_iter))][(int)*l_mer_iter];
			++l_mer_iter;
		}
		++l_mers_iter;
	}
	normalize(profile, pseudocount);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.normalize:
..summary:Determines the normalized frequencies.
..cat:Motif Search
..signature:normalize(container)
..signature:normalize(profile,pseudocount_mode)
..param.container:The @Class.FrequencyDistribution@ or @Class.StringSet@ (of @Class.FrequencyDistribution|frequency distributions@) object.
...type:Class.FrequencyDistribution
...type:Class.StringSet
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.pseudocount_mode:The @Class.Pseudocount@ object for determining the pseudocount method.
...type:Class.Pseudocount
..remarks:If necessary, pseudocounts are first added to the frequency values before normalizing them 
          when the parameter $container$ is a StringSet of @Class.FrequencyDistribution|frequency distributions@ (profile).
*/

template<typename TProfile>
void 
normalize(TProfile & profile)
{
	typename Iterator<TProfile>::Type iter = begin(profile);
	typename Iterator<TProfile>::Type iter_end = end(profile);
	while(iter!=iter_end)
	{
		normalize(*iter);
		++iter;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.completeProfile:
..summary:Concatenates the background frequency with the profile for the motif component.
..cat:Motif Search
..signature:completeProfile(profile,background_distribution)
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.background_distribution:The @Class.FrequencyDistribution@ object which represents the backround distribution.
...type:Class.FrequencyDistribution
..remarks:The first row of the final profile (probability matrix) represents the @Class.FrequencyDistribution|background distribution@.
*/

template<typename TProfile>
void 
completeProfile(TProfile & profile,
				typename Value<TProfile>::Type & background_distribution)
{
	TProfile copy(profile);
	resize(profile, length(copy)+1);
	profile[0] = background_distribution;

	typename Iterator<TProfile>::Type iter = begin(copy);
	int counter = 1;
	for(; !atEnd(iter, copy); goNext(iter))
	{
		profile[counter] = *iter;
		++counter;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.display:
..summary:Displays a given set of strings.
..cat:Motif Search
..signature:display(strings)
..param.strings:A StringSet.
...type:Class.StringSet
..remarks:This function can also be used to display a profile (probability matrix) 
          which is a set of @Class.FrequencyDistribution|frequency distributions@.
*/


template<typename TStrings>
void 
display(TStrings & strings)
{
	if(length(strings)!=0)
	{
		typename Iterator<TStrings>::Type iter = begin(strings);
		int counter = 0;
		for(; !atEnd(iter, strings); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "EMPTY STRINGS !!!\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
