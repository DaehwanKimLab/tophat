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
  $Id: seedHandlingTags.h 3420 2009-02-12 12:10:09Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/


//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file


#ifndef SEQAN_SEEDHANDLINGTAGS_H
#define SEQAN_SEEDHANDLINGTAGS_H

namespace SEQAN_NAMESPACE_MAIN
{
/**
.Tag.Seed Adding
..summary:The algorithm used to add a seed to a SeedSet.
..see:Function.addSeed
..see:Function.addSeeds
..cat:Seed Handling
..tag.Single:
	Simple adding of a seed to the set. No chaining or merging.
..tag.Blat:
	Chaining of seeds. Gap is filled with smaller matching segments.
..tag.Chaos:
	Chaining of seeds. One gap is introduced at the best position.
..tag.SimpleChain:
	Chaining of seeds.
*/
struct _addSeeding_Single;
typedef Tag<_addSeeding_Single> const Single;

//also defined in blast_base.h (103)
//moved to basic_tag
/*
struct _Chain_Blat;
typedef Tag<_Chain_Blat> const Blat;
*/

struct _Chain_Chaos;
typedef Tag<_Chain_Chaos> const Chaos;

struct _Chain_Simple;
typedef Tag<_Chain_Simple> const SimpleChain;


template <typename T>
struct BLOCK_SIZE
{
};

/**
.Tag.SeedSet
..cat:Seed Handling
..summary: Tags for the behaviour of a SeedSet
..tag.DefaultScore: 
	Enables scoring of seeds in a SeedSet.
..tag.DefaultNoScore:
	Disables scoring of seeds in a SeedSet.
*/



struct _Gap_Cost_Manhatten;
typedef Tag<_Gap_Cost_Manhatten> const Manhattan;

struct _Gap_Cost_QueryDistance;
typedef Tag<_Gap_Cost_QueryDistance> const QueryDistance;

struct _Gap_Cost_DatabaseDistance;
typedef Tag<_Gap_Cost_DatabaseDistance> const DatabaseDistance;

struct _Gap_Cost_NoGapCost;
typedef Tag<_Gap_Cost_NoGapCost> const NoGapCost;


struct _Good_Seed;
typedef Tag<_Good_Seed> const SeedScore;

struct _Good_Seed2;
typedef Tag<_Good_Seed2> const SeedLength;


template<typename T1, typename T2, typename T3>
struct Scoring_Scheme;

typedef Tag<Scoring_Scheme<SeedScore, Manhattan, int> > const DefaultScore;


typedef Tag<Scoring_Scheme<SeedLength, Manhattan, void> > const DefaultNoScore;


template<typename T>
struct GapCosts
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct GapCosts<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef TGapCosts Type;
};

template<typename T>
struct QualityFactor
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct QualityFactor<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef  TQualityFactor Type;
};

template<typename T>
struct ScoreType
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct ScoreType<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef  TScore Type;
};

/* see find_pattern_base.h
template <typename T>
struct ScoringScheme
{
	typedef T Type;
};
*/

} //namespace Seqan

#endif //#ifndef SEQAN_HEADER_
