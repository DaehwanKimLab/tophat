// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_BASIC_TAG_H
#define SEQAN_HEADER_BASIC_TAG_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.DotDrawing
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
..include:seqan/basic.h
*/

struct DotDrawing_;
typedef Tag<DotDrawing_> const DotDrawing;


/**
.Tag.HammingDistance
..summary:Switch to trigger Hamming distance, which is a measure of character substitutions.
..include:seqan/basic.h
*/

/**
.Tag.LevenshteinDistance
..summary:Switch to trigger Levenshtein distance, which is a measure of edit operations (character substitutions, deletions or insertions).
..remarks:$EditDistance$ is a synonym for $LevenshteinDistance$.
..see:Spec.EditDistance
..include:seqan/basic.h
*/

struct HammingDistance_;
struct LevenshteinDistance_;

typedef Tag<HammingDistance_>		HammingDistance;
typedef Tag<LevenshteinDistance_>	LevenshteinDistance;
typedef Tag<LevenshteinDistance_>	EditDistance; 


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Alignment: Tags
//////////////////////////////////////////////////////////////////////////////
//Sollte eigentlich nach align/, aber da jetzt ja so viele
//alignment algorithmen in graph/ gelandet sind...

/**
.Tag.Global Alignment Algorithms:
..summary:Global alignment algorithm used by globalAlignment.
..see:Function.globalAlignment
..see:Tag.Local Alignment Algorithms
..include:seqan/basic.h
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.NeedlemanWunsch:
	Dynamic programming algorithm for alignments by Needleman and Wunsch.
..include:seqan/basic.h
*/

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> NeedlemanWunsch;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.BandedNeedlemanWunsch:
	The Needleman-Wunsch alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct BandedNeedlemanWunsch_;
typedef Tag<BandedNeedlemanWunsch_> BandedNeedlemanWunsch;


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Gotoh:
	Gotoh's affine gap cost alignment algorithm.
..include:seqan/basic.h
*/
struct Gotoh_;
typedef Tag<Gotoh_> Gotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.BandedGotoh:
	Gotoh's affine gap cost alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct BandedGotoh_;
typedef Tag<BandedGotoh_> BandedGotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersBitVector:
	Myers' bit vector alignment algorithm for edit distance.
	Note that this algorithm does not returns the alignment itself, but only computes the score.
..include:seqan/basic.h
*/
struct MyersBitVector_;
typedef Tag<MyersBitVector_> const MyersBitVector;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersHirschberg:
	Myers' bit vector algorithm for edit distance combined with Hirschberg's linear space alignment algorithm.
..include:seqan/basic.h
*/
struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> const MyersHirschberg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Hirschberg:
	Hirschberg's linear space global alignment algorithm.
..include:seqan/basic.h
*/
struct Hirschberg_;
typedef Tag<Hirschberg_> const Hirschberg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Lcs:
	Longest common subsequence algorithm.
..include:seqan/basic.h
*/
struct Lcs_;
typedef Tag<Lcs_> const Lcs;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms:
..summary:Local alignment algorithm used by localAlignment.
..see:Function.localAlignment
..include:seqan/basic.h
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.SmithWaterman:
	Triggers a Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct SmithWaterman_;
typedef Tag<SmithWaterman_> const SmithWaterman;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.BandedSmithWaterman:
	Triggers a banded version of the Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct BandedSmithWaterman_;
typedef Tag<BandedSmithWaterman_> const BandedSmithWaterman;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.WatermanEggert:
	Local alignment algorithm by Waterman and Eggert with "declumping" (i.e. only non-overlapping local alignments are computed).
.Tag.Local Alignment Algorithms.value.SmithWatermanClump:
	Same as $WatermanEggert$.
..include:seqan/basic.h
*/
struct SmithWatermanClump_;
typedef Tag<SmithWatermanClump_> const SmithWatermanClump;
typedef Tag<SmithWatermanClump_> const WatermanEggert;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.BandedWatermanEggert:
	Triggers a banded version of the local alignment algorithm by Waterman and Eggert with "declumping".
.Tag.Local Alignment Algorithms.value.BandedSmithWatermanClump:
	Same as $BandedWatermanEggert$.
..include:seqan/basic.h
*/
struct BandedWatermanEggert_;
typedef Tag<BandedWatermanEggert_> const BandedSmithWatermanClump;
typedef Tag<BandedWatermanEggert_> const BandedWatermanEggert;

//////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Tag.RNA Folding Algorithms.value.Nussinov:
	Nussinov style RNA folding algorithm
..include:seqan/basic.h
*/
struct Nussinov_;
typedef Tag<Nussinov_> const Nussinov;

//////////////////////////////////////////////////////////////////////////////

struct Blat_;
typedef Tag<Blat_> const Blat;


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
