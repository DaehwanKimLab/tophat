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

#ifndef SEQAN_HEADER_ALIGN_ALGORITHMS_H
#define SEQAN_HEADER_ALIGN_ALGORITHMS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// _alignTracePrint: this function is called by various alignment algorithm
// to build up the alignment during traceback

template <typename TSize, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(AlignTraceback<TSize> & tb,
				   TStringSet const &,
				   TId const,
				   TPos const,
				   TId const,
				   TPos const,
				   TPos const segLen,
				   TTraceValue const tv)
{
SEQAN_CHECKPOINT
	appendValue(tb.sizes, segLen);
	appendValue(tb.tvs, tv);
}

//////////////////////////////////////////////////////////////////////////////
// _pumpTraceToAlign: build alignment accoring to the traceback stored in trace
// note that the traceback in trace is "reverse" (from back to front)
template <typename TSource, typename TSpec, typename TTrace> 
void
_pumpTraceToAlign(Align<TSource, TSpec> & align_,
					TTrace trace)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;

	//pump trace into align_ (note: this is relatively slow code here. it could be improved if specialized to the Align Specs)
	clearGaps(align_);
	clearClipping(align_);

	TSize i = length(trace.sizes); //scan trace backwards
	TRowIterator it0 = begin(row(align_, 0));
	TRowIterator it1 = begin(row(align_, 1));
	while (i > 0)
	{
		--i;
		TSize siz = trace.sizes[i];
		switch ((int) trace.tvs[i])
		{
		case 1: //horizontal:
			insertGaps(it1, siz);
			break;

		case 2: //vertical:
			insertGaps(it0, siz);
			break;
		}
		goFurther(it0, siz);
		goFurther(it1, siz);
	}
}

/**
.Function.integrateAlign:
..summary:Integrates an alignment into another by copying the gaps.
..cat:Alignments
...type:Class.Align
..signature:integrateAlign(align1, align2[, positions])
..param.align1:Alignment object into which align2 is to be integrated.
...type:Class.Align
..param.align2:Alignment object that is to be integrated into align1.
...type:Class.Align
..param.positions:The integration positions in align1 for all rows (view positions).
...type:Class.String
..remarks:If the integration positions are not specified, the sources of align2 have to be @Metafunction.Infix@es of the sources of align1.
..include:seqan/align.h
 */
template <typename TSource1, typename TSpec1, typename TSource2, typename TSpec2, typename TPos> 
void
integrateAlign(Align<TSource1, TSpec1> & align,
			   Align<TSource2, TSpec2> const & infixAlign,
			   String<TPos> viewPos) {
SEQAN_CHECKPOINT
	typedef Align<TSource1, TSpec1> TAlign;
	typedef Align<TSource2, TSpec2> TInfixAlign;
	typedef typename Size<TAlign>::Type TSize;
	TSize maxLen = 0;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Row<TInfixAlign>::Type TInfixRow;
	TInfixRow infixRow;

	// iterators on align and infixAlign
	typename Iterator<TRow>::Type it;
	typename Iterator<TInfixRow>::Type infixIt, infixEnd;
	
	for (TSize i = 0; i < length(rows(align)); ++i) {
		infixRow = row(infixAlign, i);

		// init iterators
		it = iter(row(align, i), value(viewPos, i));
		infixIt = begin(infixRow);
		infixEnd = end(infixRow);
		
		// insert leading gaps
		if (beginPosition(infixRow) != 0) {
			insertGaps(it, beginPosition(infixRow));
			goFurther(it, beginPosition(infixRow));
		}

		// walk through Gaps containers and copy gaps
		while (infixIt != infixEnd) {
			TSize gapSize = countGaps(infixIt);
			if (gapSize != 0) {
				insertGaps(it, gapSize);
			}
			goFurther(it, gapSize+1);
			goFurther(infixIt, gapSize+1);
		}

		// find longest row
		if (maxLen < endPosition(infixRow)) {
			maxLen = endPosition(infixRow);
		}
	}

	// insert trailing gaps
	for (TSize i = 0; i < length(rows(align)); ++i) {
		infixRow = row(infixAlign, i);
		TSize diffLen = maxLen - endPosition(infixRow);
		if (diffLen > 0) {
			insertGaps(iter(row(align, i), value(viewPos, i) + endPosition(infixRow)), diffLen);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec1, typename TSpec2> 
void
integrateAlign(Align<TSource, TSpec1> & align,
			   Align<typename Infix<TSource>::Type, TSpec2> const & infixAlign) {
SEQAN_CHECKPOINT
	typedef typename Size<TSource>::Type TSize;
	typedef typename Position<typename Row<Align<TSource, TSpec1> >::Type>::Type TPos;

	String<TPos> viewPos;
	TPos pos;
	for (TSize i = 0; i < length(rows(infixAlign)); ++i) {
		pos = beginPosition(source(row(infixAlign, i))) + clippedBeginPosition(row(infixAlign, i));
		appendValue(viewPos, toViewPosition(row(align, i), pos));
	}

	integrateAlign(align, infixAlign, viewPos);
}


//////////////////////////////////////////////////////////////////////////////
// globalAlignment Interface

//____________________________________________________________________________
//notational sugar

template <typename TAlign, typename TScoreValue, typename TScoreSpec>
TScoreValue
globalAlignment(TAlign & align_,
				Score<TScoreValue, TScoreSpec> const & score_)
{
SEQAN_CHECKPOINT
	if (scoreGapOpen(score_)==scoreGapExtend(score_))
	{//linear gap costs
		return globalAlignment(align_, score_, NeedlemanWunsch());
	}
	else
	{//affine gap costs
		return globalAlignment(align_, score_, Gotoh());
	}
}


//____________________________________________________________________________
// for Align

template <typename TSource, typename TSpec, typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TAlgorithm>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, TScoreSpec> const & score_,
				TAlignConfig tag_align_config,
				TAlgorithm tag_algorithm)
{
SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;

	AlignTraceback<TSize> trace;

	TScoreValue ret_score =  _globalAlignment(trace, stringSet(align_), score_, tag_align_config, tag_algorithm);

	_pumpTraceToAlign(align_, trace);

	return ret_score;
}

template <typename TSource, typename TSpec, typename TScoreValue, typename TScoreSpec, typename TAlgorithm>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, TScoreSpec> const & score_,
				TAlgorithm tag_algorithm)
{
SEQAN_CHECKPOINT

	return globalAlignment(align_, score_, AlignConfig<>(), tag_algorithm);
}



//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
