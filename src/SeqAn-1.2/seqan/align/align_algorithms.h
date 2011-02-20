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
  $Id: align_dynprog.h 1432 2007-12-19 15:11:24Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_ALGORITHMS_H
#define SEQAN_HEADER_ALIGN_ALGORITHMS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// helper class to store traceback

template <typename TSize>
struct _Align_Traceback
{
	String<TSize> sizes;
	String<TraceBack> tvs; //trace values: 0 = diagonal, 1 = horizontal, 2 = vertical
};

//////////////////////////////////////////////////////////////////////////////
// _align_trace_print: this function is called by various alignment algorithm
// to build up the alignment during traceback

template <typename TSize, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(_Align_Traceback<TSize> & tb,
				   TStringSet const &,
				   TId const,
				   TPos const,
				   TId const,
				   TPos const,
				   TPos const segLen,
				   TTraceValue const tv)
{
	appendValue(tb.sizes, segLen);
	appendValue(tb.tvs, tv);
}

//////////////////////////////////////////////////////////////////////////////
// _pump_trace_2_Align: build alignment accoring to the traceback stored in trace
// note that the traceback in trace is "reverse" (from back to front)

template <typename TSource, typename TSpec, typename TTrace> 
void
_pump_trace_2_Align(Align<TSource, TSpec> & align_,
					TTrace trace)
{
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;

	//pump trace into align_ (note: this is relatively slow code here. it could be improved if specialized to the Align Specs)
	clearGaps(align_);

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

//////////////////////////////////////////////////////////////////////////////
// globalAlignment Interface

//____________________________________________________________________________
//notational sugar

template <typename TAlign, typename TScoreValue, typename TScoreSpec>
TScoreValue
globalAlignment(TAlign & align_,
				Score<TScoreValue, TScoreSpec> const & score_)
{
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

	_Align_Traceback<TSize> trace;

	TScoreValue ret_score =  _globalAlignment(trace, stringSet(align_), score_, tag_align_config, tag_algorithm);

	_pump_trace_2_Align(align_, trace);

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
