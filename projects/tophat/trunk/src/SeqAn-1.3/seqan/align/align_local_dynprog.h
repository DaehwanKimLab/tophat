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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_ALIGN_LOCAL_DYNPROG_H
#define SEQAN_HEADER_ALIGN_LOCAL_DYNPROG_H

#include <functional>


namespace SEQAN_NAMESPACE_MAIN
{


//simple class that stores a value with an ID
template <typename TValue, typename TID>
class ScoreAndID 
{

public:
	TValue value_;
	TID id_;

	ScoreAndID()
	{
SEQAN_CHECKPOINT

	}

	ScoreAndID(TValue score, TID id_pos)
	{
SEQAN_CHECKPOINT
		value_ = score;
		id_ = id_pos;
	}
}; 



template <typename TValue, typename TID>
bool operator >(ScoreAndID<TValue,TID> & a,
				ScoreAndID<TValue,TID> & b)
{
	return (a.value_>b.value_);
}

template <typename TValue, typename TID>
bool operator >(const ScoreAndID<TValue,TID> & a,
				const ScoreAndID<TValue,TID> & b)
{
	return (a.value_>b.value_);
}



template <typename TValue, typename TID>
bool operator <(ScoreAndID<TValue,TID> & a,
				ScoreAndID<TValue,TID> & b)
{
	return (a.value_<b.value_);
}

template <typename TValue, typename TID>
bool operator <(const ScoreAndID<TValue,TID> & a,
				const ScoreAndID<TValue,TID> & b)
{
	return (a.value_<b.value_);
}



//////////////////////////////////////////////////////////////////////
/**
.Class.LocalAlignmentFinder:
..cat:Miscellaneous
..summary:Stores the information necessary for local alignment dynamic programming.
..signature:LocalAlignmentFinder<TScoreValue>
..param.TScoreValue:The value type that is used for scoring the alignments.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..see:Function.localAlignment
.Memfunc.LocalAlignmentFinder#LocalAlignmentFinder
..class:Class.LocalAlignmentFinder
..summary:Constructor
..signature:LocalAlignmentFinder(align)
..param.align:An @Class.Align@ object that is already initialized with the sequences.
..include:seqan/align.h
*/
template<typename TScoreValue = int>
class LocalAlignmentFinder{

public:
//____________________________________________________________________________

	typedef Matrix<TScoreValue> TMatrix;
	typedef typename Position<TMatrix>::Type TMatrixPosition;
    typedef typename Size<TMatrix>::Type TSize;
	typedef ScoreAndID<TScoreValue,TMatrixPosition> TPQEntry;

	typedef Iter<TMatrix,PositionIterator> TMatrixIterator;
	typedef PriorityType<TPQEntry> TPriorityQ;
	typedef String<bool> TBoolMatrix;

//____________________________________________________________________________
	
	//DP-matrix
	TMatrix matrix;
	//matrix that memorizes the cells from which not to go diagonal
	TBoolMatrix forbidden;
	//priority queue for quickly finding the maximum score in the DP-matrix
	TPriorityQ pQ;
	//position of maximum score (where traceback is started from) 
	TMatrixPosition bestEndPos;
	//position where traceback ended and where declumping begins
	TMatrixPosition bestBeginPos;
    //traceback path that is set to forbidden while declumping
    AlignTraceback<TSize> trace;

	bool needReinit; //true: call "smithWaterman", false: call "smithWatermanGetNext" 


//____________________________________________________________________________

	//LocalAlignmentFinder()
	//	: pQ(), matrix(), forbidden()
	//{
	//
	//}

	LocalAlignmentFinder():
		needReinit(true)
	{
SEQAN_CHECKPOINT
	}

    // TODO(holtgrew): Remove and replace all occurences with default constructor.
    template<typename TAlign>
	LocalAlignmentFinder(TAlign const &)
	    : needReinit(true)
	{
SEQAN_CHECKPOINT
    }

	~LocalAlignmentFinder()
	{
SEQAN_CHECKPOINT
	}


};
//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TTag>
void
_initLocalAlignmentFinder(TStringSet & str,
                          LocalAlignmentFinder<TScoreValue> & finder,
                          TTag) {
SEQAN_CHECKPOINT
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    typedef typename TFinder::TMatrix TMatrix;
    typedef typename Size<TMatrix>::Type TSize;

    TSize len0 = length(value(str, 0));
    TSize len1 = length(value(str, 1));

    setDimension(finder.matrix, 2);
    setLength(finder.matrix, 0, len0 + 1);
    setLength(finder.matrix, 1, len1 + 1);
    resize(finder.matrix);

    resize(finder.forbidden, (len0 + 1) * (len1 + 1), false);

	finder.bestEndPos = minValue<typename TFinder::TMatrixPosition>();
	finder.bestBeginPos = minValue<typename TFinder::TMatrixPosition>();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue>
void clear(LocalAlignmentFinder<TScoreValue> & sw_finder)
{
	sw_finder.needReinit = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue>
TScoreValue getScore(LocalAlignmentFinder<TScoreValue> & sw)
{
	typedef LocalAlignmentFinder<TScoreValue> TFinder;
	if(sw.bestEndPos !=  minValue<typename TFinder::TMatrixPosition>())
		return getValue(sw.matrix,sw.bestEndPos);
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//Smith-Waterman algorithm
template <typename TScoreValue, typename TString>
TScoreValue
_smithWatermanGetMatrix(LocalAlignmentFinder<TScoreValue> & sw,
						  TString const & str1_,
						  TString const & str2_,
						  Score<TScoreValue, Simple> const & score_,
						  TScoreValue cutoff)
{
SEQAN_CHECKPOINT

	// typedefs
	typedef Matrix<TScoreValue> TMatrix;
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Position<TMatrix>::Type TMatrixPosition;
	typedef Iter<TMatrix,PositionIterator> TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TValue;

	//-------------------------------------------------------------------------
	//define some variables


//	TSize str1_length = length(str1_);
//	TSize str2_length = length(str2_);
	TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

	TStringIterator x = x_end;
	TStringIterator y;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TScoreValue h = 0;
	TScoreValue v = 0;

	TMatrixIterator col_ = end(sw.matrix) - 1;
	TMatrixIterator finger1;
	TMatrixIterator finger2;

	//-------------------------------------------------------------------------
	// init

	finger1 = col_;
	*finger1 = 0;
	//std::cout <<"  ";
	for (x = x_end; x != x_begin; --x)
	{
		goPrevious(finger1, 0);
		*finger1 = 0;
	}

	//-------------------------------------------------------------------------
	//fill matrix
	for (y = y_end; y != y_begin; --y)
	{
		TValue cy = *y;

		h = 0;
		v = 0;

		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			
			if (*x == cy)
			{
				v = h + score_match;
				h = *finger2;
			}
			else
			{
				TScoreValue s1 = h + score_mismatch;
				h = *finger2;
				TScoreValue s2 = score_gap + ((h > v) ? h : v);
				v = (s1 > s2) ? s1 : s2;
				if (v < 0) v = 0;
			
			}
			*finger1 = v;
			if (v >= cutoff)
			{
				push(sw.pQ,ScoreAndID<TScoreValue,TMatrixPosition>(v,position(finger1)));
			}
		}
	}

	// check if any scores >= cutoff were found
	if(!empty(sw.pQ))
	{
        ScoreAndID<TScoreValue,TMatrixPosition> best = top(sw.pQ);
		v = getValue(sw.matrix,best.id_);
		sw.bestEndPos = best.id_;
	}
	else 
		v=0;

	return v;
}


///////////////////////////////////////////////////////////////////////////////
// declumping
template <typename TScoreValue, typename TSource, typename TSpec>
void
_smithWatermanDeclump(LocalAlignmentFinder<TScoreValue> & sw ,
				Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT

//typedefs
	typedef typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition TMatrixPosition;
	typedef Iter<typename LocalAlignmentFinder<TScoreValue>::TMatrix, PositionIterator> TMatrixIterator;
	
	typedef TSource TString;
	typedef typename Iterator<TString, Rooted>::Type TStringIterator;
	typedef typename Value<TString>::Type TValue;

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

//-------------------------------------------------------------------------
//variables
	TRow row0 = row(align_,0);
	TRow row1 = row(align_,1);

	TAlignIterator ali_it0_stop = iter(row0,beginPosition(row0));
	TAlignIterator ali_it1_stop = iter(row1,beginPosition(row1));

	SEQAN_ASSERT_TRUE( endPosition(row0)- beginPosition(row0) == endPosition(row1)- beginPosition(row1) );

	TAlignIterator ali_it0 = iter(row0,endPosition(row0));
	TAlignIterator ali_it1 = iter(row1,endPosition(row1));

	TStringIterator x_begin = begin(source(row0))-1; 
	TStringIterator y_begin = begin(source(row1))-1; 
	TStringIterator x_end = iter(source(row0),clippedEndPosition(row0))-1;
	TStringIterator y_end = iter(source(row1),clippedEndPosition(row1))-1;

	TStringIterator x = x_end;
	TStringIterator y = y_end;
	TStringIterator x_stop = x_end;
	
	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);
	TScoreValue h,v;

	TMatrixIterator finger0 = iter(sw.matrix,sw.bestBeginPos);
	TMatrixIterator end_col = finger0;
	TMatrixIterator finger1 = finger0;
	TMatrixIterator forbidden = finger0;

	bool different = true;
	bool forbidden_reached = true;
	bool end_reached = false;
	bool skip_row = false;
	

/*	int str0_length = length(source(row(align_,0)))+1;
	int str1_length = length(source(row(align_,1)))+1;
	for(int i = 0; i <str1_length; ++i){
 		for(int j=0;j<str0_length;++j){
 			std::cout << getValue(sw.matrix,(str0_length*i)+j);
 			if(sw.forbidden[(str0_length*i)+j]==true)
 				std::cout <<"(1) ";
 			else
 				std::cout <<"(0) ";
 		}
 		std::cout <<"\n";
 	}*/
	
	setClippedBeginPosition(row(align_, 0),0);
	setClippedBeginPosition(row(align_, 1),0);
	

	for (y = y_end; (y != y_begin) ; --y)
	{
		different = true;
		//compute next "forbidden" cell (where you are not allowed to go diagonal)
		if(forbidden_reached && !end_reached)
		{

			if(ali_it0==ali_it0_stop && ali_it1==ali_it1_stop)
			{
				end_reached = true;
			}

			if(!end_reached)
			{
				--ali_it0;
				goPrevious(forbidden,1);
				while(isGap(ali_it0)&& ali_it0!=ali_it0_stop)
				{
					skip_row = true;			
					--ali_it0;
					--ali_it1;
					goPrevious(forbidden,1);
				}

				--ali_it1;
				goPrevious(forbidden,0);
				while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
				{
					--ali_it1;
					--ali_it0;
					goPrevious(forbidden,0);
				}
				// mark the forbidden cell 
				sw.forbidden[position(forbidden)]=true;
				forbidden_reached = false;
			}

		}

		TValue cy = *y;

		h = *end_col;
		
		finger1 = end_col;		//points to last column
		goPrevious(end_col, 1);	//points to this column
		finger0 = end_col;

		v = *finger0;
		x = x_end;

		// declump the current row until the cells in the declumped and 
		// the cells in the original matrix are the same (or the border is reached)
		// indicated by bool different
		while (x != x_begin && different)
		{
			goPrevious(finger0, 0);
			goPrevious(finger1, 0);
			if (*x == cy && !(sw.forbidden[position(finger0)]))
			{
				v = h + score_match;
				h = *finger1;
			}
			else
			{
				TScoreValue s1;
			
				if(finger0 == forbidden)
				{
						skip_row = false;
						forbidden_reached = true;
						s1 = 0;
				}
				else
				{
					if(sw.forbidden[position(finger0)]) s1 = 0;
					else s1 = h + score_mismatch;
				}

				h = *finger1;
				TScoreValue s2 = score_gap + ((h > v) ? h : v);
				v = (s1 > s2) ? s1 : s2;
				if (v < 0) v = 0;
			
			}

			// value is the same as in the original matrix
			if(*finger0==v)
			{
				//x_stop is as far as we have to go at least
				if(x<x_stop)
				{
					different=false;
			//		x_stop=x;
				}
				else
				{
					// x_end indicates where to start in the next row
					if(x==x_end && ((!forbidden_reached && !skip_row)||end_reached))
					{
						--x_end;
						goPrevious(end_col, 0);	
					}
				}
			} 
			if(x<x_stop)// && different)
			{
				x_stop=x;
			}
			*finger0 = v;
			--x;
		}
		if(x_end==x_begin)
			break;
	}


//	cout <<"...declumped.\n";



}




//////////////////////////////////////////////////////////////////////////////
//traceback
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
typename Iterator<Matrix<TScoreValue, DIMENSION>, Standard >::Type
_smithWatermanTrace(Align<TTargetSource, TTargetSpec> & target_,
					 typename LocalAlignmentFinder<TScoreValue>::TBoolMatrix & fb_matrix, 
					 Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
					 Score<TScoreValue, Simple> const & scoring_) {
SEQAN_CHECKPOINT
	//typedefs
	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

//	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;
	typedef typename Iterator<TTargetSource, Standard>::Type TStringIterator;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TTargetIterator;

	//-------------------------------------------------------------------------
	//variables
	TPosition pos_0 = coordinate(source_, 0);
	TPosition pos_1 = coordinate(source_, 1);

	TTargetSource str_0 = source(row(target_, 0));
	TTargetSource str_1 = source(row(target_, 1));
	
	TTargetIterator target_0 = iter(row(target_, 0), pos_0);
	TTargetIterator target_1 = iter(row(target_, 1), pos_1);

	TStringIterator it_0 = iter(str_0, pos_0, Standard());
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, pos_1, Standard());
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_mismatch = scoreMismatch(scoring_);
	TScoreValue score_gap = scoreGapExtend(scoring_);

	//-------------------------------------------------------------------------
	//follow the trace until 0 is reached
	while ((*source_!=0) && (it_0 != it_0_end) && (it_1 != it_1_end))
	{
		bool gv;
		bool gh;
		bool forbidden = fb_matrix[position(source_)];

		if (*it_0 == *it_1 && !forbidden)
		{
			gv = gh = true;
		}
		else
		{
			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue v = *it_ + score_gap;

			TScoreValue d;
			if(forbidden)
				d = 0;
			else{
				goNext(it_, 1);
				d = *it_ + score_mismatch;
			}

			it_ = source_;
			goNext(it_, 1);
			TScoreValue h = *it_ + score_gap;

			gv = (v >= h) | (d >= h);
			gh = (h >  v) | (d >= v);
		}

		if (gv)
		{
			++it_0;
			goNext(source_, 0);
		}
		else
		{
			insertGap(target_0);
		}

		if (gh) 
		{
			++it_1;
			goNext(source_, 1);
		}
		else
		{
			insertGap(target_1);
		}
		++target_0;
		++target_1;
	}

	setClippedBeginPosition(row(target_, 0),pos_0);
	setClippedBeginPosition(row(target_, 1),pos_1);
	setBeginPosition(row(target_, 0),0);
	setBeginPosition(row(target_, 1),0);
	setClippedEndPosition(row(target_, 0),position(it_0, str_0));
	setClippedEndPosition(row(target_, 1),position(it_1, str_1));
	
	return source_;

}

/////////////////////////////////////////////////////////////////////
//adjust the priority queue of scores until the true maximum is found
template <typename TScoreValue>
typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition
_getNextBestEndPosition(LocalAlignmentFinder<TScoreValue> & sw ,
                            TScoreValue cutoff) {
SEQAN_CHECKPOINT
    // get maximal score from priority queue
	TScoreValue topScore = 0;
    if (!empty(sw.pQ)) {
        topScore = getValue(sw.matrix, top(sw.pQ).id_);
    }

    // check if matrix entry of topScore did not change while declumping
	if (!empty(sw.pQ)) {
		while (top(sw.pQ).value_ != topScore) {
			if (topScore >= cutoff) {
				((sw.pQ).heap[0]).value_ = topScore;
				adjustTop(sw.pQ); 
			} else {
				pop(sw.pQ);
			}
			if (!empty(sw.pQ)) topScore = getValue(sw.matrix, top(sw.pQ).id_);
			else break;
		}
	}

    // priority queue with top scores is empty
    if(empty(sw.pQ)) {//||top(sw.pQ).value_<cutoff) {
		sw.needReinit = true;
		return 0;
	}

	typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition ret_pos = top(sw.pQ).id_;
	sw.bestEndPos = ret_pos;
	pop(sw.pQ);
	
	return ret_pos;
}

/*DISABLED
.Function.smithWaterman:
..summary:Computes the best local alignment of the (two) sequences given in align.
..cat:Alignments
..signature:smithWaterman(align, sw_finder, score, cutoff)
..param.align:The alignment object having the sequences to be aligned as sources.
...type:Class.Align
..param.sw_finder:The local alignment finder object.
...type:Class.LocalAlignmentFinder
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.cutoff:A score limit.
...remarks:Alignments with scores < cutoff will be discarded (starTCGCGCCTGGAC
CAAGACA
ts being useful when
	not only the best, but also sub-optimal local aligments are of interest). 
	Low cutoff scores cause a higher runtime.
..returns:The score value of the best scoring local alignment or 0 if there was no alignment with score >= cutoff.
...param.align:The corresponding alignment.
..remarks:So far, only linear gap costs are allowed.
..see:Function.smithWatermanGetNext
..see:Function.localAlignment
..include:seqan/align.h
*/
///////////////////////////////////////////////////////////////////////////////////
//wrapper that computes the matrix and does the backtracking for the best alignment
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
_smithWaterman(Align<TSource, TSpec> & align_,
			    LocalAlignmentFinder<TScoreValue> & sw_finder,
			    Score<TScoreValue, Simple> const & score_, 
			    TScoreValue cutoff)
{
SEQAN_CHECKPOINT
	typedef LocalAlignmentFinder<TScoreValue> TFinder;
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

    StringSet<TSource> str;
    for (unsigned i = 0; i < length(rows(align_)); ++i) {
        appendValue(str, sourceSegment(row(align_, i)));
    }
    _initLocalAlignmentFinder(str, sw_finder, SmithWaterman());
	
	TScoreValue ret = _smithWatermanGetMatrix(sw_finder, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_,cutoff);
	
	if(ret==0)
		return ret;
	sw_finder.needReinit = false;

	typedef Iter<typename LocalAlignmentFinder<TScoreValue>::TMatrix,PositionIterator > TMatrixIterator;
	TMatrixIterator best_begin;
	
	// TODO: sw_finder statt kram
	best_begin = _smithWatermanTrace(align_,sw_finder.forbidden,iter(sw_finder.matrix,(top(sw_finder.pQ)).id_), score_);

	sw_finder.bestBeginPos = position(best_begin);
	
	pop(sw_finder.pQ);

	return ret;
}

/*DISABLED
.Function.smithWatermanGetNext:
..summary:Computes next best local alignment.
..description:Declumps the matrix and computes the next best local alignment of the (two) sequences given in align 
according to the score values given in score. Declumping means forbidding all matches and mismatches that were used
in previously found alignment to be used again.
..cat:Alignments
..signature:smithWatermanGetNext(align, sw_finder, score, cutoff)
..param.align:The alignment object having the sequences to be aligned as sources.
...type:Class.Align
..param.sw_finder:The local alignment finder object (that has been passed to @Function.smithWaterman@ before).
...type:Class.LocalAlignmentFinder
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.cutoff:Alignments with scores < cutoff will be discarded. 
...remarks:Only use a cut off score that is greater or equal to the one that was used when calling @Function.smithWaterman@.
..returns:The score value of the next best local alignment or 0 if there was no alignment with score >= cutoff.
..returns:The corresponding alignment can be found in align.
..see:Function.smithWaterman
..see:Function.localAlignment
..include:seqan/align.h
*/
///////////////////////////////////////////////////////////////////////////
// wrapper that declumps the matrix and traces back the next best alignment
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
_smithWatermanGetNext(Align<TSource, TSpec> & align_,
					 LocalAlignmentFinder<TScoreValue> & sw_finder ,
					 Score<TScoreValue, Simple> const & score_, 
					 TScoreValue cutoff)
{	
SEQAN_CHECKPOINT

	_smithWatermanDeclump(sw_finder, align_, score_);

	clearGaps(row(align_,0));
	clearGaps(row(align_,1));
	setClippedEndPosition(row(align_, 0),endPosition(source(row(align_,0))));
	setClippedEndPosition(row(align_, 1),endPosition(source(row(align_,1))));

	typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition next_best_end;
	next_best_end = _getNextBestEndPosition(sw_finder,cutoff);
	if(next_best_end==0)
		return 0;
	typename LocalAlignmentFinder<TScoreValue>::TMatrixIterator next_best_begin;
	next_best_begin= _smithWatermanTrace(align_,sw_finder.forbidden,iter(sw_finder.matrix,next_best_end), score_);
	sw_finder.bestBeginPos = position(next_best_begin);
	
	return getValue(sw_finder.matrix,next_best_end);
}

//////////////////////////////////////////////////////////////////////////////
//interface for Function.localAlignment

/**
.Function.localAlignment:
..cat:Alignments
...type:Class.Align
..signature:localAlignment(align, score, tag)
..signature:localAlignment(align, score)
..signature:localAlignment(align, sw_finder, score, cutoff, tag)
..signature:localAlignment(align, sw_finder, score, cutoff, diagLow, diagHigh, tag)
..param.align:Alignment object to use.
...type:Class.Align
..param.cutoff:Alignments with scores < cutoff will be discarded.
..param.diagLow:The lowest diagonal of the alignment matrix that will be computed for banded alignment.
..param.diagHigh:The highest diagonal of the alignment matrix that will be computed for banded alignment.
..param.tag:
...type:Tag.Local Alignment Algorithms.value.SmithWaterman
...type:Tag.Local Alignment Algorithms.value.WatermanEggert
...type:Tag.Local Alignment Algorithms.value.BandedWatermanEggert
.remarks:TODO
..include:seqan/align.h
 */	

//1. only Align object
template <typename TSource, typename TSpec, typename TScoreValue>
inline TScoreValue
localAlignment(Align<TSource, TSpec> & align_,
			   Score<TScoreValue, Simple> const & score_, 
			   SmithWaterman)
{
	LocalAlignmentFinder<TScoreValue> sw_finder(align_);

	return _smithWaterman(align_, sw_finder, score_, 0);
}

template <typename TSource, typename TSpec, typename TScoreValue>
inline TScoreValue
localAlignment(Align<TSource, TSpec> & align_,
			   Score<TScoreValue, Simple> const & score_)
{
	return localAlignment(align_, score_, SmithWaterman());
}


//2. Align, LocalAlignmentFinder, and cutoff arguments
template <typename TSource, typename TSpec, typename TScoreValue1, typename TScoreValue2, typename TScoreValue3>
inline TScoreValue1
localAlignment(Align<TSource, TSpec> & align_,
			   LocalAlignmentFinder<TScoreValue1> & sw_finder,
			   Score<TScoreValue2, Simple> const & score_, 
			   TScoreValue3 cutoff,
			   WatermanEggert)
{
	if (sw_finder.needReinit)
	{
		return _smithWaterman(align_, sw_finder, score_, cutoff);
	}
	else
	{
		return _smithWatermanGetNext(align_, sw_finder, score_, cutoff);
	}
}
template <typename TSource, typename TSpec, typename TScoreValue1, typename TScoreValue2, typename TScoreValue3>
inline TScoreValue1
localAlignment(Align<TSource, TSpec> & align_,
			   LocalAlignmentFinder<TScoreValue1> & sw_finder,
			   Score<TScoreValue2, Simple> const & score_, 
			   TScoreValue3 cutoff)
{
	return localAlignment(align_, sw_finder, score_, cutoff, WatermanEggert());
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
