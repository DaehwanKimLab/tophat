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
  $Id: align_local_dynprog.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

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
..see:Function.smithWaterman
*/
template<typename TScoreValue = int>
class LocalAlignmentFinder{

public:
//____________________________________________________________________________

	typedef Matrix<TScoreValue> TMatrix;
	typedef typename Position<TMatrix>::Type TMatrixPosition;
	typedef ScoreAndID<TScoreValue,TMatrixPosition> TPQEntry;

	typedef Iter<TMatrix,PositionIterator> TMatrixIterator;
	typedef PriorityType<TPQEntry> TPriorityQ;
	typedef String<bool> TBoolMatrix;

//____________________________________________________________________________
	
	//DP-matrix
	TMatrix matrix_;
	//matrix that memorizes the cell from which not to go diagonal
	TBoolMatrix forbidden_;
	//priority queue for quickly finding the maximum score in the DP-matrix
	TPriorityQ pq_;
	//position of maximum score (where traceback is started from) 
	TMatrixPosition best_end_pos_;
	//position where traceback ended and where declumping begins
	TMatrixPosition best_begin_pos_;

	bool _needReinit; //true: call "smithWaterman", false: call "smithWatermanGetNext" 


//____________________________________________________________________________

	//LocalAlignmentFinder()
	//	: pq_(), matrix_(), forbidden_()
	//{
	//
	//}

	template <typename TSource,typename TSpec>
	LocalAlignmentFinder(Align<TSource,TSpec> & align_):
		_needReinit(true)
	{
SEQAN_CHECKPOINT

		typedef typename Size<TMatrix>::Type TSize;

		TSize str0_length = length(sourceSegment(row(align_,0)));
		TSize str1_length = length(sourceSegment(row(align_,1)));

		setDimension(matrix_, 2);
		setLength(matrix_, 0, str0_length + 1);
		setLength(matrix_, 1, str1_length + 1);
		resize(matrix_);

		fill(forbidden_,(str0_length + 1)*(str1_length + 1),false);
	}

	~LocalAlignmentFinder()
	{
SEQAN_CHECKPOINT
	}


};
//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue>
void clear(LocalAlignmentFinder<TScoreValue> & sw_finder)
{
	sw_finder._needReinit = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue>
TScoreValue getScore(LocalAlignmentFinder<TScoreValue> & sw)
{
	if(!empty(sw.pq_))
	{
		return getValue(sw.matrix_, sw.best_end_pos_);
	}
	else 
	{
		return 0;
	}
}

//////////////////////////////////////////////////////////////////////////////
//Smith-Waterman algorithm
template <typename TScoreValue, typename TString>
TScoreValue
smith_waterman_get_matrix(LocalAlignmentFinder<TScoreValue> & sw,
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

	TMatrixIterator col_ = end(sw.matrix_) - 1;
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
			if (v > cutoff)
			{
				push(sw.pq_,ScoreAndID<TScoreValue,TMatrixPosition>(v,position(finger1)));
			}
		}
	}

	// check if any scores > cutoff were found
	if(!empty(sw.pq_))
	{
        ScoreAndID<TScoreValue,TMatrixPosition> best = top(sw.pq_);
		v = getValue(sw.matrix_,best.id_);
		sw.best_end_pos_ = best.id_;
	}
	else 
		v=0;

	return v;
}


///////////////////////////////////////////////////////////////////////////////
// declumping
template <typename TScoreValue, typename TSource, typename TSpec>
void
smith_waterman_declump(LocalAlignmentFinder<TScoreValue> & sw ,
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

	SEQAN_TASSERT( endPosition(row0)- beginPosition(row0) == endPosition(row1)- beginPosition(row1) )

	TAlignIterator ali_it0 = iter(row0,endPosition(row0));
	TAlignIterator ali_it1 = iter(row1,endPosition(row1));

	TStringIterator x_begin = begin(source(row0))-1; 
	TStringIterator y_begin = begin(source(row1))-1; 
	TStringIterator x_end = iter(source(row0),sourceEndPosition(row0))-1;
	TStringIterator y_end = iter(source(row1),sourceEndPosition(row1))-1;

	TStringIterator x = x_end;
	TStringIterator y = y_end;
	TStringIterator x_stop = x_end;
	
	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);
	TScoreValue h,v;

	TMatrixIterator finger0 = iter(sw.matrix_,sw.best_begin_pos_);
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
 			std::cout << getValue(sw.matrix_,(str0_length*i)+j);
 			if(sw.forbidden_[(str0_length*i)+j]==true)
 				std::cout <<"(1) ";
 			else
 				std::cout <<"(0) ";
 		}
 		std::cout <<"\n";
 	}*/
	
	setSourceBeginPosition(row(align_, 0),0);
	setSourceBeginPosition(row(align_, 1),0);
	

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
				sw.forbidden_[position(forbidden)]=true;
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
			if (*x == cy && !(sw.forbidden_[position(finger0)]))
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
					if(sw.forbidden_[position(finger0)]) s1 = 0;
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
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TSourceSpec>
typename Iterator<Matrix<TScoreValue, TSourceSpec>, Standard >::Type
smith_waterman_trace(Align<TTargetSource, TTargetSpec> & target_,
					 typename LocalAlignmentFinder<TScoreValue>::TBoolMatrix & fb_matrix, 
					 Iter< Matrix<TScoreValue, TSourceSpec>, PositionIterator > source_,
					 Score<TScoreValue, Simple> const &)
{
SEQAN_CHECKPOINT

	//typedefs
	typedef Iter<Matrix<TScoreValue, TSourceSpec>, PositionIterator > TMatrixIterator;
	typedef typename Position<Matrix<TScoreValue, TSourceSpec> >::Type TPosition;

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
			TScoreValue v = *it_;

			TScoreValue d;
			if(forbidden)
				d = 0;
			else{
				goNext(it_, 1);
				d = *it_;
			}

			it_ = source_;
			goNext(it_, 1);
			TScoreValue h = *it_;

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

	setSourceBeginPosition(row(target_, 0),pos_0);
	setSourceBeginPosition(row(target_, 1),pos_1);
	setBeginPosition(row(target_, 0),0);
	setBeginPosition(row(target_, 1),0);
	setSourceEndPosition(row(target_, 0),position(it_0, str_0));
	setSourceEndPosition(row(target_, 1),position(it_1, str_1));
	
	return source_;

}

/////////////////////////////////////////////////////////////////////
//adjust the priority queue of scores until the true maximum is found
template <typename TScoreValue>
typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition
get_next_best_end_position(LocalAlignmentFinder<TScoreValue> & sw ,
				 TScoreValue cutoff)
{
SEQAN_CHECKPOINT

	TScoreValue top_score = 0;
	if(!empty(sw.pq_))
	        top_score = getValue(sw.matrix_,top(sw.pq_).id_);

	while (!empty(sw.pq_) && (top(sw.pq_).value_ != top_score))
	{
		if (top_score > cutoff)
		{
			((sw.pq_).heap[0]).value_ = top_score;
			adjustTop(sw.pq_); 
		}
		else
			pop(sw.pq_);

		top_score = getValue(sw.matrix_,top(sw.pq_).id_);
	}

	if(empty(sw.pq_))//||top(sw.pq_).value_<cutoff)
	{
		sw._needReinit = true;
		return 0;
	}

	typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition ret_pos = top(sw.pq_).id_;
	sw.best_end_pos_ = ret_pos;
	pop(sw.pq_);
	
	return ret_pos;

}

/**
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
..returns:The score value of the best scoring local alignment or 0 if there was no alignment with score > cutoff.
...param.align:The corresponding alignment.
..remarks:So far, only linear gap costs are allowed.
..see:Function.smithWatermanGetNext
..see:Function.localAlignment
*/
///////////////////////////////////////////////////////////////////////////////////
//wrapper that computes the matrix and does the backtracking for the best alignment
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
smithWaterman(Align<TSource, TSpec> & align_,
			  LocalAlignmentFinder<TScoreValue> & sw_finder ,
			  Score<TScoreValue, Simple> const & score_, 
			  TScoreValue cutoff)
{
SEQAN_CHECKPOINT
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	TScoreValue ret = smith_waterman_get_matrix(sw_finder, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_,cutoff);
	
	if(ret==0)
		return ret;

	sw_finder._needReinit = false;

	typedef Iter<typename LocalAlignmentFinder<TScoreValue>::TMatrix,PositionIterator > TMatrixIterator;
	TMatrixIterator best_begin;
	
	//sw_finder statt kram
	best_begin = smith_waterman_trace(align_,sw_finder.forbidden_,iter(sw_finder.matrix_,(top(sw_finder.pq_)).id_), score_);
 					
	sw_finder.best_begin_pos_ = position(best_begin);
	
	pop(sw_finder.pq_);

	return ret;
}

///////////////////////////////////////////////////////////////////////////
// wrapper that declumps the matrix and traces back the next best alignment
/**
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
..returns:The score value of the next best local alignment or 0 if there was no alignment with score > cutoff.
..returns:The corresponding alignment can be found in align.
..see:Function.smithWaterman
..see:Function.localAlignment
*/
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
smithWatermanGetNext(Align<TSource, TSpec> & align_,
					 LocalAlignmentFinder<TScoreValue> & sw_finder ,
					 Score<TScoreValue, Simple> const & score_, 
					 TScoreValue cutoff)
{	
SEQAN_CHECKPOINT

	smith_waterman_declump(sw_finder, align_, score_);

	clearGaps(row(align_,0));
	clearGaps(row(align_,1));
	setSourceEndPosition(row(align_, 0),endPosition(source(row(align_,0))));
	setSourceEndPosition(row(align_, 1),endPosition(source(row(align_,1))));

	typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition next_best_end;
	next_best_end = get_next_best_end_position(sw_finder,cutoff);
	if(next_best_end==0)
		return 0;
	typename LocalAlignmentFinder<TScoreValue>::TMatrixIterator next_best_begin;
	next_best_begin= smith_waterman_trace(align_,sw_finder.forbidden_,iter(sw_finder.matrix_,next_best_end), score_);
	sw_finder.best_begin_pos_ = position(next_best_begin);
	
	return getValue(sw_finder.matrix_,next_best_end);
}

//////////////////////////////////////////////////////////////////////////////
//interface for Function.localAlignment

//1. only Align object

template <typename TSource, typename TSpec, typename TScoreValue>
inline TScoreValue
localAlignment(Align<TSource, TSpec> & align_,
			   Score<TScoreValue, Simple> const & score_, 
			   SmithWaterman)
{
	LocalAlignmentFinder<TScoreValue> sw_finder(align_);

	return smithWaterman(align_, sw_finder, score_, 0);
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
	if (sw_finder._needReinit)
	{
		return smithWaterman(align_, sw_finder, score_, cutoff);
	}
	else
	{
		return smithWatermanGetNext(align_, sw_finder, score_, cutoff);
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
