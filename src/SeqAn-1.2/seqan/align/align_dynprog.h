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

#ifndef SEQAN_HEADER_ALIGN_DYNPROG_H
#define SEQAN_HEADER_ALIGN_DYNPROG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//needleman wunsch alignment
template <typename TScoreValue, typename TMatrixSpec, typename TString>
TScoreValue
_needleman_wunsch(Matrix<TScoreValue, TMatrixSpec> & matrix_,
				  TString const & str1_,
				  TString const & str2_,
				  Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT

	typedef Matrix<TScoreValue, TMatrixSpec> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TValue;

	//-------------------------------------------------------------------------
	//define some variables
	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
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
	TScoreValue border_ = score_gap;
	TScoreValue v = border_;

	setDimension(matrix_, 2);
	setLength(matrix_, 0, str1_length + 1);
	setLength(matrix_, 1, str2_length + 1);
	resize(matrix_);

	TMatrixIterator col_ = end(matrix_) - 1;
	TMatrixIterator finger1;
	TMatrixIterator finger2;

	//-------------------------------------------------------------------------
	// init

	finger1 = col_;
	*finger1 = 0;
	for (x = x_end; x != x_begin; --x)
	{
		goPrevious(finger1, 0);
		*finger1 = border_;
		border_ += score_gap;
	}

	//-------------------------------------------------------------------------
	//fill matrix

	border_ = 0;
	for (y = y_end; y != y_begin; --y)
	{
		TValue cy = *y;

		h = border_;
		border_ += score_gap;
		v = border_;

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
			}
			*finger1 = v;
		}
	}

	return v;
}

//////////////////////////////////////////////////////////////////////////////
//traceback through needleman wunsch matrix
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TSourceSpec>
void
_needleman_wunsch_trace(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, TSourceSpec>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT
	typedef Iter<Matrix<TScoreValue, TSourceSpec>, PositionIterator > TMatrixIterator;
	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));

	typedef typename Position<Matrix<TScoreValue, TSourceSpec> >::Type TPosition;
	TPosition pos_0 = coordinate(source_, 0);
	TPosition pos_1 = coordinate(source_, 1);

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), pos_0);
	TTargetIterator target_1 = iter(row(target_, 1), pos_1);

	typedef typename Iterator<TTargetSourceSegment, Standard>::Type TStringIterator;
	TStringIterator it_0 = iter(str_0, pos_0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, pos_1);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		bool gv;
		bool gh;

		if (*it_0 == *it_1)
		{
			gv = gh = true;
		}
		else
		{

			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue v = *it_;

			goNext(it_, 1);
			TScoreValue d = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue h = *it_;

			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);

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
}

///////////////////////////////////////////////////////////////////////////////////////
//Gotoh
//Global alignment with affine gap costs
template <typename TScoreValue, typename TMatrixSpec, typename TString>
TScoreValue
_gotoh(Matrix<TScoreValue, TMatrixSpec> & diag_matrix_,
	   Matrix<TScoreValue, TMatrixSpec> & vert_matrix_,
	   Matrix<TScoreValue, TMatrixSpec> & hori_matrix_,
	   TString const & str1_,
	   TString const & str2_,
	   Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT


	typedef Matrix<TScoreValue, TMatrixSpec> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TValue;

	//-------------------------------------------------------------------------
	//define some variables
	TSize str1_length = length(str1_);
	TSize str2_length = length(str2_);
	TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

	TStringIterator x = x_end;
	TStringIterator y;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);

	TScoreValue border_ = score_gap_open;
	TScoreValue v;

	setDimension(diag_matrix_, 2);
	setLength(diag_matrix_, 0, str1_length + 1);
	setLength(diag_matrix_, 1, str2_length + 1);
	resize(diag_matrix_);
	setDimension(vert_matrix_, 2);
	setLength(vert_matrix_, 0, str1_length + 1);
	setLength(vert_matrix_, 1, str2_length + 1);
	resize(vert_matrix_);
	setDimension(hori_matrix_, 2);
	setLength(hori_matrix_, 0, str1_length + 1);
	setLength(hori_matrix_, 1, str2_length + 1);
	resize(hori_matrix_);

	TMatrixIterator diag_col_ = end(diag_matrix_) - 1;
	TMatrixIterator diag_finger1;
	TMatrixIterator diag_finger2;
	TMatrixIterator vert_col_ = end(vert_matrix_) - 1;
	TMatrixIterator vert_finger1;
	TMatrixIterator vert_finger2;
	TMatrixIterator hori_col_ = end(hori_matrix_) - 1;
	TMatrixIterator hori_finger1;
	TMatrixIterator hori_finger2;
	
	//-------------------------------------------------------------------------
	// init

	diag_finger1 = diag_col_;
	*diag_finger1 = 0;
	vert_finger1 = vert_col_;
	*vert_finger1 = -1000000;
	hori_finger1 = hori_col_;
	*hori_finger1 = -1000000;
	for (x = x_end; x != x_begin; --x)
	{
		goPrevious(diag_finger1, 0);
		*diag_finger1 = border_;
		goPrevious(hori_finger1, 0);
		*hori_finger1 = border_;
		goPrevious(vert_finger1, 0);
		*vert_finger1 = -1000000;//-inf
		border_ += score_gap_extend;
	}

	//------------------------------------------------------------------------
	//fill matrix

	border_ = score_gap_open;
	for (y = y_end; y != y_begin; --y)
	{
		TValue cy = *y;
		v = border_;

		vert_finger2 = vert_col_;	//points to last column
		goPrevious(vert_col_, 1);	//points to this column
		vert_finger1 = vert_col_;
		*vert_finger1 = v;          //initialize first column

		diag_finger2 = diag_col_;	
		goPrevious(diag_col_, 1);	
		diag_finger1 = diag_col_;
		*diag_finger1 = v;			//initialize first column

		hori_finger2 = hori_col_;	
		goPrevious(hori_col_, 1);	
		hori_finger1 = hori_col_;
		*hori_finger1 = -1000000;	//initialize first column		

		for (x = x_end; x != x_begin; --x)
		{
			//compute entry in diag_matrix
			goPrevious(diag_finger1, 0);
			v = (*diag_finger2 > *vert_finger2) ? *diag_finger2 : *vert_finger2;
			v = (*hori_finger2 > v) ? *hori_finger2 : v;
			if (*x == cy) *diag_finger1 = v + score_match;
			else *diag_finger1 = v + score_mismatch;

			//compute entry in hori_matrix
			v = *hori_finger1; 
			goPrevious(hori_finger1, 0);
			v += score_gap_extend;
			goPrevious(diag_finger2, 1);
			*hori_finger1 = (v > (*diag_finger2 + score_gap_open)) ? v : (*diag_finger2 + score_gap_open);
			goPrevious(hori_finger2, 0);

			//compute entry in vert_matrix
			goPrevious(vert_finger2, 0);
			goPrevious(vert_finger1, 0);
			v = *vert_finger2 + score_gap_extend;
			goNext(diag_finger2, 1);
			goPrevious(diag_finger2, 0);
			*vert_finger1 = (v > (*diag_finger2 + score_gap_open)) ? v : (*diag_finger2 + score_gap_open);
			
		}
		border_ += score_gap_extend;
	}

	v = (*vert_finger1 > *hori_finger1) ? *vert_finger1 : *hori_finger1;
	v = (*diag_finger1 > v) ? *diag_finger1 : v;

	return v;
}


//////////////////////////////////////////////////////////////////////////////
//gotoh trace
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TMatrixSpec>
void
_gotoh_trace(Align<TTargetSource, TTargetSpec> & target_,
			 Matrix<TScoreValue, TMatrixSpec> & diag_matrix_,
			 Matrix<TScoreValue, TMatrixSpec> & vert_matrix_,
			 Matrix<TScoreValue, TMatrixSpec> & hori_matrix_,
			 Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT
	typedef Iter<Matrix<TScoreValue, TMatrixSpec>, PositionIterator > TMatrixIterator;
	typedef typename Position<Matrix<TScoreValue, TMatrixSpec> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;

	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;
	typedef typename Iterator<TTargetSource, Standard>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(str_0) + 1;
	//typename Size<TTargetSourceSegment>::Type dim_1_len = length(str_1) + 1;
	
	TScoreValue score_gap_open = scoreGapOpen(score_);
	TScoreValue score_gap_diff = score_gap_open;
	//TScoreValue score_gap_extend = scoreGapExtend(score_);
	
	TMatrixIterator diag_source_ = begin(diag_matrix_);
	TMatrixIterator hori_source_ = begin(hori_matrix_);
	TMatrixIterator vert_source_ = begin(vert_matrix_);
	TPosition pos_0 = coordinate(diag_source_, 0);
	TPosition pos_1 = coordinate(diag_source_, 1);
	TPosition pos = position(diag_source_);

	TTargetIterator target_0 = iter(row(target_, 0), pos_0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), pos_1, Standard());

	TStringIterator it_0 = iter(str_0, pos_0, Standard());
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, pos_1, Standard());
	TStringIterator it_1_end = end(str_1);

	// indicate which matrix we are in
	bool hori = false, vert = false, diag = false;
	if (*diag_source_ > *hori_source_)
		if (*diag_source_ > *vert_source_) diag = true;
		else vert = true;
	else
		if (*hori_source_ > *vert_source_) hori = true;
		else vert = true;

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if(diag)
		{
			++it_0;
			++it_1;
			pos += dim_0_len + 1;

			if (getValue(diag_matrix_,pos) >= getValue(hori_matrix_,pos))
			{
				if (getValue(diag_matrix_,pos) < getValue(vert_matrix_,pos))
				{
					vert = true;
					diag = false;
				}
			}
			else
			{
				diag = false;
				if (getValue(hori_matrix_,pos) >= getValue(vert_matrix_,pos)) hori = true;
				else vert = true;
			}
		}
		else
		{
			if(vert)
			{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len;
				if (getValue(diag_matrix_,pos) + score_gap_diff >= getValue(vert_matrix_,pos))
				{
					diag = true; 
					vert = false;
				}
			}
			else
			{
				if(hori)
				{
					++it_0;
					insertGap(target_1);
					++pos;
					if (getValue(diag_matrix_,pos) + score_gap_diff >= getValue(hori_matrix_,pos))
					{
						diag = true; 
						hori = false;
					}
				}
			}
		}

		++target_0;
		++target_1;
	}





}



///////////////////////////////////////////////////////////////////////////
//if gap open == 0 regular needleman wunsch alignment, else gotoh alignment
/*DISABLED
.Function.needlemanWunsch:
..summary:Computes the best global alignment of the (two) sequences given in align according to the score values given in score.
..cat:Alignments
..signature:needlemanWunsch(align, score)
..param.align:The alignment object having the sequences to be aligned as sources.
...type:Class.Align
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..returns:The score value of the best scoring global alignment.
..returns:The corresponding alignment can be found in align.
..remarks:Depending on the Score object either the regular Needleman Wunsch algorithm (gap open = 0) or the Gotoh algorithm (gap open != 0) is applied.
..see:Function.smithWaterman
*/
/*
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT

	if(scoreGapOpen(score_)==scoreGapExtend(score_))
	{//linear gap costs
		return globalAlignment(align_, score_, NeedlemanWunsch());
	}
	else
	{//affine gap costs
		return globalAlignment(align_, score_, Gotoh());
	}
}


template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const & score_,
				NeedlemanWunsch)
{
SEQAN_CHECKPOINT
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	TScoreValue ret;

	Matrix<TScoreValue> matr;
	ret = _needleman_wunsch(matr, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_);
	_needleman_wunsch_trace(align_, begin(matr), score_);

	return ret;
}

template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const & score_,
				Gotoh)
{
SEQAN_CHECKPOINT
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	TScoreValue ret;

	Matrix<TScoreValue> d_matr;
	Matrix<TScoreValue> v_matr;
	Matrix<TScoreValue> h_matr;
	ret = _gotoh(d_matr, v_matr, h_matr, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_);
	_gotoh_trace(align_, d_matr, v_matr, h_matr, score_);	

	return ret;
}
*/

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
