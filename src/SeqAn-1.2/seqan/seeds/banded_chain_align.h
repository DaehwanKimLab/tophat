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
  $Id: banded_chain_align.h 3361 2009-02-05 16:38:14Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BANDED_CHAIN_ALIGN_H
#define SEQAN_HEADER_BANDED_CHAIN_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Function.bandedChainAlignment:
..summary: Calculates a banded alignment around a chain of seeds. 
..cat:Seed Handling
..signature:bandedChainAlignment(seedChain, k, alignment, scoreMatrix);
..param.seedChain:A chain of seeds.
..param.k:Half of the width of the band.
..param.alignment:The alignment where the result is stored.
...type:Class.Align
..param.scoreMatrix: The score matrix.
...type:Spec.Simple Score
...remarks: Depending on the score matrix the Needleman-Wunsch or the Gotoh algorithm is used.
..returns: The score of the alignment.
*/

template<typename TContainer, typename TValue, typename TScore, typename TAlign>
TScore
bandedChainAlignment(TContainer const &seedChain, 
					 TValue k,
					 TAlign &whole_alignment, 
					 Score< TScore, Simple> const &scoreMatrix)
{

	//teste die Reihenfolge der seeds und erstelle gegebenenfalls eine die rueckwaerts laeuft
	typedef typename Iterator<const TContainer, Standard>::Type TIterator;	//Iterator for the seed chain
	TIterator it = begin(seedChain);
	if (it != end(seedChain))
	{
		TIterator it2 = it; ++it2;
		if (it2 != end(seedChain))
		{
			if (leftDim0(*it) < leftDim0(*it2))
			{
				TContainer copyChain(seedChain);
				std::reverse(begin(copyChain), end(copyChain));

				if(scoreGapOpen(scoreMatrix)==scoreGapExtend(scoreMatrix))
					return chain_to_alignment_needlemanwunsch(copyChain, k, whole_alignment, scoreMatrix);
				else
					return chain_to_alignment_gotoh(copyChain, k, whole_alignment, scoreMatrix);
			}
		}
	}

	if(scoreGapOpen(scoreMatrix)==scoreGapExtend(scoreMatrix))
		return chain_to_alignment_needlemanwunsch(seedChain, k, whole_alignment, scoreMatrix);
	else
		return chain_to_alignment_gotoh(seedChain, k, whole_alignment, scoreMatrix);
}

template<typename TContainer, typename TValue, typename TScore, typename TAlign>
TScore
chain_to_alignment_needlemanwunsch(TContainer const &seedChain, 
						TValue k,
						TAlign & whole_alignment, 
						Score< TScore, Simple> const &scoreMatrix)
{
	typedef typename Value<TAlign>::Type TChar;
	typedef  String<TChar> TString;									//Sequence in an align object
	typedef typename Size<TString>::Type TSize;								//Size of the string
	typedef typename Iterator<const TContainer, Standard>::Type TIterator;	//Iterator for the seed chain
	typedef typename Value<TContainer>::Type TSeed;							//Type of Seed
	typedef typename Value<TSeed>::Type TPosition;
	typedef String<TScore> TScoreString;
	typedef Matrix<TScore> TMatrix;
	typedef Iter<TMatrix, PositionIterator > TMatrixIterator;
	typedef ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > TAlignVector;
	typedef typename ::std::map<TValue, Pair<TValue, TAlign> >::iterator TMapIterator;
	
	TScoreString score_str;
	TAlignVector alignmentVector;
	TString *p_seq1 = &source(row(whole_alignment,0));
	TString *p_seq2 = &source(row(whole_alignment,1));
	
	TValue score_length = 0;
	TMatrix matrix_;
	TIterator it = begin(seedChain);

	//calculation begins at the end
	//rectangle between last seed and end
	TValue k_begin = k;
	TValue k_end = k;
	if ((length(*it) <= k) ||((rightDim1(*it)-leftDim1(*it)) <=k))
 		k_end = length(*it) -1;

	_calculateLastRectangle(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoreMatrix);
	_calculateBandedSeed(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoreMatrix);

	TIterator it_begin = end(seedChain);
	TIterator it2 = it;
	++it2;
	while (it2 != it_begin)
	{
		k_begin = k_end;
		if ((length(*it2) <= k) || ((rightDim1(*it2)-leftDim1(*it2)) <=k))
			k_end = length(*it2) -2;
		else 
			k_end = k;
		
		_calculateRectangle(*it, *it2, k_begin, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoreMatrix);	
		_calculateBandedSeed(*it2, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoreMatrix);
		++it;
		++it2;
	}
	_calculateFirstRectangle(*it, k_end, matrix_, p_seq1, p_seq2, score_str, score_length, alignmentVector, scoreMatrix);

	_constructAlignment(alignmentVector, whole_alignment);
	return score_str[0];
}

//"Glues" single alignments together
template<typename TValue, typename TAlign, typename TAlign2>
void
_constructAlignment(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > >const &me,
					TAlign2 &wholeAlignment)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::const_iterator TIterator;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	typedef typename Row<TAlign2>::Type TRow2;
	typedef typename Iterator<TRow2, Standard>::Type TTargetIterator2;

	TTargetIterator2 align_it0 = iter(row(wholeAlignment, 0), 0);
	TTargetIterator2 align_it1 = iter(row(wholeAlignment, 1), 0);
	TValue position = 0;
	int length_ = me.size()-1;
	TIterator it;

	typedef typename Position<TAlign>::Type TPosition;
	for (int i = length_; i != -1; --i)
	{
		//cout << "LENGHT: " << me[i].size() << endl;
		it = me[i].begin();//find(position);
		position = it->second.i1;
		TRow tmp_row0 = row(it->second.i2, 0);
		TRow tmp_row1 = row(it->second.i2, 1);
		TPosition end_ = endPosition(cols(it->second.i2));
		unsigned int j = 0;
		for (TPosition it_ = beginPosition(cols(it->second.i2)); it_ != end_; ++it_)
		{
			if(isGap(tmp_row0, j))
				insertGap(align_it0);
			if(isGap(tmp_row1, j))
				insertGap(align_it1);
			++align_it0;
			++align_it1;
			++j;
		}
	}
}

template <typename TScoreValue, typename TMatrixSpec, typename TString, typename TValue, typename TValue2, typename TSpecSeed>
TScoreValue
_banded_needleman_wunsch(Matrix<TScoreValue, TMatrixSpec> & matrix_,
						 Seed<TValue, TSpecSeed> const &seed,
						 TValue2 k,
						 TString const & str1_,
						 TString const & str2_,
						 Score<TScoreValue, Simple> const & score_,
						 String<TScoreValue> init)
{
	SEQAN_CHECKPOINT
	typedef Matrix<TScoreValue, TMatrixSpec> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables
	TSize height = leftDiagonal(seed) - rightDiagonal(seed)+1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TStringIterator x_begin = begin(str1_);
	TStringIterator x_end = end(str1_) -1;
	TStringIterator y_begin = begin(str2_);
	--y_begin;
	TStringIterator y_end = end(str2_)-1;
	
	TValue up_height = leftDiagonal(seed) - startDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TValue down_height = endDiagonal(seed) - rightDiagonal(seed) + k; //equals length of empty lower triangle
	TValue down_width = str1_length - down_height;
	TSize length_right_diag = str2_length - down_height;
	TSize length_left_diag = str2_length - up_height;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TScoreValue horizontalValue = 0;
	TScoreValue border_ = score_gap;
	TScoreValue verticalValue = border_;

	setDimension(matrix_, 2);
	setLength(matrix_, 0, str1_length + 2);
	setLength(matrix_, 1, str2_length + 1);
	resize(matrix_);
	
	TMatrixIterator col_ = begin(matrix_);
	TMatrixIterator finger1;
	TMatrixIterator finger2;
	TMatrixIterator finger3;

	TStringIterator x = x_begin;
	TStringIterator y = y_begin;

	//-------------------------------------------------------------------------
	// init

	TValue inf = -1000000;
	finger2 = col_;

	setPosition(finger2, length(matrix_, 0)-1);
	for (unsigned int i = 1; i != length_right_diag; ++i){
		*finger2 = inf;
		goNext(finger2, 1);
	}

	TValue pos = -1;

	*finger2 = inf;
	for (int i = -1; i != down_height; ++i){
		goPrevious(finger2, 0);
		goNext(finger2,1);
		*finger2 = init[++pos];
	}
		
	finger3 = finger2;

	for (int i = 0; i != down_width; ++i){
		goPrevious(finger2, 0);
		*finger2 = init[++pos];
	}

	*finger2 = inf;
	TScoreValue tmp;
	//-------------------------------------------------------------------------
	//first section
	
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

	for (int i = -1; i != down_height; ++i){	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		goNext(finger3,0);
		horizontalValue = *finger3;
		x = x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure < length_left_diag)
			++run_length;
		--y;
	}
	--run_length;

	goPrevious(finger3);

	while (y!= y_begin)
	{	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		x = --x_end;
		horizontalValue = inf;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure >= length_left_diag)
			--run_length;
		--y;
	}
	return *(++finger2);

}


//////////////////////////////////////////////////////////////////////////////
//traceback through needleman wunsch matrix


//Position berechnen!
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TSourceSpec>
typename Size<Matrix<TScoreValue, TSourceSpec> >::Type
_banded_needleman_wunsch_trace2(Align<TTargetSource, TTargetSpec> & target_,
								Matrix<TScoreValue, TSourceSpec> & matrix,
								Iter< Matrix<TScoreValue, TSourceSpec>, PositionIterator > source_,
								Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT
	typedef Iter<Matrix<TScoreValue, TSourceSpec>, PositionIterator > TMatrixIterator;
	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));


	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);



	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;
	TStringIterator it_0 = begin(str_0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = begin(str_1);
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
			TScoreValue h = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue d = *it_;
		
			goPrevious(it_, 0);
			TScoreValue v = *it_;

			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);

		}

	
		if (gv){
			++it_1;
			goNext(source_, 1);
			goPrevious(source_, 0);
		}else{
			insertGap(target_1);
		}

		if (gh) {
			++it_0;
			goNext(source_, 0);
	
		}else{
			insertGap(target_0);
		}
		++target_0;
		++target_1;
	}

	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));

	return length(matrix,0) - coordinate(source_,0) -2;
}


//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TSourceSpec>
int
_needleman_wunsch_trace_lastRectangle(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, TSourceSpec>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT
	typedef Iter<Matrix<TScoreValue, TSourceSpec>, PositionIterator > TMatrixIterator;
	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);

	typedef typename Iterator<TTargetSourceSegment, Standard>::Type TStringIterator;
	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
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
	return 0;
}

template<typename TAlign, typename TValue>
void
_rec_delete(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &vec,	//alignment vector
		   TValue index,										//position im vector
		   TValue position)										//alignment to delete
{
	if (position != -1)
	{
		typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator it, it1, it2;
		it = vec[index].find(position);
		bool x = true; //true = successor of it can be deleted
		if (it != vec[index].begin())
		{
			it1 = it;
			--it1;
			x = (it1->second.i1 != it->second.i1);
		}
		it2 = it;
		++it2;
		if (x && (it2 != vec[index].end()))
		{
			x = (it2->second.i1 != it->second.i1);
		}
		if (x)
			_rec_delete(vec, index-1, it->second.i1);
		vec[index].erase(position);
	}
}

template<typename TValue, typename TAlign, typename TSize>
void
_deleteAlignment(::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &me,
				TSize old_end,
				TSize new_end)
{
	//cout << "me: " << me[0].size() << endl;
	int length = me.size()-2;
	typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator it, it1, it2;
	for(TSize i = old_end+1; i < new_end; ++i)
	{
	
		it = me[length].find(i);
		bool x = true; //true = successor of it can be deleted
		if (it != me[length].begin())
		{
			it1 = it;
			--it1;
			x = (it1->second.i1 != it->second.i1);
		}
		it2 = it;
		++it2;
		if (x && (it2 != me[length].end()))
		{
			x = (it2->second.i1 != it->second.i1);
		}
		if (x)
			_rec_delete(me, length-1, it->second.i1);
		me[length].erase(i);
		//cout << "me2: " << me[0].size() << endl;
	}
	
}

//calculation and backtracking of the banded alignment of a seed
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateBandedSeed(TSeed const &seed,
					 TDiff k,
					 TMatrix &matrix_,
					 TString *p_seq1,
					 TString *p_seq2,
					 TScoreString &score_str,
					 TValue &score_length,
					 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					 TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
	Segment<TString, InfixSegment> seg1_align(*p_seq1, leftDim0(seed), rightDim0(seed)+1);
	Segment<TString, InfixSegment> seg2_align(*p_seq2, leftDim1(seed), rightDim1(seed)+1);
	_banded_needleman_wunsch(matrix_, seed, k, seg1_align, seg2_align, scoreMatrix, score_str);

	TValue height_diag = leftDiagonal(seed)-startDiagonal(seed)+k;
	TValue width_diag = startDiagonal(seed)-rightDiagonal(seed)+k;
	TValue overall = height_diag + width_diag + 1;
	
	resize(score_str,overall);
	TMatrixIterator matr_it = begin(matrix_);
	setPosition(matr_it, length(matrix_,0)-2);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	TValue width_align = leftDim0(seed) + width_diag;
	TValue height_align = leftDim1(seed);

	//diagonal back_track
	TValue new_connect;
	TValue old_connect = -1;
	for(TValue j = 0; j<=width_diag; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, rightDim0(seed)+1);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align,rightDim1(seed)+1);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		new_connect = _banded_needleman_wunsch_trace2(mapIt->second.i2, matrix_,  matr_it, scoreMatrix);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		--width_align;
		if (old_connect != new_connect)
		_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
			
	++width_align;
	++height_align;

	for(TValue j = width_diag+1; j < overall; ++j)
	{
		goNext(matr_it,1);
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, rightDim0(seed)+1);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align,rightDim1(seed)+1);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		
		new_connect = _banded_needleman_wunsch_trace2(mapIt->second.i2, matrix_,  matr_it, scoreMatrix);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		++height_align;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str);
}

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateFirstRectangle(TSeed const &seed,
						 TDiff k,
						 TMatrix &matrix_,
						 TString *p_seq1,
						 TString *p_seq2,
						 TScoreString &score_str,
						 TValue &score_length,
						 ::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						 TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
	TValue new_connect;
	Segment<TString, InfixSegment> seg1b_align(*p_seq1, 0, leftDim0(seed) + startDiagonal(seed) - rightDiagonal(seed) + k);
	Segment<TString, InfixSegment> seg2b_align(*p_seq2, 0, leftDim1(seed) + leftDiagonal(seed)  - startDiagonal(seed) + k);

	_banded_needleman_wunsch_rectangle_first(matrix_, seed, k, seg1b_align, seg2b_align, scoreMatrix,score_str);

	TValue w_d2 = startDiagonal(seed) - rightDiagonal(seed) + k;
	TValue h_d2 = leftDiagonal(seed) - startDiagonal(seed) + k;

	resize(score_str,1);

	TMatrixIterator matr_it = begin(matrix_);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());

	TValue width_stop = leftDim0(seed);
	TValue height_stop = leftDim1(seed);

	alignmentVector.back().insert(std::make_pair(0,Pair<TValue, TAlign> ()));
	TMapIterator mapIt = --(alignmentVector.back().end());
	Segment<TString, InfixSegment> seg1_align(*p_seq1, 0, leftDim0(seed) + w_d2);
	Segment<TString, InfixSegment> seg2_align(*p_seq2, 0, leftDim1(seed)+ h_d2);
	
	resize(rows(mapIt->second.i2),2);
	assignSource(row(mapIt->second.i2,0), seg1_align);
	assignSource(row(mapIt->second.i2,1), seg2_align);

	new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoreMatrix, matrix_, width_stop, height_stop);
	score_str[0] = *matr_it;
	mapIt->second.i1 = new_connect;
	_deleteAlignment(alignmentVector, -1, new_connect);
	_deleteAlignment(alignmentVector, new_connect, score_length);
}

template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateLastRectangle(TSeed const &seed,
						TDiff k,
						TMatrix &matrix_,
						TString *p_seq1,
						TString *p_seq2,
						TScoreString &score_str,
						TValue &score_length,
						::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
						TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;

	TValue seq1_length = length(*p_seq1);
	TValue seq2_length = length(*p_seq2);
	TValue width_diag=leftDiagonal(seed) - endDiagonal(seed)+k;
	TValue width =  width_diag + seq1_length - rightDim0(seed);
	TValue height_diag = endDiagonal(seed) - rightDiagonal(seed)+k;
	TValue height = height_diag + seq2_length - rightDim1(seed);
	
	resize(score_str, width_diag+height_diag+1);
	
	Segment<TString, InfixSegment> seg1(*p_seq1,seq1_length-width+1,seq1_length);
	Segment<TString, InfixSegment> seg2(*p_seq2,seq2_length-height+1,seq2_length);

	_needleman_wunsch(matrix_, seg1, seg2, scoreMatrix);

	TValue width_align = seq1_length - width  + width_diag +1;
	TValue height_align = seq2_length - height+1;

	TMatrixIterator iter_ = begin(matrix_);

	TValue x = width_diag;
	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	
	//last rectangle
	for(TValue i = 0; i<height_diag; ++i)
	{
		alignmentVector[0].insert(std::make_pair(i,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --alignmentVector[0].end();
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, seq1_length);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align, seq2_length);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0),seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		setPosition(iter_, x);
		score_str[i] = *iter_;
		mapIt->second.i1 = -1;
		_needleman_wunsch_trace_lastRectangle(mapIt->second.i2, iter_, scoreMatrix);
		x+=width;
		++height_align;
	}
	TValue a = height_diag + width_diag+1;
	
	for(TValue i = height_diag; i<a; ++i)
	{
		alignmentVector[0].insert(std::make_pair(i,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --alignmentVector[0].end();
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, seq1_length);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align, seq2_length);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0),seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		setPosition(iter_, x);
		score_str[i] = *iter_;
		mapIt->second.i1 = -1;
		_needleman_wunsch_trace_lastRectangle(mapIt->second.i2, iter_, scoreMatrix);
		--x;
		--width_align;
	}
	score_length = length(score_str);
}

//calculation and backtracking of the edit matrix between two seeds
template<typename TSeed, typename TString, typename TDiff, typename TMatrix, typename TScoreString, typename TValue, typename TAlign, typename TScoreMatrix>
void
_calculateRectangle(TSeed const &seed,
					TSeed const &seed2,
					TDiff k_begin,
					TDiff k_end,
					TMatrix &matrix_,
					TString *p_seq1,
					TString *p_seq2,
					TScoreString &score_str,
					TValue &score_length,
					::std::vector< ::std::map<TValue,Pair<TValue, TAlign> > > &alignmentVector,
					TScoreMatrix const &scoreMatrix)
{
	typedef typename ::std::map<TValue,Pair<TValue, TAlign> >::iterator TMapIterator;
	typedef Iter<TMatrix, PositionIterator> TMatrixIterator;
	Segment<TString, InfixSegment> seg1b_align(*p_seq1, rightDim0(seed2)-(leftDiagonal(seed2) - endDiagonal(seed2)   + k_end) + 1, leftDim0(seed) + startDiagonal(seed) - rightDiagonal(seed) + k_begin);
	Segment<TString, InfixSegment> seg2b_align(*p_seq2, rightDim1(seed2)-(endDiagonal(seed2) -  rightDiagonal(seed2) + k_end) + 1, leftDim1(seed) + leftDiagonal(seed)  - startDiagonal(seed) + k_begin);

	_needleman_wunsch_rectangle(matrix_, seed, seed2, k_begin, k_end, seg1b_align, seg2b_align, scoreMatrix,score_str);

	TValue width_diag = leftDiagonal(seed2) - endDiagonal(seed2)+k_end;
	TValue height_diag = endDiagonal(seed2) - rightDiagonal(seed2)+k_end;
	
	TValue w_d2 = startDiagonal(seed) - rightDiagonal(seed) + k_begin;
	TValue h_d2 = leftDiagonal(seed) - startDiagonal(seed) + k_begin;

	TValue overall = height_diag + width_diag + 1;

	resize(score_str, overall);
	TMatrixIterator matr_it = begin(matrix_);
	setPosition(matr_it, width_diag);

	alignmentVector.push_back(::std::map<TValue,Pair<TValue, TAlign> >());
	TValue width_align = rightDim0(seed2)+1;
	TValue height_align = rightDim1(seed2) - height_diag + 1;

	TValue width_stop = leftDim0(seed) - rightDim0(seed2)-1;
	TValue height_stop = leftDim1(seed) - (rightDim1(seed2) - height_diag)-1;

	TValue old_connect = -1;
	TValue new_connect;
	for(TValue j = 0; j<height_diag; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, leftDim0(seed) + w_d2);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align, leftDim1(seed)+ h_d2);

		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);

		//cout << "ALL: " << length(seg1_align) << " " << length(seg2_align) << " " << width_stop << " " <<height_stop << endl;
		new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoreMatrix, matrix_, width_stop, height_stop);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goNext(matr_it, 1);
		++height_align;
		--height_stop;
		if (old_connect != new_connect){
			//cout << "ups "<< old_connect << " " << new_connect << endl;
			_deleteAlignment(alignmentVector, old_connect, new_connect);
			//cout << "ups_ENDE" << endl;
		}
		old_connect = new_connect;
	}
	for(TValue j = height_diag; j < overall; ++j)
	{
		alignmentVector.back().insert(std::make_pair(j,Pair<TValue, TAlign> ()));
		TMapIterator mapIt = --(alignmentVector.back().end());
		Segment<TString, InfixSegment> seg1_align(*p_seq1, width_align, leftDim0(seed) + w_d2);
		Segment<TString, InfixSegment> seg2_align(*p_seq2, height_align, leftDim1(seed)+ h_d2);
		resize(rows(mapIt->second.i2),2);
		assignSource(row(mapIt->second.i2,0), seg1_align);
		assignSource(row(mapIt->second.i2,1), seg2_align);
		new_connect = _needleman_wunsch_trace_rectangle(mapIt->second.i2,  matr_it, scoreMatrix, matrix_, width_stop, height_stop);
		score_str[j] = *matr_it;
		mapIt->second.i1 = new_connect;
		goPrevious(matr_it,0);
		--width_align;
		++width_stop;
		if (old_connect != new_connect)
			_deleteAlignment(alignmentVector, old_connect, new_connect);
		old_connect = new_connect;
	}
	_deleteAlignment(alignmentVector, old_connect, score_length);
	score_length = length(score_str);
}



//Rectangle calculation between two seeds
template <typename TScoreValue, typename TMatrixSpec, typename TString, typename TValue, typename TValue2, typename TSpecSeed>
void
_needleman_wunsch_rectangle(Matrix<TScoreValue, TMatrixSpec> & matrix_,			//edit matrix
							Seed<TValue, TSpecSeed> const &seed1,		//Seed nearer to the end
							Seed<TValue, TSpecSeed> const &seed2,		//Seed nearer to the start
							TValue2 k_begin,							//upper diagonal extension
							TValue2 k_end,								//lower diagonal extension
							TString const & str1_,						//first sequence
							TString const & str2_,						//secondSequence
							Score<TScoreValue, Simple> const & score_,	//score matrix
							String<TScoreValue> init)//Values for initialisation
{
SEQAN_CHECKPOINT
	typedef Matrix<TScoreValue, TMatrixSpec> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, PositionIterator>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TValue inf = -1000000;

	TValue diag_height1 = leftDiagonal(seed1) - startDiagonal(seed1) + k_begin;
	TValue diag_width1 = startDiagonal(seed1) - rightDiagonal(seed1) + k_begin;

	TValue diag_height2 = endDiagonal(seed2) - rightDiagonal(seed2) + k_end;
	TValue diag_width2 = leftDiagonal(seed2) - endDiagonal(seed2) + k_end;


	TValue rectangle_height = leftDim1(seed1) - rightDim1(seed2) -1;
	TValue rectangle_width = length(str1_);

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TValue str1_length = length(str1_);
	TValue str2_length = length(str2_);
	setLength(matrix_, 0, str1_length + 1);
	setLength(matrix_, 1, str2_length + 1);
	
	resize(matrix_);

	TMatrixIterator col_ = begin(matrix_);
	//Initialisierung

	setPosition(col_, (diag_height2 + rectangle_height+1)*(str1_length+1)-1);

	for (int i = 0; i != diag_width1; ++i)
	{
		
		*col_ = init[i];
		--col_;
	}

	TValue len = length(init);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_=init[i];
		goNext(col_,1);
	}

	goPrevious(col_,1);
	TMatrixIterator finger2 =col_;
	
	--col_;



	TValue width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_ = inf;
		--col_;
	}


	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	setPosition(x_end,width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;

	col_ = finger2;
	TMatrixIterator finger1;
	
	TScoreValue h, v;
	TScoreValue s1;
	for (int i = 0; i < diag_height1; ++i)
	{
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;
		v = *finger1;
		h = *finger2;
		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
		}
		--y_end;
	}

	setPosition(col_,(diag_height2 + rectangle_height+1)*(rectangle_width+1)-1);

	h =*col_;
	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			//cout << position(finger1,0) << endl;
			*finger1 = v;
		}
		h = inf;
		--y_end;
	}

	setPosition(x_begin,diag_width2-1);

	for (int i = 0; i < diag_height2; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
			
		}
		h = inf;
		--y_end;
	}
}

//changed version of usual trace-back
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, typename TSourceSpec, typename TValue, typename TMatrix>
TValue
_needleman_wunsch_trace_rectangle(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, TSourceSpec>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_,
						TMatrix matrix, 
						TValue width_stop, 
						TValue height_stop)
{
SEQAN_CHECKPOINT
	typedef Iter<Matrix<TScoreValue, TSourceSpec>, PositionIterator > TMatrixIterator;
	typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);

	typedef typename Iterator<TTargetSourceSegment, Rooted>::Type TStringIterator;
	TStringIterator it_0 = iter(str_0, 0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);
	//TScoreValue score_match = scoreMatch(score_);

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached

	while ((static_cast<TValue>(position(it_0)) < width_stop) || (static_cast<TValue>(position(it_1)) < height_stop))
	{
		
		bool gv;
		bool gh;
		TMatrixIterator it3_ = source_;
		goNext(it3_,1);
		goNext(it3_,0);
		if ((*it_0 == *it_1)&&(*it3_ != -1000000))
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

	setSourceEndPosition(row(target_,1),position(it_1));
	setSourceEndPosition(row(target_,0),position(it_0));

	return length(matrix,0) - coordinate(source_,0) + position(it_1) - height_stop -1;
}


//Rectangle calculation between two seeds
template <typename TScoreValue, typename TMatrixSpec, typename TString, typename TValue, typename TValue2, typename TSpecSeed>
void
_banded_needleman_wunsch_rectangle_first(Matrix<TScoreValue, TMatrixSpec> & matrix_,	//edit matrix
								   Seed<TValue, TSpecSeed> const &seed,		//Seed
								   TValue2 k,									//diagonal extension
								   TString const & str1_,						//first sequence
								   TString const & str2_,						//secondSequence
								   Score<TScoreValue, Simple> const & score_,	//score matrix
								   String<TScoreValue> init)					//Values for initialisation
{
SEQAN_CHECKPOINT

	typedef Matrix<TScoreValue, TMatrixSpec> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TValue inf = -1000000;

	TValue diag_height1 = leftDiagonal(seed) - startDiagonal(seed) + k;
	TValue diag_width1 = startDiagonal(seed) - rightDiagonal(seed) + k;

	TValue rectangle_height = leftDim1(seed);
	TValue rectangle_width = length(str1_);
	


	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TValue str1_length = length(str1_);
	TValue str2_length = length(str2_);
	setLength(matrix_, 0, str1_length + 1);
	setLength(matrix_, 1, str2_length + 1);
	
	resize(matrix_);
	TMatrixIterator col_ = begin(matrix_);

	//Initialisierung
	setPosition(col_, (rectangle_height+1)*(rectangle_width+1)-1);

	for (int i = 0; i != diag_width1; ++i)
	{
		*col_ = init[i];
		--col_;
	}

	TValue len = length(init);
	for (int i = diag_width1; i != len; ++i)
	{
		*col_=init[i];
		goNext(col_,1);
	}

	goPrevious(col_,1);
	TMatrixIterator finger2 =col_;
	--col_;

	TValue width_align = str1_length - diag_width1;
	for (int i = 0; i!= width_align; ++i)
	{
		*col_ = inf;
		--col_;
	}

	//calculation of matrix

	TStringIterator x_begin = begin(str1_);
	--x_begin;
	TStringIterator x_end = end(str1_) -1;
	
	setPosition(x_end, width_align-1);
	TStringIterator y_begin = begin(str2_);
	TStringIterator y_end = end(str2_)-1;
	TStringIterator x,y;
	col_ = finger2;
	TMatrixIterator finger1;
	
	TScoreValue h, v;
	TScoreValue s1;
	for (int i = 0; i < diag_height1; ++i)
	{
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;
		v = *finger1;
		h = *finger2;
		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
		}
		--y_end;
	}

	setPosition(col_, (rectangle_height+1)*(rectangle_width+1)-1);


	h =*col_;
	x_end = end(str1_) -1;
	for (int i = 0; i < rectangle_height; ++i)
	{
		v = inf;
		finger2 = col_;		//points to last column
		goPrevious(col_, 1);	//points to this column
		finger1 = col_;

		*finger1 = v;

		for (x = x_end; x != x_begin; --x)
		{
			goPrevious(finger1, 0);
			goPrevious(finger2, 0);
			s1 = h + ((*x == *y_end) ? score_match : score_mismatch);
			h = *finger2;
			TScoreValue s2 = score_gap + ((h > v) ? h : v);
			v = (s1 > s2) ? s1 : s2;
			*finger1 = v;
			
		}
		h = inf;
		--y_end;
	}	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
