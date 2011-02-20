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
  $Id: align_hirschberg.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_ALIGN_HIRSCHBERG_H
#define SEQAN_HEADER_ALIGN_HIRSCHBERG_H

#include <stack>
#include <seqan/align/hirschberg_set.h>

//#define SEQAN_HIRSCHBERG_DEBUG_CUT

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
	template<typename TSource>
	void write_debug_matrix(TSource s1,TSource s2)
	{
		int l1 = length(s1);
		int l2 = length(s2);
	    
		int i,j,sg,sd;

		String<String<int> > fMatrix,rMatrix,tMatrix;

		resize(fMatrix,l1 + 1);
		resize(rMatrix,l1 + 1);
		resize(tMatrix,l1 + 1);

		for(i = 0;i <= l1;++i)
		{
			resize(fMatrix[i],l2 + 1);
			resize(rMatrix[i],l2 + 1);
			resize(tMatrix[i],l2 + 1);
		}

		for(i = 0;i <= l1;++i)
			fMatrix[i][0] = i * (-1);

		for(i = l1;i >= 0;--i)
			rMatrix[i][l2] = (l1 - i) * (-1);

		// calculate forward matrix
		for(j = 1;j <= l2;++j)
		{
			fMatrix[0][j] = j*(-1);
			for(i = 1;i <= l1;++i)
			{
				sg = -1 + ((fMatrix[i-1][j] > fMatrix[i][j-1]) ? fMatrix[i-1][j] : fMatrix[i][j-1]);
				sd = fMatrix[i-1][j-1] + ((s1[i - 1] == s2[j-1]) ? 0 : -1 );
		
				fMatrix[i][j] = ((sg > sd) ? sg : sd);
			}
		}

		// calculate reverse matrix
		for(j = l2 - 1;j >= 0;--j)
		{	
			rMatrix[l1][j] = (l2 - j)*(-1);
			for(i = l1 - 1;i >= 0;--i)
			{
				sg = -1 + ((rMatrix[i+1][j] > rMatrix[i][j+1]) ? rMatrix[i+1][j] : rMatrix[i][j+1]);
				sd = rMatrix[i+1][j+1] + ((s1[i] == s2[j]) ? 0 : -1 );
		
				rMatrix[i][j] = ((sg > sd) ? sg : sd);
			}
		}

		// print fMatrix
		std::cout << ";-;";
		for(i = 0;i < l1;++i)
			std::cout << s1[i] << ";";

		std::cout << std::endl << "-;";
		for(j = 0;j <= l2;++j)
		{	
			if(j != 0) std::cout << s2[j-1] << ";";
			for(i = 0;i <= l1;++i)
			{
				std::cout << fMatrix[i][j] << ";";
			}
			std::cout << std::endl;
		}
		// print rMatrix
		std::cout << ";";
		for(i = 0;i < l1;++i)
			std::cout << s1[i] << ";";
		std::cout << "-;" << std::endl;

		for(j = 0;j <= l2;++j)
		{	
			if(j != l2) std::cout << s2[j] << ";";
			else std::cout << "-;";
			for(i = 0;i <= l1;++i)
			{
				std::cout << rMatrix[i][j] << ";";
			}
			std::cout << std::endl;
		}

		// fill and print target matrix
		std::cout << ";-;";
		for(i = 0;i < l1;++i)
			std::cout << s1[i] << ";";

		std::cout << std::endl << "-;";
		for(j = 0;j <= l2;++j)
		{	
			if(j != 0) std::cout << s2[j-1] << ";";
			for(i = 0;i <= l1;++i)
			{
				tMatrix[i][j] = fMatrix[i][j] + rMatrix[i][j];
				std::cout << tMatrix[i][j] << ";";
			}
			std::cout << std::endl;
		}
	}

#endif

// debug flag .. define to see where Hirschberg cuts the sequences	
//#define SEQAN_HIRSCHBERG_DEBUG_CUT

/*DISABLED
.Function.hirschberg:
..cat:Alignment
..summary:Computes a global Alignment for the passed Alignment-Container with the specified scoring scheme
..signature:hirschberg(Align<TSource, TSpec> & align,Score<TScoreValue, Simple> const & score)
..param.align: Reference to the Alignment-Object
..param.score: Const Reference to the Scoring Scheme
..remarks: The alignment is based on the algorithm proposed by Hirschberg. The general idea is to divide the DP (dynamic programming) matrix,
to compute a global alignment in linear space. Instead of computing half of the
DP matrix in forward direction and the other half in reverse, a pointer to the cell of the DP matrix, were the actual, optimal alignment
passes the mid column ist saved, during the computation of the second part of the Matrix.
*/

template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
                Score<TScoreValue, Simple> const & score_,
                Hirschberg)
{
SEQAN_CHECKPOINT
	TSource s1 = sourceSegment(row(align_, 0));
	TSource s2 = sourceSegment(row(align_, 1));
	
	TScoreValue total_score;

	typedef typename Value<TSource>::Type TStringValue;
	typedef typename Size<TSource>::Type TStringSize;
	
	typedef typename Iterator<TSource>::Type TStringIterator;

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TTargetIterator;

	TTargetIterator target_0 = iter(row(align_, 0), 0);
	TTargetIterator target_1 = iter(row(align_, 1), 0);

	typedef typename Size<Matrix<TScoreValue> >::Type TSize;
	typedef typename Iterator<Matrix<TScoreValue> >::Type TMatrixIterator;

	TStringValue v;

	TStringSize len1 = length(s1);
	TStringSize len2 = length(s2);

	// string to store the score values for the currently active cell
	String<TScoreValue> c_score;
	resize(c_score,len2 + 1);
	// string to strore the backpointers
	String<int> pointer;
	resize(pointer,len2 + 1);

	// scoring-scheme specific score values
	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TScoreValue border,s,sg,sd,sg1,sg2;
	int dp;

	std::stack<_HirschbergSet> to_process;
	_HirschbergSet target;
	
	int i,j;
	
	_HirschbergSet hs_complete(0,len1,0,len2,0);
	to_process.push(hs_complete);

	while(!to_process.empty())
	{
SEQAN_CHECKPOINT
		target = to_process.top();
		to_process.pop();

		if(_begin2(target) == _end2(target))
		{
			for(i = 0;i < (_end1(target) - _begin1(target));++i)
			{
				insertGap(target_1);
				++target_0;
				++target_1;
			}
		}
		else if(_begin1(target) + 1 == _end1(target) || _begin2(target) + 1 == _end2(target))
		{
			/* ALIGN */			
SEQAN_CHECKPOINT
#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
			std::cout << "align s1 " << _begin1(target) << " to " << _end1(target) << " and s2 " << _begin2(target) << " to " << _end2(target) << std::endl;
			std::cout << "align " << infix(s1,_begin1(target),_end1(target)) << " and " << infix(s2,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif

			TStringSize len_1 = _end1(target) - _begin1(target);
			TStringSize len_2 = _end2(target) - _begin2(target);

			Matrix<TScoreValue> matrix_;

			setDimension(matrix_, 2);
			setLength(matrix_, 0, len_1 + 1);
			setLength(matrix_, 1, len_2 + 1);
			resize(matrix_);

			/* init matrix */
			TStringIterator x_begin = iter(s1,_begin1(target)) - 1;
			TStringIterator x_end = iter(s1,_end1(target)) - 1;
			TStringIterator y_begin = iter(s2,_begin2(target)) - 1;
			TStringIterator y_end = iter(s2,_end2(target)) - 1;

			TStringIterator x = x_end;
			TStringIterator y;

			TMatrixIterator col_ = end(matrix_) - 1;
			TMatrixIterator finger1;
			TMatrixIterator finger2;


			TScoreValue h = 0;
			TScoreValue border_ = score_gap;
			TScoreValue v = border_;


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
				TStringValue cy = *y;
				h = border_;
				border_ += score_gap;
				v = border_;

				finger2 = col_;	
				goPrevious(col_, 1);
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

			/* TRACE BACK */
			finger1 = begin(matrix_);
			x = iter(s1,_begin1(target));
			y = iter(s2,_begin2(target));
			x_end = iter(s1,_end1(target));
			y_end = iter(s2,_end2(target));

			while ((x != x_end) && (y != y_end))
			{
				bool gv;
				bool gh;

				if (*x == *y)
				{
					gv = gh = true;
				}
				else
				{
					TMatrixIterator it_ = finger1;

					goNext(it_, 0);
					TScoreValue v = *it_;

					goNext(it_, 1);
					TScoreValue d = *it_;

					it_ = finger1;
					goNext(it_, 1);
					TScoreValue h = *it_;

					gv = (v >= h) | (d >= h);
					gh = (h >= v) | (d >= v);
				}

				if (gv)
				{
					++x;
					goNext(finger1, 0);
				}
				else
				{
					insertGap(target_0);
				}

				if (gh) 
				{
					++y;
					goNext(finger1, 1);
				}
				else
				{
					insertGap(target_1);
				}

				++target_0;
				++target_1;
			}

			// if x or y did not reached there end position, fill the rest with gaps
			while(x != x_end)
			{
				insertGap(target_1);
				++target_0;
				++target_1;
				++x;
			}

			while(y != y_end)
			{
				insertGap(target_0);
				++target_0;
				++target_1;
				++y;
			}
			/* END ALIGN */
		}
		else
		{	
			/* 
				Calculate cut using the algorithm as proposed in the lecture of Clemens Gröpl 
				using a backpointer to remember the position where the optimal alignment passes
				the mid column
			*/
SEQAN_CHECKPOINT
			int mid = static_cast<int>(floor( static_cast<double>((_begin1(target) + _end1(target))/2) ));

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
			std::cout << "calculate cut for s1 " << _begin1(target) << " to " << _end1(target) << " and s2 " << _begin2(target) << " to " << _end2(target) << std::endl;
			std::cout << "calculate cut for " << infix(s1,_begin1(target),_end1(target)) << " and " << infix(s2,_begin2(target),_end2(target)) << std::endl;
			std::cout << "cut is in row " << mid << " symbol is " << getValue(s1,mid-1) << std::endl << std::endl;


			write_debug_matrix(infix(s1,_begin1(target),_end1(target)),infix(s2,_begin2(target),_end2(target)));
#endif

			border = 0;
			for(i = _begin2(target);i <= _end2(target);++i)
			{
				c_score[i] = border;
				border += score_gap; 
				pointer[i] = i;
			}

			// iterate over s1 until the mid column is reached
			border = score_gap;
			for(i = _begin1(target) + 1;i <= mid;++i)
			{
				s = c_score[_begin2(target)];
				c_score[_begin2(target)] = border;
				border += score_gap;
				v = getValue(s1,i-1);
				for(j = _begin2(target) + 1;j <= _end2(target);++j)
				{
					sg = score_gap + ((c_score[j] > c_score[j - 1]) ? c_score[j] : c_score[j - 1]);
					sd = s + ((v == getValue(s2,j-1)) ? score_match : score_mismatch);

					s = c_score[j];
					c_score[j] = (sg > sd) ? sg : sd;	
				}
			}

			// from here, rememeber the cell of mid-column, where optimal alignment passed
			for(i = mid + 1;i <= _end1(target);++i)
			{
				s = c_score[_begin2(target)];
				c_score[_begin2(target)] = border;
				border += score_gap;
				v = getValue(s1,i-1);
			
				dp = _begin2(target);
			
				for(j = _begin2(target) + 1;j <= _end2(target);++j)
				{
					sg1 = score_gap + c_score[j];
					sg2 = score_gap + c_score[j - 1];

					sd = s + ((v == getValue(s2,j-1)) ? score_match : score_mismatch);

					s = c_score[j];
					sg = pointer[j];
					if(sd >= max(sg1,sg2))
					{
						c_score[j] = sd;
						pointer[j] = dp;
					}
					else
					{
						if(sg2 > sg1)
						{
							c_score[j] = sg2;
							pointer[j] = pointer[j-1];
						}
						else
						{
							// gap introduced from left
							// no update for the pointer
							c_score[j] = sg1;
						}
					}
					dp = sg;
				}
			}

			// if computed the whole matrix max = alignment score
			if(hs_complete == target)
				total_score = c_score[_end2(target)];

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
			std::cout << "hirschberg calculates cut in column " << mid << " and row " << pointer[_end2(target)] << std::endl;
			std::cout << "requested position in c_score and pointer is " << _end2(target) << std::endl;
			std::cout << "alignment score is " << c_score[_end2(target)] << std::endl << std::endl;
#endif
			to_process.push(_HirschbergSet(mid,_end1(target),pointer[_end2(target)],_end2(target),0));
			to_process.push(_HirschbergSet(_begin1(target),mid,_begin2(target),pointer[_end2(target)],0));
		}
		/* END CUT */
	}
	return total_score;
}


}
#endif
