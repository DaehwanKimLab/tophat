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
  $Id: align_myers.h 1770 2008-03-12 13:10:57Z aiche@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_MYERS_H
#define SEQAN_HEADER_ALIGN_MYERS_H

// is necessary to include extended alphabet 
// can possibly be omitted by extending the alphabet defenition
#ifdef BAC_ALIGNER 
#include <../apps/bac_aligner/extended_iupac_alphabet.h>
#endif

namespace SEQAN_NAMESPACE_MAIN
{
//#define MYERS_HIRSCHBERG_VERBOSE

#ifdef MYERS_HIRSCHBERG_VERBOSE
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
/*DISABLED
.Function.align_myers:
..cat:Alignment
..summary:(Name of the Function is going to be changed in future)Computes the global Alignmentscore for the passed Alignment-Container with an Edit-Distance Scoring Scheme
..signature:hirschberg(Align<TSource, TSpec> & align,Score<TScoreValue, Simple> const & score)
..param.align:The alignment object having the sequences to be aligned as sources.
...type:Class.Align
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..remarks: The Computation of the Alignment-Score is based on Myers-Bitvektor-Algorihm for Approximate Stringmatching. The Algorithm was customized,
as proposed by Hyrroe to compute edit distance. To compute a complete Alignment use @Function.hirschberg_myers@ or @Function.hirschberg@ for smaller Instances of the Problem.
*/
template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const &,
				MyersBitVector)
{
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	// use size of unsigned int as blocksize for bit-vectors
	const unsigned int BLOCK_SIZE = BitsPerValue<unsigned int>::VALUE;

	typedef typename Value<TSource>::Type TAlphabet;
	typedef typename Size<TSource>::Type TSourceSize;

	// switch x and y .. y should be the shorter one, to allocate less memory for the bitMasks
	TSource x,y;
	if(length(sourceSegment(row(align_, 0))) < length(sourceSegment(row(align_, 1))))
	{
		y = sourceSegment(row(align_, 0));
		x = sourceSegment(row(align_, 1));
	}
	else
	{
		x = sourceSegment(row(align_, 0));
		y = sourceSegment(row(align_, 1));
	}
	
	TSourceSize len_x = length(x);
	unsigned int pos = 0;

	// init variables 
	unsigned int len_y = length(y);
	unsigned int score = (-1)*len_y;
	unsigned int alphabetSize = ValueSize<TAlphabet>::VALUE;
	unsigned int blockCount = (len_y + BLOCK_SIZE - 1) / BLOCK_SIZE;

	unsigned int scoreMask = 1 << ((len_y % BLOCK_SIZE) - 1);	// the mask with a bit set at the position of the last active cell

	unsigned int * VP;
	unsigned int * VN;
	unsigned int * bitMask;

	allocate (align_, VP, blockCount);
	arrayFill (VP, VP + blockCount, ~0);

	allocate (align_, VN, blockCount);
	arrayFill (VN, VN + blockCount, 0);

	// first bitMask will be constructed from the shorter sequence
	allocate (align_, bitMask, alphabetSize * blockCount);
	arrayFill(bitMask, bitMask + alphabetSize * blockCount, 0);

	// encoding the letters as bit-vectors
    for (unsigned int j = 0; j < len_y; j++)
		bitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] = bitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);

#ifdef BAC_ALIGNER // see definition of BAC_ALIGNER
	//extend the bitMasks for ambigous alphabets
	//	possibly intergrate Tag for Alphabet-Class
	if(_ClassIdentifier<TAlphabet>::getID() == _ClassIdentifier<EIupac>::getID())
	{
		unsigned int i,j,m;
		unsigned int * copyMask;

		// copy of bitmask
		allocate (align_, copyMask, alphabetSize * blockCount);
		arrayCopy(bitMask,bitMask + alphabetSize * blockCount,copyMask);
		
		for(i = 0;i < alphabetSize;++i)
		{
			for(j = 0;j < alphabetSize;++j)
			{
				if(static_cast<TAlphabet>(i) == static_cast<TAlphabet>(j))
				{
					unsigned int char_ind_i = blockCount * ordValue(static_cast<TAlphabet>(i));
					unsigned int char_ind_j = blockCount * ordValue(static_cast<TAlphabet>(j));
					for(m = 0;m < blockCount;++m)
					{
						bitMask[char_ind_i] |= copyMask[char_ind_j];
						++char_ind_i;
						++char_ind_j;
					}
				}
			}

		}

		deallocate(align_, copyMask, alphabetSize * blockCount);
	}
#endif
	// compute score
	unsigned int X, D0, HN, HP;
	if(blockCount == 1)
	{
		while (pos < len_x) {
			X = bitMask[ordValue(getValue(x,pos))] | VN[0];

			D0 = ((VP[0] + (X & VP[0])) ^ VP[0]) | X;
			HN = VP[0] & D0;
			HP = VN[0] | ~(VP[0] | D0);

			// customized to compute edit distance
			X = (HP << 1) | 1;
			VN[0] = X & D0;
			VP[0] = (HN << 1) | ~(X | D0);

			if (HP & scoreMask)
				score--;
			else if (HN & scoreMask)
				score++;

			++pos;
		}
	} // end compute score - short pattern
	else
	{
		unsigned int temp, shift, currentBlock;
		unsigned int carryD0, carryHP, carryHN;

		while (pos < len_x) 
		{
			// set vars
			carryD0 = carryHP = carryHN = 0;
			shift = blockCount * ordValue(getValue(x,pos));

			// computing first the top most block
			X = bitMask[shift] | VN[0];
	
			temp = VP[0] + (X & VP[0]);
			carryD0 = temp < VP[0];
			
			D0 = (temp ^ VP[0]) | X;
			HN = VP[0] & D0;
			HP = VN[0] | ~(VP[0] | D0);
			
			// customized to compute edit distance
			X = (HP << 1) | 1;
			carryHP = HP >> (BLOCK_SIZE - 1);
			
			VN[0] = X & D0;

			temp = (HN << 1);
			carryHN = HN >> (BLOCK_SIZE - 1);
								
		 	VP[0] = temp | ~(X | D0);

			// computing the necessary blocks, carries between blocks following one another are stored
			for (currentBlock = 1; currentBlock < blockCount; currentBlock++) {
				X = bitMask[shift + currentBlock] | VN[currentBlock];
		
				temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;
				
				carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);
			
				D0 = (temp ^ VP[currentBlock]) | X;
				HN = VP[currentBlock] & D0;
				HP = VN[currentBlock] | ~(VP[currentBlock] | D0);
				
				X = (HP << 1) | carryHP;
				carryHP = HP >> (BLOCK_SIZE-1);
				
				VN[currentBlock] = X & D0;

				temp = (HN << 1) | carryHN;
				carryHN = HN >> (BLOCK_SIZE - 1);
									
		 		VP[currentBlock] = temp | ~(X | D0);
			}

			// update score with the HP and HN values of the last block the last block
			if (HP & scoreMask)
				score--;
			else if (HN & scoreMask)
				score++;
			++pos;
		}

	} // end compute score - long pattern

	// clean up
	deallocate(align_, bitMask, alphabetSize * blockCount);
	deallocate(align_, VP, blockCount);
	deallocate(align_, VN, blockCount);

	return score;

} // end align_myers_score

/*DISABLED
.Function.hirschberg_myers:
..cat:Alignment
..summary:Computes the global Alignment for the passed Alignment-Container with an Edit-Distance Scoring Scheme.
..signature:hirschberg_myers(Align<TSource, TSpec> & align,Score<TScoreValue, Simple> const & score)
..param.align:The alignment object having the sequences to be aligned as sources.
...type:Class.Align
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..remarks: The Computation of the Alignment-Score is based on a combination of Myers-Bitvektor-Algorihm for Approximate Stringmatching and
the Algorithm proposed by Hirschberg to compute Sequence Alignments with linear space.
*/

template <typename TSource, typename TSpec, typename TScoreValue>
TScoreValue
globalAlignment(Align<TSource, TSpec> & align_,
				Score<TScoreValue, Simple> const &,
				MyersHirschberg)
{
SEQAN_CHECKPOINT
	
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	// use size of unsigned int as blocksize for bit-vectors
	const unsigned int BLOCK_SIZE = BitsPerValue<unsigned int>::VALUE;

	// saves the score value that will be returned
	TScoreValue score,total_score = 0;

	// switch x and y .. y should be the shorter one, to allocate less memory for the bitMasks
	TSource x,y;

	typedef typename Value<TSource>::Type TAlphabet;
	typedef typename Size<TSource>::Type TStringSize;
	
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

	typedef Align<TSource, TSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TTargetIterator;

	typedef typename Iterator<Matrix<TScoreValue>, Rooted>::Type TMatrixIterator;

	TTargetIterator target_0,target_1;

	if(length(sourceSegment(row(align_, 0))) < length(sourceSegment(row(align_, 1))))
	{
		x = sourceSegment(row(align_, 1));
		y = sourceSegment(row(align_, 0));
		
		target_0 = iter(row(align_, 1), 0);
		target_1 = iter(row(align_, 0), 0);
	}
	else
	{
		x = sourceSegment(row(align_, 0));
		y = sourceSegment(row(align_, 1));

		target_0 = iter(row(align_, 0), 0);
		target_1 = iter(row(align_, 1), 0);
	}


	TStringSize len_x = length(x);
	TStringSize len_y = length(y);

	// string to store the score values for the currently active cell
	String<TScoreValue> c_score;
	resize(c_score,len_x + 1);
	
	// scoring-scheme specific score values
	TScoreValue score_match = 0;
	TScoreValue score_mismatch = -1;
	TScoreValue score_gap = -1;

	// additional vars
	int i;

	// stack with parts of matrix that have to be processed
	std::stack<_HirschbergSet> to_process;
	_HirschbergSet target;

	// myers specific vars and preprocessing
	unsigned int alphabetSize = ValueSize<TAlphabet>::VALUE;
	unsigned int blockCount = (len_y + BLOCK_SIZE - 1) / BLOCK_SIZE; // maximal count of blocks

	unsigned int * VP;
	unsigned int * VN;
	
	unsigned int * forwardBitMask;			// encoding the alphabet as bit-masks
	unsigned int * reverseBitMask;
	
	allocate (align_, VP, blockCount);
	arrayFill (VP, VP + blockCount, ~0);

	allocate (align_, VN, blockCount);
	arrayFill (VN, VN + blockCount, 0);

	// first bitMask will be constructed from the shorter sequence
	allocate (align_, forwardBitMask, alphabetSize * blockCount);
	arrayFill(forwardBitMask, forwardBitMask + alphabetSize * blockCount, 0);

	allocate (align_, reverseBitMask, alphabetSize * blockCount);
	arrayFill(reverseBitMask, reverseBitMask + alphabetSize * blockCount, 0);

	// encoding the letters as bit-vectors
    for (unsigned int j = 0; j < len_y; j++){
SEQAN_CHECKPOINT
		forwardBitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] = forwardBitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);
		reverseBitMask[blockCount * ordValue(getValue(y,len_y - j - 1)) + j/BLOCK_SIZE] = reverseBitMask[blockCount * ordValue(getValue(y,len_y - j - 1)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);
	}

#ifdef BAC_ALIGNER // see definition of BAC_ALIGNER

	// extend the bitMasks for ambigous alphabets
	//		possibly intergrate Tag for Alphabet-Class
	if(_ClassIdentifier<TAlphabet>::getID() == _ClassIdentifier<EIupac>::getID())
	{
		unsigned int * fCopyMask;
		unsigned int * rCopyMask;
		
		// allocate memory for temporary copy of bitMasks
		allocate (align_, fCopyMask, alphabetSize * blockCount);
		allocate (align_, rCopyMask, alphabetSize * blockCount);

		arrayCopy(forwardBitMask, forwardBitMask + alphabetSize * blockCount, fCopyMask);
		arrayCopy(reverseBitMask, reverseBitMask + alphabetSize * blockCount, rCopyMask);

		unsigned int i,j,m;
		for(i=0;i < alphabetSize;++i)
		{
			// iterate over the whole alphabet
			// 1. all letters that were allready processed
			for(j = 0;j < i;++j)
			{
				if(static_cast<TAlphabet>(i) == static_cast<TAlphabet>(j))
				{
					unsigned int char_ind_i = blockCount*ordValue(static_cast<TAlphabet>(i));
					unsigned int char_ind_j = blockCount*ordValue(static_cast<TAlphabet>(j));

					for(m = 0;m < blockCount;++m)
					{
						forwardBitMask[char_ind_i] |= fCopyMask[char_ind_j];
						reverseBitMask[char_ind_i] |= rCopyMask[char_ind_j];
						++char_ind_i;
						++char_ind_j;
					}
				}
			}
			// 2. all unprocessed letters
			for(j = i+1;j < alphabetSize;++j)
			{
				if(static_cast<TAlphabet>(i) == static_cast<TAlphabet>(j))
				{
					unsigned int char_ind_i = blockCount*ordValue(static_cast<TAlphabet>(i));
					unsigned int char_ind_j = blockCount*ordValue(static_cast<TAlphabet>(j));
					
					for(m = 0;m < blockCount;++m)
					{
						forwardBitMask[char_ind_i] |= fCopyMask[char_ind_j];
						reverseBitMask[char_ind_i] |= rCopyMask[char_ind_j];
						++char_ind_i;
						++char_ind_j;
					}
				}
			}
		}

		/* deallocate space used for the temporary bitMasks */
		deallocate(align_, fCopyMask, alphabetSize * blockCount);
		deallocate(align_, rCopyMask, alphabetSize * blockCount);
	}
#endif
	_HirschbergSet hs_complete(0,len_x,0,len_y,1);
	to_process.push(hs_complete);

	while(!to_process.empty())
	{
SEQAN_CHECKPOINT
		target = to_process.top();
		to_process.pop();
		/* if score is zero, the whole part of the sequence can be simply skipped */
		if(_score(target) == 0)
		{
SEQAN_CHECKPOINT
			/* coukd work faster */
			for(i = 0;i < (_end1(target) - _begin1(target));++i)
			{
				++target_0;
				++target_1;
			}

#ifdef MYERS_HIRSCHBERG_VERBOSE
			printf("skipped %i to %i in first sequence\n",_begin1(target),_end1(target));
#endif
		}
		else if(_begin1(target) == _end1(target))
		{
SEQAN_CHECKPOINT

#ifdef MYERS_HIRSCHBERG_VERBOSE
			std::cout << "align y " << _begin2(target) << " to " << _end2(target) << std::endl;
			std::cout << "align " << infix(y,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif	
			for(i = 0;i < (_end2(target) - _begin2(target));++i)
			{
				insertGap(target_0);
				++target_0;
				++target_1;
			}
		}
		else if(_begin2(target) + 1 == _end2(target))
		{
			/* ALIGN */			
SEQAN_CHECKPOINT
#ifdef MYERS_HIRSCHBERG_VERBOSE
			std::cout << "align x " << _begin1(target) << " to " << _end1(target) << " and y " << _begin2(target) << " to " << _end2(target) << std::endl;
			std::cout << "align " << infix(x,_begin1(target),_end1(target)) << " and " << infix(y,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif

			TStringSize len_1 = _end1(target) - _begin1(target);
			TStringSize len_2 = _end2(target) - _begin2(target);

			Matrix<TScoreValue> matrix_;

			setDimension(matrix_, 2);
			setLength(matrix_, 0, len_1 + 1);
			setLength(matrix_, 1, len_2 + 1);
			resize(matrix_);

			/* init matrix */
			TStringIterator xs_begin = iter(x,_begin1(target)) - 1;
			TStringIterator xs_end = iter(x,_end1(target)) - 1;
			TStringIterator ys_begin = iter(y,_begin2(target)) - 1;
			TStringIterator ys_end = iter(y,_end2(target)) - 1;

			TStringIterator xs = xs_end;
			TStringIterator ys;

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
			for (xs = xs_end; xs != xs_begin; --xs)
			{
				goPrevious(finger1, 0);
				*finger1 = border_;
				border_ += score_gap;
			}

			//-------------------------------------------------------------------------
			//fill matrix

			border_ = 0;
			for (ys = ys_end; ys != ys_begin; --ys)
			{
				TAlphabet cy = *ys;
				h = border_;
				border_ += score_gap;
				v = border_;

				finger2 = col_;		
				goPrevious(col_, 1);	
				finger1 = col_;

				*finger1 = v;

				for (xs = xs_end; xs != xs_begin; --xs)
				{
					goPrevious(finger1, 0);
					goPrevious(finger2, 0);
					if (*xs == cy)
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

			// if computed the whole matrix last value of v = alignment score
			if(target == hs_complete)   total_score = v;

			/* TRACE BACK */
			finger1 = begin(matrix_);
			xs = iter(x,_begin1(target));
			ys = iter(y,_begin2(target));
			xs_end = iter(x,_end1(target));
			ys_end = iter(y,_end2(target));

			while ((xs != xs_end) && (ys != ys_end))
			{
				bool gv;
				bool gh;

				if (*xs == *ys)
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
					++xs;
					goNext(finger1, 0);
				}
				else
				{
					insertGap(target_0);
				}

				if (gh) 
				{
					++ys;
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
			while(xs != xs_end)
			{
				insertGap(target_1);
				++target_0;
				++target_1;
				++xs;
			}

			while(ys != ys_end)
			{
				insertGap(target_0);
				++target_0;
				++target_1;
				++ys;
			}
			/* END ALIGN */


#ifdef MYERS_HIRSCHBERG_VERBOSE
			std::cout << std::endl << align_ << std::endl << std::endl;
#endif

		}
		else
		{
SEQAN_CHECKPOINT
			/*
				---------------------------------------------------------------
				Calculate cut position using extended Myers-Bitvector-Algorithm
			    --------------------------------------------------------------- 
			*/

			/* declare variables */
			unsigned int X, D0, HN, HP;

			/* compute cut position */
			int mid = static_cast<int>(floor( static_cast<double>((_begin2(target) + _end2(target))/2) ));

			/* debug infos */
#ifdef MYERS_HIRSCHBERG_VERBOSE
			std::cout << "calculate cut for x " << _begin1(target) << " to " << _end1(target) << " and y " << _begin2(target) << " to " << _end2(target) << std::endl;
			std::cout << "calculate cut for " << infix(x,_begin1(target),_end1(target)) << " and " << infix(y,_begin2(target),_end2(target)) << std::endl;
			std::cout << "cut is in row " << mid << " symbol is " << getValue(x,mid-1) << std::endl << std::endl;

			std::cout << std::endl;
			write_debug_matrix(infix(x,_begin1(target),_end1(target)),infix(y,_begin2(target),_end2(target)));
			std::cout << std::endl;
#endif
			/* compute blocks and score masks */
			int fStartBlock = _begin2(target) / BLOCK_SIZE;
			int fEndBlock = (mid - 1) / BLOCK_SIZE;
			int fSpannedBlocks = (fEndBlock - fStartBlock) + 1;

			unsigned int fScoreMask = 1 << ((mid  - 1) % BLOCK_SIZE);
			
			unsigned int fOffSet = _begin2(target) % BLOCK_SIZE;
			unsigned int fSilencer = ~0;
			fSilencer <<= fOffSet;

			/* reset v-bitvectors */
			arrayFill (VP + fStartBlock, VP + fEndBlock + 1, ~0);
			arrayFill (VN + fStartBlock, VN + fEndBlock + 1, 0);


			/* determine start-position and start-score */
			int pos = _begin1(target);			
			score = (mid - _begin2(target)) * score_gap;
			c_score[pos] = score;

			/* compute with myers - forward - begin */
			if(fSpannedBlocks == 1)
			{
SEQAN_CHECKPOINT
				while (pos < _end1(target)) {
					X = (fSilencer & forwardBitMask[(blockCount * ordValue(getValue(x,pos))) + fStartBlock]) | VN[fStartBlock];

					D0 = ((VP[fStartBlock] + (X & VP[fStartBlock])) ^ VP[fStartBlock]) | X;
					HN = VP[fStartBlock] & D0;
					HP = VN[fStartBlock] | ~(VP[fStartBlock] | D0);

					X = (HP << 1) | (1 << fOffSet);
					VN[fStartBlock] = X & D0;
					VP[fStartBlock] = (HN << 1) | ~(X | D0);

					if (HP & fScoreMask)
						score--;
					else if (HN & fScoreMask)
						score++;

					c_score[pos + 1] = score;

					++pos;
				}
			} /* end - short patten */
			else
			{
SEQAN_CHECKPOINT

				int shift, currentBlock;
				unsigned int temp, carryD0, carryHP, carryHN;

				while (pos < _end1(target))
				{
					carryD0 = carryHP = carryHN = 0;
					shift = blockCount * ordValue(getValue(x,pos));

					// computing first the top most block
					X = (fSilencer & forwardBitMask[shift + fStartBlock]) | VN[fStartBlock];
			
					temp = VP[fStartBlock] + (X & VP[fStartBlock]);
					carryD0 = temp < VP[fStartBlock];
					
					D0 = (temp ^ VP[fStartBlock]) | X;
					HN = VP[fStartBlock] & D0;
					HP = VN[fStartBlock] | ~(VP[fStartBlock] | D0);
					
					X = (HP << 1) | (1 << fOffSet);
					carryHP = HP >> (BLOCK_SIZE - 1);
					
					VN[fStartBlock] = X & D0;

					temp = (HN << 1);
					carryHN = HN >> (BLOCK_SIZE - 1);
										
		 			VP[fStartBlock] = temp | ~(X | D0);

					// compute the remaining blocks
					for (currentBlock = fStartBlock + 1; currentBlock <= fEndBlock; currentBlock++) {
						X = forwardBitMask[shift + currentBlock] | VN[currentBlock];
				
						temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;
						
						carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);
					
						D0 = (temp ^ VP[currentBlock]) | X;
						HN = VP[currentBlock] & D0;
						HP = VN[currentBlock] | ~(VP[currentBlock] | D0);
						
						X = (HP << 1) | carryHP;
						carryHP = HP >> (BLOCK_SIZE-1);
						
						VN[currentBlock] = X & D0;

						temp = (HN << 1) | carryHN;
						carryHN = HN >> (BLOCK_SIZE - 1);
											
		 				VP[currentBlock] = temp | ~(X | D0);
					}
					
					/* update score */
					if (HP & fScoreMask)
						score--;
					else if (HN & fScoreMask)
						score++;

					c_score[pos + 1] = score;

					++pos;
				}

			} /* end - long patten */
			/* compute with myers - forward - end */
			
			/* compute blocks and score masks */
			int rStartBlock = (len_y - _end2(target)) / BLOCK_SIZE;
			int rEndBlock = (len_y - mid - 1) / BLOCK_SIZE;
			int rSpannedBlocks = (rEndBlock - rStartBlock) + 1;

			unsigned int rScoreMask = 1 <<  ((len_y - mid - 1) % BLOCK_SIZE);
			unsigned int rOffSet = (len_y - _end2(target)) % BLOCK_SIZE;
			unsigned int rSilencer = ~0;
			rSilencer <<= rOffSet;

			/* reset v-bitvectors */
			arrayFill (VP + rStartBlock, VP + rEndBlock + 1, ~0);
			arrayFill (VN + rStartBlock, VN + rEndBlock + 1, 0);

			/* determine start-position and start-score */
			pos = _end1(target)-1;			
			score = (_end2(target) - mid) * score_gap;

			/* set start score */
			c_score[_end1(target)] += score;

			/* determine optimal cut position -- score extension */
			TScoreValue max = c_score[_end1(target)];
			TScoreValue rmax = score;
			unsigned int pos_max = _end1(target);

			/* compute with myers - reverse - begin */
			if(rSpannedBlocks == 1)
			{
SEQAN_CHECKPOINT

				while (pos >= _begin1(target)) {
					X = (rSilencer & reverseBitMask[(blockCount * ordValue(getValue(x,pos))) + rStartBlock]) | VN[rStartBlock];

					D0 = ((VP[rStartBlock] + (X & VP[rStartBlock])) ^ VP[rStartBlock]) | X;
					HN = VP[rStartBlock] & D0;
					HP = VN[rStartBlock] | ~(VP[rStartBlock] | D0);

					X = (HP << 1) | (1 << rOffSet);
					VN[rStartBlock] = X & D0;
					VP[rStartBlock] = (HN << 1) | ~(X | D0);

					if (HP & rScoreMask)
						--score;
					else if (HN & rScoreMask)
						++score;

					c_score[pos] += score;

					/* check for optimality -- score extension */
					if(c_score[pos]> max)
					{
						pos_max = pos;
						max = c_score[pos];
						rmax =  score;
					}

					--pos;
				}
			} /* end - short pattern */
			else
			{
SEQAN_CHECKPOINT
				int shift, currentBlock;
				unsigned int temp, carryD0, carryHP, carryHN;

				while (pos >= _begin1(target))
				{
					carryD0 = carryHP = carryHN = 0;
					shift = blockCount * ordValue(getValue(x,pos));

					// compute first the top most block
					X = (rSilencer & reverseBitMask[shift + rStartBlock]) | VN[rStartBlock];
			
					temp = VP[rStartBlock] + (X & VP[rStartBlock]);
					carryD0 = temp < VP[rStartBlock];
					
					D0 = (temp ^ VP[rStartBlock]) | X;
					HN = VP[rStartBlock] & D0;
					HP = VN[rStartBlock] | ~(VP[rStartBlock] | D0);
					
					X = (HP << 1) | (1 << rOffSet);
					carryHP = HP >> (BLOCK_SIZE - 1);
					
					VN[rStartBlock] = X & D0;

					temp = (HN << 1);
					carryHN = HN >> (BLOCK_SIZE - 1);
										
		 			VP[rStartBlock] = temp | ~(X | D0);

					// compute the remaining blocks
					for (currentBlock = rStartBlock + 1; currentBlock <= rEndBlock; currentBlock++) {
						X = reverseBitMask[shift + currentBlock] | VN[currentBlock];
				
						temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;
						
						carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);
					
						D0 = (temp ^ VP[currentBlock]) | X;
						HN = VP[currentBlock] & D0;
						HP = VN[currentBlock] | ~(VP[currentBlock] | D0);
						
						X = (HP << 1) | carryHP;
						carryHP = HP >> (BLOCK_SIZE-1);
						
						VN[currentBlock] = X & D0;

						temp = (HN << 1) | carryHN;
						carryHN = HN >> (BLOCK_SIZE - 1);
											
		 				VP[currentBlock] = temp | ~(X | D0);
					}

					if (HP & rScoreMask)
						--score;
					else if (HN & rScoreMask)
						++score;

					c_score[pos] += score;
					
					/* check for optimality -- score extension*/
					if(c_score[pos] > max)
					{
						pos_max = pos;
						max = c_score[pos];
						rmax = score;
					}
					
					--pos;
				}

			}  /* end - long pattern */			
			/* compute with myers - reverse - end */

			// if computed the whole matrix max = alignment score
			if(target == hs_complete)
				total_score = max;

#ifdef MYERS_HIRSCHBERG_VERBOSE
			printf("Optimal cut is at %i and %i with forward score %i and reverse score %i\n\n",mid,pos_max,(max - rmax),rmax);
#endif
			/* push the two computed parts of the dp-matrix on process stack */
			to_process.push(_HirschbergSet(pos_max,_end1(target),mid,_end2(target),rmax));
			to_process.push(_HirschbergSet(_begin1(target),pos_max,_begin2(target),mid,max - rmax));

		}
		/* END CUT */
	}
	
	/* clean up */
	deallocate(align_, forwardBitMask, alphabetSize * blockCount);
	deallocate(align_, reverseBitMask, alphabetSize * blockCount);

	deallocate(align_, VP, blockCount);
	deallocate(align_, VN, blockCount);

	return total_score;
}

} // end namespace
#endif // end ifndef
