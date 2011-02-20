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
  $Id: graph_align_hirschberg.h 1812 2008-03-31 15:54:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H
#define SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Hirschberg Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


// Note: The hirschberg trace is different from all the other traces because it simply contains the trace points!!!
template <typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TTrace>
void
_align_hirschberg_trace(TAlign& align,
						TStringSet const& str,
						Score<TScoreValue, TSpec> const& sc,
						TTrace& trace)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef unsigned char TTraceValue;   
	
	// TraceBack values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
		
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = length(str[0]);
	// Current trace point (currentPointer,len1)
	TSize currentPointer = getValue(trace, len1);
	// Next trace point (movePointer, len1 - 1);
	TSize movePointer = getValue(trace, len1 - 1);
	do {
		// Diagonal walk
		if (((currentPointer - movePointer) == 1)  &&
			(scoreGapOpen(sc) + scoreGapOpen(sc) <= score(const_cast<Score<TScoreValue, TSpec>&>(sc), len1-1, movePointer, str[0], str[1]))) {
			TSize segLen = 0;
			while (((currentPointer - movePointer) == 1)   &&
					(scoreGapOpen(sc) + scoreGapOpen(sc) <= score(const_cast<Score<TScoreValue, TSpec>&>(sc), len1-1, movePointer, str[0], str[1]))) {
				++segLen;
				--len1;
				--currentPointer;
				if (len1>0) movePointer = getValue(trace, len1 - 1);
				else break;
			}
			_align_trace_print(align, str, id1, len1, id2, currentPointer, segLen, Diagonal);
		} // Diagonal walk, but gapOpen and gapOpen yields a higher score
		// Note: The case gapExtension + gapOpen cannot occur because then we would have gotten a different trace point!!!
		else if (((currentPointer - movePointer) == 1)  &&
			(scoreGapOpen(sc) + scoreGapOpen(sc) > score(const_cast<Score<TScoreValue, TSpec>&>(sc), len1-1, movePointer, str[0], str[1]))) {
				_align_trace_print(align, str, id1, (TSize) 0, id2, --currentPointer, (TSize) 1, Vertical);
				_align_trace_print(align, str, id1, --len1, id2,(TSize) 0, (TSize) 1, Horizontal);
				if (len1>0) movePointer = getValue(trace, len1 - 1);
		} // Horizontal walk 
		else if ((currentPointer - movePointer) == 0) {
			TSize segLen = 0;
			while ((currentPointer - movePointer) == 0) {
				++segLen;
				--len1;
				if (len1>0) movePointer = getValue(trace, len1 - 1);
				else break;
			}
			_align_trace_print(align, str, id1, len1, id2, (TSize) 0, segLen, Horizontal);
		} // Vertical walk 
		else {
			// Vertical walk finishes with a horizontal vs. vertical gap
			if ((len1 > 2) && 
				(movePointer - getValue(trace, len1 - 2) == 0) &&
				(scoreGapExtend(sc) + scoreGapExtend(sc) > score(const_cast<Score<TScoreValue, TSpec>&>(sc), len1-1, movePointer, str[0], str[1]))) {
					_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) movePointer, (TSize) currentPointer - movePointer, Vertical);
					currentPointer = movePointer;
			} // Vertical walk that finishes with a diagonal, but gapOpen + gapOpen is better
			else if ((scoreGapOpen(sc) + scoreGapOpen(sc) > score(const_cast<Score<TScoreValue, TSpec>&>(sc), len1-1, movePointer, str[0], str[1]))) {
					_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) movePointer, (TSize) currentPointer - movePointer, Vertical);
					currentPointer = movePointer;
			} // Normal vertical walk that continues with a diagonal walk
			else {
				_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) movePointer + 1, (TSize) currentPointer - (movePointer + 1), Vertical);
				currentPointer = movePointer + 1;
			}
		}
	} while (len1 != 0);
	// Maybe we have to align the last vertical piece
	if (getValue(trace, 0) != 0) {
		_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) 0, (TSize) getValue(trace, 0), Vertical);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue, typename TSpec>
TScoreValue
_align_hirschberg(TTrace& trace,
				  TStringSet const& str,
				  Score<TScoreValue, TSpec> const& sc,
				  Hirschberg) 
{
	SEQAN_CHECKPOINT
	typedef String<TScoreValue> TColumn;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef typename Value<TStringSet const>::Type TString;
	typedef typename Iterator<TString const, Standard>::Type TStringIter;
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TCharacter;
	typedef unsigned char TTraceValue;   

	// TraceBack values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue maxVal = 0; // The final score

	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	resize(trace, len1+1);  // The set of trace positions
	typedef std::pair<TSize, TSize> TPoint;   // Trace Point...
	typedef std::pair<TPoint, TTraceValue> TListElement;  // ... and how we traversed it
	typedef std::list<TListElement> TMidPointQueue;  // List of tracepoints -> Iterative, not recursive
	TMidPointQueue midpoints;
	
	// Walk multiple time through the whole DP Matrix
	do {
		// Initialization
		TSize x1 = 0;	// Current alignment window
		TSize x2 = 0;
		TSize y1 = 0;
		TSize y2 = 0;
		bool firstRun = true;		
		bool wasHorizontal = false;		// How we entered the alignment window
		bool wasVertical = false;
		TScoreValue upperLeft = 0;		// The old score for the upper left initialization of the DP matrix
		TTraceValue oldWay = Diagonal;
		
		
	  
		// Now iterate through all tracepoints found so far
		typename TMidPointQueue::iterator it = midpoints.begin();
		do {
			TTraceValue tvMat, tvHorizontal, tvVertical;  // Where did the maximum for the diagonal, horizontal, and vertical matrix came from

			// Step1: Get the alignment window
			if (firstRun) {
				x1 = 0;
				y1 = 0;
				firstRun = false;
			} else {
				x1 = x2;
				y1 = y2;
				++it;
			}

			if (it == midpoints.end()) {
				x2 = len1;
				y2 = len2;
				// In this case we don't care how we aligned len1, len2 before
			} else {
				TListElement p = *it;
				x2 = p.first.first;
				y2 = p.first.second;
				oldWay = p.second;  // Old alignment direction is important, we have to restore it!!!
			}

			// Step2: Do the alignment
			TSize middle = x1 + ((x2-x1)/2);		// New cut-point in global coordinates
			TSize inf_middle = (TSize) ((x2-x1)/2);   // New cut-point in "infix" coordinates
			TSize inf_len1 = x2-x1;
			TSize inf_len2 = y2-y1;
			TScoreValue tmp = 0;		

			TScoreValue inf = infimumValue<TScoreValue>() / 2;
			// The DP Matrix for diagonal walks
			TColumn mat;
			resize(mat, (inf_len2+1));   // One column for the diagonal matrix
			
			// The DP Matrix for gaps from the left
			TColumn horizontal;
			fill(horizontal, (inf_len2+1), inf); 
			
			// The DP Matrix for gaps from the top, a single value is enough
			TScoreValue vert;
			vert = inf;



			// The pointers
			// Every back pointer stores the row it came from, so (middle, row) determines the trace point
			// Additionally, every pointer has to remember how its trace point was traversed
			// We do need diagonal, vertical, and horizontal pointers!!!
			// The vertical and horizontal once are sometimes fixed depending on the way how we leave the trace point!!!
			typedef std::pair<TSize, TTraceValue> TBackPointer;
			String<TBackPointer> verticalPointer;
			String<TBackPointer> horizontalPointer;
			String<TBackPointer> diagonalPointer;
			TBackPointer diagonalPointerOld;
			TBackPointer tmpPointer;

			// Initialization
			resize(horizontalPointer, (inf_len2+1));  
			resize(verticalPointer, (inf_len2+1));  
			resize(diagonalPointer, (inf_len2+1)); 
				
			// Classical DP
			assignValue(mat, 0, upperLeft);
			// Handle the base cases
			if ((wasHorizontal) || (inf_len2 == 0)) {
				assignValue(horizontalPointer, 0, std::make_pair(0, Horizontal));
				assignValue(verticalPointer, 0, std::make_pair(0, Horizontal));
				assignValue(diagonalPointer, 0, std::make_pair(0, Horizontal));
			} else if (wasVertical) {
				assignValue(horizontalPointer, 0, std::make_pair(0, Vertical));
				assignValue(verticalPointer, 0, std::make_pair(0, Vertical));
				assignValue(diagonalPointer, 0, std::make_pair(0, Vertical));
			} else {
				assignValue(horizontalPointer, 0, std::make_pair(0, Horizontal));
				assignValue(verticalPointer, 0, std::make_pair(0, Vertical));
				assignValue(diagonalPointer, 0, std::make_pair(0, Diagonal));
			}
			// Initialize the first column vector
			for(TSize row = 1; row <= inf_len2; ++row) {
				if (wasVertical) {
					assignValue(mat, row, upperLeft + gap + (row - 1) * gap);
					assignValue(horizontal, row, getValue(mat, row) + gapOpen - gap);
				} else {
					assignValue(mat, row, upperLeft + gapOpen + (row - 1) * gap);
					assignValue(horizontal, row, getValue(mat, row) + gapOpen - gap);
				}
				// All pointers in the first column point back to 0,0
				assignValue(verticalPointer, row, std::make_pair(0, Vertical));
				assignValue(horizontalPointer, row, std::make_pair(0, Vertical));
				assignValue(diagonalPointer, row, std::make_pair(0, Vertical));
			}
			// Now traverse all columns
//			TStringIter strIter0 = begin(str[0]);
//			strIter0 += x1;
			for(TSize col = 1, posX = x1; col <= inf_len1; ++col, ++posX) {
				//if ((x1+col) % 10 == 0) std::cout << x1+col << '-' << std::flush;

				TScoreValue diagValMat = getValue(mat, 0);
				// Initialize the first row (just the first element of each column vector)
				if (wasHorizontal) {
					assignValue(mat, 0, upperLeft + gap + (col - 1) * gap);
					vert = getValue(mat, 0) + gapOpen - gap;
				} else {			    
					assignValue(mat, 0, upperLeft + gapOpen + (col - 1) * gap);
					vert = getValue(mat, 0) + gapOpen - gap;
				}
				// If there are no rows we step horizontally !!!
				if (inf_len2 < 1) {
					tvHorizontal = Horizontal;
					tvMat = Horizontal;
					assignValue(horizontal,0,getValue(mat,0));
				}
				// Initialize the pointers
				if (col >= inf_middle) {
					diagonalPointerOld = getValue(diagonalPointer,0);  // Remember the old diagonal pointer
					// All points in the first row point back to 0,0
					assignValue(horizontalPointer, 0, std::make_pair(0, Horizontal));
					assignValue(diagonalPointer, 0, std::make_pair(0, Horizontal));
					assignValue(verticalPointer, 0, std::make_pair(0, Horizontal));
				}
//				TCharacter charTop = getValue(strIter0);
//				goNext(strIter0);
//				TStringIter strIter1 = begin(str[1]);
//				strIter1 += y1;
				for(TSize row = 1, posY = y1; row <= inf_len2; ++row, ++posY) {
					// Get the new maximum for vertical
					if ((tmp = getValue(mat, row - 1) + gapOpen) > vert + gap) {
						vert = tmp;
						tvVertical = Diagonal;
					} else {
						vert = vert + gap;
						tvVertical = Vertical;
					}
					
					// Get the new maximum for horizontal
					if ((tmp = getValue(mat, row) + gapOpen) > getValue(horizontal, row) + gap) {
						assignValue(horizontal, row, tmp);
						tvHorizontal = Diagonal;
					} else {
						assignValue(horizontal, row, getValue(horizontal, row) + gap);
						tvHorizontal = Horizontal;
					}
					
					// Get the new maximum for diagonal
					//TScoreValue sc_ = score(const_cast<Score<TScoreValue, Simple>&>(sc), charTop, getValue(str[1],y1+row-1));
					//tmp = diagValMat + sc_;
//					tmp = diagValMat + score(const_cast<Score<TScoreValue, TSpec>&>(sc), charTop, getValue(strIter1));	// deprecated
//					goNext(strIter1);
					tmp = diagValMat + score(const_cast<Score<TScoreValue, TSpec>&>(sc), posX, posY, str[0], str[1]);
					tvMat = Diagonal;
					if (vert > tmp) {
						tmp = vert;
						tvMat = Vertical;
					}
					if (getValue(horizontal, row) > tmp) {
						tmp = getValue(horizontal,row);
						tvMat = Horizontal;
					}
					
					// Assign the new diagonal value
					diagValMat = getValue(mat, row);
					assignValue(mat, row, tmp);

					// //Debug
					//std::cout << "Diagonal: " << row << ',' << col << ':' << getValue(mat, row) << ';';
					//std::cout << "Horizontal: " << row << ',' << col << ':' << getValue(horizontal, row) << ';';
					//std::cout << "Vertical: " << row << ',' << col << ':' << vert << std::endl;
					
					if (col >= inf_middle) {
						
						// If we are at the middle, initialize the pointers
						if (col == inf_middle) {
							assignValue(diagonalPointer, row, std::make_pair(row, tvMat));
							assignValue(verticalPointer, row, std::make_pair(row, tvMat));
							// If the horizontal value is as good as the vertical one, take horizontal because a gap extension is cheaper then a gap opening
							if ((tvMat == Vertical) && 
								(getValue(mat, row) == getValue(horizontal, row))) {
								assignValue(horizontalPointer, row, std::make_pair(row, Horizontal));
							} else assignValue(horizontalPointer, row, std::make_pair(row, tvMat));
							// Fix the vertical pointer
							if (tvVertical== Vertical) {
								// Overwrite the diagonal value
								// If it is needed it is still in the diagonal pointer
								if ((getValue(verticalPointer, row-1)).second == Diagonal) value(verticalPointer, row-1).second = Vertical;
							
								// If I have to pay a gapOpen anyway, I better pay it right away because more gap extensions might follow
								else if ((getValue(verticalPointer, row-1)).second == Horizontal) value(verticalPointer, row-1).second = Vertical;
							}
						}

						// Fix the horizontal pointer
						if (col == inf_middle + 1) {
							if (tvHorizontal == Horizontal) {
								// Overwrite the diagonal value
								// If it is needed it is still in the diagonal pointer
								if ((getValue(horizontalPointer, row)).second == Diagonal) value(horizontalPointer, row).second = Horizontal;

								// If I have to pay a gapOpen anyway, I better pay it right away because more gap extensions might follow
								else if ((getValue(horizontalPointer, row)).second == Vertical) value(horizontalPointer, row).second = Horizontal;
							}
						}

						// If we passed the middle we need to take care of the pointers
						if (col > inf_middle) {
							if (tvVertical == Diagonal) assignValue(verticalPointer, row, getValue(diagonalPointer, row - 1));
							else assignValue(verticalPointer, row, getValue(verticalPointer, row-1));
							if (tvHorizontal == Diagonal) assignValue(horizontalPointer, row, getValue(diagonalPointer, row));
						
							tmpPointer = getValue(diagonalPointer, row);
							if (tvMat == Diagonal) assignValue(diagonalPointer, row, diagonalPointerOld);
							else if (tvMat == Horizontal) assignValue(diagonalPointer, row, getValue(horizontalPointer, row));
							else if (tvMat == Vertical) assignValue(diagonalPointer, row, getValue(verticalPointer, row));
							diagonalPointerOld = tmpPointer;
						}
					}
				}
			}
			// Alignment is done, what is the best trace point?
			// If this is the first alignment iteration...
			if (it == midpoints.end()) {
				maxVal = getValue(mat, inf_len2);
				tmpPointer = getValue(diagonalPointer, inf_len2);
				if (getValue(horizontal, inf_len2) ==  maxVal) {
					tmpPointer = getValue(horizontalPointer, inf_len2);
				}
				if (vert == maxVal) {
					tmpPointer = getValue(verticalPointer, inf_len2);
				}
			} // ... and if it is not
			else {
				if (oldWay == Diagonal) {
					maxVal = getValue(mat, inf_len2);
					tmpPointer = getValue(diagonalPointer, inf_len2);
					wasVertical = false;
					wasHorizontal = false;
				} else if (oldWay == Horizontal) {
					maxVal = getValue(horizontal,inf_len2);
					tmpPointer = getValue(horizontalPointer, inf_len2);
					wasHorizontal = true;
					wasVertical = false;
				} else {
					maxVal = vert;
					tmpPointer = getValue(verticalPointer, inf_len2);;
					wasHorizontal = false;
					wasVertical = true;
				}
			}



			TSize cut = y1 + tmpPointer.first;
			upperLeft = maxVal;

			// Step3: Insert the new midpoint
			if ((x2 - x1) > 1) {
				it = midpoints.insert(it, std::make_pair(std::make_pair(middle,cut), tmpPointer.second));
				++it;
			}


			////Debug code
			//typedef typename Infix<TString const>::Type TInfix;
			//TInfix infix1 = infix(str[0], x1, x2);
			//TInfix infix2 = infix(str[1], y1, y2);
			//std::cout << infix1 << std::endl;
			//std::cout << infix2 << std::endl;
			//std::cout << "Trace-Point " << middle << ',' << cut << std::endl;
			//std::cout << "Max-Score " << maxVal << std::endl;
			//std::cout << "Direction " << (TSize) tmpPointer.second << std::endl;
			//std::cout << "--" << std::endl;

		} while (it != midpoints.end());

		//std::cout << midpoints.size() << ',';
		// Exit if all positions have a trace point
	} while ((midpoints.size() < len1 - 1));
	//std::cout << std::endl;

	// Transform the midpoints into simple trace positions
	TTraceIter traceIter = begin(trace, Standard() );
	typename TMidPointQueue::iterator it = midpoints.begin();
	if ((!midpoints.empty()) && (it->second == Diagonal)) assignValue(traceIter, it->first.second - 1); 
	else if ((!midpoints.empty()) && (it->second == Horizontal)) assignValue(traceIter, it->first.second); 
	else if ((!midpoints.empty()) && (it->second == Vertical)) {
		// Where is the best exit
		TScoreValue best = score(const_cast<Score<TScoreValue, TSpec>&>(sc), 0, 0, str[0], str[1]);
		TSize index = 0;
		if (it->first.second > 1) best += gapOpen;
		if (it->first.second > 2) best += gap * (it->first.second - 2);
		for(TSize row = 2; row < it->first.second;++row) {
			TScoreValue newBest = gapOpen + score(const_cast<Score<TScoreValue, TSpec>&>(sc), 0, row-1, str[0], str[1]) + gapOpen;
			newBest += gap * (it->first.second - 3);
			if (newBest > best) {
				best = newBest;
				index = row-1;
			}
		}
		assignValue(traceIter, index); 
	} else assignValue(traceIter, 0); 
	goNext(traceIter);
	for(;it!=midpoints.end();++it) {
		assignValue(traceIter, (it->first.second));
		goNext(traceIter);
	}
	assignValue(traceIter, len2); // End point of the trace
	
	return maxVal;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScore, typename TSpec>
inline typename Value<TScore>::Type
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 TScore const& sc,
				 AlignConfig<false, false, false, false, TSpec> const,
				 Hirschberg)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	
	// Trace
	String<TSize> trace;
	TScoreValue maxScore = _align_hirschberg(trace, str, sc, Hirschberg());

	// Follow the trace and create the graph
	_align_hirschberg_trace(align, str, sc, trace);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScore, typename TSpec>
inline typename Value<TScore>::Type
_globalAlignment(TStringSet const& str,
				 TScore const& sc,
				 AlignConfig<false, false, false, false, TSpec> const,
				 Hirschberg)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	String<TSize> trace;
	return _align_hirschberg(trace, str, sc, Hirschberg());	
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
