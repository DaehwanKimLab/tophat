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
  $Id: graph_algorithm_refine_fragment.h 1911 2008-05-02 09:28:04Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_REFINE_FRAG_H
#define SEQAN_HEADER_GRAPH_REFINE_FRAG_H



namespace SEQAN_NAMESPACE_MAIN
{



	
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Fragments
//project onto other sequence for Graph<Alignment>
template<typename TFragSize, typename TFragSpec,typename TValue, typename TMap>
void
_getOtherSequenceAndProject(Fragment<TFragSize,TFragSpec> & segment,
						   TMap &,
						   TValue seq_i_id,
						   TValue pos_i,
						   TValue & seq_j_id,
						   TValue & pos_j)
{
SEQAN_CHECKPOINT
	getProjectedPosition(segment,seq_i_id, pos_i,seq_j_id,pos_j);
	
	//if(seq_i_id == sequenceId(segment,0))
	//	seq_j_id = sequenceId(segment,1);
	//else
	//	seq_j_id = sequenceId(segment,0);
}


//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TFragSize, typename TFragSpec, typename TValue>
void
_getSeqBeginAndEnd(Fragment<TFragSize,TFragSpec> & segment,
				  std::map<const void * ,int> &, 
				  TValue & seq_i_id, 
				  TValue & begin_i, 
				  TValue & end_i,
				  TValue seq)
{
SEQAN_CHECKPOINT
	seq_i_id = sequenceId(segment,seq);
	begin_i = fragmentBegin(segment,seq_i_id);
	end_i = begin_i + fragmentLength(segment,seq_i_id);
}


////////////////////////////////////////////////////////////////////////////////////////
// 50000 getScore Functions
////////////////////////////////////////////////////////////////////////////////////////
//get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
//template<typename TScore,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec,typename TValue>
//typename Value<TScore>::Type
//getScore(TScore & score_type, 
//		 TStringSet & seqs,
//		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
//		 TValue pos_i, 
//		 TValue pos_j,
//		 TValue len1, 
//		 TValue len2)
//{
//SEQAN_CHECKPOINT
//	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
//	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));
//	int i = 0;
//	typename Value<TScore>::Type ret_score = 0;
//	TValue len = (len1 < len2) ? len1 : len2;
//	while(i < len)
//	{
//		ret_score += score(score_type,label0[i],label1[i]);
//		++i;
//	}
//	len = (len1 > len2) ? len1 : len2;
//	ret_score += (len - i) * scoreGapExtend(score_type);
//	return ret_score;
//}				
//
//
///////////////////////////////////////////////////////////////////////////////////////////////
////die nächsten beiden funktionen: für Fragmente und Score vom Typ Simple
////für den fall dass es keine mismatches innerhalb der segmente gibt und Score vom typ Simple ist
////TODO: müsste für einen bestimmten TFragSpec sein (Exact oder noMismatches)
////get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
////and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
////process was stopped (the cut is not exact)
//template<typename TScoreValue,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
//TScoreValue
//getScore(Score<TScoreValue, Simple> & score_type,
//		 TStringSet & seqs, 
//		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
//		 TFragPos pos_i, 
//		 TFragPos pos_j, 
//		 TFragSize len1, 
//		 TFragSize len2)
//{
//SEQAN_CHECKPOINT
//	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
//	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));
//	TScoreValue ret_score = 0;
//	TFragSize len;
//	if (len1 < len2) len = len1;
//	else len = len2;
//	if(len1 <= len2)
//	{
//		ret_score += len1 * scoreMatch(score_type);
//		ret_score += (len2 - len1) * scoreGapExtend(score_type);
//	}
//	else{
//		ret_score += len2 * scoreMatch(score_type);
//		ret_score += (len1 - len2) * scoreGapExtend(score_type);
//	}
//	return ret_score;
//}				


//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScoreValue,typename TScoreSpec,typename TStringSet,typename TFragment,typename TFragPos,typename TFragSize>
TScoreValue
getScore(Score<TScoreValue,TScoreSpec> & score_type,
		 TStringSet & seqs,
		 TFragment& segment,
		 TFragPos pos_i,
		 TFragPos pos_j,
		 TFragSize len,
		 TFragSize)
{
SEQAN_CHECKPOINT
	typedef typename Infix<typename Value<TStringSet>::Type>::Type TSegmentLabel;
	TSegmentLabel label0 = label(segment,seqs, sequenceId(segment, 0));
	TSegmentLabel label1 = label(segment,seqs, sequenceId(segment, 1));
	typename Iterator<TSegmentLabel, Rooted>::Type label_it0 = begin(label0) + (pos_i - fragmentBegin(segment,sequenceId(segment,0)));
	typename Iterator<TSegmentLabel, Rooted>::Type label_it1 = begin(label1) + (pos_j - fragmentBegin(segment,sequenceId(segment,1)));
	int i = 0;
	TScoreValue ret_score = 0;
	while(i < (int) len)
	{
		ret_score += score(score_type,*label_it0,*label_it1);
		++label_it0;
		++label_it1;
		++i;
	}
	return ret_score;
}				


//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScoreValue,typename TStringSet,typename TFragPos,typename TFragSize, typename TSpec>
TScoreValue
getScore(Score<TScoreValue, Simple> & score_type,
		 TStringSet &,
		 Fragment<TFragSize,ExactFragment<TSpec> > &,
		 TFragPos,
		 TFragPos,
		 TFragSize len,
		 TFragSize)
{
SEQAN_CHECKPOINT
	return len*scoreMatch(score_type);
}				


}
#endif //#ifndef SEQAN_HEADER_...
