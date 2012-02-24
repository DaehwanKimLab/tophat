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

#ifndef SEQAN_HEADER_GRAPH_REFINE_ALIGN_H
#define SEQAN_HEADER_GRAPH_REFINE_ALIGN_H



namespace SEQAN_NAMESPACE_MAIN
{


///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Align<TSource,TSpec>
//project onto other sequence 
template<typename TSource,typename TSpec,typename TValue,typename TMap>
void
_getOtherSequenceAndProject(Align<TSource,TSpec> & segment, 
				TValue seg_num,
							TMap & seq_map, 
						   TValue , 
						   TValue node_i, 
						   TValue & seq_j_id, 
						   TValue & node_j)
{
SEQAN_CHECKPOINT

	if(seg_num == 0)
	{
		seq_j_id = seq_map[id(source(row(segment,1)))];
		node_j = toSourcePosition(row(segment,1),toViewPosition(row(segment,0),node_i));
	}
	else
	{
		seq_j_id  = seq_map[id(source(row(segment,0)))];
		node_j = toSourcePosition(row(segment,0),toViewPosition(row(segment,1),node_i));
	}
}


//unspektakul�re funktion, die die int ID zur�ckgibt (braucht man damit es f�r alle alignment typen geht)
//template<typename TSource,typename TSpec, typename TValue, typename TSeqMap>					
//int 
//_getSeqMapId(TSeqMap & seq_map,
//			Align<TSource,TSpec> & segment,
//			TValue seq_i)
//{
//SEQAN_CHECKPOINT
//	return seq_map[id(source(row(segment,seq_i)))];
//}
//
//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TAliSource,typename TAliSpec, typename TValue>
void
_getSeqBeginAndEnd(Align<TAliSource,TAliSpec> & segment,
				  std::map<const void * ,int> & seq_map, 
				  TValue & seq_i_id, 
				  TValue & begin_i, 
				  TValue & end_i,
				  TValue seq)
{
	seq_i_id = seq_map[id(source(row(segment,seq)))];
	begin_i = clippedBeginPosition(row(segment,seq));
	end_i = clippedEndPosition(row(segment,seq));
}



////////////////////////////////////////////////////////////////////////////////////////
// 50000 _getRefinedMatchScore Functions
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////
//for Align<TAliSource,TAliSpec>
//get score for alignment of length len starting at pos_i on first sequence
//and pos_j on second sequence 
template<typename TScoreValue,typename TScoreSpec,typename TStringSet,typename TAliSource,typename TAliSpec,typename TValue>
TScoreValue
_getRefinedMatchScore(Score<TScoreValue,TScoreSpec> & score_type,
		 TStringSet &,
		 Align<TAliSource,TAliSpec> & segment,
		 TValue pos_i,
		 TValue pos_j,
		 TValue len,
		 TValue)
{
SEQAN_CHECKPOINT     
	typedef Align<TAliSource,TAliSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
//	typedef typename Iterator<TRow,GapsIterator<ArrayGaps> >::Type TIterator;	
	typedef typename Iterator<TRow, Rooted>::Type TIterator;
	TIterator row0_it, row1_it;
	row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_i));
	row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_j));
	len = toViewPosition(row(segment,0),pos_i + len) - toViewPosition(row(segment,0),pos_i);
	TValue i = 0;
	TScoreValue ret_score = 0;
	while(i < len)
	{
		if(isGap(row1_it)||isGap(row0_it))
			ret_score += scoreGapExtend(score_type);
		else
			ret_score += score(score_type,getValue(row0_it),getValue(row1_it));
		++i;
		++row0_it; 
		++row1_it; 
	}
	return ret_score;
}				
					

//get score for alignment starting at pos_i on first sequence 
//and pos_j on second sequence, if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
//template<typename TScore,typename TStringSet, typename TAliSource,typename TAliSpec,typename TValue>
//typename Value<TScore>::Type
//_getRefinedMatchScore(TScore & score_type,
//		 TStringSet &, 
//		 Align<TAliSource,TAliSpec> & segment,
//		 TValue pos_i,
//		 TValue pos_j,
//		 TValue len1,
//		 TValue len2)
//{
//SEQAN_CHECKPOINT
//	typedef Align<TAliSource,TAliSpec> TAlign;
//	typedef typename Row<TAlign>::Type TRow;
//	typedef typename Iterator<TRow>::Type TIterator;	
//	TIterator row0_it, row1_it;
//	TValue len;
//	row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_i));
//	row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_j));
//	len1 = toViewPosition(row(segment,0),pos_i + len1) - toViewPosition(row(segment,0),pos_i);
//	len2 = toViewPosition(row(segment,1),pos_j + len2) - toViewPosition(row(segment,1),pos_j);
//	len = (len1 < len2) ? len1 : len2;
//	int i = 0;
//	typename Value<TScore>::Type ret_score = 0;
//	
//	//calculate score for aligned region
//	while(i < len)
//	{
//		ret_score += score(score_type,getValue(row0_it),getValue(row1_it));
//		++i;
//		++row0_it;
//		++row1_it;
//	}
//	//fill up with gaps if one sequence is longer than the other
//	len = (len1 > len2) ? len1 : len2;
//	ret_score += (len - i) * scoreGapExtend(score_type);
//	
//	return ret_score;
//}				


}
#endif //#ifndef SEQAN_HEADER_...
