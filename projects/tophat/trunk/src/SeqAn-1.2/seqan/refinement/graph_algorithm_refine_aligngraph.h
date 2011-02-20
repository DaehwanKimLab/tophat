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
  $Id: graph_algorithm_refine_aligngraph.h 1757 2008-02-27 16:26:20Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_REFINE_ALIGRAPH_H
#define SEQAN_HEADER_GRAPH_REFINE_ALIGRAPH_H

//SEQAN_NO_DDDOC: do not generate documentation for this file


namespace SEQAN_NAMESPACE_MAIN
{

	
/////////////////////// NOCH NICHT EINSATZF�HIG!!!!! /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functios for Align Graphs
//project onto other sequence for Graph<Alignment>
template<typename TAlignment,typename TValue, typename TMap>
void
_getOtherSequenceAndProject(Graph<TAlignment> & segment, 
                            TMap &, 
 							TValue seq_i_id, 
							TValue pos_i, 
							TValue & seq_j_id, 
							TValue & pos_j)
{
	getProjectedPosition(segment,seq_i_id, pos_i,seq_j_id,pos_j);

}




//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TValue,typename TAlign,typename TId>
void
_getSeqBeginAndEnd(Graph<TAlign> & segment,
				  std::map<const void * ,int> &, 
				  TId & seq_i_id, 
				  TValue & begin_i, 
				  TValue & end_i,
				  TValue seq)
{
SEQAN_CHECKPOINT
	//walk through edges, take first edge, target, source,
	//define: seq == 0   ==> seq_i_id = id of source of first edge
	//        seq == 1   ==> seq_i_id = id of target of first edge
	typedef Graph<TAlign> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	
	TEdgeIterator ed_it(segment);
	//goBegin(ed_it);
	if(seq==0)
	{
		TVertexDescriptor src = sourceVertex(ed_it);
		seq_i_id = sequenceId(segment,src);
		
	}
	else
	{
		TVertexDescriptor trg = targetVertex(ed_it);
		seq_i_id = sequenceId(segment,trg);
	}
	begin_i = getFirstCoveredPosition(segment,seq_i_id);
	end_i = getLastCoveredPosition(segment,seq_i_id);


}




//////////////////////////
//for Graph<TAlign>
//vorsichtig! noch nicht richtig, bis jetzt nur ungapped exact matches...
//template<typename TScore,typename TStringSet,typename TAlignment,typename TValue>
//typename Value<TScore>::Type
//getScore(TScore & score_type,
//		 TStringSet & seqs,
//		 Graph<TAlignment> & segment,
//		 TValue pos_i,
//		 TValue pos_j,
//		 TValue len)
//{
//SEQAN_CHECKPOINT
//	int pseudo_map = 0;
//	TValue pos_j_check,seq_j_id;
//	_getOtherSequenceAndProject(segment,pseudo_map,seq_i_id,pos_i,seq_j_id,pos_j_check);
//	SEQAN_TASSERT(pos_j_check==pos_j)
//	TValue last_pos_i = pos_i + len;
//	TValue last_pos_j;
//	_getOtherSequenceAndProject(segment,pseudo_map,seq_i_id,last_pos_i,seq_j_id,last_pos_j);
//
////	typename Infix<typename Value<TStringSet>::Type>::Type label0 = infix(getValueById(seqs,seq_i_id),pos_i,last_pos_i);
////	typename Infix<typename Value<TStringSet>::Type>::Type label1 = infix(getValueById(seqs,seq_j_id),pos_j,last_pos_j);
//	
////b����������������������
//	typename Value<TScore>::Type score = 0;
//	TValue i = 0;
//	while (i < len)
//	{
//		next_pos_j = getProjectedPosition(segment,seq_i_id,pos_i);
//		//gaps
//		if(pos_j+1 != next_pos_j)
//		{
//			if(pos_j == next_pos_j)
//			{
//				score += scoreGapExtend(score_type);
//				++pos_i;
//				++i;
//				continue;
//			}
//		}		
//		pos_j = next_pos_j;
//		score += score(score_type,getValueById(seqs,seq_i_id)[pos_i],getValueById(seqs,seq_j_id)[pos_j]);
//		++i;
//		++pos_i;
//	}
//
//
//
//
//
//
//
//
//	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,stringSet(segment)[0]);
//	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,stringSet(segment)[1]);
//
//
//	//typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,0);
//	//typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,1);
//	int i = 0;
//	typename Value<TScore>::Type ret_score = 0;
//	while(i < len)
//	{
//		ret_score += score(score_type,label0[i],label1[i]);
//		++i;
//	}
//	//ret_score = scoreMatch(score_type);
//	//ret_score *= len;
//	return ret_score;
//}				


//////////////////////////
//for Graph<TAlign>
//only for complete refinement 
template<typename TScoreValue,typename TScoreSpec,typename TStringSet,typename TAlignment,typename TValue>
TScoreValue
getScore(Score<TScoreValue,TScoreSpec> & score_type,
		 TStringSet & seqs,
		 Graph<TAlignment> & segment,
		 TValue pos_i,
		 TValue pos_j,
		 TValue len,
		 TValue len_j)
{
SEQAN_CHECKPOINT

	typedef Graph<TAlignment> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename TGraph::TPosToVertexMap TPosToVertexMap;
	typedef typename TPosToVertexMap::const_iterator TVertexMapIter;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	TValue seq_i_id;
	TEdgeIterator ed(segment);
	//goBegin(ed);
	seq_i_id = sequenceId(segment,sourceVertex(ed));
		
	int pseudo_map = 0;
	TValue pos_j_check,seq_j_id;
	_getOtherSequenceAndProject(segment,pseudo_map,seq_i_id,pos_i,seq_j_id,pos_j_check);
	SEQAN_TASSERT(pos_j_check==pos_j)
	TValue last_pos_j = pos_j + len_j;
		
	TScoreValue ret_score = 0;
	bool last_one_was_aligned = false;

	while(len != 0)
	{
		TValue rest = 0;
		TVertexMapIter it = segment.data_pvMap.upper_bound(std::make_pair(seq_i_id, pos_i));
		// vd definitely exists (unless the input parameters for getScore are wrong)	
		// no! falsch gedacht... nilVertex, falls pos in einem gap liegt
		if(it->second == nilVertex)
		{	
			++it;
			if(it != segment.data_pvMap.end() && it->first.first == seq_i_id)
			{
				rest = fragmentBegin(segment,it->second)-pos_i;
				last_one_was_aligned = false;
				if(rest < len)//add rest many gaps
						ret_score += rest * scoreGapExtend(score_type);
				else
				{//add len many gaps
					ret_score += len * scoreGapExtend(score_type);
					return ret_score;	//and done!
				}
			}
			else
			{//add len many gaps
				ret_score += len * scoreGapExtend(score_type);
				return ret_score;	//and done!
			}
		}
		else{
			TVertexDescriptor vd = it->second;
			rest = fragmentBegin(segment,vd)+fragmentLength(segment,vd)-pos_i;
			TOutEdgeIterator ed_it(segment,vd);
			if(!atEnd(ed_it)) //aligned stretch
			{
				if(last_one_was_aligned)
				{
					TValue next_pos_j,temp;
					getProjectedPosition(segment,seq_i_id,pos_i,temp,next_pos_j);
					ret_score += (next_pos_j-pos_j) * scoreGapExtend(score_type);
					pos_j = next_pos_j; 
				}
				last_one_was_aligned = true;
				TValue i = 0;
				while(i < rest && i < len)
				{
					ret_score += score(score_type,getValueById(seqs,seq_i_id)[pos_i++],getValueById(seqs,seq_j_id)[pos_j++]);
					++i;
				}
				if(rest>len) return ret_score; //done (last_pos_i is somewhere inside the current node)
				else 
				{//last_pos_i is the first position of the next node
					if(rest == len) //check if there is an unalgined stretch on seq_j_id that needs to be included
					{
						if(pos_j != last_pos_j)
							ret_score += (last_pos_j-pos_j) * scoreGapExtend(score_type);
					}
				}
			}
			else //gap
			{
				last_one_was_aligned = false;
				if(rest < len)//add rest many gaps
						ret_score += rest * scoreGapExtend(score_type);
				else
				{//add len many gaps
					ret_score += len * scoreGapExtend(score_type);
					return ret_score;	//and done!
				}
			}
		}
		len -= rest;	
	}


	return ret_score;
}				


}
#endif //#ifndef SEQAN_HEADER_...
