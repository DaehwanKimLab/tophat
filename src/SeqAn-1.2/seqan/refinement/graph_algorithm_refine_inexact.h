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
  $Id: graph_algorithm_refine_inexact.h 1757 2008-02-27 16:26:20Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_REFINE_INEXACT_H
#define SEQAN_HEADER_GRAPH_REFINE_INEXACT_H


//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{


struct TagInexactRefinement_;
typedef Tag<TagInexactRefinement_> const InexactRefinement;


///inexact refinement (cuts that would produce segments shorter than min_len are not made)
template<typename TValue>
inline bool
cutIsOk(String<std::set<TValue> > & all_nodes,
		TValue seq_i_pos,
		TValue pos_i,
		typename std::set<TValue>::iterator iter,
		TValue min_len,
		Tag<TagInexactRefinement_> const)
{
SEQAN_CHECKPOINT
	
	//cut already exists
	if(iter != all_nodes[seq_i_pos].end())
		return false;
	typename std::set<TValue>::iterator tmp_iter = all_nodes[seq_i_pos].upper_bound(pos_i);
	if(tmp_iter != all_nodes[seq_i_pos].end())
		if((*tmp_iter - pos_i) < min_len)
			return false;
	if(tmp_iter != all_nodes[seq_i_pos].begin())
	{
		--tmp_iter;
		if((pos_i - *tmp_iter) < min_len)
			return false;
	}
	return true;
}






//////////////////////////////////////////////////////////////////////////////////////////////
//step 2 of constructing the refined alignment graph: add all edges    
//version for inexact refinement
template<typename TAliGraph, typename TVertexDescriptor,typename TValue>
TValue  
_getClosestRefinedNeighbor(TAliGraph & ali_g,
						   TVertexDescriptor & vd,
						   TValue seq,
						   TValue pos)
{
SEQAN_CHECKPOINT
	if(pos-fragmentBegin(ali_g,vd) < fragmentBegin(ali_g,vd)+fragmentLength(ali_g,vd)-pos)
		return fragmentBegin(ali_g,vd);
	else
		return fragmentBegin(ali_g,vd) + fragmentLength(ali_g,vd);
}



template<typename TAliGraph,typename TValue>
void
_getCutEndPos(TAliGraph & ali_g, 
			  typename VertexDescriptor<TAliGraph>::Type & end_knot,
			  TValue seq,
			  TValue end_pos,
			  TValue & cut_end_pos)
{
SEQAN_CHECKPOINT
	end_knot = findVertex(ali_g,seq,end_pos-1);//end_pos1 is the first position of the next node
	if(end_pos == fragmentBegin(ali_g,end_knot) + fragmentBegin(ali_g,end_knot))
		cut_end_pos = end_pos;
	else
	{
		cut_end_pos = _getClosestRefinedNeighbor(ali_g,end_knot,seq,end_pos);
		end_knot =  findVertex(ali_g,seq,cut_end_pos-1);
		SEQAN_TASSERT(cut_end_pos == fragmentBegin(ali_g,end_knot)+fragmentLength(ali_g,end_knot))
	}
		
}


template<typename TAliGraph,typename TValue>
void
_getCutBeginPos(TAliGraph & ali_g, 
			  typename VertexDescriptor<TAliGraph>::Type & act_knot,
			  TValue seq,
			  TValue act_pos,
			  TValue & cut_act_pos)
{
SEQAN_CHECKPOINT
	
	act_knot = findVertex(ali_g,seq,act_pos);
	//if completely refined
	if(act_pos == fragmentBegin(ali_g,act_knot))
		cut_act_pos = act_pos;
	else //if incompletely refined
	{
		cut_act_pos = _getClosestRefinedNeighbor(ali_g,act_knot,seq,act_pos);
		act_knot =  findVertex(ali_g,seq,cut_act_pos);
		SEQAN_TASSERT(cut_act_pos == fragmentBegin(ali_g,act_knot))
	}
}


//step 2 of constructing the refined alignment graph: add all edges    
//version for inexact refinement
template<typename TAlignmentString,typename TPropertyMap,typename TStringSet,typename TSeqMap, typename TScore,typename TAliGraph>
void
_makeRefinedGraphEdges(TAlignmentString & alis,
					   TPropertyMap & pm,
					  TStringSet & seqs,
				      TSeqMap & seq_map,
				      TScore & score_type,
					  TAliGraph & ali_g,
					  Tag<TagInexactRefinement_> const)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Value<TAlign>::Type TValue;
	typedef typename Iterator<TAlignmentString, Rooted>::Type TAliIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAliGraph>::Type TEdgeDescriptor;
	typedef typename Cargo<TAliGraph>::Type TCargo;
	//make edges
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);
	//for each segment/fragment/alignment
	while(ali_it != ali_end)
	{
		//get first sequence that takes part in the alignment + boundaries of the ali
		TValue seq1,begin_pos1,end_pos1;
		_getSeqBeginAndEnd(*ali_it,seq_map,seq1,begin_pos1,end_pos1,(TValue)0);
		//get the last node that is within the current ali
		TVertexDescriptor end_knot1;
		TValue cut_end_pos1;
		_getCutEndPos(ali_g,end_knot1,seq1,end_pos1,cut_end_pos1);
	
		//get the node that represents the current interval (begin_pos until next_cut_pos or end_pos)
		TVertexDescriptor act_knot1;
		TValue cut_act_pos1,act_pos1;
		act_pos1 = begin_pos1;
		_getCutBeginPos(ali_g,act_knot1,seq1,act_pos1,cut_act_pos1);
		TValue act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
		//walk through cuts on the first sequence
//		while (act_end_pos1 <= cut_end_pos1)
		while (true)
		{
			//get other sequence and projected position
			TValue seq2,act_pos2;
			_getOtherSequenceAndProject(*ali_it,seq_map,seq1,act_pos1,seq2,act_pos2);
		
			//get node that corresponds to that position
			TVertexDescriptor act_knot2;
			TValue cut_act_pos2;
			_getCutBeginPos(ali_g,act_knot2,seq2,act_pos2,cut_act_pos2);
			//corresponding end on seq2 (there might be more than one node on seq2 that corresponds
			//to the same interval (=node) on seq1)
			TValue act_end_pos2;
			_getOtherSequenceAndProject(*ali_it,seq_map,seq1,act_end_pos1-1,seq2,act_end_pos2);
			++act_end_pos2;
			TVertexDescriptor act_end_knot2;
			TValue cut_act_end_pos2;
			_getCutEndPos(ali_g,act_end_knot2,seq2,act_end_pos2,cut_act_end_pos2);
			
			if(cut_act_pos2 == cut_act_end_pos2)
				break;
			while(true)
			{
				//should at the moment return score for:
				//
				//seq1 = ....cr...rc....
				//            ||||||
				//seq2 = ...c.r...rc....
				//bzw
				//seq1 = ..cr.....x....   man will aber nur    ..cr......x....
				//          |||||||-							 ---||||||  
				//seq2 = ...r.c...rc... 					   ...r.c...rc....
				typename Value<TScore>::Type score = 0;
				score = getScore(score_type,seqs,*ali_it,act_pos1,act_pos2,act_end_pos1-act_pos1,cut_act_end_pos2);
				//score *= getAnnoScore(ali_g,pm,act_knot1,act_knot2,score_type);
				//add score for
				//
				//seq1 = ...-cr....x....
				//          ||
				//seq2 = ...c.r...rc....
//					score += getLeftRestScore(score_type,seqs,seq1,seq2,act_pos1,cut_act_pos1,act_pos2,cut_act_pos2);
				if(score > 0)
					if(findEdge(ali_g,act_knot1,act_knot2)==0)
						addEdge(ali_g,act_knot1,act_knot2,score);
				
				if(act_knot2==act_end_knot2)
					break;
				act_pos2 = cut_act_pos2 + fragmentLength(ali_g,act_knot2);
				_getCutBeginPos(ali_g,act_knot2,seq2,act_pos2,cut_act_pos2);
			}
			if(act_knot1 == end_knot1)
				break;
			act_pos1 = act_end_pos1;
			act_knot1 = findVertex(ali_g,seq1,act_pos1);
			cut_act_pos1 = act_pos1;
			act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
		}
		++ali_it;
	}
}




/**
.Function.matchRefinement:
..signature:matchRefinement(matches,stringSet,scoringScheme,refinedGraph,minFragmentLen)
..param.minFragmentLen:The minimal segment length allowed (unsigned int).
...remarks:If in the refinement process a cut would result in
a segment shorter than minFragmentLen, then the cut is not made and a heuristic is applied to refine this short overlap.
...remarks:If no minFragmentLen is given, then all cuts are made. This corresponds to a minFragmentLen of 1.
*/
//score type given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				Score<TScoreValue,TScoreSpec> & score_type,
				TOutGraph & ali_graph,
				unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
	bool anno = false;
	if(min_frag_len > 1)
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,InexactRefinement());
	else
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,ExactRefinement());
}


/**
.Function.matchRefinement:
..signature:matchRefinement(matches,stringSet,refinedGraph,minFragmentLen)
*/
//score type not given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TOutGraph & ali_graph,
				unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
//	Score<int,FakeScore > fake_score;
	typename Cargo<TOutGraph>::Type fake_score = 1;
	bool anno = false;
	if(min_frag_len > 1)
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,InexactRefinement());
	else
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,ExactRefinement());
}
	
}
#endif //#ifndef SEQAN_HEADER_...
