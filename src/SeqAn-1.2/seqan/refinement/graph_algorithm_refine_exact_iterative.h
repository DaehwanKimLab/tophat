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
  $Id: graph_algorithm_refine_exact.h 1919 2008-05-02 15:54:46Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_REFINE_EXACT_H
#define SEQAN_HEADER_GRAPH_REFINE_EXACT_H


namespace SEQAN_NAMESPACE_MAIN
{


	
struct TagExactRefinement_;
typedef Tag<TagExactRefinement_> const ExactRefinement;

//exact method, every cut is made (unless it already exists)
template<typename TValue>
inline bool
cutIsOk(String<std::set<TValue> > & all_nodes,
		TValue seq_i_pos,
		TValue,
		typename std::set<TValue>::iterator iter,
		TValue,
		Tag<TagExactRefinement_> const)
{
SEQAN_CHECKPOINT
	//cut already exists
	if(iter != all_nodes[seq_i_pos].end())
		return false;
	return true;
}


template<typename TSize, typename TSpec,typename TPos>
inline void
updateCutPosition(Fragment<TSize, ExactReversableFragment<TSpec> > & f, TPos & pos_j)
{
	if(f.reversed)
		++pos_j;
}
//template<typename TSize, typename TSpec,typename TPos>
//inline void
//updateCutPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f, TPos & pos_j)
//{
//	if(f.reversed)
//		++pos_j;
//}


template<typename TFrag,typename TPos>
inline void
updateCutPosition(TFrag &, TPos &)
{
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Recursive Refinement
//refine position node_i on sequence seq_i
template<typename TValue, typename TAlignmentString, typename TStringSet,typename TGraph, typename TPropertyMap,typename TSeqMap, typename TTagSpec>
inline void
_refine(TValue node_i, 
	 TValue seq_i_id, 
	 TStringSet & seqs,
	 TSeqMap & seq_map,
	 TAlignmentString & alis, 
	 String<TGraph> & gs, 
	 String<TPropertyMap> & pms, 
     String<std::set<TValue> > & all_nodes, 
	 TValue min_len,
	 Tag<TTagSpec> tag)
{
SEQAN_CHECKPOINT
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer>, Rooted>::Type TSegmentIterator;
	//find all segment matches that contain the current position (node_i)
	String<TAlignmentPointer> relevant_segments;
	TValue seq_i_pos = idToPosition(seqs,seq_i_id);
	findIntervalsExcludeTouching(gs[seq_i_pos],pms[seq_i_pos],node_i,relevant_segments);

	
	TSegmentIterator segment_it = begin(relevant_segments);
	TSegmentIterator segment_end = end(relevant_segments);
	//foreach of those segments
	while(segment_it != segment_end)
	{
		//get the sequence that node_i needs to be projected onto (seq_j)
		//and get the projected position (pos_j)
		TValue seq_j_id, node_j;
		_getOtherSequenceAndProject(alis[*segment_it],seq_map,seq_i_id,node_i,seq_j_id,node_j);
		TValue seq_j_pos = idToPosition(seqs,seq_j_id);
		updateCutPosition(alis[*segment_it],node_j);

		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j_pos].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(cutIsOk(all_nodes,seq_j_pos,node_j,iter,min_len,tag))
		{
			all_nodes[seq_j_pos].insert(node_j);
			_refine(node_j,seq_j_id,seqs,seq_map,alis,gs,pms,all_nodes,min_len,tag);
			//TODO: else //verschmelzen, abschneiden und �bergehen, erst sp�ter... 	
			//do nothing or resolve problems  
		}
	
		++segment_it;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Construct interval trees 
////////////////////////////////////////////////////////////////////////////////////////////////////
//construct intervals from allignments for each sequence (other Alignment types)
template<typename TInterval, typename TStringSet, typename TAlignmentString, typename TSeqMap>
void
_buildIntervalsForAllSequences(TAlignmentString & alis, 
							   String<String<TInterval> > & intervals, 
	   						   TStringSet & seqs,
							   TSeqMap & seq_map)
{
SEQAN_CHECKPOINT
	
	typedef typename Value<TInterval>::Type TValue;
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Iterator<TAlignmentString,Standard>::Type TAliIterator;
	TAliIterator ali_it = begin(alis,Standard());
	TAliIterator ali_end = end(alis,Standard());
	TValue ali_counter = 0;
	//foreach alignment
	while(ali_it != ali_end)
	{
		TValue seq_i_id,begin_,end_;
	
		//get the first sequence (and its begin and end) that takes part in the alignment (seq_i)
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_,end_,0);
		TValue seq_i_pos = idToPosition(seqs, seq_i_id);
		//and append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i_pos],IntervalAndCargo<TValue,TCargo>(begin_,end_,ali_counter)); 
	
		//get the second sequence (and its begin and end) that takes part in the alignment (seq_i)
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_,end_,1);
		seq_i_pos = idToPosition(seqs, seq_i_id);
		//and again append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i_pos],IntervalAndCargo<TValue,TCargo>(begin_,end_,ali_counter)); 
	
		++ali_counter;
		++ali_it;
	}
}


//get all intervals from the alignments and construct an interval tree for each sequence
template<typename TGraph, typename TPropertyMap, typename TAlignmentString, typename TSequence, typename TSetSpec, typename TValue, typename TSeqMap>
void
createTreesForAllSequences(String<TGraph> & gs, 
						   String<TPropertyMap> & pms, 
						   TAlignmentString & alis, 
						   StringSet<TSequence,TSetSpec> & seqs,
                           TSeqMap & seq_map,
						   TValue numSequences)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlignment;
//	typedef TAlignment* TCargo;
	typedef TValue TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	//std::cout <<"create interval trees...";
	clock_t start, finish1;
	double duration;
	start = clock();
	//one tree for each sequence
	resize(gs,numSequences);
	resize(pms,numSequences);
	
	//and one string of intervals for each sequence
	String<String<TInterval> > intervals;
	resize(intervals,numSequences);
	//fill intervals
	_buildIntervalsForAllSequences(alis,intervals,seqs,seq_map);
	
	TValue i = 0;
	
	while(i < numSequences)
	{
		//std::cout << (numSequences-i) <<" more ("<<length(intervals[i])<<" intervals)... "<<std::flush;
		//vllt zum speicher sparen: numSequences mal alle alis durchgehen
		//und jedes mal nur buildIntervalsForJustOneSequence(); 
		TValue center = length(seqs[i])/2; // center raus, hat hier nix zu suchen
		//create interval tree!
		createIntervalTree(gs[i],pms[i],intervals[i],center);
		
		//intervals for sequence i are not needed anymore
		clear(intervals[i]);
		++i;
	}
	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
}


//step 1 of constructing the refined alignment graph: create all the nodes
template<typename TStringSet,typename TValue,typename TAliGraph>
void
_makeRefinedGraphNodes(String<std::set<TValue> > & all_nodes,
					  TStringSet & seqs,
					  TAliGraph & ali_g)
{
SEQAN_CHECKPOINT
	typedef typename std::set<TValue>::iterator TSetIterator;
	//for each sequence look at all cut positions and create nodes between them
	for(unsigned int seq_i_pos = 0; seq_i_pos < length(seqs); ++seq_i_pos)
	{
		TValue seq_i_id = positionToId(stringSet(ali_g), seq_i_pos);
		TSetIterator it = all_nodes[seq_i_pos].begin();
		TSetIterator end_it = all_nodes[seq_i_pos].end();
		TSetIterator next_it = it;
		if(next_it != end_it)
			++next_it;
		else
			addVertex(ali_g, seq_i_id, 0, length(seqs[seq_i_pos]));
		
		//first unaligned node
		if(it != end_it && *it != 0)
			addVertex(ali_g, seq_i_id, 0, *it);
		//a new node for each interval
		while(next_it != end_it)
		{
			TValue pos_i = *it;
			addVertex(ali_g, seq_i_id, pos_i, *next_it - pos_i); 
			++it;
			++next_it;
		}
		//last unaligned node
		if(it !=end_it && *it<length(seqs[seq_i_pos]))
			addVertex(ali_g, seq_i_id, *it, (length(seqs[seq_i_pos])) - *it);
		all_nodes[seq_i_pos].clear();
	}
}


//step 2 of constructing the refined alignment graph: add all edges    
//version for exact refinement
template<typename TAlignmentString,typename TStringSet,typename TSeqMap, typename TPropertyMap,typename TScore,typename TAliGraph>
void
_makeRefinedGraphEdges(TAlignmentString & alis,
					   TPropertyMap & pm,
					  TStringSet & seqs,
				      TSeqMap & seq_map,
				      TScore & score_type,
					  TAliGraph & ali_g,
					  Tag<TagExactRefinement_> const)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Size<TAlign>::Type TValue;
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
		//get sequence, begin position and end position
		TValue seq_id,begin_pos,end_pos;
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_id,begin_pos,end_pos,(TValue)0);
		
		//get the node represents the current interval (begin_pos until next_cut_pos or end_pos)
		TVertexDescriptor act_knot = findVertex(ali_g,seq_id,begin_pos);
		TValue act_pos = begin_pos;
	
		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{
			//get other sequence and projected position
			TValue seq_j_id,pos_j;
			_getOtherSequenceAndProject(*ali_it,seq_map,seq_id,act_pos,seq_j_id,pos_j);
			//find node that contains the projected position (pos_j)
			TVertexDescriptor vd = findVertex(ali_g, seq_j_id, pos_j);
		
			SEQAN_TASSERT(fragmentBegin(ali_g,vd)==pos_j)
			typename Value<TScore>::Type score = getScore(score_type,seqs,*ali_it,act_pos,pos_j,fragmentLength(ali_g,act_knot),fragmentLength(ali_g,vd));//,fragmentLength(ali_g,vd));
	//		typename Value<TScore>::Type score = fragmentLength(ali_g,vd);
			score *= getAnnoScore(ali_g,pm,vd,act_knot,score_type);
		//this needs to be generalized (makes sense for positive scores only)
			if(score <= 0) score = 1;
			if(score > 0)
			{
				if (findEdge(ali_g, act_knot, vd) == 0) addEdge(ali_g,act_knot,vd,(TCargo)score);
				else {
					TEdgeDescriptor ed = findEdge(ali_g, act_knot, vd);
					//if((TCargo)score > getCargo(ed))
						//assignCargo(ed, score);
					assignCargo(ed, getCargo(ed)+score);
				}
			}
			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
			act_knot = findVertex(ali_g,seq_id,act_pos);
		
		}
		++ali_it;
	}
}





////////////////////////////////////////////////////////////////////////////////////////
//build refined alignment graph 
////////////////////////////////////////////////////////////////////////////////////////
//nodes are numbered ascendingly:
//seq1   0  1  2  3  4 
//seq2   5  6  7  8  9 10
//seq3  11 12 13 14 15 
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TAliGraph,typename TSeqMap,typename TTagSpec>
void
_makeAlignmentGraphFromRefinedSegments(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   TAliGraph & ali_g,
			   	   Tag<TTagSpec> const tag, 
				   bool)
{
SEQAN_CHECKPOINT
	//std::cout << "making refined alignment graph...";
	//clock_t start, finish1;
	//double duration;
	//start = clock();
	
	//make nodes (same function for inexact and exact refinement)
	_makeRefinedGraphNodes(all_nodes,seqs,ali_g);

	bool pm = false;
	//add edges (different functions depending on exact/inexact refinement)
	_makeRefinedGraphEdges(alis,pm,seqs,seq_map,score_type,ali_g,tag);
	
	//std::cout << "check\n";
	//finish1 = clock();
	//duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
}


      
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TAliGraph,typename TSeqMap,typename TAnnoString,typename TTagSpec>
void
_makeAlignmentGraphFromRefinedSegments(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   TAliGraph & ali_g,
			   	   Tag<TTagSpec> const tag,
				   TAnnoString & annotation)
{
SEQAN_CHECKPOINT
	//std::cout << "making refined alignment graph...";
	//clock_t start, finish1;
	//double duration;
	//start = clock();
	
	//make nodes (same function for inexact and exact refinement)
	_makeRefinedGraphNodes(all_nodes,seqs,ali_g);

	//add annotation to nodes
	typedef typename Value<TAnnoString>::Type TAnnotation;
	//typedef typename Value<TAnnotation>::Type TLabel;
	typedef char TLabel;
	String<String<TLabel> > pm;
	_addNodeAnnotation(seqs,seq_map,annotation,pm,ali_g,tag);

	//add edges (different functions depending on exact/inexact refinement)
	_makeRefinedGraphEdges(alis,pm,seqs,seq_map,score_type,ali_g,tag);
	
	//std::cout << "check\n";
	//finish1 = clock();
	//duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
}





////////////////////////////////////////////////////////////////////////////////////////
//The big matchRefinement function that does everything: build interval trees, do the 
//refinement and construct a refined alignment graph
////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlignmentString, typename TAnnotation, typename TOutGraph, typename TSequence, typename TSetSpec, typename TScore,typename TTagSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TScore & score_type,
				TOutGraph & ali_graph,
				typename Size<typename Value<TAlignmentString>::Type>::Type min_fragment_len,
				TAnnotation & annotation,
				Tag<TTagSpec> const tag)
{
SEQAN_CHECKPOINT
	////////////////////////////////////////////////////////////////
	//typedefs
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString, Rooted>::Type TAliIterator;
	typedef typename Size<TAlign>::Type TValue;
	typedef TValue TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TCargo> TList;
	typedef typename std::set<TValue>::iterator TSetIterator;
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer>, Rooted>::Type TSegmentIterator;
	
	////////////////////////////////////////////////////////////////
	TValue numSequences = length(seq);
	//weird ID --> good ID map
	std::map<const void * ,int> seq_map;
	for(int i = 0; i < (int) numSequences; ++i)
		seq_map[id(seq[i])] = i;
	////////////////////////////////////////////////////////////////
	//build interval trees
	String<TGraph> gs;
	String<TPropertyMap> pms;
	createTreesForAllSequences(gs, pms, alis, seq, seq_map, numSequences);
	
	////////////////////////////////////////////////////////////////
	//do refinement
	//std::cout <<"refining..."<<std::flush;
	clock_t start, finish1;
	double duration;
	start = clock();
	
	//all_nodes = set of all cut positions
	String<std::set<TValue> > all_nodes;
	resize(all_nodes,numSequences);

	//all_nodes that need to be processed set of all cut positions
	String<std::set<TValue> > all_node_queues;
	resize(all_node_queues,numSequences);

	//call function _refine for each startknoten
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);
	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//for each of the two sequences
		for(TValue i = 0; i < 2; ++i)
		{
			TValue seq_i_id,begin_i,end_i;
			_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_i,end_i,i);
			TValue seq_i_pos = idToPosition(seq,seq_i_id);
			
			all_node_queues[seq_i_pos].insert(begin_i);
			all_node_queues[seq_i_pos].insert(end_i);
		}	
		++ali_it;
	}

	TSetIterator queueIt;
	bool done = false;
	while(!done)
	{
		for(unsigned seq_i_pos = 0; seq_i_pos < numSequences; ++seq_i_pos)
		{
			queueIt = all_node_queues[seq_i_pos].begin();
			while (queueIt != all_node_queues[seq_i_pos].end())
			{
				TValue node_i = *queueIt;
				TSetIterator iter = all_nodes[seq_i_pos].find(node_i);		
				if(iter == all_nodes[seq_i_pos].end())
				{
					TValue seq_i_id = positionToId(seq, seq_i_pos);
					all_nodes[seq_i_pos].insert(node_i);
					String<TAlignmentPointer> relevant_segments;
					findIntervalsExcludeTouching(gs[seq_i_pos],pms[seq_i_pos],node_i,relevant_segments);
					
					TSegmentIterator segment_it = begin(relevant_segments);
					TSegmentIterator segment_end = end(relevant_segments);
					//foreach of those segments
					while(segment_it != segment_end)
					{
						//get the sequence that node_i needs to be projected onto (seq_j)
						//and get the projected position (pos_j)
						TValue seq_j_id, node_j;
						_getOtherSequenceAndProject(alis[*segment_it],seq_map,seq_i_id,node_i,seq_j_id,node_j);
						TValue seq_j_pos = idToPosition(seq,seq_j_id);
						updateCutPosition(alis[*segment_it],node_j);

						typename std::set<TValue>::iterator iter_j;
						iter_j = all_nodes[seq_j_pos].find(node_j);
						
						//if node does not exist yet ---> insert and continue cutting
						if(iter_j == all_nodes[seq_j_pos].end())
						{
							all_node_queues[seq_j_pos].insert(node_j);
						}
					
						++segment_it;
					}
						
				}
				++queueIt;
			}
			all_node_queues[seq_i_pos].clear();
		}
		unsigned i;
		for(i = 0; i < numSequences; ++i)
		{
			queueIt = all_node_queues[i].begin();
			if (queueIt != all_node_queues[i].end())
				break;
		}
		if(i==numSequences)
			done=true;
	}
	_addAnnotationCuts(all_nodes,alis,gs,pms,seq,seq_map,annotation,min_fragment_len,tag);

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
	//for(int seq_i = 0; seq_i < length(seq); ++seq_i)
	//{
	//	typename std::set<TValue>::iterator it = all_nodes[seq_i].begin();
	//	typename std::set<TValue>::iterator end_it = all_nodes[seq_i].end();
	//
	//	while(it != end_it)
	//	{
	//		std::cout << *it << ",";
	//		++it;
	//	}
	//	std::cout << "\n";
	//}
	//std::cout <<"building tree..."<<std::flush;
	
	////////////////////////////////////////////////////////////////
	//build refined alignment graph
	_makeAlignmentGraphFromRefinedSegments(all_nodes,alis,score_type,seq,seq_map,ali_graph,tag,annotation);
}


///////WRAPPERS

/**
.Function.matchRefinement:
..signature:matchRefinement(matches,stringSet,scoringScheme,refinedGraph)
..param.matches:The set of matches.
..param.scoringScheme:The scoring scheme used to score the refined matches (scores are attached to 
edges in the refined Alignment Graph).
...remarks:If no scoring scheme is given, all edges get weight 1.
...type:Class.Score
*/
//exact refinement, score type given
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec,typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				Score<TScoreValue,TScoreSpec> & score_type,
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT
	//min_fragment_len = 1   ==> Exact cutting
	bool anno = false;
	matchRefinement(alis,seq,score_type,ali_graph,1,anno,ExactRefinement());
}



/**
.Function.matchRefinement:
..cat:Alignments
..summary:Refines (i.e. cuts into smaller parts) a set of pairwise segment 
matches in such a way that none of the segments partly overlap. They are either 
identical (fully overlapping) or non-overlapping.
..signature:matchRefinement(matches,stringSet,refinedGraph)
..param.matches:The set of matches.
..param.stringSet:The StringSet containing the sequences which the matches lie on.
...type:Class.StringSet
..param.refinedGraph:The resulting refined set of matches stored in a graph.
...type:Spec.Alignment Graph
*/
//exact refinement, score type not given
template<typename TFragmentString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TFragmentString & matches,
				StringSet<TSequence, TSetSpec> & strSet, 
				TOutGraph & ali_graph)
{
	SEQAN_CHECKPOINT
	typename Cargo<TOutGraph>::Type fake_score = 1;
	bool anno = false;
	matchRefinement(matches,strSet,fake_score,ali_graph,1,anno,ExactRefinement());
}


}
#endif //#ifndef SEQAN_HEADER_...
