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


#ifndef SEQAN_HEADER_GRAPH_REFINE_ANNOTATION_H
#define SEQAN_HEADER_GRAPH_REFINE_ANNOTATION_H

//SEQAN_NO_DDDOC: do not generate documentation for this file


namespace SEQAN_NAMESPACE_MAIN
{


template<typename TSequence, typename TValue, typename TSpec = Simple>
class Annotation;



/**
.Class.Annotation:
..cat:Sequences
..summary:Class for annotating sequences. 
..signature:Annotation<TSequence, TLabel, TSpec>  
..param.TSequence:The sequence that annotation is available for.
..param.TLabel:The label type (e.g. int or String<char>)
..param.TSpec:The specializing type.
...default:Simple
..include:refinement.h
*/
template<typename TSequence,typename TValue>
class Annotation<TSequence,TValue,Simple>{

public:
	typedef typename Id<TSequence>::Type TId_;
	typedef typename Position<TSequence>::Type TPos_;
	typedef typename Size<TSequence>::Type TSize_;

	TId_ data_seq_id;
	TPos_ data_begin;
	TSize_ data_length;
	//String<char> data_label;
	TValue data_label;

	Annotation()
	{
	}

	Annotation(TId_ seqId, TPos_ begin, TSize_ len, TValue label) :
			data_seq_id(seqId),
			data_begin(begin), 
			data_length(len),
			data_label(label)
	{
	}

	~Annotation()
	{
	}

};



template<typename TSequence,typename TValue,typename TSpec>	
typename Id<TSequence>::Type&
sequenceId(Annotation<TSequence,TValue,TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_seq_id;
}

template<typename TSequence,typename TValue,typename TSpec>	
typename Position<TSequence>::Type&
fragmentBegin(Annotation<TSequence,TValue,TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_begin;
}

template<typename TSequence,typename TValue,typename TSpec>	
typename Size<TSequence>::Type&
fragmentLength(Annotation<TSequence,TValue,TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_length;
}

template<typename TSequence,typename TValue,typename TSpec>	
TValue
label(Annotation<TSequence,TValue,TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_label;
}


template<typename TSequence,typename TValue,typename TSpec>	
struct Value<Annotation<TSequence,TValue,TSpec> >
{
	typedef TValue Type;
};


template<typename TValue, typename TAliString, typename TGraph, typename TPropertyMap, typename TStringSet, typename TMap, typename TTagSpec>
inline void
_addAnnotationCuts(String<std::set<TValue> > &,
				   TAliString &, 
				   String<TGraph> &, 
				   String<TPropertyMap> &, 
				   TStringSet &,
				   TMap &,
				   bool,
				   TValue,
				   Tag<TTagSpec>)
{
SEQAN_CHECKPOINT
	return;
}




//refine positions where annotation changes (borders of annotated stretches)
template<typename TValue, typename TAliString, typename TGraph, typename TPropertyMap, typename TStringSet, typename TMap, typename TAnnoString, typename TTagSpec>
inline void
_addAnnotationCuts(String<std::set<TValue> > & all_nodes,
				   TAliString & alis, 
				   String<TGraph> & gs, 
				   String<TPropertyMap> & pms, 
				   TStringSet & seq,
				   TMap & seq_map,
				   TAnnoString & annotation,
				   TValue min_fragment_len,
				   Tag<TTagSpec> tag)
{
SEQAN_CHECKPOINT

	//typedef typename Value<TAnnoString>::Type TAnnotation;
	typedef typename Iterator<TAnnoString,Standard>::Type TAnnoIter;
	typedef typename std::set<TValue>::iterator TSetIterator;
	//call function _refine for each position that annotation is given for (the borders of annotated stretches)
	TAnnoIter anno_it = begin(annotation,Standard());
	TAnnoIter anno_end = end(annotation,Standard());
	//for each annotated stretch
	while(anno_it != anno_end)
	{
		
		TValue seq_i_id = sequenceId(*anno_it);
		TValue begin_i = fragmentBegin(*anno_it);
		TValue end_i = begin_i + fragmentLength(*anno_it);
		TValue seq_i_pos = idToPosition(seq,seq_i_id);
		
		//refine begin
		TSetIterator iter = all_nodes[seq_i_pos].find(begin_i);		
		if(_cutIsValid(all_nodes,seq_i_pos,begin_i,iter,min_fragment_len,tag))
		{
			all_nodes[seq_i_pos].insert(begin_i);
			_refine(begin_i, seq_i_id, seq, seq_map, alis, gs,pms,all_nodes,min_fragment_len,tag);//TStop());
		}
		//and end position
		iter = all_nodes[seq_i_pos].find(end_i);		
		if(_cutIsValid(all_nodes,seq_i_pos,end_i,iter,min_fragment_len,tag))
		{
			all_nodes[seq_i_pos].insert(end_i);
			_refine(end_i, seq_i_id, seq, seq_map, alis, gs,pms,all_nodes,min_fragment_len,tag);//TStop());
		}
		++anno_it;
	}
	      	
}
	      


template<typename TPropertyMap, typename TStringSet, typename TMap, typename TAnnoString, typename TAliGraph,typename TTagSpec>
inline void
_addNodeAnnotation(TStringSet &,
				   TMap &,
				   TAnnoString & annotation,
				   TPropertyMap & pm,
				   TAliGraph & ali_g,
				   Tag<TTagSpec>)
{
SEQAN_CHECKPOINT

	resizeVertexMap(ali_g, pm);

	typedef typename Value<TAnnoString>::Type TAnnotation;
	typedef typename Id<TAnnotation>::Type TId;
	typedef typename Position<TAnnotation>::Type TPos;
	typedef typename Value<TAnnotation>::Type TLabel;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	
	typedef typename Iterator<TAnnoString,Standard>::Type TAnnoIter;
	TAnnoIter anno_it = begin(annotation,Standard());
	TAnnoIter anno_end = end(annotation,Standard());

	//for each annotated stretch
	while(anno_it != anno_end)
	{
		TLabel label_ = label(*anno_it);
		TId seq_id = sequenceId(*anno_it);
		TPos act_pos = fragmentBegin(*anno_it);
		TPos end_pos = act_pos + fragmentLength(*anno_it);

		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{
			//get the node represents the current interval (begin_pos until next_cut_pos or end_pos)
			TVertexDescriptor act_knot = findVertex(ali_g,seq_id,act_pos);

			String<TLabel> property = getProperty(pm, act_knot);
			appendValue(property,label_);
			assignProperty(pm, act_knot, property);
			SEQAN_ASSERT_TRUE(fragmentBegin(ali_g,act_knot)==act_pos);

			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
		}
		++anno_it;

	}



}

template<typename TAliGraph,typename TScore, typename TPropertyMap>
typename Value<TScore>::Type 
_getRefinedAnnoScore(TAliGraph &,
			 TPropertyMap & pm,
			 typename VertexDescriptor<TAliGraph>::Type vd1,
			 typename VertexDescriptor<TAliGraph>::Type vd2,
			 TScore &)
{
SEQAN_CHECKPOINT
	typedef typename Value<TPropertyMap>::Type TProperty;
	//typedef typename Value<TProperty>::Type TChar;
	typedef typename Iterator<TProperty,Standard>::Type TIterator;

	TIterator prop1_it = begin(property(pm,vd1),Standard());
	TIterator prop1_end = end(property(pm,vd1),Standard());
	while(prop1_it != prop1_end)
	{
		TIterator prop2_it = begin(property(pm,vd2),Standard());
		TIterator prop2_end = end(property(pm,vd2),Standard());
		while(prop2_it != prop2_end)
		{
			if(*prop2_it==*prop1_it)
				return 2;/*scoreMatch(score_type);*/
			++prop2_it;
		}
		++prop1_it;
	}
	return 1;

}


template<typename TAliGraph,typename TScore>
typename Value<TScore>::Type 
_getRefinedAnnoScore(TAliGraph &,
			 bool,
			 typename VertexDescriptor<TAliGraph>::Type,
			 typename VertexDescriptor<TAliGraph>::Type,
			 TScore &)
{
SEQAN_CHECKPOINT
	return (typename Value<TScore>::Type) 1;
}




/**
.Function.matchRefinement:
..signature:matchRefinement(matches,annotation,stringSet,scoringScheme,refinedGraph)
..param.annotation:Sequence annotation data. 
...remarks: Additional semgent match subdivisions will be made at sequence positions at which the annotation label changes.
...type:Class.Annotation
..include:seqan/refinement.h
*/
//annotation given,exact refinement, score type given
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec,typename TAnnoString,typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				TAnnoString & anno,
				StringSet<TSequence, TSetSpec> & seq, 
				Score<TScoreValue,TScoreSpec> & score_type,
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT
	//min_fragment_len = 1   ==> Exact cutting
	matchRefinement(alis,seq,score_type,ali_graph,1,anno,ExactRefinement());
}


/**
.Function.matchRefinement:
..signature:matchRefinement(matches,annotation,stringSet,scoringScheme,refinedGraph,minFragmentLen)
..include:seqan/refinement.h
*/
//annotation given,score type given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec,typename TAnnoString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				TAnnoString & anno,
				StringSet<TSequence, TSetSpec> & seq, 
				Score<TScoreValue,TScoreSpec> & score_type,
				TOutGraph & ali_graph,
				unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
	if(min_frag_len > 1)
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,InexactRefinement());
	else
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,ExactRefinement());
}



/**
.Function.matchRefinement:
..signature:matchRefinement(matches,annotation,stringSet,refinedGraph,minFragmentLen)
..include:seqan/refinement.h
*/
//annotation given,score type not given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TOutGraph, typename TAnnoString, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				TAnnoString & anno,
				StringSet<TSequence, TSetSpec> & seq, 
				TOutGraph & ali_graph,
				unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
//	Score<int,FakeScore > fake_score;
	typename Cargo<TOutGraph>::Type fake_score = 1;
	if(min_frag_len > 1)
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,InexactRefinement());
	else
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,ExactRefinement());
}
	


/**
.Function.matchRefinement:
..signature:matchRefinement(matches,annotation,stringSet,refinedGraph)
..include:seqan/refinement.h
*/
//annotation given,exact refinement, score type not given
template<typename TAlignmentString,typename TAnnoString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				TAnnoString & anno,
				StringSet<TSequence, TSetSpec> & seq, 
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT
//	Score<int,FakeScore > fake_score;
	typename Cargo<TOutGraph>::Type fake_score = 1;
	matchRefinement(alis,seq,fake_score,ali_graph,1,anno,ExactRefinement());
}

}
#endif //#ifndef SEQAN_HEADER_...
