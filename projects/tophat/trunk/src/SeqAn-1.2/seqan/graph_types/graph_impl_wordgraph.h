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
  $Id: graph_impl_wordgraph.h 1757 2008-02-27 16:26:20Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H
#define SEQAN_HEADER_GRAPH_IMPL_WORDGRAPH_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSpec = Default>
struct WordGraph;

//////////////////////////////////////////////////////////////////////////////
// Graph - WordGraph
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Word Graph:
..cat:Graph
..general:Spec.Automaton
..summary:A special automaton that stores words instead of single characters along its edges.
..signature:Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >
..param.TAlphabet:The alphabet type that is used for the words.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the value type of the words.
...default:$char$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
*/
template<typename TAlphabet, typename TSpec>
class Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > 
{
	public:
		typedef typename VertexIdHandler<Graph>::Type TVertexIdManager;
		typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename EdgeType<Graph>::Type TEdge;	

		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		TVertexIdManager data_id_managerV;
		TEdgeIdManager data_id_managerE;
		TVertexDescriptor data_root;


//____________________________________________________________________________


		Graph() : data_root(0) {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			_copyGraph(_other, *this);
			return *this;
		}
};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		String<TAlphabet> const& label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	
	TAlphabet firstChar = getValue(label, 0);
	TEdgeDescriptor e = findEdge(g, source, firstChar);
	TId id = obtainId(g.data_id_managerE);
	_assignId(e, id);
	assignTarget(e, target);
	String<TAlphabet> suf(suffix(label,1));
	assignCargo(e, suf);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TChars>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TChars const* chars) 
{
	SEQAN_CHECKPOINT
	return addEdge(g,source,target,String<TAlphabet>(chars));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
addEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label,
		TEdgeCargo const cargo)
{
	// No additional cargo allowed. Cargo is used for the words in the graph.
	// Use external property map.
	SEQAN_ASSERT(false)
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		String<TAlphabet> const& label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	TAlphabet firstChar = getValue(label, 0);
	removeEdge(g, findEdge(g,source, firstChar));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
	  TIDString const &,
	  Raw)
{
	typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

	_streamWrite(target,"WordGraph - Directed:\n");
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Rooted>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
		for(TSize i=0;i<table_length;++i) {
			TEdge const* ed = &g.data_vertex[sourceVertex].data_edge[i];
			if (getTarget(ed) ==  nilVal) continue;
			_streamPutInt(target, sourceVertex);
			_streamWrite(target,"->");
			_streamPutInt(target, getTarget(ed));
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "Label: ");
			_streamPut(target, TAlphabet(i));
			_streamWrite(target, getCargo(ed));
			_streamPut(target, '\n');
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
getSuccessor(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
			 TVertexDescriptor vertex,
			 TCharacters const& chars)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TAlphabet>::Type TSize;
	TEdgeStump* ed = findEdge(g, vertex, getValue(chars, 0));
	if (getCargo(ed) == suffix(chars, 1)) {
		return getTarget(ed);
	} else {
		return getNil<TVertexDescriptor>();
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
getSuccessor(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
			 TVertexDescriptor vertex,
			 TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return getSuccessor(g,vertex,String<TAlphabet>(chars));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > >::Type 
parseString(Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > const& g,
			TVertexDescriptor const vertex,
			TIterator beginIt,
			TIterator endIt)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<TSpec> > > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	while (beginIt!=endIt) {
		String<TAlphabet> label(*beginIt);
		TSize range = 1;
		TVertexDescriptor tmp = getSuccessor(g,succ,label);
		while ((tmp == nilVal) &&
				(beginIt+range != endIt))
		{
			appendValue(label, *(beginIt + range));
			tmp = getSuccessor(g,succ,label);
			++range;
		}
		if (tmp == nilVal) break;
		succ = tmp;
		beginIt = beginIt+range;
	}
	return succ;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
