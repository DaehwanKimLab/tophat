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
  $Id: graph_impl_automaton.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H
#define SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Automaton
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TEdge, typename TAlphabet>
class AutomatonEdgeArray {
public:
	TEdge data_edge[ValueSize<TAlphabet>::VALUE];

	AutomatonEdgeArray() 
	{
		typedef typename VertexDescriptor<TEdge>::Type TVertexDescriptor;
		typedef typename Size<TAlphabet>::Type TSize;
		TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
		for(TSize i=0;i < (TSize) ValueSize<TAlphabet>::VALUE;++i) 
			assignTarget(&data_edge[i], nilVal);
	}
};

/**
.Spec.Automaton:
..cat:Graph
..general:Class.Graph
..summary:An Automaton has directed edges, labeled with input symbols, and a distinct start state, called root. The
input symbols require the use of a third parameter: The alphabet of the input symbols.
..description:
...image:automatonGraph|An automaton, where $0$ is the start state.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
*/
template<typename TAlphabet, typename TCargo, typename TSpec>
class Graph<Automaton<TAlphabet, TCargo, TSpec> > 
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
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline String<AutomatonEdgeArray<typename EdgeType<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type, TAlphabet> >&
_getVertexString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	return const_cast<String<AutomatonEdgeArray<TEdge, TAlphabet> >&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
	return const_cast<TVertexIdManager&>(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
	return const_cast<TEdgeIdManager&>(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		   Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest,
		   bool transpose)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Cargo<TEdge>::Type TCargoType;
	typedef typename Size<TAlphabet>::Type TSize;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Standard>::Type TIterConst;
	
	clear(dest);
	resize(dest.data_vertex, length(_getVertexString(source)));
	dest.data_root = source.data_root;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TIterConst it = begin(source.data_vertex, Standard());
	TIterConst itEnd = end(source.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor sourceVertex = pos;
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor const edSource = (TEdgeDescriptor) &source.data_vertex[sourceVertex].data_edge[i];
			TVertexDescriptor targetVertex = getTarget(edSource);
			if (targetVertex == nilVal) continue;
			TEdgeDescriptor edTarget;
			if (!transpose) {
				edTarget = &dest.data_vertex[sourceVertex].data_edge[i];
				assignTarget(edTarget, targetVertex);
			} else {
				edTarget = &dest.data_vertex[targetVertex].data_edge[i];
				assignTarget(edTarget, sourceVertex);
			}
			_assignId(edTarget, _getId(edSource));
			assignCargo(edTarget, getCargo(edSource));
		}
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		   Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest)
{
	_copyGraph(source,dest,false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		  Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	_copyGraph(source, dest, true);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	Graph<Automaton<TAlphabet, TCargo, TSpec> > dest;
	_copyGraph(g, dest, true);
	g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
numEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
numVertices(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline bool 
empty(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	clear(g.data_vertex);
	releaseAll(g.data_id_managerE);
	resize(g.data_vertex, getIdUpperBound(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
clear(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
outDegree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		  TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE; ++i) 
		if ( (TVertexDescriptor) getTarget(&g.data_vertex[vertex].data_edge[i])!=nilVal) ++count;
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
inDegree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		 TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd; ++it, ++pos) {
		if (idInUse(g.data_id_managerV, pos)) {
			for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE;++i) 
				if ( (TVertexDescriptor) getTarget(&(*it).data_edge[i]) == vertex) ++count;			
		}
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
degree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g, 
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (inDegree(g,vertex)+outDegree(g,vertex));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;

	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) appendValue(g.data_vertex, AutomatonEdgeArray<TEdge, TAlphabet>()); 
	else g.data_vertex[vd] =  AutomatonEdgeArray<TEdge, TAlphabet>();
	return vd;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TEdgeDescriptor e = findEdge(g, source, (TAlphabet) label);
	_assignId(e, obtainId(g.data_id_managerE));
	assignTarget(e, target);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
				  TVertexDescriptor const source, 
				  TVertexDescriptor const target,
				  TLabel const label,
				  TEdgeCargo const cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TEdgeDescriptor e = addEdge(g,source,target, (TAlphabet) label);
	assignCargo(e,cargo);
	TId id = obtainId(g.data_id_managerE);
	_assignId(e, id);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
	       TVertexDescriptor const source,
	       TVertexDescriptor const target,
		   TLabel const label)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	removeEdge(g, &g.data_vertex[source].data_edge[ordValue((TAlphabet) label)]);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, _getId(edge)) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	releaseId(g.data_id_managerE, _getId(edge));
	assignTarget(edge, getNil<TVertexDescriptor>());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			   TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE;++i) {
		TEdgeDescriptor ed = &g.data_vertex[vertex].data_edge[i];
		if ( (TVertexDescriptor) getTarget(ed) == nilVal) continue;
		assignTarget(ed, nilVal);
		releaseId(g.data_id_managerE, _getId(ed));
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			  TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> >, Standard>::Type TIter;
	TIter it = begin(g.data_vertex, Standard());
	TIter itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		if (idInUse(g.data_id_managerV, pos)) {
			for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE;++i) {
				TEdgeDescriptor ed = &(value(it)).data_edge[i];
				if (( (TVertexDescriptor) getTarget(ed) == nilVal) ||  ( (TVertexDescriptor) getTarget(ed) != vertex)) continue;
				assignTarget(ed, nilVal);
				releaseId(g.data_id_managerE, _getId(ed));	
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
targetVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return (getTarget(edge));
}

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
sourceVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
		
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> >, Standard>::Type TIter;
	TIter it = begin(g.data_vertex, Standard());
	TIter itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd; ++it, ++pos) {
		if (idInUse(g.data_id_managerV, pos)) {
			for(TSize i=0;i<table_length;++i) {
				TEdgeDescriptor ed = &(*it).data_edge[i];
				if (getTarget(ed) == nilVal) continue;
				if (ed==edge) return pos;
			}
		}
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Size<TMatrix>::Type TMatrixSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TMatrixSize len = getIdUpperBound(g.data_id_managerV);
	fill(mat, len*len, 0);
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		if (!idInUse(g.data_id_managerV, pos)) continue;
		for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE;++i) {
			if (((*it).data_edge[i].data_target!=nilVal))
			{
				TVertexDescriptor const source = pos;
				TVertexDescriptor const target = (*it).data_edge[i].data_target;
				++mat[source*len+target];
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
findEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		 TVertexDescriptor const v,
		 TLabel const c)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	return &g.data_vertex[v].data_edge[ordValue((TAlphabet) c)];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
findEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		 TVertexDescriptor const v,
		 TLabel const c)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	TGraph* graph = const_cast<TGraph*>(&g);
	return findEdge(*graph, v, c);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

	_streamWrite(target,"Automaton - State: (Input / NextState)\n");
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		if (!idInUse(g.data_id_managerV, pos)) continue;
		TVertexDescriptor sourceVertex = pos;
		_streamPutInt(target, sourceVertex);
		_streamWrite(target,": ");
		for(TSize i=0;i< (TSize) ValueSize<TAlphabet>::VALUE;++i) {
			_streamPut(target, ' ');
			_streamPut(target, '(');
			_streamPut(target, TAlphabet(i));
			_streamPut(target, ' ');
			_streamPut(target, '/');
			_streamPut(target, ' ');
			if (g.data_vertex[sourceVertex].data_edge[i].data_target ==  nilVal) _streamWrite(target,"nil");
			else _streamPutInt(target, g.data_vertex[sourceVertex].data_edge[i].data_target);
			_streamPut(target, ')');
			_streamPut(target, ' ');
		}
		_streamPut(target, '\n');
	}
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.createRoot:
..cat:Graph
..summary:Creates the root in a tree or an automaton.
..signature:createRoot(g)
..param.g:A tree or an automaton.
...type:Spec.Automaton
...type:Spec.Tree
..returns:void
*/

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
createRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	g.data_root = addVertex(g);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignRoot:
..cat:Graph
..summary:Assigns a new root vertex to the graph.
..signature:assignRoot(g, vertex)
..param.g:A tree or automaton.
...type:Spec.Automaton
...type:Spec.Tree
..param.vertex:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.getRoot
..see:Function.root
..see:Function.Graph#isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		   TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	g.data_root = vertex;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.root:
..cat:Graph
..summary:Gets a reference to the root of the tree.
..signature:root(g)
..param.g:A tree or an automaton.
...type:Spec.Automaton
...type:Spec.Tree
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.getRoot
..see:Function.assignRoot
..see:Function.Graph#isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type&
root(Graph<Automaton<TAlphabet, TCargo, TSpec> > & g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getRoot:
..cat:Graph
..summary:Get method for the root of a tree or an automaton.
..signature:getRoot(g)
..param.g:An automaton or a tree.
...type:Spec.Automaton
...type:Spec.Tree
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.assignRoot
..see:Function.root
..see:Function.Graph#isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
getRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return g.data_root;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.Graph#isRoot:
..cat:Graph
..summary:Tests whether a given vertex is the root or not.
..signature:isRoot(g, v)
..param.g:An automaton or a tree.
...type:Spec.Automaton
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:True if vertex is the root.
..see:Function.getRoot
..see:Function.assignRoot
..see:Function.root
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.data_root == v);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.getSuccessor:
..cat:Graph
..summary:Gets the successor for a given vertex and an edge label.
For an automaton a single character is required whereas for a word graph getSuccessor takes a string.
..signature:getSuccessor(g, v, c)
..signature:getSuccessor(g, v, str)
..param.g:An automaton or a word graph.
...type:Spec.Automaton
...type:Spec.Word Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.c:A label.
...type:Metafunction.Alphabet
..param.str:A string.
...type:Class.String
..returns:A vertex descriptor or nil if successor is not defined.
...type:Metafunction.VertexDescriptor
..see:Function.parseString
..see:Function.getNil
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			 TVertexDescriptor vertex,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	return getTarget(findEdge(g, vertex, c));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.parseString:
..cat:Graph
..summary:Parses a string one character at a time and moves accordingly in the automaton.
..signature:parseString(g, v, beginIt, endIt)
..signature:parseString(g, v, str)
..param.g:An automaton or a word graph.
...type:Spec.Automaton
...type:Spec.Word Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.beginIt:An iterator to the beginning of the string. (Input and Output)
...type:Metafunction.Iterator.
..param.endIt:An iterator to the end of the string.
...type:Metafunction.Iterator.
..param.str:A string.
...type:Class.String
..returns:The vertex descriptor of the state that was reached after parsing.
...type:Metafunction.VertexDescriptor
..returns.param.beginIt:Iterator to the first character that could not be parsed, 
	or $endIt$ if the complete string was parsed.
..remarks:The parsing stops before @Function.getSuccessor@ reaches the $nil$ state, 
	or if the the complete sequences is read.
..see:Function.getSuccessor
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TIterator & beginIt,
			TIterator const & endIt)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	while (beginIt!=endIt) {
		TVertexDescriptor tmp = getSuccessor(g,succ,*beginIt);
		if (tmp == nilVal) break;
		succ = tmp;
		++beginIt;
	}
	return succ;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TIterator const & beginIt,
			TIterator const & endIt)
{
	SEQAN_CHECKPOINT
	TIterator beginIt2 = beginIt;
	return parseString(g, vertex, beginIt2, endIt);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo,  typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TCharacters const& chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,begin(chars),end(chars));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,chars,chars+length(chars));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.canParseString:
..cat:Graph
..summary:Test whether an automaton can parse a string completely.
..signature:canParseString(g [, v], str)
..param.g:An automaton or a word graph.
...type:Spec.Automaton
...type:Spec.Word Graph
..param.v:A vertex descriptor. (optional)
...type:Metafunction.VertexDescriptor
...default:The @Function.root@ of $g$.
..param.str:A string.
...type:Class.String
..returns:$true$ if $g$ parses $str$ starting at $v$ completely, $false$ otherwise.
..see:Function.getSuccessor
..see:Function.parseString
*/

template<typename TAlphabet, typename TCargo,  typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline bool 
canParseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const & g,
			   TVertexDescriptor const vertex,
			   TCharacters const & chars)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TCharacters const, Standard>::Type TIterator;
	TIterator it = begin(chars, Standard());
	TIterator it_end = end(chars, Standard());
	parseString(g, vertex, it, it_end);
	return (it == it_end);
}
template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline bool
canParseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const & g,
			   TVertexDescriptor const vertex,
			   TCharacters const * chars)
{
	SEQAN_CHECKPOINT
	typedef TCharacters const * TIterator;
	TIterator it = begin(chars, Standard());
	TIterator it_end = end(chars, Standard());
	parseString(g, vertex, it, it_end);
	return (it == it_end);
}

template<typename TAlphabet, typename TCargo,  typename TSpec, typename TCharacters>
inline bool 
canParseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > & g,
			   TCharacters const & chars)
{
	SEQAN_CHECKPOINT
	return canParseString(g, root(g), chars);
}

template<typename TAlphabet, typename TCargo,  typename TSpec, typename TCharacters>
inline bool 
canParseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > & g,
			   TCharacters const * chars)
{
	SEQAN_CHECKPOINT
	return canParseString(g, root(g), chars);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
