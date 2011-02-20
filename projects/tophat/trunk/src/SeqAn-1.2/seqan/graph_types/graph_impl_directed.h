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
  $Id: graph_impl_directed.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_DIRECTED_H
#define SEQAN_HEADER_GRAPH_IMPL_DIRECTED_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Directed
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Directed graph:
..cat:Graph
..general:Class.Graph
..summary:A directed graph that stores the edges in an adjacency list.
..description:
...image:directedGraph|A directed graph.
..signature:Graph<Directed<TCargo, TSpec> >
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of a directed graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
*/
template<typename TCargo, typename TSpec>
class Graph<Directed<TCargo, TSpec> > 
{
	public:
		typedef typename VertexIdHandler<Graph>::Type TVertexIdManager;
		typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager;
		typedef typename EdgeType<Graph>::Type TEdgeStump;	
		typedef Allocator<SinglePool<sizeof(TEdgeStump)> > TAllocator;
		
		String<TEdgeStump*> data_vertex;			// Pointers to EdgeStump lists
		TVertexIdManager data_id_managerV;
		TEdgeIdManager data_id_managerE;		
		TAllocator data_allocator;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}

		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) :
			data_allocator(_other.data_allocator)
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);		
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			clear(*this);
			data_allocator = _other.data_allocator;
			_copyGraph(_other, *this);
			return *this;
		}
};


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Directed<TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Directed<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	return const_cast<String<TEdgeStump*>&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Directed<TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Directed<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
	return const_cast<TVertexIdManager&>(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Directed<TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Directed<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
	return const_cast<TEdgeIdManager&>(g.data_id_managerE);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Directed<TCargo, TSpec> > const& source,
		   Graph<Directed<TCargo, TSpec> >& dest,
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(_getVertexString(source)));
	TIter itInit = begin(dest.data_vertex, Standard());
	TIter itInitEnd = end(dest.data_vertex, Standard());
	for(;itInit != itInitEnd; ++itInit) *itInit = (TEdgeStump*) 0;
	TIterConst it = begin(source.data_vertex, Standard());
	TIterConst itEnd = end(source.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it != itEnd; ++it, ++pos) {
		TEdgeStump* current = *it;
		TVertexDescriptor sourceVertex = pos;
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor targetVertex = current->data_target;
			// Create missing vertices
			if (sourceVertex>targetVertex) _createVertices(dest,sourceVertex);
			else _createVertices(dest,targetVertex);
			// Add edge
			TEdgeDescriptor e;
			if (!transpose) e = addEdge(dest, sourceVertex, targetVertex);
			else e = addEdge(dest, targetVertex, sourceVertex);
			_assignId(e, _getId(current));
			assignCargo(e, getCargo(current));
			current = getNextT(current);
		}
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Directed<TCargo, TSpec> > const& source,
		   Graph<Directed<TCargo, TSpec> >& dest) 
{
	_copyGraph(source, dest, false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.transpose:
..cat:Graph
..summary:Transposes a graph, either in-place or from source to dest.
..signature:transpose(source [, dest])
..param.source:Source graph.
...type:Class.Graph
..param.dest:Destination graph.
...type:Class.Graph
..returns:void
*/

template<typename TCargo, typename TSpec>
inline void 
transpose(Graph<Directed<TCargo, TSpec> > const& source,
		  Graph<Directed<TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	_copyGraph(source, dest, true);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void 
transpose(Graph<Directed<TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	Graph<Directed<TCargo, TSpec> > dest;
	_copyGraph(g, dest, true);
	g = dest;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.numEdges:
..cat:Graph
..summary:Number of edges in a graph.
..signature:numEdges(g)
..param.g:A graph.
...type:Class.Graph
..returns:Number of edges.
..see:Function.numVertices
*/

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type 
numEdges(Graph<Directed<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.numVertices:
..cat:Graph
..summary:Number of vertices in a graph.
..signature:numVertices(g)
..param.g:A graph.
...type:Class.Graph
..returns:Number of vertices.
..see:Function.numEdges
*/

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type 
numVertices(Graph<Directed<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.empty:
..cat:Graph
..param.object.type:Class.Graph
*/

template<typename TCargo, typename TSpec>
inline bool 
empty(Graph<Directed<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.clearEdges:
..cat:Graph
..summary:Removes all edges in a graph.
..signature:clearEdges(g)
..param.g:A graph.
...type:Class.Graph
..returns:void
..see:Function.clearVertices
..see:Function.clear
*/

template<typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Directed<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
	TIter it = begin(g.data_vertex, Standard());
	TIter itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it != itEnd; ++it, ++pos) 
		if (*it != (TEdgeStump*) 0) removeOutEdges(g, pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.clearVertices:
..cat:Graph
..summary:Removes all vertices in a graph.
..signature:clearVertices(g)
..param.g:A graph.
...type:Class.Graph
..returns:void
..see:Function.clearEdges
..see:Function.clear
*/

template<typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Directed<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}


//////////////////////////////////////////////////////////////////////////////


/**
.Function.clear:
..cat:Graph
..param.object.type:Class.Graph
..remarks:If $object$ is a @Class.Graph.graph@, then all vertices and all edges are removed.
*/

template<typename TCargo, typename TSpec>
inline void 
clear(Graph<Directed<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.outDegree:
..cat:Graph
..summary:Number of outgoing edges for a given vertex.
..signature:outDegree(g, vertex)
..param.g:A graph.
...type:Class.Graph
..param.g:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:Number of out-edges.
..see:Function.inDegree
..see:Function.degree
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type 
outDegree(Graph<Directed<TCargo, TSpec> > const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStump* current = getValue(g.data_vertex, vertex);
	while(current!=0) {
		current = getNextT(current);
		++count;
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.inDegree:
..cat:Graph
..summary:Number of incoming edges for a given vertex.
..signature:inDegree(g, vertex)
..param.g:A graph.
...type:Class.Graph
..param.g:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:Number of in-edges.
..see:Function.outDegree
..see:Function.degree
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type 
inDegree(Graph<Directed<TCargo, TSpec> > const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());

	TSize count=0;
	for(;it!=itEnd;++it) {
		TEdgeStump* current = *it;
		while(current!=0) {
			if ( (TVertexDescriptor) getTarget(current) == vertex) ++count;
			current = getNextT(current);
		}
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.degree:
..cat:Graph
..summary:Number of incident edges for a given vertex.
..signature:degree(g, vertex)
..param.g:A graph.
...type:Class.Graph
..param.g:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:Number of incident edges.
..see:Function.outDegree
..see:Function.inDegree
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Directed<TCargo, TSpec> > >::Type 
degree(Graph<Directed<TCargo, TSpec> > const& g, 
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (inDegree(g,vertex)+outDegree(g,vertex));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addVertex:
..cat:Graph
..summary:Adds a new vertex to the graph.
..signature:addVertex(g)
..param.g:A graph.
...type:Class.Graph
..returns:A new vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.removeVertex
*/

template<typename TCargo, typename TSpec> 
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
addVertex(Graph<Directed<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT	
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) appendValue(g.data_vertex, (TEdgeStump*) 0); 
	else g.data_vertex[vd] = (TEdgeStump*) 0;
	return vd;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeVertex:
..cat:Graph
..summary:Removes a vertex.
..signature:removeVertex(g, v)
..param.g:A graph.
...type:Class.Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.addVertex
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Directed<TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addEdge:
..cat:Graph
..summary:Adds a new edge to the graph, either with or without cargo.
..remarks:For automatons a label is required.
..signature:addEdge(g, source, target [,cargo | ,label])
..signature:addEdge(g, source, target [,label ,cargo])
..param.g:A graph.
...type:Class.Graph
..param.source:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.target:A second vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.label:A label for the edge.
...type:Metafunction.Alphabet
...remarks:Can only be used with @Spec.Automaton@.
..param.cargo:A cargo object.
...type:Metafunction.Cargo
..returns:A new edge descriptor.
...type:Metafunction.EdgeDescriptor
..see:Function.removeEdge
..see:Function.addEdges
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
addEdge(Graph<Directed<TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Id<TGraph>::Type TId;

	TEdgeStump* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	assignTarget(edge_ptr, target);
	assignNextT(edge_ptr, (TEdgeStump*) 0);
	TId id = obtainId(g.data_id_managerE);
	_assignId(edge_ptr, id);
	if (g.data_vertex[source]!=0) assignNextT(edge_ptr, getValue(g.data_vertex, source));
	value(g.data_vertex, source)=edge_ptr;
	return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
addEdge(Graph<Directed<TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo const cargo) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,source,target);
	assignCargo(e,cargo);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeEdge:
..cat:Graph
..summary:Removes an edge from the graph. For automatons a label is required.
..signature:removeEdge(g, source, target [, label])
..signature:removeEdge(g, e)
..param.g:A graph.
...type:Class.Graph
..param.source:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.target:A second vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.label:A label for the edge.
...type:Metafunction.Alphabet
...remarks:Required to remove an edge in an @Spec.Automaton@.
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..returns:void
..see:Function.addEdge
..see:Function.addEdges
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeEdge(Graph<Directed<TCargo, TSpec> >& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	// Find edge and predecessor
	TEdgeStump* pred = 0;
	TEdgeStump* current = g.data_vertex[source];
	while(current != (TEdgeStump*) 0) {
		if ( (TVertexDescriptor) getTarget(current) == target) break;
		pred = current;
		current = getNextT(current);
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;
	
	// Relink the next pointer of predecessor
	if (pred != (TEdgeStump*) 0) assignNextT(pred, getNextT(current));
	else g.data_vertex[source] = getNextT(current);
	
	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void 
removeEdge(Graph<Directed<TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, sourceVertex(g,edge)) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, targetVertex(g,edge)) == true)
	
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	// Find edge and predecessor
	TEdgeStump* pred = 0;
	TEdgeStump* current = g.data_vertex[sourceVertex(g,edge)];
	while(current != (TEdgeStump*) 0) {
		if (current == edge) break;
		pred = current;
		current = getNextT(current);
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;
	
	// Relink the next pointer of predecessor
	if (pred != (TEdgeStump*) 0) assignNextT(pred, getNextT(current));
	else g.data_vertex[sourceVertex(g,edge)] = getNextT(current);
	
	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeOutEdges:
..cat:Graph
..summary:Removes the outgoing edges of a given vertex.
..signature:removeOutEdges(g, v)
..param.g:A graph.
...type:Class.Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.removeInEdges
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<Directed<TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	while(g.data_vertex[v] != (TEdgeStump*) 0) {
		TVertexDescriptor target = targetVertex(g, g.data_vertex[v]);
		removeEdge(g,v,target);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeInEdges:
..cat:Graph
..summary:Removes the incoming edges of a given vertex.
..signature:removeInEdges(g, v)
..param.g:A graph.
...type:Class.Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.removeOutEdges
*/

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<Directed<TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
	TIter it = begin(g.data_vertex, Standard());
	TIter itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		TEdgeStump* current = *it;
		TVertexDescriptor const sourceVertex = pos;
		while(current!=0) {
			if ( (TVertexDescriptor) current->data_target==v) {
				removeEdge(g, sourceVertex, v);
				current = g.data_vertex[sourceVertex];
			} else {
				current = getNextT(current);
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.targetVertex:
..cat:Graph
..summary:Returns the target vertex of an edge.
..remarks:In a tree the target vertex is always the child. 
In an undirected graph the larger vertex descriptor of the two endpoints is the target.
For an out-edge iterator the target is always the vertex the out-edge iterator has not been initialized with.
..signature:targetVertex(g, e)
..signature:targetVertex(it)
..param.g:A graph.
...type:Class.Graph
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..param.it:An edge iterator.
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.sourceVertex
*/

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
targetVertex(Graph<Directed<TCargo, TSpec> > const&,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.sourceVertex:
..cat:Graph
..summary:Returns the source vertex of an edge.
..remarks:In a tree the source vertex is always the parent. 
In an undirected graph the smaller vertex descriptor is the source.
Note: If source vertices are not stored in the EdgeStump this operation is expensive.
Consider using sourceVertex directly on an edge iterator where this operation is fast!
..signature:sourceVertex(g, e)
..signature:sourceVertex(it)
..param.g:A graph.
...type:Class.Graph
..param.e:An edge descriptor.
...type:Metafunction.EdgeDescriptor
..param.it:An edge iterator.
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.targetVertex
*/

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
sourceVertex(Graph<Directed<TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		TEdgeDescriptor current = *it;
		while(current!=(TEdgeDescriptor) 0) {
			if (current == edge) return pos;
			current=getNextT(current);
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getAdjacencyMatrix:
..cat:Graph
..summary:Returns an adjacency matrix representation of the graph.
..signature:getAdjacencyMatrix(g, mat)
..param.g:In-parameter: A graph.
...type:Class.Graph
..param.mat:Out-parameter: A matrix.
...type:Class.Matrix
..returns:void
*/

template<typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Directed<TCargo, TSpec> > const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	typedef typename Value<TMatrix>::Type TMatValue;
	TSize len = getIdUpperBound(g.data_id_managerV);
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	fill(mat, len * len, (TMatValue) 0);
	TVertexDescriptor pos = 0;
	for(;it!=itEnd; ++it, ++pos) {
		TEdgeStump* current = *it;
		TVertexDescriptor const source = pos;
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor target = targetVertex(g,current);
			++mat[source*len+target];
			current = getNextT(current);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.findEdge:
..cat:Graph
..summary:Finds an edge. 
..remarks:In an automaton an edge is uniquely defined by a vertex and a label.
In all other graphs two adjacent vertices uniquely define an edge.
If there are multiple edges between two vertices the behaviour is undefined.
..signature:findEdge(g, v, c)
..signature:findEdge(g, v, w)
..param.g:A graph.
...type:Class.Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.c:An edge label.
...type:Metafunction.Alphabet
..param.w:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:An edge descriptor or 0 if edge is not present. 
Note: In automatons there is always a valid edge descriptor but the target may be nil.
...type:Metafunction.EdgeDescriptor
*/
template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Directed<TCargo, TSpec> > >::Type 
findEdge(Graph<Directed<TCargo, TSpec> > const& g,
		 TVertexDescriptor const v,
		 TVertexDescriptor const w)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, w) == true)
	
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	
	TEdgeStump* current = g.data_vertex[v];
	while(current != (TEdgeStump*) 0) {
		if ( (TVertexDescriptor) getTarget(current) == w) return current;
		current = getNextT(current);
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Directed<TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	_streamWrite(target,"Adjacency list:\n");
	for(;it!=itEnd;++it, ++pos) {
		if (!idInUse(_getVertexIdManager(g), pos)) continue;
		TEdgeStump* current = getValue(it);
		_streamPutInt(target, pos);
		_streamWrite(target," -> ");
		while(current!=0) {
			_streamPutInt(target, getTarget(current));
			_streamPut(target, ',');
			current=getNextT(current);
		}
		_streamPut(target, '\n');
	}
	it = begin(g.data_vertex, Standard());
	pos = 0;
	_streamWrite(target,"Edge list:\n");
	for(;it!=itEnd;++it, ++pos) {
		TEdgeStump* current = getValue(it);
		while(current!=0) {
			_streamWrite(target,"Source: ");
			_streamPutInt(target, pos);		
			_streamPut(target, ',');
			_streamWrite(target,"Target: ");
			_streamPutInt(target, getTarget(current));
			_streamPut(target, ' ');
			_streamWrite(target,"(Id: ");
			_streamPutInt(target, _getId(current));
			_streamPut(target, ')');
			_streamPut(target, '\n');
			current=getNextT(current);
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
