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
  $Id: graph_impl_undirected.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_UNDIRECTED_H
#define SEQAN_HEADER_GRAPH_IMPL_UNDIRECTED_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Undirected
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Undirected graph:
..cat:Graph
..general:Class.Graph
..summary:An undirected graph that stores the edges in an adjacency list.
..description:
...image:undirectedGraph|An undirected graph with multiple components.
..signature:Graph<Undirected<TCargo, TSpec> >
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
template<typename TCargo,typename TSpec>
class Graph<Undirected<TCargo, TSpec> > 
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
inline String<typename EdgeType<Graph<Undirected<TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Undirected<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	return const_cast<String<TEdgeStump*>&>(g.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Undirected<TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Undirected<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
	return const_cast<TVertexIdManager&>(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Undirected<TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Undirected<TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
	return const_cast<TEdgeIdManager&>(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Undirected<TCargo, TSpec> > const& source,
		   Graph<Undirected<TCargo, TSpec> >& dest,
		   bool) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;
	clear(dest);
	resize(dest.data_vertex, length(_getVertexString(source)));
	TIter itInit = begin(dest.data_vertex, Standard());
	TIter itInitEnd = end(dest.data_vertex, Standard());
	for(;itInit!=itInitEnd; ++itInit) *itInit = (TEdgeStump*) 0;
	TIterConst it = begin(source.data_vertex, Standard());
	TIterConst itEnd = end(source.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		TEdgeStump* current = *it;
		TVertexDescriptor sourceVertex = pos;
		while(current != (TEdgeStump*) 0) {
			TVertexDescriptor targetVertex = getTarget(current);

			if (targetVertex != sourceVertex) {
				// Create missing vertices, targetVertex is always bigger than sourceVertex
				_createVertices(dest,targetVertex);
			
				// Add edge
				TEdgeDescriptor e = addEdge(dest, sourceVertex, targetVertex);
				_assignId(e, _getId(current));
				assignCargo(e, getCargo(current));
				current = getNextS(current);  // Follow the link belonging to the source id
			} else {
				// Do nothing here because we don't want to create edges twice!!!
				current = getNextT(current);
			}
		}
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}


template<typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Undirected<TCargo, TSpec> > const& source,
		   Graph<Undirected<TCargo, TSpec> >& dest) 
{
	_copyGraph(source, dest, false); // Never transpose because undirected and transposed undirected graph are equal
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Undirected<TCargo, TSpec> > const& source,
		  Graph<Undirected<TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	// Undirected graph, no transpose just copy
	_copyGraph(source, dest, false);
	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
transpose(Graph<Undirected<TCargo, TSpec> > const&)
{
	// Nothing to do
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Undirected<TCargo, TSpec> > >::Type 
numEdges(Graph<Undirected<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Size<Graph<Undirected<TCargo, TSpec> > >::Type 
numVertices(Graph<Undirected<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline bool 
empty(Graph<Undirected<TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Undirected<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*>, Standard>::Type TIter;

	// Collect all edges
	String<TEdgeDescriptor> edges;
	TIter it = begin(g.data_vertex, Standard());
	TIter itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd; ++it, ++pos) {
		TEdgeStump* current = *it;
		TVertexDescriptor sourceVertex = pos;
		while(current != (TEdgeStump*) 0) {
			if (getTarget(current) != sourceVertex) {
				appendValue(edges, current, Generous());
				current = getNextS(current);
			}
			// Do nothing here because we don't want to create edges twice!!!
			else current = getNextT(current);
		}
		*it = (TEdgeStump*) 0;
	}
	SEQAN_TASSERT(numEdges(g) == length(edges))

	// Release all edges
	typedef typename Iterator<String<TEdgeDescriptor>, Standard>::Type TStringIter;
	TStringIter edgeIt = begin(edges, Standard());
	TStringIter edgeEndIt = end(edges, Standard());
	for(; edgeIt != edgeEndIt; ++edgeIt) {
		valueDestruct(*edgeIt);
		deallocate(g.data_allocator, *edgeIt, 1);	
	}
	releaseAll(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Undirected<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline void 
clear(Graph<Undirected<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Undirected<TCargo, TSpec> > >::Type 
outDegree(Graph<Undirected<TCargo, TSpec> > const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TEdgeStump* current = g.data_vertex[vertex];
	while(current!=0) {
		if ( (TVertexDescriptor) getTarget(current)==vertex) current = getNextT(current);
		else current = getNextS(current);
		++count;
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Undirected<TCargo, TSpec> > >::Type 
inDegree(Graph<Undirected<TCargo, TSpec> > const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	// In-degree and out-degree are equal for undirected graphs
	return outDegree(g,vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Undirected<TCargo, TSpec> > >::Type 
degree(Graph<Undirected<TCargo, TSpec> > const& g,
	   TVertexDescriptor const vertex) 
{
	return outDegree(g,vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec> 
inline typename VertexDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
addVertex(Graph<Undirected<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) appendValue(g.data_vertex, (TEdgeStump*) 0, Generous()); 
	else g.data_vertex[vd] = (TEdgeStump*) 0;
	return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Undirected<TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
addEdge(Graph<Undirected<TCargo, TSpec> >& g, 
		TVertexDescriptor source, 
		TVertexDescriptor target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(source != target)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)

	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Id<TGraph>::Type TId;

	// Source must be the smaller vertex id
	if (source > target) {TVertexDescriptor tmp = target; target = source; source = tmp; }

	TEdgeStump* edge_ptr;
	allocate(g.data_allocator, edge_ptr, 1);
	valueConstruct(edge_ptr);
	assignSource(edge_ptr, source);
	assignTarget(edge_ptr, target);
	assignNextS(edge_ptr, (TEdgeStump*) 0);
	assignNextT(edge_ptr, (TEdgeStump*) 0);
	TId id = obtainId(g.data_id_managerE);
	_assignId(edge_ptr, id);
	if (g.data_vertex[source]!=0) assignNextS(edge_ptr, g.data_vertex[source]);
	if (g.data_vertex[target]!=0) assignNextT(edge_ptr, g.data_vertex[target]);
	g.data_vertex[source]=edge_ptr;
	g.data_vertex[target]=edge_ptr;
	return edge_ptr;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCargo2> 
inline typename EdgeDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
addEdge(Graph<Undirected<TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo2 const cargo) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(source != target)
	
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g,source,target);
	assignCargo(e,cargo);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeEdge(Graph<Undirected<TCargo, TSpec> >& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;

	// Find edge and source predecessor
	TEdgeStump* predSource = 0;
	TEdgeStump* current = g.data_vertex[source];
	while(current != (TEdgeStump*) 0) {
		TVertexDescriptor adjV = (TVertexDescriptor) getTarget(current);
		if (adjV != source) {
			if ( adjV == target) break;
			predSource = current;
			current = getNextS(current);
		} else {
			adjV = (TVertexDescriptor) getSource(current);
			if ( adjV == target) break;
			predSource = current;
			current = getNextT(current);
		}
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;

	// Find edge and target predecessor
	TEdgeStump* predTarget = 0;
	current = g.data_vertex[target];
	while(current != (TEdgeStump*) 0) {
		TVertexDescriptor adjV = (TVertexDescriptor) getTarget(current);
		if (adjV != target) {
			if ( adjV == source) break;
			predTarget = current;
			current = getNextS(current);
		} else {
			adjV = (TVertexDescriptor) getSource(current);
			if ( adjV == source) break;
			predTarget = current;
			current = getNextT(current);
		}
	}

	
	// Relink the next pointer of source predecessor
	if (predSource != (TEdgeStump*) 0) {
		if (source != (TVertexDescriptor) getTarget(current)) {
			if (source != (TVertexDescriptor) getTarget(predSource)) assignNextS(predSource, getNextS(current));
			else assignNextT(predSource, getNextS(current));
		} else {
			if (source != (TVertexDescriptor) getTarget(predSource)) assignNextS(predSource, getNextT(current));
			else assignNextT(predSource, getNextT(current));
		}
	}
	else {
		if (source != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, source) = getNextS(current);
		else value(g.data_vertex, source) = getNextT(current);
	}

	// Relink the next pointer of target predecessor
	if (predTarget != (TEdgeStump*) 0) {
		if (target != (TVertexDescriptor) getTarget(current)) {
			if (target != (TVertexDescriptor) getTarget(predTarget)) assignNextS(predTarget, getNextS(current));
			else assignNextT(predTarget, getNextS(current));
		} else {
			if (target != (TVertexDescriptor) getTarget(predTarget)) assignNextS(predTarget, getNextT(current));
			else assignNextT(predTarget, getNextT(current));
		}
	}
	else {
		if (target != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, target) = getNextS(current);
		else value(g.data_vertex, target) = getNextT(current);
	}

	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void 
removeEdge(Graph<Undirected<TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, sourceVertex(g,edge)) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, targetVertex(g,edge)) == true)
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor source = sourceVertex(g,edge);
	TVertexDescriptor target = targetVertex(g,edge);

	// Find edge and source predecessor
	TEdgeStump* predSource = 0;
	TEdgeStump* current = g.data_vertex[source];
	while(current != (TEdgeStump*) 0) {
		TVertexDescriptor adjV = (TVertexDescriptor) getTarget(current);
		if (adjV != source) {
			if ((adjV == target) && (current == edge)) break;
			predSource = current;
			current = getNextS(current);
		} else {
			adjV = (TVertexDescriptor) getSource(current);
			if ((adjV == target) && (current == edge)) break;
			predSource = current;
			current = getNextT(current);
		}
	}
	
	// Not found?
	if (current == (TEdgeStump*) 0) return;

	// Find edge and target predecessor
	TEdgeStump* predTarget = 0;
	current = g.data_vertex[target];
	while(current != (TEdgeStump*) 0) {
		TVertexDescriptor adjV = (TVertexDescriptor) getTarget(current);
		if (adjV != target) {
			if ((adjV == source) && (current == edge)) break;
			predTarget = current;
			current = getNextS(current);
		} else {
			adjV = (TVertexDescriptor) getSource(current);
			if ((adjV == source) && (current == edge)) break;
			predTarget = current;
			current = getNextT(current);
		}
	}

	
	// Relink the next pointer of source predecessor
	if (predSource != (TEdgeStump*) 0) {
		if (source != (TVertexDescriptor) getTarget(current)) {
			if (source != (TVertexDescriptor) getTarget(predSource)) assignNextS(predSource, getNextS(current));
			else assignNextT(predSource, getNextS(current));
		} else {
			if (source != (TVertexDescriptor) getTarget(predSource)) assignNextS(predSource, getNextT(current));
			else assignNextT(predSource, getNextT(current));
		}
	}
	else {
		if (source != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, source) = getNextS(current);
		else value(g.data_vertex, source) = getNextT(current);
	}

	// Relink the next pointer of target predecessor
	if (predTarget != (TEdgeStump*) 0) {
		if (target != (TVertexDescriptor) getTarget(current)) {
			if (target != (TVertexDescriptor) getTarget(predTarget)) assignNextS(predTarget, getNextS(current));
			else assignNextT(predTarget, getNextS(current));
		} else {
			if (target != (TVertexDescriptor) getTarget(predTarget)) assignNextS(predTarget, getNextT(current));
			else assignNextT(predTarget, getNextT(current));
		}
	}
	else {
		if (target != (TVertexDescriptor) getTarget(current)) value(g.data_vertex, target) = getNextS(current);
		else value(g.data_vertex, target) = getNextT(current);
	}

	// Deallocate
	releaseId(g.data_id_managerE, _getId(current));
	valueDestruct(current);
	deallocate(g.data_allocator, current, 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<Undirected<TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor eD = g.data_vertex[v];
	while(eD != (TEdgeStump*) 0) {
		removeEdge(g,eD);
		eD = g.data_vertex[v];
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<Undirected<TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeOutEdges(g,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
targetVertex(Graph<Undirected<TCargo, TSpec> > const&,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getTarget(edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
sourceVertex(Graph<Undirected<TCargo, TSpec> > const&,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return getSource(edge);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Undirected<TCargo, TSpec> > const& g, 
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Size<TMatrix>::Type TSize;
	TSize len = getIdUpperBound(g.data_id_managerV);
	fill(mat, len*len, 0);

	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	TVertexDescriptor pos = 0;
	for(;it!=itEnd;++it, ++pos) {
		TVertexDescriptor sourceV = pos;
		TEdgeStump* current = *it;
		while(current!=0) {
			TVertexDescriptor adjV = getTarget(current);
			if (adjV != sourceV) {
				++mat[sourceV*len+adjV];
				current=getNextS(current);
			} else {
				adjV = getSource(current);
				++mat[sourceV*len+adjV];
				current=getNextT(current);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Undirected<TCargo, TSpec> > >::Type 
findEdge(Graph<Undirected<TCargo, TSpec> > const& g,
		 TVertexDescriptor const v,
		 TVertexDescriptor const w)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, w) == true)
	
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	
	TEdgeStump* current = g.data_vertex[v];
	while(current != (TEdgeStump*) 0) {
		TVertexDescriptor adjV = getTarget(current);
		if (adjV != v) {
			if ( adjV == w) return current;
			current=getNextS(current);
		} else {
			adjV = getSource(current);
			if ( adjV == w) return current;
			current=getNextT(current);
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Undirected<TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const, Standard>::Type TIterConst;
	TIterConst it = begin(g.data_vertex, Standard());
	TIterConst itEnd = end(g.data_vertex, Standard());
	_streamWrite(target,"Adjacency list:\n");
	TVertexDescriptor pos = 0;
	for(;it != itEnd; ++it, ++pos) {
		TVertexDescriptor sourceV = pos;
		_streamPutInt(target, sourceV);
		_streamWrite(target," -> ");
		TEdgeStump* current = *it;
		while(current!=0) {
			TVertexDescriptor adjV = getTarget(current);
			if (adjV != sourceV) {
				_streamPutInt(target, adjV);
				_streamPut(target, ',');
				current=getNextS(current);
			} else {
				adjV = getSource(current);
				_streamPutInt(target, adjV);
				_streamPut(target, ',');
				current=getNextT(current);
			}
		}
		_streamPut(target, '\n');
	}
	it = begin(g.data_vertex, Standard());
	pos = 0;
	_streamWrite(target,"Edge list:\n");
	for(; it != itEnd; ++it, ++pos) {
		TVertexDescriptor sourceV = pos;
		TEdgeStump* current = *it;
		while(current!=0) {
			TVertexDescriptor targetV = getTarget(current);
			if (sourceV != targetV) {
				_streamWrite(target,"Source: ");
				_streamPutInt(target, sourceV);		
				_streamPut(target, ',');
				_streamWrite(target,"Target: ");
				_streamPutInt(target, targetV);
				_streamPut(target, ' ');
				_streamWrite(target,"(Id: ");
				_streamPutInt(target, _getId(current));
				_streamPut(target, ')');
				_streamPut(target, '\n');
				current=getNextS(current);
			} else {
				current=getNextT(current);
			}
		}
	}
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
