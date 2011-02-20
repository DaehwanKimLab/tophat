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
  $Id: graph_base.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeDescriptor:
..summary:Type of an object that represents an edge descriptor.
..signature:EdgeDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs use a pointer to an edge stump as an edge descriptor.
..returns.param.Type:EdgeDescriptor type.
..remarks.text:The edge descriptor is a unique handle to a given edge in a graph.
It is used in various graph functions, e.g., to remove edges, to assign a cargo to an edge or to get the endpoints of an edge.
It is also used to attach properties to edges.
..example.code:EdgeDescriptor<Graph<> >::Type eD; //eD is an edge descriptor
*/
template<typename T>
struct EdgeDescriptor;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Cargo
..example.code:Cargo<Graph<Directed<int> > >::Type c; //c has type int
*/
template<typename T>
struct Cargo;


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeType:
..summary:Edge type of a graph object.
..signature:EdgeType<T>::Type
..param.T:Type T must be a graph.
..returns.param.Type:Edge type.
..remarks.text:The specific edge stump type that is used in a graph.
..example.code:EdgeType<TGraph>::Type e; //e is an edge in TGraph
*/
template<typename T>
struct EdgeType;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Alphabet:
..summary:Access to the Alphabet type.
..signature:Alphabet<T>::Type
..param.T:Type T must be a type that uses some kind of alphabet internally.
..returns.param.Type:Alphabet type.
..remarks.text:Type T can be for example an automaton where the alphabet type describes the domain of the transition labels.
..example.code:Alphabet<Graph<Automaton<Dna> > >::Type alph; //alph is of type Dna
*/
template<typename T>
struct Alphabet;


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.EdgeIdHandler:
..summary:Type of an object that represents an Id Manager.
..signature:EdgeIdHandler<T>::Type
..param.T:A graph.
...type:Class.Graph
..returns.param.Type:IdManager type.
..remarks.text:The exact IdManager type depends on the edge stump.
If the edge stump is id-free the IdManager simply counts edge ids, 
otherwise it manages a list of free and used ids.
*/
template<typename T>
struct EdgeIdHandler;


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.VertexIdHandler:
..summary:Type of an object that represents an Id Manager.
..signature:VertexIdHandler<T>::Type
..param.T:A graph.
..returns.param.Type:IdManager type.
*/
template<typename T>
struct VertexIdHandler;


//////////////////////////////////////////////////////////////////////////////
// General Graph Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;

//////////////////////////////////////////////////////////////////////////////

struct TreeTag_;
typedef Tag<TreeTag_> const TreeTag;


//////////////////////////////////////////////////////////////////////////////
// Graph Iterator Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Graph Iterator:
..summary:A specification of the iterator to traverse a graph.
*/

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.VertexIterator:
	Traverses all vertices of a graph.
*/
struct VertexIterator_;
typedef Tag<VertexIterator_> const VertexIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.EdgeIterator:
	Traverses all edges of a graph.
*/
struct EdgeIterator_;
typedef Tag<EdgeIterator_> const EdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.OutEdgeIterator:
	Traverses all edges of a graph given a vertex.
*/
struct OutEdgeIterator_;
typedef Tag<OutEdgeIterator_> const OutEdgeIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.AdjacencyIterator:
	Traverses all neighbors of a graph given a vertex.
*/
struct AdjacencyIterator_;
typedef Tag<AdjacencyIterator_> const AdjacencyIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.BfsIterator:
	Traverses all vertices of a graph in Bfs order.
*/
struct BfsIterator_;
typedef Tag<BfsIterator_> const BfsIterator;

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Graph Iterator.value.DfsPreorder:
	Traverses all vertices of a graph in Dfs order.
*/
struct DfsPreorder_;
typedef Tag<DfsPreorder_> const DfsPreorder;





//////////////////////////////////////////////////////////////////////////////
// Graph - Default edge stump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo = void, bool TList = true, bool TSource = false, bool TId = true, typename TSpec = Default>
class EdgeStump;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> > 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> 
{
	typedef typename Id<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type Type;
};




//////////////////////////////////////////////////////////////////////////////
// Graph - Default Id Manager
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TIdType = unsigned int, typename TSpec = Default>
class IdManager;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeIdHandler.param.T.type:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, false, TSpec> > {
	typedef IdManager<void> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
struct EdgeIdHandler<EdgeStump<TCargo, TList, TSource, true, TSpec> > {
	typedef IdManager<typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> >::Type> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct VertexIdHandler {
	typedef IdManager<> Type;
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
