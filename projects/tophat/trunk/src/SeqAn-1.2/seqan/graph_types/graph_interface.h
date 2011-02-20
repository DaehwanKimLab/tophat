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
  $Id: graph_interface.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_INTERFACE_H
#define SEQAN_HEADER_GRAPH_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

// Default directed graph
template<typename TCargo = void, typename TSpec = Default>
struct Directed;

// Default undirected graph
template<typename TCargo = void, typename TSpec = Default>
struct Undirected;

// Default Tree
template<typename TCargo = void, typename TSpec = Default>
struct Tree;

// Default Automaton
template<typename TAlphabet = char, typename TCargo = void, typename TSpec = Default>
struct Automaton;

// Default Hmm
template<typename TAlphabet = Dna, typename TCargo = double, typename TSpec = Default>
struct Hmm;


//////////////////////////////////////////////////////////////////////////////
// Graph
//////////////////////////////////////////////////////////////////////////////

/**
.Class.Graph:
..cat:Graph
..summary:Generic graph.
..signature:Graph<TSpec>
..param.TSpec:The specializing type determines the kind of graph, e.g., directed, undirected, tree, or automaton.
...remarks:The default Graph<> corresponds to a directed graph.
...default:Directed<>
*/
template<typename TSpec = Directed<> >
class Graph;

//////////////////////////////////////////////////////////////////////////////
// General Graph Metafunction
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.Graph

template<typename TSpec>
struct Spec<Graph<TSpec> > 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
struct Spec<Graph<TSpec> const>
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeDescriptor.param.T.type:Class.Graph

template<typename TSpec>
struct EdgeDescriptor<Graph<TSpec> > 
{
	typedef typename EdgeType<Graph<TSpec> >::Type* Type;
};

template<typename TSpec>
struct EdgeDescriptor<Graph<TSpec> const>
{
	typedef typename EdgeType<Graph<TSpec> const>::Type* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.VertexDescriptor.param.T.type:Class.Graph

template<typename TSpec>
struct VertexDescriptor<Graph<TSpec> > 
{
	typedef typename Id<Graph<TSpec> >::Type Type;
};

template<typename TSpec>
struct VertexDescriptor<Graph<TSpec> const>
{
	typedef typename Id<Graph<TSpec> >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeType.param.T.type:Class.Graph

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, false, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Directed<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, false, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, true, false, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Directed<TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, true, false, false, WithoutEdgeId> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, false, false, TreeTag> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Tree<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, false, false, TreeTag> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, true, true, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
struct EdgeType<Graph<Undirected<TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, true, true, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, true, true, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
struct EdgeType<Graph<Undirected<TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, true, true, false, WithoutEdgeId> const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TSpec> > > {
	typedef EdgeStump<TCargo, false, false, true, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TSpec> > const> {
	typedef EdgeStump<TCargo, false, false, true, TSpec> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId> > > {
	typedef EdgeStump<TCargo, false, false, false, WithoutEdgeId> Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, WithoutEdgeId> > const> {
	typedef EdgeStump<TCargo, false, false, false, WithoutEdgeId> const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Hmm<TAlphabet, TCargo, TSpec> > const> {
	typedef typename EdgeType<Graph<Directed<TCargo, TSpec> > const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct EdgeType<Graph<Hmm<TAlphabet, TCargo, TSpec> > > {
	typedef typename EdgeType<Graph<Directed<TCargo, TSpec> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.Graph

template<typename TSpec>
struct Cargo<Graph<TSpec> > {
	typedef typename Cargo<typename EdgeType<Graph<TSpec> >::Type>::Type Type;
};


template<typename TSpec>
struct Cargo<Graph<TSpec> const> {
	typedef typename Cargo<typename EdgeType<Graph<TSpec> const>::Type>::Type Type;
};



//////////////////////////////////////////////////////////////////////////////

///.Metafunction.EdgeIdHandler.param.T.type:Class.Graph

template<typename TSpec>
struct EdgeIdHandler<Graph<TSpec> const> {
	typedef typename EdgeIdHandler<typename EdgeType<Graph<TSpec> const>::Type>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
struct EdgeIdHandler<Graph<TSpec> > {
	typedef typename EdgeIdHandler<typename EdgeType<Graph<TSpec> >::Type>::Type Type;
};




//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Alphabet.param.T.type:Class.Graph

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > > {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > const> {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Hmm<TAlphabet, TCargo, TSpec> > > {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
struct Alphabet<Graph<Hmm<TAlphabet, TCargo, TSpec> > const> {
	typedef TAlphabet Type;
};

//////////////////////////////////////////////////////////////////////////////
// Generic Graph Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getNil:
..cat:Graph
..summary:Utility function returning a value that represents nil.
Useful for various graph algorithms, e.g., missing predecessors, vertices that have not been visited, etc.
..signature:getNil<T>()
..returns:Pseudo nil value for type T.
*/
template <typename T>
inline T
getNil(T *)
{
	return ~0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
getNil()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return getNil(_tag);
}


//////////////////////////////////////////////////////////////////////////////
// Purely internal!!! Never compare to _getInfinity()!!!.
// Just returns a very large value.
//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
_getInfinity()
{
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////

template <>
inline double
_getInfinity()
{
	return 1000000000;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TWeightMap>
inline typename Value<TWeightMap>::Type
_getInfinityDistance(TWeightMap const&)
{
	// We need to divide by 2 because of addition in some graph algorithms: infinity + something
	return (_getInfinity<typename Value<TWeightMap>::Type>()/2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T
_getInfinityDistance()
{
	return (_getInfinity<T>() / 2);
}

//////////////////////////////////////////////////////////////////////////////

// Simple _getId function to get the id for a vertex descriptor which is the id!
template<typename TId>
inline TId
_getId(TId const id)
{
	SEQAN_CHECKPOINT
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor>
inline void
_createVertices(Graph<TSpec>& g,
				TVertexDescriptor const maxId) 
{
		// Create missing vertices
		while (maxId >= getIdUpperBound(g.data_id_managerV)) addVertex(g);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addEdges:
..cat:Graph
..summary:Shortcut to add multiple edges at once.
Creates vertices implicitly.
..signature:addEdge(g, edges, size)
..param.g:A graph.
...type:Class.Graph
..param.edges:An array of vertex descriptors. It is assumed that the
edges are stored in the following way: Source1, Target1, Source2, Target2, Source3, ...
For a tree the root must be the first vertex in this array and the enumeration is Parent, Child, Parent, Child, ...
...type:Metafunction.VertexDescriptor
..param.size:Size of the array. Must be a multiple of 2.
...type:Metafunction.Size
..returns:void
..see:Function.addEdge
*/
template<typename TSpec, typename TEdgeArray, typename TSize>
inline void
addEdges(Graph<TSpec>& dest,
		 TEdgeArray const & edges,
		 TSize const size) 
{
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	for(TSize i=0;i<size;++i) {
		TVertexDescriptor source = edges[2*i];
		TVertexDescriptor target = edges[2*i+1];
		// Create missing vertices
		if (source>target) _createVertices(dest,source);
		else _createVertices(dest,target);
		// Add edge
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, source) == true)
		SEQAN_ASSERT(idInUse(dest.data_id_managerV, target) == true)
		addEdge(dest, source, target);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Graph<TSpec> const& source)
{
	write(target, source);
	return target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
