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
  $Id: graph_iterator.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_H
#define SEQAN_HEADER_GRAPH_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph Iterators
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
struct GraphIterator;

	// Vertex iterator
	template <typename TSpec = Default>
	struct InternalVertexIterator;

	// Edge iterator
	template <typename TSpec = Default>
	struct InternalEdgeIterator;

	// OutEdge iterator
	template <typename TSpec = Default>
	struct InternalOutEdgeIterator;

	// Adjacency iterator
	template <typename TSpec = Default>
	struct InternalAdjacencyIterator;

	// Bfs iterator
	template <typename TSpec = Default>
	struct InternalBfsIterator;

	// Dfs iterator
	template <typename TSpec = Default>
	struct InternalDfsIterator;

//////////////////////////////////////////////////////////////////////////////
// Graph Iterators - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Graph

template<typename TGraph, typename TIteratorSpec>
struct Host<Iter<TGraph, GraphIterator<TIteratorSpec> > >
{	
	typedef TGraph Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Host<Iter<TGraph const, GraphIterator<TIteratorSpec> > >
{	
	typedef TGraph const Type;
};


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
