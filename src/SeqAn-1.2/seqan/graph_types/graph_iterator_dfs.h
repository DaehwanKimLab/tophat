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
  $Id: graph_iterator_dfs.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_DFS_H
#define SEQAN_HEADER_GRAPH_ITERATOR_DFS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph DfsIterator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Dfs Preorder Iterator:
..cat:Graph
..summary:Depth-first search iterator for @Class.Graph@.
..remarks:Preorder means that a vertex is enumerated before its adjacent vertices have been explored.
..signature:Iterator<TGraph, DfsPreorder>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Vertex Iterator
..see:Spec.Out-Edge Iterator
..see:Spec.Edge Iterator
..see:Spec.Adjacency Iterator
..see:Spec.Bfs Iterator
*/
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	String<bool> data_tokenMap;			// Which vertices have been visited
	String<TVertexDescriptor> data_stack;
	
	void _init() {
		resizeVertexMap(*data_host,data_tokenMap);
		typedef typename Iterator<String<bool>, Rooted>::Type TIter;
		TIter it = begin(data_tokenMap);
		for(;!atEnd(it);goNext(it)) {
			assignValue(it,false);
		}
		assignProperty(data_tokenMap, data_source, true);
		clear(data_stack);
		appendValue(data_stack, data_source, Generous());
	}

	Iter()
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph& _graph, TVertexDescriptor v) : 
		data_host(&_graph),
		data_source(v)
	{
		SEQAN_CHECKPOINT
		_init();
	}
	
	
	Iter(Iter const& _iter) :
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_tokenMap(_iter.data_tokenMap),
		data_stack(_iter.data_stack)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host=_other.data_host;
		data_source=_other.data_source;
		data_tokenMap=_other.data_tokenMap;
		data_stack=_other.data_stack;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalDfsIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph>
struct Iterator<TGraph, DfsPreorder>
{	
	typedef Iter<TGraph, GraphIterator<InternalDfsIterator<DfsPreorder> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, DfsPreorder>
{	
	typedef Iter<TGraph const, GraphIterator<InternalDfsIterator<DfsPreorder> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalDfsIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalBfsIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return getValue(it.data_stack, length(it.data_stack) - 1);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	// We don't want vertex ids to be changed
	return getValue(it);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
}

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (empty(it.data_stack)) return false;
	else return (getValue(it.data_stack, length(it.data_stack) - 1) == it.data_source);
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it._init();
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (empty(it.data_stack));
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	clear(it.data_stack);
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (empty(it.data_stack)) return;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor u = getValue(it.data_stack, length(it.data_stack) - 1);
	resize(it.data_stack, length(it.data_stack) - 1);
	typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itad(*it.data_host,u);
	for(;!atEnd(itad);goNext(itad)) {
		TVertexDescriptor v = getValue(itad);
		if (getProperty(it.data_tokenMap, v) == false) {
			assignProperty(it.data_tokenMap, v, true);
			appendValue(it.data_stack, v, Generous());
		}
	}
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_source==it2.data_source) &&
			(it1.data_tokenMap==it2.data_tokenMap) &&
			(it1.data_stack==it2.data_stack));
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalDfsIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_source!=it2.data_source) ||
			(it1.data_tokenMap!=it2.data_tokenMap) ||
			(it1.data_stack!=it2.data_stack));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
