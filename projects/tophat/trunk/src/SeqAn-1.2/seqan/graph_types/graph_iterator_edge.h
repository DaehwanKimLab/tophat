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
  $Id: graph_iterator_edge.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_EDGE_H
#define SEQAN_HEADER_GRAPH_ITERATOR_EDGE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph EdgeIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Edge Iterator:
..cat:Graph
..summary:Edge iterator for @Class.Graph@.
..signature:Iterator<TGraph, EdgeIterator>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Vertex Iterator
..see:Spec.Out-Edge Iterator
..see:Spec.Adjacency Iterator
..see:Spec.Bfs Iterator
..see:Spec.Dfs Preorder Iterator
*/

template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TVertexIterator data_vertex_it;
	TOutEdgeIterator data_edge_it;
	TVertexDescriptor data_first_slot;


	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph) : 
		data_vertex_it(_graph),
		data_edge_it(_graph, getIdLowerBound(_getVertexIdManager(_graph)))  
	{
		SEQAN_CHECKPOINT
		while((atEnd(data_edge_it)) && (!atEnd(data_vertex_it))) 
		{
				goNext(data_vertex_it);
				typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
				data_edge_it = TOutEdgeIterator(hostGraph(*this), value(data_vertex_it));			
		}
		data_first_slot = value(data_vertex_it);
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : 
		data_vertex_it(_iter.data_vertex_it),
		data_edge_it(_iter.data_edge_it),
		data_first_slot(_iter.data_first_slot)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_vertex_it = _other.data_vertex_it;
		data_edge_it = _other.data_edge_it;
		data_first_slot = _other.data_first_slot;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
struct Iterator<TGraph, EdgeIterator>
{	
	typedef Iter<TGraph, GraphIterator<InternalEdgeIterator<EdgeIterator> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, EdgeIterator>
{	
	typedef Iter<TGraph const, GraphIterator<InternalEdgeIterator<EdgeIterator> > > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph const, GraphIterator<InternalEdgeIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalEdgeIterator - Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return getValue(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return value(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return hostGraph(it.data_vertex_it);
} 

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return ((it.data_first_slot == getValue(it.data_vertex_it)) &&
			atBegin(it.data_edge_it));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	goBegin(it.data_vertex_it);
	while (it.data_first_slot != getValue(it.data_vertex_it)) ++it.data_vertex_it;
	it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return atEnd(it.data_vertex_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goEnd(it.data_edge_it);
	goEnd(it.data_vertex_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
_goNextInternal(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) {
		goNext(it.data_edge_it);
		if(!atEnd(it.data_edge_it)) return;
		else {
			while ((atEnd(it.data_edge_it)) && (!atEnd(it.data_vertex_it))) {
				goNext(it.data_vertex_it);
				if ((atEnd(it.data_vertex_it)) ||
					(!idInUse(_getVertexIdManager(hostGraph(it)), getValue(it.data_vertex_it)))) continue;
				typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
				it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));			
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	_goNextInternal(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	_goNextInternal(it);
	TVertexDescriptor sourceV = sourceVertex(it.data_edge_it);
	while((!atEnd(it)) && (targetVertex(hostGraph(it), getValue(it.data_edge_it)) == sourceV)) {
		_goNextInternal(it);
		sourceV = sourceVertex(it.data_edge_it);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
_goPreviousInternal(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) {
		if ((atEnd(it)) ||
			(atBegin(it.data_edge_it))) 
		{
			--it.data_vertex_it;
			typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
			it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));
			while((!atBegin(it.data_vertex_it)) && 
					(atEnd(it.data_edge_it))) {
				--it.data_vertex_it;
				typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
				it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));			
			}
			goEnd(it.data_edge_it);
		}
		--it.data_edge_it;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	_goPreviousInternal(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Undirected<TCargo, TGraphSpec> >, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	_goPreviousInternal(it);
	TVertexDescriptor sourceV = sourceVertex(it.data_edge_it);
	while((!atBegin(it)) && (targetVertex(hostGraph(it), getValue(it.data_edge_it)) == sourceV)) {
		_goPreviousInternal(it);
		sourceV = sourceVertex(it.data_edge_it);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_vertex_it==it2.data_vertex_it) && 
			(it1.data_edge_it==it2.data_edge_it));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_vertex_it!=it2.data_vertex_it) || 
			(it1.data_edge_it!=it2.data_edge_it));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
sourceVertex(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return sourceVertex(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
targetVertex(Iter<TGraph, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Alphabet<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> > >::Type
label(Iter<Graph<Automaton<TAlphabet, TCargo, TGraphSpec> >, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	return label(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
