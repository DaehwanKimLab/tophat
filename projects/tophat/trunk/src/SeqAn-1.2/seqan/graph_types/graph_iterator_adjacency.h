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
  $Id: graph_iterator_adjacency.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_ADJACENCY_H
#define SEQAN_HEADER_GRAPH_ITERATOR_ADJACENCY_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph AdjacencyIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Adjacency Iterator:
..cat:Graph
..summary:Adjacency iterator for @Class.Graph@.
..signature:Iterator<TGraph, AdjacencyIterator>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Vertex Iterator
..see:Spec.Out-Edge Iterator
..see:Spec.Edge Iterator
..see:Spec.Bfs Iterator
*/
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator data_edge_it;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_edge_it(_graph, v)
	{
		SEQAN_CHECKPOINT
	}
	
	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : data_edge_it(_iter.data_edge_it)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_edge_it = _other.data_edge_it;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalAdjacencyIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
struct Iterator<TGraph, AdjacencyIterator>
{	
	typedef Iter<TGraph, GraphIterator<InternalAdjacencyIterator<AdjacencyIterator> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, AdjacencyIterator>
{	
	typedef Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<AdjacencyIterator> > > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph const, GraphIterator<InternalAdjacencyIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalAdjacencyIterator - Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return getValue(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return hostGraph(it.data_edge_it);
} 

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return atBegin(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goBegin(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (atEnd(it.data_edge_it));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goEnd(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goNext(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it==it2.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalAdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it!=it2.data_edge_it);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
