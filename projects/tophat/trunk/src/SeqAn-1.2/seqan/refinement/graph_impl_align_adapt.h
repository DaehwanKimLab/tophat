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
  $Id: graph_impl_align_adapt.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_ALIGN_ADAPT_H
#define SEQAN_HEADER_GRAPH_IMPL_ALIGN_ADAPT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Adaptations so that the alignment graph works like any other graph
//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// Alignment Graph OutEdgeIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator for Alignment
//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TIteratorSpec>
class Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TIteratorSpec> > > 
{
public:
	typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	TEdgeDescriptor data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_host(&_graph),
		data_source(v)
	{
		SEQAN_CHECKPOINT
		if (empty(_graph)) data_edge = 0;
		else data_edge = getValue(_graph.data_align.data_vertex,v);
	}
	
	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
	{
	}

	~Iter() {
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, OutEdgeIterator>
{	
	typedef Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};

template<typename TStringSet, typename TCargo, typename TGraphSpec>
struct Iterator<Graph<Alignment<TStringSet, TCargo, TGraphSpec> > const, OutEdgeIterator>
{	
	typedef Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> > const, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > Type;
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalOutEdgeIterator - Functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(_getVertexString(*it.data_host), it.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(_getVertexString(*it.data_host),it.data_source);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) {
		if (it.data_source == getSource(it.data_edge)) it.data_edge = getNextS(it.data_edge);
		else it.data_edge = getNextT(it.data_edge);
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TGraphSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	TEdge* current = getValue(_getVertexString(*it.data_host), it.data_source);
	if (current == it.data_edge) return;
	while (current != 0) {
		if (it.data_source == getSource(current)) {
			if (it.data_edge == getNextS(current)) break;
			else current = getNextS(current);
		} else {
			if (it.data_edge == getNextT(current)) break;
			else current = getNextT(current);
		}
	}
	it.data_edge = current;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TGraphSpec> > >::Type
targetVertex(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalOutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor target = targetVertex(*it.data_host, it.data_edge);
	if (target != it.data_source) return target;
	else return sourceVertex(*it.data_host, it.data_edge);
}



//////////////////////////////////////////////////////////////////////////////
// Go Next function for EdgeIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Alignment<TStringSet, TCargo, TGraphSpec> >, GraphIterator<InternalEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	_goNextInternal(it);
	TVertexDescriptor sourceV = sourceVertex(it.data_edge_it);
	while((!atEnd(it)) && (targetVertex(hostGraph(it), getValue(it.data_edge_it)) == sourceV)) {
		_goNextInternal(it);
		sourceV = sourceVertex(it.data_edge_it);
	}
}




//////////////////////////////////////////////////////////////////////////////
// Graph drawing adaptations
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TNodeAttributes>
inline void
_createNodeAttributes(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
					  TNodeAttributes& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	resizeVertexMap(g, nodeMap);

	unsigned int scaling = 20 / length(value(g.data_sequence));
	if (scaling == 0) scaling = 1;

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		TIdType id = sequenceId(g, *it);
		std::ostringstream outs; 
		outs << "label = \"";
		outs << label(g, *it);
		outs << "\", group = \"";
		outs << id;
		outs << "\"";
		append(property(nodeMap, *it), outs.str().c_str());		
		//std::cout << property(nodeMap, *it) << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Alignment<TStringSet, void, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g.data_align,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TCargo c = getCargo(*itEd);
		std::ostringstream outs; 
		outs << "label = \"";
		outs << (TCargo) c;
		outs << "\"";
		append(property(edgeMap, *itEd), outs.str().c_str());
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile & file,
				  Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				  DotDrawing)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TStringSet& str = stringSet(g);
	TSize len = length(str);
	for(TSize i = 0; i<len; ++i) {
		TId seqId = positionToId(str, i);
		TSize j = 0;
		TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();	
		TVertexDescriptor previousVertex = nilVertex;
		while(j<length(str[i])) {
			TVertexDescriptor nextVertex = findVertex(const_cast<TGraph&>(g), seqId, j);
			if (nextVertex == nilVertex) {
				++j;
				continue;
			}
			if (previousVertex != nilVertex) {
				_streamPutInt(file, previousVertex);
				_streamWrite(file, " -- ");
				_streamPutInt(file, nextVertex);
				_streamWrite(file, " [");
				_streamWrite(file, "weight=100.0, len=4.0, arrowhead=vee");
				_streamWrite(file, "];\n");
			}
			previousVertex = nextVertex;
			j += fragmentLength(g, nextVertex);
		}
	}
	_streamPut(file, '\n');
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
				DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
			   DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, " -- ");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
