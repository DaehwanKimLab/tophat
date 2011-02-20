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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Drawing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// WRITING
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Directed<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Undirected<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Tree<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TNodeMap>
inline void
_createTrieNodeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
						  String<String<unsigned int> > pos,
						  TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeVertexMap(g, nodeMap);
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		String<char> tmp;
		std::stringstream s;
		s << *it;
		String<unsigned int> endPositions = getProperty(pos,*it);
		if (!empty(endPositions)) {
			s <<  " {";
			append(tmp, "shape = box, ");
			typename Iterator<String<unsigned int>, Rooted>::Type itP = begin(endPositions);
			typename Iterator<String<unsigned int>, Rooted>::Type beginP = itP;
			for(;!atEnd(itP);goNext(itP)) {
				if (beginP != itP) s << ", ";
				s << *itP;
			}
			s << "}";
		}
		
		append(tmp, "label = \"");
		append(tmp, s.str().c_str());
		append(tmp, "\"");
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap)
{
	SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << *it;
		outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes, typename TNameMap>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap,
					  TNameMap const& nameMap)
{
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << getProperty(nameMap,*it);
		outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TEdgeAttributes>
inline void
_createEmptyEdgeAttributes(Graph<TSpec> const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		assignProperty(edgeMap, *itEd, String<char>(""));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Directed<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Undirected<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Tree<void, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void 
_createEdgeAttributes(Graph<Tree<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << (TCargo) getCargo(*itEd);
		outs << "\"";
		append(property(edgeMap, *itEd), outs.str().c_str());		
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<char> tmp("label = \"");
		append(tmp, label(itEd));
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<TAlphabet> labelTmp = getCargo(*itEd);
		String<char> str;
		resize(str,length(labelTmp)+1);
		value(str,0) = label(itEd);
		typename Iterator<String<TAlphabet>, Rooted>::Type it = begin(labelTmp);
		for(;!atEnd(it);++it) {
			char c = convert<char>(getValue(it));
			value(str,position(it) + 1) = c;
		}
		String<char> tmp("label = \"");
		append(tmp, str);
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Directed<TCargo, TSpec> > const&,
				  DotDrawing)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Undirected<TCargo, TSpec> > const&,
				  DotDrawing)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Tree<TCargo, TSpec> > const&,
				  DotDrawing)
{
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				  DotDrawing)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Directed<TCargo, TSpec> > const&,
				DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Undirected<TCargo, TSpec> > const&,
				DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Tree<TCargo, TSpec> > const&,
				DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
			   DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Directed<TCargo, TSpec> > const&,
			   DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Undirected<TCargo, TSpec> > const&,
			   DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, " -- ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Tree<TCargo, TSpec> > const&,
			   DotDrawing)
{
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void 
write(TFile & file, 
	  Graph<TSpec> const& g,
	  TNodeAttributes const& nodeMap,
	  TEdgeAttributes const& edgeMap,
	  DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	_writeGraphType(file,g,DotDrawing());
	_streamWrite(file, " G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = circle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		_streamPutInt(file, *it);
		_streamWrite(file, " [");
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamWrite(file, "];\n");
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		_streamPutInt(file, sc);
		_writeEdgeType(file, g, DotDrawing());
		_streamPutInt(file, tr);
		_streamWrite(file, " [");
		_streamWrite(file, getProperty(edgeMap, *itEd));
		_streamWrite(file, "];\n");
	}
	_streamPut(file, '\n');

	_writeGraphFooter(file,g,DotDrawing());

	_streamWrite(file, "}\n");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNodeAttributes>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  TNodeAttributes const& nodeMap,
	  DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeAttributes(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}



/*

//////////////////////////////////////////////////////////////////////////////
// READING
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addNode(Graph<TSpec>& g,
		 TStatement& node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes&,			  
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	if (nodeIdMap.find(node_id) == nodeIdMap.end()) {
		TVertexDescriptor id = addVertex(g);
		nodeIdMap.insert(std::make_pair(node_id, id));
		resizeVertexMap(g, nodeMap);
		assignProperty(nodeMap, id, attr_list);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Directed<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Undirected<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline typename Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, TSpec> >&,
				  TString& str)
{
	return str[0];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline String<TAlphabet>
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > >&,
				  TString& str)
{
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// We need the label
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	typedef typename Position<TIter>::Type TPos;
	
	String<TValue> label;
	TIter it = begin(attr_list);
	bool found = false;
	for(;!atEnd(it);goNext(it)) {
		TPos pos = position(it);
		if (*it == ',') {
			found = false;
		} else if (found) {
			append(label, *it);
		} else if ((pos + 5 < length(attr_list)) &&
			(infix(attr_list, it, it + 5) == "label")) 
		{
				found = true;
				it += 5;
		}
	}
	TEdgeDescriptor e = addEdge(g, sourceV, targetV, _getInternalLabel(g, label));
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addEdge(Graph<TSpec>& g,
		 TStatement& left_node_id,
		 TStatement& right_node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes& edgeMap,
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TStatement>::Type TValue;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;

	TVertexDescriptor sourceV;
	TVertexDescriptor targetV;

	typename TMap::iterator pos;
	pos = nodeIdMap.find(left_node_id);
	if (pos == nodeIdMap.end()) return;
	else sourceV = pos->second;

	pos = nodeIdMap.find(right_node_id);
	if (pos == nodeIdMap.end()) return;
	else targetV = pos->second;

	_addEdge(g, sourceV, targetV, nodeMap, edgeMap, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processNodeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	for(;!atEnd(it);goNext(it)) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else {
			append(node_id, *it);
		}
	}
	_addNode(g, node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TPosition, typename TNodeIdMap>
inline void
_processEdgeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TPosition pos,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> left_node_id;
	String<TValue> right_node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	for(;!atEnd(it);goNext(it)) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else if (position(it) < pos) {
			append(left_node_id, *it);
		} else if (position(it) > pos+1) {
			append(right_node_id, *it);
		}
	}
	_addEdge(g, left_node_id, right_node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processStatement(Graph<TSpec>& g,
				  TStatement& stmt,
				  TNodeAttributes& nodeMap,
				  TEdgeAttributes& edgeMap,
				  TNodeIdMap& nodeIdMap) 
{
	// Clear everything up to the last line
	Finder<TStatement> finder(stmt);
	TStatement needle("\n");
	Pattern<TStatement, ShiftOr> pattern(needle);
	String<unsigned int> pos;
	while (find(finder, pattern))
		append(pos,position(finder));
	if (!empty(pos)) {
		stmt = suffix(stmt, pos[length(pos) - 1] + 1);
	}

	// Ignore all statements about attributes or subgraphs
	if ((prefix(stmt, 5) == "graph") ||
		(prefix(stmt, 4) == "node") ||
		(prefix(stmt, 4) == "edge") ||
		(prefix(stmt, 8) == "subgraph") ||
		(prefix(stmt, 1) == "\n") ||
		(prefix(stmt, 1) == "\r")) {
			clear(stmt);
			return;
	}

	// Node or Edge statement ?
	Finder<TStatement> finder2(stmt);
	needle = "--";
	setHost(pattern, needle);
	clear(pos);
	while (find(finder2, pattern))
		append(pos,position(finder2));
	if (!empty(pos)) {
		// Undirected Edge
		_processEdgeStatement(g, stmt, nodeMap, edgeMap, pos[0], nodeIdMap);
	} else {
		Finder<TStatement> finder3(stmt);
		needle = "->";
		setHost(pattern, needle);
		clear(pos);
		while (find(finder3, pattern))
			append(pos,position(finder3));
		if (!empty(pos)) {
			// Directed edge
			_processEdgeStatement(g, stmt, nodeMap, edgeMap, pos[0], nodeIdMap);
		} else {
			// process node statement
			_processNodeStatement(g, stmt, nodeMap, edgeMap, nodeIdMap);
		}
	}	
	clear(stmt);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void read(TFile & file,
		  Graph<TSpec>& g,
		  TNodeAttributes& nodeMap,
		  TEdgeAttributes& edgeMap,
		  DotDrawing) 
{
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;
	TMap nodeIdMap;

	TValue c;
	String<TValue> stmt;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		
		if (c == ';') _processStatement(g,stmt, nodeMap, edgeMap, nodeIdMap);
		else if ((c == '\n') ||
				(c == '\r')) {
					clear(stmt);
		}
		else append(stmt,c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec>
void read(TFile & file,
		  Graph<TSpec>& g,
		  DotDrawing) 
{
	String<String<char> > nodeMap;
	String<String<char> > edgeMap;
	read(file,g,nodeMap,edgeMap,DotDrawing());
}

*/

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
