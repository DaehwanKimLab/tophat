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
  $Id: graph_impl_oracle.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_ORACLE_H
#define SEQAN_HEADER_GRAPH_IMPL_ORACLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Oracle
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Oracle:
..cat:Graph
..general:Class.Graph
..summary:A factor oracle.
..remarks:A factor oracle is a special automaton and thus, it is not implemented in its own class.
It solely provides create functions where based upon a string an oracle is created.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
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

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TPropertyMap, typename TChar>
inline void
_addLetterToOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				   TPropertyMap& supplyState,
				   TChar const c)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor newState = addVertex(g);
	TVertexDescriptor pred = newState - 1;
	addEdge(g, pred, newState, c);
	TVertexDescriptor k = getProperty(supplyState, pred);
	while ((k!=nilVal) &&
			(getTarget(&g.data_vertex[k].data_edge[ordValue(TAlphabet(c))])==nilVal))
	{
		addEdge(g,k,newState,c);
		k = getProperty(supplyState, k);
	}
	TVertexDescriptor s;
	if (k==nilVal) s=0;
	else s = getTarget(&g.data_vertex[k].data_edge[ordValue(TAlphabet(c))]);
	assignProperty(supplyState, newState, s);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.createOracle:
..cat:Graph
..summary:Creates a factor oracle.
..signature:createOracle(g, text)
..param.g:Out-parameter: An oracle.
...type:Spec.Oracle
..param.text:In-parameter: A string.
...type:Class.String
..returns:void
..see:Function.createOracleOnReverse
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			 TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignRoot(g,v1);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = 0; i<len; ++i) _addLetterToOracle(g, supplyState, getValue(text,i));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.createOracleOnReverse:
..cat:Graph
..summary:Creates a factor oracle for the reversed string.
..signature:createOracleOnReverse(g, text)
..param.g:Out-parameter: An oracle.
...type:Spec.Oracle
..param.text:In-parameter: A string.
...type:Class.String
..returns:void
..see:Function.createOracle
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracleOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
					  TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignRoot(g,v1);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = len-1; i>0; --i) _addLetterToOracle(g, supplyState, getValue(text,i));
	_addLetterToOracle(g, supplyState, getValue(text,0));
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createSetOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				TTerminalStateMap& terminalStateMap,
				TKeywords const& keywords)
{
SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Position<TKeywords>::Type TPos;
	typedef typename Value<TKeywords>::Type TKeyword;
	typedef typename Iterator<TKeyword, Standard>::Type TIterator;
	typedef typename Value<TKeywords>::Type TValue;

	createTrie(g, terminalStateMap, keywords);

	String<TVertexDescriptor> supplyState;
	resizeVertexMap(g, supplyState);
	String<bool> visited;
	resizeVertexMap(g, visited);
	arrayFill(begin(visited), end(visited), false);

	TVertexDescriptor nil = getNil<TVertexDescriptor>();
	assignProperty(supplyState, root(g), nil);

	TVertexDescriptor _root = getRoot(g);

	TPos len = length(keywords);
	String<TVertexDescriptor> _here_v;
	fill(_here_v, len, _root);
	String<TIterator> _here_it;
	resize(_here_it, len);
	for (TPos i = 0; i < len; ++i)
	{
		_here_it[i] = begin(keywords[i], Standard());
	}
	TPos _active_count = len;
	while (_active_count)
	{
		for (TPos i = 0; i < len; ++i)
		{
			TIterator & it = _here_it[i];
			TIterator it_end = end(keywords[i], Standard());
			TVertexDescriptor & _parent = _here_v[i];
			
			if (it != it_end)
			{
				TVertexDescriptor _current = getSuccessor(g, _parent, *it);

				if (!getProperty(visited, _current))
				{
					assignProperty(visited, _current, true);

					TVertexDescriptor _down = getProperty(supplyState, _parent);
					TVertexDescriptor _supply = _root;
					while (_down != nil)
					{
						TVertexDescriptor _next = getSuccessor(g, _down, *it);
						if (_next != nil)
						{
							_supply = _next;
							break;
						}

						addEdge(g, _down, _current, *it);
						_down = getProperty(supplyState, _down);
					}
					assignProperty(supplyState, _current, _supply);
				}
				_parent = _current;
				++it;
				if (it == it_end)
				{
					--_active_count;
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
