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
  $Id: graph_impl_trie.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Trie
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Trie:
..cat:Graph
..general:Class.Graph
..summary:A keyword trie.
..description:
...image:trieGraph|A trie for the words announce, annual, and annually.
..remarks:A keyword trie is a special automaton and thus, it is not implemented in its own class.
It solely provides create functions where based upon a set of strings a keyword trie is created.
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

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeyword, typename TPos>
inline void
_addStringToTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				 TTerminalStateMap& terminalStateMap,
				 TKeyword const& str,
				 TPos const& keywordIndex)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TKeyword>::Type TSize;

	TVertexDescriptor current = getRoot(g);
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	typename Iterator<TKeyword const, Rooted>::Type sIt = begin(str);
	for(;!atEnd(sIt);goNext(sIt)) {
		if (getSuccessor(g, current, *sIt) == nilVal) break;
		current = getSuccessor(g, current, *sIt);
	}
	for(;!atEnd(sIt);goNext(sIt)) {
		TVertexDescriptor newState = addVertex(g);
		resize(terminalStateMap, numVertices(g), Generous());
		assignProperty(terminalStateMap,newState,String<TPos>());
		addEdge(g,current,newState,*sIt);
		current = newState;
	}
	String<TPos> tmp = getProperty(terminalStateMap,current);
	appendValue(tmp, keywordIndex);
	assignProperty(terminalStateMap,current,tmp);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.createTrie:
..cat:Graph
..summary:Creates a trie.
..signature:createTrie(g, terminalStateMap, keywords)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. This is the case in the Aho-Corasick
algorithm if one pattern is a suffix of another pattern! Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.keywords:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrieOnReverse
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		   TTerminalStateMap& terminalStateMap,
		   TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Position<TKeywords>::Type TPos;
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPos>());
	typename Iterator<TKeywords const, Rooted>::Type it = begin(keywords);
	for(;!atEnd(it);goNext(it)) _addStringToTrie(g,terminalStateMap,*it,position(it));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.createTrieOnReverse:
..cat:Graph
..summary:Creates a trie for all reversed keywords.
..signature:createTrieOnReverse(g, terminalStateMap, keywords)
..returns.param.g:Out-parameter: An automaton.
...type:Spec.Trie
..returns.param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. This is the case in the Aho-Corasick
algorithm if one pattern is a suffix of another pattern! Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.keywords:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrie
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrieOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
					TTerminalStateMap& terminalStateMap,
					TKeywords const& keywords)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Position<TKeywords>::Type TPos;
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPos>());
	typename Iterator<TKeywords const, Rooted>::Type it = begin(keywords);
	for(;!atEnd(it);goNext(it)) {
		typedef typename Value<TKeywords>::Type TKeyword;
		TKeyword tmp;
		typename Iterator<TKeyword const, Rooted>::Type sIt = end(*it);
		while(!atBegin(sIt)) {
			goPrevious(sIt);
			appendValue(tmp,getValue(sIt));
		}
		_addStringToTrie(g,terminalStateMap,tmp,position(it));
	}
}



/**
.Function.createSuffixTrie:
..cat:Graph
..summary:Creates a trie of all suffixes of a text.
..signature:createSuffixTrie(g, terminalStateMap, text)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> >.
..param.text:In-parameter: A text.
...type:Class.String
..returns:void
..see:Function.createTrie
..see:Function.createTrieOnReverse
*/
template <typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TText>
inline void
createSuffixTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				 TTerminalStateMap& terminalStateMap,
				 TText const& text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TText const>::Type TPosition;
	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPosition>());

	for (TPosition i = 0; i < length(text); ++i)
	{
		_addStringToTrie(g,terminalStateMap,suffix(text, i),i);
	}
}

//////////////////////////////////////////////////////////////////////////////


template <typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TTexts>
inline void
createSetSuffixTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
					TTerminalStateMap& terminalStateMap,
					TTexts const& texts)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TTexts const>::Type TTextsPosition;
	typedef typename Value<TTexts const>::Type TText;
	typedef typename Position<TText const>::Type TPosition;

	TVertexDescriptor root = addVertex(g);
	assignRoot(g,root);
	resize(terminalStateMap, numVertices(g), Generous());
	assignProperty(terminalStateMap,root,String<TPosition>());

	for (TTextsPosition j = 0; j < length(texts); ++j)
	{
		TText const & text = texts[j];
		for (TPosition i = 0; i < length(text); ++i)
		{
			_addStringToTrie(g,terminalStateMap,suffix(text, i),j);
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
