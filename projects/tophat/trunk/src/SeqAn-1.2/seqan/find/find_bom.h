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
  $Id: find_bom.h 2653 2008-09-05 10:34:13Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BOM_H
#define SEQAN_HEADER_FIND_BOM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BomAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BFAM:
..summary:Backward Factor Automaton Matching algorithm.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, BFAM<TAutomaton> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TAutomaton:A tag that specifies the used automaton.
...default:@Spec.BFAM<Oracle>@
..remarks.text:To be used in combination with the default specialization of @Class.Finder@.
*/

/**
.Spec.BFAM<Oracle>:
..summary:Backward Oracle Matching algorithm.
..general:Class.BFAM
..cat:Searching
..signature:Pattern<TNeedle, BFAM<Oracle> >
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:To be used in combination with the default specialization of @Class.Finder@.
*/
/**
.Spec.BFAM<Trie>:
..summary:Backward Suffix Trie Matching algorithm.
..general:Class.BFAM
..cat:Searching
..signature:Pattern<TNeedle, BFAM<Trie> >
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:To be used in combination with the default specialization of @Class.Finder@.
*/

///.Class.Pattern.param.TSpec.type:Spec.BFAM

struct Oracle; //Oracle Tag => "BOM"
struct Trie; //Trie Tag => "BTM"

template <typename TSpec = Oracle>
struct BFAM; //backward factor automaton searching

typedef BFAM<Oracle> BomAlgo; //deprecated, still there for compatibility reasons

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
class Pattern<TNeedle, BFAM<TSpec> > {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef typename Size<TNeedle>::Type TSize;
	Holder<TNeedle> data_host;
	TSize needleLength;		
	TSize haystackLength;
	TSize step;
	Graph<Automaton<TAlphabet, void, WithoutEdgeId> > automaton;

//____________________________________________________________________________

	Pattern() {	
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
SEQAN_CHECKPOINT
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

//BFAM<Oracle>: BOM Algorithm
template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BFAM<Oracle> > & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	clear(me.automaton);
	createOracleOnReverse(me.automaton,needle);
	setValue(me.data_host, needle);
}
//BFAM<Trie>: BTM Algorithm (the same as BOM, but with an trie)
template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BFAM<Trie> > & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	clear(me.automaton);

	String<String<unsigned int> > terminal_state_map; //dummy
	typedef typename Value<TNeedle2 const>::Type TValue;
	String<TValue> reverse_string = needle;
	reverseInPlace(reverse_string);

	createSuffixTrie(me.automaton, terminal_state_map, reverse_string);

	setValue(me.data_host, needle);
}

template <typename TNeedle, typename TNeedle2, typename TSpec>
inline void 
setHost (Pattern<TNeedle, BFAM<TSpec> > & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle, typename TSpec>
inline void _patternInit (Pattern<TNeedle, BFAM<TSpec> > & me) 
{
SEQAN_CHECKPOINT
	me.step = 0;
}


//____________________________________________________________________________


template <typename TFinder, typename TNeedle, typename TSpec>
inline bool 
find(TFinder & finder, Pattern<TNeedle, BFAM<TSpec> > & me) 
{
	SEQAN_CHECKPOINT
	
	if (empty(finder)) {
		_patternInit(me);
		_setFinderLength(finder, length(needle(me)));
		_finderSetNonEmpty(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.step;

	if (me.haystackLength < me.needleLength) return false;
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TOracle;
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename VertexDescriptor<TOracle>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TOracle>::Type TEdgeDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	while (position(finder) <= me.haystackLength - me.needleLength) {
		TVertexDescriptor current = getRoot(me.automaton);
		TSize j = me.needleLength;
		while ((j>0) &&	(current != nilVal))
		{
			TAlphabet c = *(finder+(j-1));
			current = targetVertex(me.automaton, findEdge(me.automaton, current, c));
			--j;
		}
		if (current != nilVal) {
			me.step = j + 1;
			_setFinderEnd(finder, position(finder) + me.needleLength);
			return true;
		}
		finder += j + 1;
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
