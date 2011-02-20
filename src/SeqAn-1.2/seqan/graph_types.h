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
  $Id: graph_types.h 1901 2008-04-28 13:07:56Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_TYPES_H
#define SEQAN_HEADER_GRAPH_TYPES_H

// External / STL
#include <deque>


// Seqan
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_types/graph_types_generated_forwards.h>
#endif


// Basic graph stuff
#include <seqan/graph_types/graph_base.h>
#include <seqan/graph_types/graph_idmanager.h>	// Id manager
#include <seqan/graph_types/graph_edgestump.h>	// EdgeStumps
#include <seqan/graph_types/graph_interface.h>	// Graph metafunctions

// Graph types
#include <seqan/graph_types/graph_impl_directed.h>		// Directed Graph
#include <seqan/graph_types/graph_impl_undirected.h>	// Undirected graph
#include <seqan/graph_types/graph_impl_automaton.h>		// Automaton
#include <seqan/graph_types/graph_impl_wordgraph.h>		// Specialized automaton: Word graph
#include <seqan/graph_types/graph_impl_tree.h>			// Tree
#include <seqan/graph_types/graph_impl_fragment.h>		// Fragment
#include <seqan/graph_types/graph_impl_hmm.h>			// HMM

// Graph iterators
#include <seqan/graph_types/graph_iterator.h>
#include <seqan/graph_types/graph_iterator_vertex.h>
#include <seqan/graph_types/graph_iterator_outedge.h>
#include <seqan/graph_types/graph_iterator_adjacency.h>
#include <seqan/graph_types/graph_iterator_edge.h>

// Graph property maps
#include <seqan/graph_types/graph_property.h>

// Specializations
#include <seqan/graph_types/graph_impl_oracle.h>	// Oracle
#include <seqan/graph_types/graph_impl_trie.h>		// Trie

// Specialized iterators
#include <seqan/graph_types/graph_iterator_bfs.h>
#include <seqan/graph_types/graph_iterator_dfs.h>

// Graph drawing and some file parsing
#include <seqan/graph_types/graph_drawing.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/graph_types/graph_utility_parsing.h>

#endif //#ifndef SEQAN_HEADER_...
