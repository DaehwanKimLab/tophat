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
  $Id: refinement.h 1901 2008-04-28 13:07:56Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_REFINEMENT_H
#define SEQAN_HEADER_REFINEMENT_H

// External STL
#include <map>

// Seqan
#include <seqan/score.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/refinement/refinement_generated_forwards.h>
#endif

// Alignment graph
#include <seqan/refinement/graph_impl_align.h>
#include <seqan/refinement/graph_impl_align_adapt.h>

// Intervall trees
#include <seqan/refinement/graph_impl_interval_types.h>
#include <seqan/refinement/graph_impl_interval_tree.h>

// Refinement
//#include <seqan/refinement/graph_algorithm_refine.h>
#include <seqan/refinement/graph_algorithm_refine_scoring.h>
#include <seqan/refinement/graph_algorithm_refine_fragment.h>
#include <seqan/refinement/graph_algorithm_refine_aligngraph.h>
#include <seqan/refinement/graph_algorithm_refine_align.h>
#include <seqan/refinement/graph_algorithm_refine_annotation.h>
//#include <seqan/refinement/graph_algorithm_refine_exact.h>
#include <seqan/refinement/graph_algorithm_refine_exact_iterative.h>
#include <seqan/refinement/graph_algorithm_refine_inexact.h>

#endif //#ifndef SEQAN_HEADER_...
