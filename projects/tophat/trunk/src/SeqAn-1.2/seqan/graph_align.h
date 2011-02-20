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
  $Id: graph_align.h 1901 2008-04-28 13:07:56Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_H
#define SEQAN_HEADER_GRAPH_ALIGN_H

// Seqan
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/refinement.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_align/graph_align_generated_forwards.h>
#endif

// Alignment
#include <seqan/graph_align/graph_align_base.h>
#include <seqan/graph_align/graph_align_config.h>
#include <seqan/graph_align/graph_align_interface.h>
#include <seqan/graph_align/graph_align_needleman_wunsch.h>
#include <seqan/graph_align/graph_align_banded_needleman_wunsch.h>
#include <seqan/graph_align/graph_align_gotoh.h>
#include <seqan/graph_align/graph_align_banded_gotoh.h>
#include <seqan/graph_align/graph_align_hirschberg.h>
#include <seqan/graph_align/graph_align_smith_waterman.h>
#include <seqan/graph_align/graph_align_smith_waterman_clump.h>


#endif //#ifndef SEQAN_HEADER_...
