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
  $Id: consensus.h 1901 2008-04-28 13:07:56Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_H
#define SEQAN_HEADER_CONSENSUS_H

// Seqan
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/refinement.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>
#include <seqan/store.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/consensus/consensus_generated_forwards.h>
#endif

// Consensus tool
#include <seqan/consensus/consensus_base.h>
#include <seqan/consensus/consensus_score.h>
#include <seqan/consensus/consensus_realign.h>
#include <seqan/consensus/consensus_library.h>

#endif //#ifndef SEQAN_HEADER_...
