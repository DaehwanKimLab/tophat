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
  $Id: chaining.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CHAINING_H
#define SEQAN_HEADER_CHAINING_H

#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/seeds.h>

#include <seqan/misc/misc_random.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/chaining/chaining_generated_forwards.h>
#endif

#include <seqan/chaining/score_chain.h>
#include <seqan/chaining/score_manhattan.h>
#include <seqan/chaining/score_zero.h>
#include <seqan/chaining/score_chain_sop.h>

#include <algorithm>
#include <vector>

#include <seqan/chaining/range_max_tree.h>
#include <seqan/chaining/chain_base.h>
#include <seqan/chaining/fragment.h>
#include <seqan/chaining/chain_point.h>
#include <seqan/chaining/chain_wrapper_point.h>
#include <seqan/chaining/chain_meta_fragment.h>
#include <seqan/chaining/tree_chain_utils.h>
#include <seqan/chaining/tree_chain.h>
#include <seqan/chaining/tree_chain_sop.h>


#include <seqan/chaining/chain_generic.h>



namespace seqan{

}

#endif

