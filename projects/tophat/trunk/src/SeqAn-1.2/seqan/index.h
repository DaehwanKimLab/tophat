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
  $Id: index.h 2653 2008-09-05 10:34:13Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_H
#define SEQAN_HEADER_INDEX_H

//____________________________________________________________________________
// prerequisites

#include <seqan/sequence.h>
#include <seqan/pipe.h>
#include <seqan/modifier.h>

#include <seqan/find.h>
#include <seqan/misc/misc_set.h>

#include <climits>
#include <functional>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <string.h> // memset


//////////////////////////////////////////////////////////////////////////////
// INDEX CONSTRUCTION
//////////////////////////////////////////////////////////////////////////////


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/index/index_manual_forwards.h>
#include <seqan/index/index_generated_forwards.h>
#endif

#include <seqan/index/index_base.h>

//____________________________________________________________________________
// suffix array creators

#include <seqan/index/radix.h>
#include <seqan/index/index_sa_btree.h>
#include <seqan/index/index_sa_lss.h>
#include <seqan/index/index_sa_mm.h>
#include <seqan/index/index_sa_qsort.h>

#include <seqan/index/pump_extender3.h>
#include <seqan/index/pipe_merger3.h>
#include <seqan/index/index_skew3.h>

#include <seqan/index/pump_extender7.h>
#include <seqan/index/pipe_merger7.h>
#include <seqan/index/index_skew7.h>

#include <seqan/index/pump_separator7.h>
#include <seqan/index/index_skew7_multi.h>

//____________________________________________________________________________
// enhanced table creators

#include <seqan/index/pump_lcp_core.h>
#include <seqan/index/index_lcp.h>
#include <seqan/index/index_lcp_tree.h>

#include <seqan/index/index_childtab.h>
#include <seqan/index/index_bwt.h>

//____________________________________________________________________________
// q-gram index creator

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>
#include <seqan/index/shape_onegapped.h>
#include <seqan/index/shape_predefined.h>
#include <seqan/index/index_qgram.h>
//#include <seqan/index/index_qgram_nested.h>


//////////////////////////////////////////////////////////////////////////////
// INDEX USAGE
//////////////////////////////////////////////////////////////////////////////

#include <seqan/index/index_shims.h>
#include <seqan/index/index_esa_find.h>
#include <seqan/index/index_qgram_find.h>
#include <seqan/index/index_find.h>

//____________________________________________________________________________
// (virtual) suffix trees

#include <seqan/index/index_esa_base.h>
#include <seqan/index/index_esa_stree.h>
#include <seqan/index/index_wotd.h>
#include <seqan/index/index_dfi.h>

//____________________________________________________________________________
// suffix tree algorithms

#include <seqan/index/index_esa_algs.h>
#include <seqan/index/index_esa_algs_multi.h>
#include <seqan/index/index_esa_drawing.h>
#include <seqan/index/repeat_base.h>

//____________________________________________________________________________
// Pizza & Chili interface (compressed indices)

#include <seqan/index/index_pizzachili.h>
#include <seqan/index/index_pizzachili_find.h>

//____________________________________________________________________________
// Shawarma interface (suffix array creators)

#include <seqan/index/index_shawarma.h>

#endif //#ifndef SEQAN_HEADER_...
