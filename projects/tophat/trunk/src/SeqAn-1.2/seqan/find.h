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
  $Id: find.h 2619 2008-08-30 20:00:10Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_H
#define SEQAN_HEADER_FIND_H

//____________________________________________________________________________
// prerequisites

#include <cmath>

#include <deque>

#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/score.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/map.h>

//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/find/find_generated_forwards.h>
#endif

#include <seqan/find/find_base.h>
#include <seqan/find/find_pattern_base.h>

//____________________________________________________________________________
// exact pattern matching
#include <seqan/find/find_simple.h>
#include <seqan/find/find_horspool.h>
#include <seqan/find/find_shiftand.h>
#include <seqan/find/find_shiftor.h>
#include <seqan/find/find_bndm.h>
#include <seqan/find/find_bom.h>

//____________________________________________________________________________
//complex pattern matching
#include <seqan/find/find_wild_shiftand.h>

//____________________________________________________________________________
//multiple pattern search
#include <seqan/find/find_ahocorasick.h>
#include <seqan/find/find_multiple_shiftand.h>
#include <seqan/find/find_set_horspool.h>

//#include <seqan/find/find_multi.h> //wegwerfen
#include <seqan/find/find_wumanber.h> 
#include <seqan/find/find_multiple_bfam.h>

//____________________________________________________________________________
// approximate pattern matching
#include <seqan/find/find_begin.h>

#include <seqan/find/find_score.h>
#include <seqan/find/find_myers_ukkonen.h>
#include <seqan/find/find_abndm.h>
#include <seqan/find/find_pex.h>

//#include <seqan/find/find_hamming_horspool.h>

//____________________________________________________________________________
//filtering
#include <seqan/find/find_quasar.h>

//____________________________________________________________________________

#endif //#ifndef SEQAN_HEADER_...
