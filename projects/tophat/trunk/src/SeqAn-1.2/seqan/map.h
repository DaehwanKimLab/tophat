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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MAP_H
#define SEQAN_HEADER_MAP_H

//____________________________________________________________________________
// prerequisites

#include <set>
#include <map>

#include <seqan/sequence.h>
#include <seqan/misc/misc_random.h>


//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/map/map_generated_forwards.h>
#endif

//____________________________________________________________________________

#include <seqan/map/map_base.h>
#include <seqan/map/map_vector.h>
#include <seqan/map/map_skiplist.h>
#include <seqan/map/map_chooser.h>
#include <seqan/map/map_adapter_stl.h>

//____________________________________________________________________________

#include <seqan/map/sumlist.h>
#include <seqan/map/sumlist_mini.h>
#include <seqan/map/sumlist_skip.h>

//____________________________________________________________________________

#endif //#ifndef SEQAN_HEADER_...
