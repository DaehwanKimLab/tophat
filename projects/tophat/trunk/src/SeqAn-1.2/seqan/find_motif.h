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
  $Id:
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_MOTIF_H
#define SEQAN_HEADER_FIND_MOTIF_H

//____________________________________________________________________________
// prerequisites

#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <cfloat>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>

//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/find_motif/find_motif_generated_forwards.h>
#endif

//____________________________________________________________________________

#include <seqan/find_motif/sequence_model_types.h>

#include <seqan/find_motif/pseudocount_base.h>
#include <seqan/find_motif/pseudocount_mode_c.h>
#include <seqan/find_motif/pseudocount_mode_p.h>

#include <seqan/find_motif/frequency_distribution.h>
#include <seqan/find_motif/profile.h>

#include <seqan/find_motif/find_motif_base.h>
#include <seqan/find_motif/find_motif_pms1.h>
#include <seqan/find_motif/find_motif_pmsp.h>
#include <seqan/find_motif/find_motif_projection.h>
#include <seqan/find_motif/find_motif_epatternbranching.h>
#include <seqan/find_motif/em_algorithm.h>


#endif //#ifndef SEQAN_HEADER_...
