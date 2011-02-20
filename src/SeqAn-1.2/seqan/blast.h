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
  $Id: blast.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BLAST_H
#define SEQAN_HEADER_BLAST_H

//____________________________________________________________________________
// prerequisites

#include <iostream>
#include <fstream>

#include <seqan/align.h>
#include <seqan/graph_types.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/graph_algorithms.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/blast/blast_generated_forwards.h>
#endif


#include <seqan/blast/blast_base.h>
#include <seqan/blast/blast_hsp.h>
#include <seqan/blast/blast_hit.h>
#include <seqan/blast/blast_report.h>
#include <seqan/blast/blast_parsing.h>
#include <seqan/blast/blast_iterator.h>
#include <seqan/blast/blast_hit_iterator.h>
#include <seqan/blast/blast_hsp_iterator.h>


#include <seqan/blast/blast_stream_report.h>
#include <seqan/blast/blast_stream_hit.h>
#include <seqan/blast/blast_stream_hit_iterator.h>
#include <seqan/blast/blast_stream_hsp_iterator.h>


#include <seqan/blast/blast_run.h>


#endif //#ifndef SEQAN_HEADER_...
