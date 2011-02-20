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
  $Id: file.h 4762 2009-09-03 14:09:36Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

//____________________________________________________________________________
// prerequisites

#include <iostream>
#include <climits>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <cmath>

#include <seqan/sequence.h>
#include <seqan/modifier.h>


//____________________________________________________________________________

#include <seqan/file/file_forwards.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/file/file_generated_forwards.h>
#endif

#include <seqan/file/cstream.h>
#include <seqan/file/stream.h>

#include <seqan/file/chunk_collector.h>
#include <seqan/file/meta.h>

//____________________________________________________________________________
// files

#include <seqan/file/file_base.h>
#include <seqan/file/file_cstyle.h>
#include <seqan/file/file_array.h>

#include <seqan/system.h>	// async file (default file type of File<>)
/*#include <seqan/system/file_sync.h>
#include <seqan/system/system_event.h>
#include <seqan/system/file_async.h>
*/

//____________________________________________________________________________
// file formats

#include <seqan/file/file_filereaderiterator.h>
#include <seqan/file/file_filereader.h>

#include <seqan/file/file_format.h>

#include <seqan/file/stream_algorithms.h>

//file formats for sequences
#include <seqan/file/file_format_raw.h>
#include <seqan/file/file_format_fasta.h>
#include <seqan/file/file_format_embl.h>
#include <seqan/file/file_format_genbank.h>

//file formats for alignments
#include <seqan/file/file_format_fasta_align.h>
// TODO: include SAM file format specification
#include <seqan/file/file_format_sam.h>

//others
#include <seqan/file/file_format_cgviz.h>

//____________________________________________________________________________

//#include <seqan/file/file_format_guess.h>

//____________________________________________________________________________
// external strings

#include <seqan/file/file_page.h>
#include <seqan/file/file_page_raid0.h>
#include <seqan/file/string_external.h>
#include <seqan/file/string_mmap.h>
#include <seqan/file/file_format_mmap.h>

#endif //#ifndef SEQAN_HEADER_...
