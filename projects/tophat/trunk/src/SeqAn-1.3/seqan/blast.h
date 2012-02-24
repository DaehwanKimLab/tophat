// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

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
