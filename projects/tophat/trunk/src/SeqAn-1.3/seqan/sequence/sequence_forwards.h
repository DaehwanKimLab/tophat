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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Manual forwards for the sequence module.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_FORWARDS_H
#define SEQAN_HEADER_SEQUENCE_FORWARDS_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

// TODO(holtgrew): Still required since we dropped support for VS2003?
// Workaround (copied from generated forwards) for VS 2003.
#if defined(_MSC_VER) && (_MSC_VER < 1400)
template <unsigned int SPACE > struct Block;        // "projects\library\seqan/sequence\string_stack.h"(48)
template <typename THostspec > struct Packed;           // "projects\library\seqan/sequence\string_packed.h"(33)
template <typename TValue, typename TSpec > class String;           // "projects\library\seqan/sequence\string_base.h"(54)
template <typename TString, typename TSpec > class StringSet;           // "projects\library\seqan/sequence\sequence_multiple.h"(98)

template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > & me, Tag<TTag> const tag_);           // "projects\library\seqan/sequence\string_packed.h"(470)
template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > const & me, Tag<TTag> const tag_);           // "projects\library\seqan/sequence\string_packed.h"(478)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type end(String<TValue, Block<SPACE> > & me, Tag<TSpec> const);           // "projects\library\seqan/sequence\string_stack.h"(209)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type end(String<TValue, Block<SPACE> > const & me, Tag<TSpec> const);        // "projects\library\seqan/sequence\string_stack.h"(217)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec >, Tag<TTag> const>::Type end(StringSet< TString, TSpec > & me, Tag<TTag> const tag);        // "projects\library\seqan/sequence\sequence_multiple.h"(1398)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec > const, Tag<TTag> const>::Type end(StringSet< TString, TSpec > const & me, Tag<TTag> const tag);        // "projects\library\seqan/sequence\sequence_multiple.h"(1405)
#endif  // defined(_MSC_VER) && (_MSC_VER < 1400)

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
void read(TFile & file, TData & data);          // "projects/library/seqan/file/file_format_raw.h"(307)

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
void write(TFile & file, TData & data);         // "projects/library/seqan/file/file_format_raw.h"(327)

template <typename TFile, typename TData>
void write(TFile & file, TData const & data);   // "projects/library/seqan/file/file_format_raw.h"(335)


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif

