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

#ifndef SEQAN_HEADER_BASIC_H
#define SEQAN_HEADER_BASIC_H

//____________________________________________________________________________
// prerequisites

#include <seqan/platform.h>

//#include <cstring>
#ifdef PLATFORM_WINDOWS
#include <limits>	// limits include file exists only for g++ >= 3.0
#endif

#include <cstddef>	// size_t
#include <cstdio>	// FILE, basic_debug
#include <cstdlib>	// posix_memalign
#include <ctime>
#include <iterator>
#include <algorithm>
#include <cstring>  // memset, memcpy
#include <string>	// basic_profile
#ifdef PLATFORM_WINDOWS
#include <malloc.h>	// _aligned_malloc
#endif  // PLATFORM_WINDOWS

#define SEQAN_NAMESPACE_MAIN seqan

//____________________________________________________________________________

#include <seqan/basic/basic_forwards.h>
#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/basic/basic_generated_forwards.h>
#endif

#include <seqan/basic/basic_testing.h>
#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_testing.h>  // new, better debug
#include <seqan/basic/basic_profile.h>
#include <seqan/basic/basic_parallelism.h>  // include after basic_testing.h!
#include <seqan/basic/basic_metaprogramming.h>
#include <seqan/basic/basic_definition.h>
#include <seqan/basic/basic_type.h>
#include <seqan/basic/basic_tag.h>

//____________________________________________________________________________
// allocators

#include <seqan/basic/basic_allocator_interface.h>
#include <seqan/basic/basic_allocator_to_std.h>

#include <seqan/basic/basic_holder.h>

#include <seqan/basic/basic_allocator_simple.h>
#include <seqan/basic/basic_allocator_singlepool.h>
#include <seqan/basic/basic_allocator_multipool.h>
#include <seqan/basic/basic_allocator_chunkpool.h>

//____________________________________________________________________________

#include <seqan/basic/basic_converter.h>
#include <seqan/basic/basic_compare.h>
#include <seqan/basic/basic_operator.h>

#include <seqan/basic/basic_host.h>

//____________________________________________________________________________
// iterators

#include <seqan/basic/basic_iterator.h>
#include <seqan/basic/basic_iterator_base.h>

#include <seqan/basic/basic_transport.h>

#include <seqan/basic/basic_iterator_simple.h>
#include <seqan/basic/basic_iterator_adaptor.h>
#include <seqan/basic/basic_iterator_position.h>
#include <seqan/basic/basic_iterator_adapt_std.h>
//#include <seqan/basic_identifier.h>

#include <seqan/basic/basic_proxy.h>

#include <seqan/basic/basic_pointer.h>

//____________________________________________________________________________
// alphabets

#include <seqan/basic/basic_alphabet_interface.h>
#include <seqan/basic/basic_alphabet_trait_basic.h>

#include <seqan/basic/basic_alphabet_interface2.h>

#include <seqan/basic/basic_alphabet_simple_tabs.h>
#include <seqan/basic/basic_alphabet_simple.h>

#include <seqan/basic/basic_sse2.h>

#include <seqan/basic/basic_profchar.h>

//____________________________________________________________________________

#include <seqan/basic/basic_holder_dynamic.h>

//____________________________________________________________________________

//#include <seqan/basic/basic_counted_ptr>
#include <seqan/basic/basic_volatile_ptr.h>

#include <seqan/basic/basic_aggregates.h>

#endif //#ifndef SEQAN_HEADER_...
