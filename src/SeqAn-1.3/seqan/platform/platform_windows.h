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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#define PLATFORM_WINDOWS
#define PLATFORM_WINDOWS_VS

// Disable warning "'function' : resolved overload was found by
// argument-dependent lookup".  Visual Studio warns because Koenig
// lookup was introduced in later version and behaviour has changed at some
// point.
#pragma warning( disable : 4675 )
// Disable warning for identifer name truncation.
#pragma warning( disable : 4503 )

// Disabling warning 4267 assigning variables with different size on 32 and 64 bit.
#pragma warning( disable : 4267 )
// Disabling warning 4244, loss of data when values with different domain sizes.
#pragma warning( disable : 4244 )

#define finline __forceinline

typedef unsigned __int64 __uint64;
typedef unsigned __int32 __uint32;
typedef unsigned __int16 __uint16;
typedef unsigned __int8 __uint8;

// The symbols SEQAN_IS_64_BIT and SEQAN_IS_32_BIT can be used to check
// whether we are on a 32 bit or on a 64 bit machine.
#if defined(_WIN64)
#define SEQAN_IS_64_BIT 1
#define SEQAN_IS_32_BIT 0
#else
#define SEQAN_IS_64_BIT 0
#define SEQAN_IS_32_BIT 1
#endif  // #if defined(_WIN64)

// Workaround for missing round() from C99 in Visual Studio.
template <typename T>
inline T round(T const & x)
{
	return floor(x + 0.5);
}

// Rename some underscore-functions in Windows.
#ifndef snprintf
#define snprintf _snprintf
#endif  // #ifndef snprintf

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
//#define SEQAN_SWITCH_USE_FORWARDS
