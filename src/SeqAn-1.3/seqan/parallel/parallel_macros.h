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
// Utility macros for parallelism.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_MACROS_H_
#define SEQAN_PARALLEL_PARALLEL_MACROS_H_

// Macro: SEQAN_OMP_PRAGMA(x)
//
// Example: SEQAN_OMP_PRAGMA(omp parallel for)
//
// If OpenMP is enabled, this macro expands to "#pragma + $x".  This is useful
// for disabling OpenMP pragmas on compilers that do not support OpenMP to
// suppress warnings.
#ifdef _OPENMP
  #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
    // GCC _Pragma operator
    #define SEQAN_OMP_PRAGMA(x) _Pragma (#x)
  #else  // #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
    // MSVC __pragma-operator
    #define SEQAN_OMP_PRAGMA(x) __pragma (x)
  #endif // #if defined(PLATFORM_WINDOWS_MINGW) || defined(PLATFORM_GCC)
#else  // #ifdef _OPENMP
  #define SEQAN_OMP_PRAGMA(x)
#endif  // #ifdef _OPENMP

#endif  // SEQAN_PARALLEL_PARALLEL_MACROS_H_
