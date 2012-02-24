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

#ifndef SEQAN_HEADER_BASIC_POINTER_H
#define SEQAN_HEADER_BASIC_POINTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Adaption.char array

template <typename TValue>
struct Value< TValue * >
{
	typedef TValue Type;
};
template <typename TValue>
struct Value< TValue * const>
{
	typedef TValue Type;
};

//The next two metafunctions dont work in VC++ due to a compiler bug.
//(the default implementation in common_type.h is called instead)
//work-around: convert arrays to pointers.
template <typename TValue, size_t SIZE>
struct Value< TValue [SIZE] >
{
	typedef TValue Type;
};
template <typename TValue, size_t SIZE>
struct Value< TValue const [SIZE] >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Adaption.char array

template <typename TValue>
struct Iterator< TValue *, Standard>
{
	typedef TValue * Type;
};
template <typename TValue>
struct Iterator< TValue * const, Standard>
{
	typedef TValue * Type;
};

//____________________________________________________________________________

template <typename TValue, size_t SIZE>
struct Iterator< TValue [SIZE], Standard>:
	Iterator<TValue *, Standard>
{
};
template <typename TValue, size_t SIZE>
struct Iterator< TValue const [SIZE], Standard>:
	Iterator<TValue const *, Standard>
{
};

template <typename TValue, size_t SIZE>
struct Iterator< TValue [SIZE], Rooted>:
	Iterator<TValue *, Rooted>
{
};
template <typename TValue, size_t SIZE>
struct Iterator< TValue const [SIZE], Rooted>:
	Iterator<TValue const *, Rooted>
{
};




//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
