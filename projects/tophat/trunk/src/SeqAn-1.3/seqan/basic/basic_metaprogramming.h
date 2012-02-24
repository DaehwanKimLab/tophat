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

#ifndef SEQAN_BASIC_METAPROGRAMMING_H
#define SEQAN_BASIC_METAPROGRAMMING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
/**
.Tag.Logical Values:
..summary:Tag that represents true and false.
..tag.True:The logical value "true".
..tag.False:The logical value "false".
..include:seqan/basic.h
*/

	struct True { enum { VALUE = true }; };
	struct False { enum { VALUE = false }; };

	template <bool b>
	struct Eval
	{
		typedef False Type;
	};

	template <>
	struct Eval<true>
	{
		typedef True Type;
	};

	//////////////////////////////////////////////////////////////////////////////
	// generic "or" (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template <typename TBool1, typename TBool2>
	struct Or
	{
		typedef True Type;
	};

	template <>
	struct Or<False, False>
	{
		typedef False Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic "and" (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template <typename TBool1, typename TBool2>
	struct And
	{
		typedef False Type;
	};

	template <>
	struct And<True, True>
	{
		typedef True Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic "if" (using meta-programming)
	// if Flag is true,  the resulting type is Type1
	// if Flag is false, the resulting type is Type2
	//////////////////////////////////////////////////////////////////////////////

	template <bool Flag,class Type1, class Type2>
	struct If
	{
		typedef Type1 Type;
	};

	template <class Type1, class Type2>
	struct If<false,Type1,Type2>
	{
		typedef Type2 Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic type comparison (using meta-programming)
	// if Type1 equals Type2,		VALUE is true
	// if Type1 differs from Type2, VALUE is false
	//////////////////////////////////////////////////////////////////////////////

	template <class Type1, class Type2>
	struct IsSameType
	{
		typedef False Type;
		enum { VALUE = false };
	};

	template <class Type1>
	struct IsSameType<Type1, Type1>
	{
		typedef True Type;
		enum { VALUE = true };
	};

	//////////////////////////////////////////////////////////////////////////////
	// generic "switch" (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	const int DEFAULT = ~(~0u >> 1); // initialize with the smallest int

	struct NilCase {};
	  
	template <int tag_,class Type_,class Next_ = NilCase>
	struct Case
	{
		enum { tag = tag_ };
		typedef Type_ Type;
		typedef Next_ Next;
	};

	template <int tag,class Case>
	class Switch
	{
		typedef typename Case::Next NextCase;
		enum
		{
			caseTag = Case::tag,
			found   = (caseTag == tag || caseTag == DEFAULT)
		};
	public:
		typedef typename 
			If<
				found,
				typename Case::Type,
				typename Switch<tag,NextCase>::Type
			>::Type Type;
	};

	template <int tag>
	class Switch<tag,NilCase>
	{
	public:
		typedef NilCase Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// generic loops (using meta-programming)
	// corresponds to for(i=1; i<=I; ++i) ...
	//////////////////////////////////////////////////////////////////////////////

	// example of a loop Worker class
	struct WorkerNothing
	{
		template <typename Arg>
		static inline void body(Arg &/*arg*/, int /*I*/) {}
	};

	template <typename Worker, int I>
	class Loop {
	public:
		template <typename Arg>
		static inline void run(Arg &arg) {
			Loop<Worker, I - 1>::run(arg);
			Worker::body(arg, I);
		}
	};

	template <typename Worker>
	class Loop<Worker, 0> {
	public:
		// end of loop
		template <typename Arg>
		static inline void run(Arg &) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// generic reverse loops (using meta-programming)
	// corresponds to for(i=I; i>0; --i) ...
	//////////////////////////////////////////////////////////////////////////////

	template <typename Worker, int I>
	class LoopReverse {
	public:
		template <typename Arg>
		static inline void run(Arg &arg) {
			Worker::body(arg, I);
			LoopReverse<Worker, I - 1>::run(arg);
		}
	};

	template <typename Worker>
	class LoopReverse<Worker, 0> {
	public:
		// end of loop
		template <typename Arg>
		static inline void run(Arg &) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// logarithmus dualis (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template < __int64 numerus >
	struct Log2 {
		static const __uint64 VALUE = Log2<(numerus + 1) / 2>::VALUE + 1;		// ceil(log_2(n))
	};

	template < __int64 numerus >
	struct Log2Floor {
		static const __uint64 VALUE = Log2Floor<numerus / 2>::VALUE + 1;		// floor(log_2(n))
	};

	template <> struct Log2<1> { static const __uint64 VALUE = 0; };
	template <> struct Log2<0> { static const __uint64 VALUE = 0; };
	template <> struct Log2Floor<1> { static const __uint64 VALUE = 0; };
	template <> struct Log2Floor<0> { static const __uint64 VALUE = 0; };


	//////////////////////////////////////////////////////////////////////////////
	// exponentiation (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template < __int64 base, __int64 exponent >
	struct Power {
		static const __uint64 VALUE =
				Power<base, exponent / 2>::VALUE * 
				Power<base, exponent - (exponent / 2)>::VALUE;
	};

	template < __int64 base > struct Power<base, 1> { static const __uint64 VALUE = base; };
	template < __int64 base > struct Power<base, 0> { static const __uint64 VALUE = 1; };


	//////////////////////////////////////////////////////////////////////////////
	// memset with fill size (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	using ::std::memset;

	template <unsigned SIZE, bool direct>
	struct MemsetWorker {
		finline static void run(unsigned char* ptr, unsigned char c) { ::std::memset(ptr, c, SIZE); }
	};

	template <unsigned  SIZE>
	struct MemsetWorker<SIZE, true> {
		finline static void run(unsigned char* ptr, unsigned char c) {
			*((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
			MemsetWorker<SIZE - 4, true>::run(ptr + 4, c);
		}
	};

	template <>
	struct MemsetWorker<0, true> {
		finline static void run(unsigned char*, unsigned char) {}
	};

	template <>
	struct MemsetWorker<1, true> {
		finline static void run(unsigned char* ptr, unsigned char c) { *ptr = c; }
	};

	template <>
	struct MemsetWorker<2, true> {
		finline static void run(unsigned char* ptr, unsigned char c) { *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c; }
	};

	template <>
	struct MemsetWorker<3, true> {
		finline static void run(unsigned char* ptr, unsigned char c) {
			MemsetWorker<2, true>::run(ptr, c);
			MemsetWorker<1, true>::run(ptr + 2, c);
		}
	};

	template <unsigned SIZE>
	finline void memset(void* ptr, unsigned char c) {
		MemsetWorker<SIZE, SIZE <= 32>::run((unsigned char*)ptr, c);
	}


	//////////////////////////////////////////////////////////////////////////////
	// memset with fill value (using meta-programming)
	//////////////////////////////////////////////////////////////////////////////

	template <unsigned SIZE, bool direct, unsigned char c>
	struct MemsetConstValueWorker {
		finline static void run(unsigned char* ptr) { ::std::memset(ptr, c, SIZE); }
	};

	template <unsigned  SIZE, unsigned char c>
	struct MemsetConstValueWorker<SIZE, true, c> {
		finline static void run(unsigned char* ptr) {
			*((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
			MemsetConstValueWorker<SIZE - 4, true, c>::run(ptr + 4);
		}
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<0, true, c> {
		finline static void run(unsigned char*) {}
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<1, true, c> {
		finline static void run(unsigned char* ptr) { *ptr = c; }
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<2, true, c> {
		finline static void run(unsigned char* ptr) { *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c; }
	};

	template <unsigned char c>
	struct MemsetConstValueWorker<3, true, c> {
		finline static void run(unsigned char* ptr) {
			MemsetConstValueWorker<2, true, c>::run(ptr);
			MemsetConstValueWorker<1, true, c>::run(ptr + 2);
		}
	};

	template <unsigned SIZE, unsigned char c>
	finline void memset(void* ptr) {
		MemsetConstValueWorker<SIZE, SIZE <= 32, c>::run((unsigned char*)ptr);
	}

}

#endif
