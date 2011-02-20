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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_SSE2_H
#define SEQAN_HEADER_BASIC_SSE2_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifdef SEQAN_USE_SSE2_WORDS
#ifdef __SSE2__

#include <emmintrin.h>

namespace SEQAN_NAMESPACE_MAIN
{
	
//	#define __int128 SSE2_int128
	
	#ifndef SEQAN_SSE2_INT128
	#define SEQAN_SSE2_INT128
	#endif
		
	// ATTENTION:
	// The SSE2_int128 struct must be 16-byte aligned. Some allocators
	// don't ensure the correct alignment, e.g. the STL allocators.
	// Either avoid STL classes holding SSE2_int128 structs (maps in find_pex.h)
	// or avoid the SEQAN_USE_SSE2_WORDS define above.

	// may become obsolete when 128-bit integers will be introduced
	struct SSE2_int128
	{
	public:
		union {
			__m128i	 v;
			__uint64 v64[2];
			unsigned v32[4];
		}						data;

//____________________________________________________________________________

	public:
		SSE2_int128();
		SSE2_int128(SSE2_int128 const &);
		SSE2_int128(__m128i const &);
		SSE2_int128(__int64, __int64);
		SSE2_int128(int, int, int, int);
		SSE2_int128(short, short, short, short, short, short, short, short);

		template <typename TValue>
		SSE2_int128(TValue const &);
		
		//____________________________________________________________________________
		
		template <typename TValue>
		SSE2_int128 & operator = (TValue const &);
		
		//____________________________________________________________________________

		operator __m128i () const;
		operator __int64 () const;
/*		operator bool () const;
		operator __uint64 ();
		operator int ();
		operator unsigned int ();
		operator short ();
		operator unsigned short ();
		operator char ();
		operator signed char ();
		operator unsigned char ();
*/	};
	
	template <> struct _IsSimple<__m128i> { typedef True Type; };
	template <> struct _IsSimple<SSE2_int128> { typedef True Type; };

//____________________________________________________________________________
// clear

inline void
clear(SSE2_int128 &me)
{
	me.data.v = _mm_setzero_si128();
}

//____________________________________________________________________________
// assign

inline void
assign(SSE2_int128 &me, SSE2_int128 const &other)
{
	me.data = other.data;
}

// 1x 128bit
inline void
assign(SSE2_int128 &me, __m128i const &other)
{
	me.data.v = other;
}
inline SSE2_int128::operator __m128i () const
{
	return data.v;
}
	
// 64bit => 128bit
inline void
assign(SSE2_int128 &me, __int64 other)
{
	me.data.v = _mm_set_epi32(0, 0, other >> 32, other);
}
inline void
assign(SSE2_int128 &me, __uint64 other)
{
	me.data.v = _mm_set_epi32(0, 0, other >> 32, other);
}
// 64bit <= 128bit
inline SSE2_int128::operator __int64 () const
{
	return data.v64[0];
}
/*inline SSE2_int128::operator __uint64 ()
{
	return data.v64[0];
}
*/
// 32bit => 128bit
inline void
assign(SSE2_int128 &me, int other)
{
	me.data.v = _mm_set_epi32(0, 0, 0, other);
}
inline void
assign(SSE2_int128 &me, unsigned int other)
{
	me.data.v = _mm_set_epi32(0, 0, 0, other);
}
// 32bit <= 128bit
/*inline SSE2_int128::operator int ()
{
	return data.v32[0];
}
inline SSE2_int128::operator unsigned int ()
{
	return data.v32[0];
}
*/
// 16bit => 128bit
inline void
assign(SSE2_int128 &me, short other)
{
	me.data.v = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(SSE2_int128 &me, unsigned short other)
{
	me.data.v = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, other);
}
// 16bit <= 128bit
/*inline SSE2_int128::operator short ()
{
	return _mm_extract_epi16(data.v, 0);
}
inline SSE2_int128::operator unsigned short ()
{
	return _mm_extract_epi16(data.v, 0);
}
*/
// 8bit => 128bit
inline void
assign(SSE2_int128 &me, char other)
{
	me.data.v = _mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(SSE2_int128 &me, signed char other)
{
	me.data.v = _mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(SSE2_int128 &me, unsigned char other)
{
	me.data.v = _mm_set_epi8(
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, other);
}
// 8bit <= 128bit
/*inline SSE2_int128::operator char ()
{
	return _mm_extract_epi16(data.v, 0);
}
inline SSE2_int128::operator signed char ()
{
	return _mm_extract_epi16(data.v, 0);
}
inline SSE2_int128::operator unsigned char ()
{
	return _mm_extract_epi16(data.v, 0);
}
	
inline SSE2_int128::operator bool () const
{
	return (__int64)(SSE2_int128)_mm_or_si128(data.v, _mm_unpackhi_epi64(data.v, data.v));
}
*/
//____________________________________________________________________________
// constructors

inline SSE2_int128::SSE2_int128()
{
	clear(*this);
}
	
inline SSE2_int128::SSE2_int128(SSE2_int128 const &other)
{
	assign(*this, other);
}

inline SSE2_int128::SSE2_int128(__m128i const &other)
{
	assign(*this, other);
}

// 2x 64bit
inline SSE2_int128::SSE2_int128(__int64 q1, __int64 q0)
{
	data.v = _mm_set_epi32(q1 >> 32, q1, q0 >> 32, q0);
}

// 4x 32bit
inline SSE2_int128::SSE2_int128(int q3, int q2, int q1, int q0)
{
	data.v = _mm_set_epi32(q3, q2, q1, q0);
}

// 8x 16bit
inline SSE2_int128::SSE2_int128(
				   short q7, short q6, short q5, short q4,
				   short q3, short q2, short q1, short q0)
{
	data.v = _mm_set_epi16(q7, q6, q5, q4, q3, q2, q1, q0);
}
	
template <typename TValue>
inline SSE2_int128::SSE2_int128(TValue const &other)
{
	assign(*this, other);
}

//____________________________________________________________________________
// operator =

template <typename TValue>
inline SSE2_int128 &
SSE2_int128::operator = (TValue const &other)
{
	assign(*this, other);
	return *this;
}
	
//____________________________________________________________________________
// logical operators

inline SSE2_int128
operator & (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return _mm_and_si128((__m128i)a, (__m128i)b);
}
inline SSE2_int128
operator &= (SSE2_int128 &a, SSE2_int128 const &b)
{
	a.data.v = _mm_and_si128((__m128i)a, (__m128i)b);
	return a;
}

inline SSE2_int128
operator | (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return _mm_or_si128((__m128i)a, (__m128i)b);
}
inline SSE2_int128
operator |= (SSE2_int128 &a, SSE2_int128 const &b)
{
	a.data.v = _mm_or_si128((__m128i)a, (__m128i)b);
	return a;
}

inline SSE2_int128
operator ^ (SSE2_int128 const &a, SSE2_int128 const &b)
{
	return _mm_xor_si128((__m128i)a, (__m128i)b);
}
inline SSE2_int128
operator ^= (SSE2_int128 &a, SSE2_int128 const &b)
{
	a.data.v = _mm_xor_si128((__m128i)a, (__m128i)b);
	return a;
}

inline SSE2_int128
operator ~ (SSE2_int128 const &a)
{
	__m128i _a = (__m128i)a;
	return _mm_xor_si128(_a, _mm_cmpeq_epi32(_a, _a));
}

//____________________________________________________________________________
// shift operators

inline SSE2_int128
operator << (SSE2_int128 const &a, int n)
{
	return _mm_or_si128(
		// n <= 64
		_mm_or_si128(
			_mm_sll_epi64((__m128i)a, (SSE2_int128)n),
			_mm_srl_epi64(
				_mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
				(SSE2_int128)(64-n)
			)
		),
		// n >= 64
		_mm_sll_epi64(
			_mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
			(SSE2_int128)(n-64)
		)
	);
}
inline SSE2_int128
operator <<= (SSE2_int128 &a, int n)
{
	a.data.v = _mm_or_si128(
		// n <= 64
		_mm_or_si128(
			_mm_sll_epi64((__m128i)a, (SSE2_int128)n),
			_mm_srl_epi64(
				_mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
				(SSE2_int128)(64-n)
			)
		),
		// n >= 64
		_mm_sll_epi64(
			_mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
			(SSE2_int128)(n-64)
		)
	);
	return a;
}
inline SSE2_int128
operator << (SSE2_int128 const &a, unsigned int n)
{
	return a << (int)n;
}
inline SSE2_int128
operator <<= (SSE2_int128 &a, unsigned int n)
{
	return a <<= (int)n;
}

inline SSE2_int128
operator >> (SSE2_int128 const &a, int n)
{
	return _mm_or_si128(
		// n <= 64
		_mm_or_si128(
			_mm_srl_epi64((__m128i)a, (SSE2_int128)n),
			_mm_sll_epi64(
				_mm_unpackhi_epi64(a.data.v, _mm_setzero_si128()),
				(SSE2_int128)(64-n)
			)
		),
		// n >= 64
		_mm_srl_epi64(
			_mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
			(SSE2_int128)(n-64)
		)
	);
}
inline SSE2_int128
operator >>= (SSE2_int128 &a, int n)
{
	a.data.v = _mm_or_si128(
		// n <= 64
		_mm_or_si128(
			_mm_srl_epi64((__m128i)a, (SSE2_int128)n),
			_mm_sll_epi64(
				_mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
				(SSE2_int128)(64-n)
			)
		),
		// n >= 64
		_mm_srl_epi64(
			_mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
			(SSE2_int128)(n-64)
		)
	);
	return a;
}
inline SSE2_int128
operator >> (SSE2_int128 const &a, unsigned int n)
{
	return a >> (int)n;
}
inline SSE2_int128
operator >>= (SSE2_int128 &a, unsigned int n)
{
	return a >>= (int)n;
}

//____________________________________________________________________________
// artihmetic operators

//template <typename T>
inline SSE2_int128
operator + (SSE2_int128 const &a, SSE2_int128 const &b)
{
	static const __m128i carry = SSE2_int128(0, 1, 0, 0);

	union {
		__uint64 v64[2];
		__m128i v;
	} _sum;
	
	_sum.v = _mm_add_epi64((__m128i)a, (__m128i)b);
	if (_sum.v64[0] >= a.data.v64[0])
		return _sum.v;
	else
		return _mm_add_epi64(_sum.v, carry);
}

//template <typename T>
inline SSE2_int128
operator - (SSE2_int128 const &a, SSE2_int128 const &b)
{
	static const __m128i carry = SSE2_int128(0, 1, 0, 0);

	union {
		__uint64 v64[2];
		__m128i v;
	} _diff;
	
	_diff.v = _mm_sub_epi64((__m128i)a, (__m128i)b);
	if (_diff.v64[0] <= a.data.v64[0])
		return _diff.v;
	else
		return _mm_sub_epi64(_diff.v, carry);
}

// TODO: operator *, /

//____________________________________________________________________________
// compares

inline bool
operator == (SSE2_int128 const &a, SSE2_int128 const &b)
{
	__m128i _e = _mm_cmpeq_epi32((__m128i)a, (__m128i)b);
	return ((__int64)(SSE2_int128)_mm_and_si128(_e, _mm_unpackhi_epi64(_e, _e))) == ~(__int64)0;
}

inline bool
operator != (SSE2_int128 const &a, SSE2_int128 const &b)
{
	__m128i _e = _mm_cmpeq_epi32((__m128i)a, (__m128i)b);
	return ((__int64)(SSE2_int128)_mm_and_si128(_e, _mm_unpackhi_epi64(_e, _e))) != ~(__int64)0;
}
	
// TODO: operator <, <=, >, >=

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifdef __SSE2__
#endif //#ifdef SEQAN_USE_SSE2_WORDS

#endif //#ifndef SEQAN_HEADER_...
