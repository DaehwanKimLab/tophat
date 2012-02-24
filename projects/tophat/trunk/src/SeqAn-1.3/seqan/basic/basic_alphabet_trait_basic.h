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

#ifndef SEQAN_HEADER_BASIC_ALPHABET_TRAIT_BASIC_H
#define SEQAN_HEADER_BASIC_ALPHABET_TRAIT_BASIC_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//arrayConstruct
//////////////////////////////////////////////////////////////////////////////

template<typename TIterator>
inline void 
_arrayConstructPointer(TIterator, 
						TIterator,
						True)
{
SEQAN_CHECKPOINT
	//nothing to do
}
template<typename TIterator>
inline void 
_arrayConstructPointer(TIterator begin_, 
						TIterator end_,
						False)
{
SEQAN_CHECKPOINT
	_arrayConstructDefault(begin_, end_);
}
template<typename TValue>
inline void 
arrayConstruct(TValue * begin_, 
			   TValue * end_)
{
SEQAN_CHECKPOINT
	_arrayConstructPointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

//____________________________________________________________________________

template<typename TIterator, typename TParam>
inline void 
_arrayConstructPointer(TIterator begin_, 
						TIterator end_, 
						TParam const & param_,
						True)
{
SEQAN_CHECKPOINT
	arrayFill(begin_, end_, param_);
}
template<typename TIterator, typename TParam>
inline void 
_arrayConstructPointer(TIterator begin_, 
						TIterator end_, 
						TParam const & param_,
						False)
{
SEQAN_CHECKPOINT
	_arrayConstructDefault(begin_, end_, param_);
}
template<typename TValue, typename TParam>
inline void 
arrayConstruct(TValue * begin_, 
			   TValue * end_, 
			   TParam const & param_)
{
SEQAN_CHECKPOINT
	_arrayConstructPointer(begin_, end_, param_, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructCopy
//////////////////////////////////////////////////////////////////////////////

template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
							TValueSource * source_end, 
							TValueTarget * target_begin,
							True)
{
SEQAN_CHECKPOINT
	arrayCopyForward(source_begin, source_end, target_begin);
}
template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
							TValueSource * source_end, 
							TValueTarget const* target_begin,
							True)
{
SEQAN_CHECKPOINT
	arrayCopyForward(source_begin, source_end, const_cast<TValueTarget *>(target_begin));
}

template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
							TValueSource * source_end, 
							TValueTarget * target_begin,
							False)
{
SEQAN_CHECKPOINT
	_arrayConstructCopyDefault(source_begin, source_end, target_begin);
}
template<typename TValueSource, typename TValueTarget>
inline void 
arrayConstructCopy(TValueSource * source_begin, 
				   TValueSource * source_end, 
				   TValueTarget * target_begin)
{
SEQAN_CHECKPOINT
	_arrayConstructCopyPointer(source_begin, source_end, target_begin, typename IsSimple<TValueTarget>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayConstructMove
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayConstructMovePointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							True)
{
SEQAN_CHECKPOINT
	arrayMoveForward(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
_arrayConstructMovePointer(TValue * source_begin, 
							TValue * source_end, 
							TValue * target_begin,
							False)
{
SEQAN_CHECKPOINT
	_arrayConstructMoveDefault(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayConstructMove(TValue * source_begin, 
				   TValue * source_end, 
				   TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayConstructMovePointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayDestruct
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayDestructPointer(TValue * /*begin_*/, 
					   TValue * /*end_*/,
					   True)
{
SEQAN_CHECKPOINT
	//do nothing
}
template<typename TValue>
inline void 
_arrayDestructPointer(TValue * begin_, 
					   TValue * end_,
					   False)
{
SEQAN_CHECKPOINT
	_arrayDestructDefault(begin_, end_);
}
template<typename TValue>
inline void 
arrayDestruct(TValue * begin_, 
			  TValue * end_)
{
SEQAN_CHECKPOINT
	_arrayDestructPointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayFill
//////////////////////////////////////////////////////////////////////////////

//no specializiation for pointer to simple

//////////////////////////////////////////////////////////////////////////////
//arrayCopyForward
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayCopyForwardPointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template<typename TValue>
inline void 
_arrayCopyForwardPointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  False)
{
SEQAN_CHECKPOINT
	_arrayCopyForwardDefault(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayCopyForward(TValue * source_begin, 
				 TValue * source_end, 
				 TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayCopyBackward
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayCopyBackwardPointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void 
_arrayCopyBackwardPointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   False)
{
SEQAN_CHECKPOINT
	_arrayCopyBackwardDefault(source_begin, source_end, target_begin); 
}
template<typename TValue>
inline void 
arrayCopyBackward(TValue * source_begin, 
				  TValue * source_end, 
				  TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayCopyBackwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveForward
//////////////////////////////////////////////////////////////////////////////

template<typename TValue>
inline void 
_arrayMoveForwardPointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template<typename TValue>
inline void 
_arrayMoveForwardPointer(TValue * source_begin, 
						  TValue * source_end, 
						  TValue * target_begin,
						  False)
{
SEQAN_CHECKPOINT
	_arrayMoveForwardDefault(source_begin, source_end, target_begin);
}
template<typename TValue>
inline void 
arrayMoveForward(TValue * source_begin, 
				 TValue * source_end, 
				 TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayMoveBackward
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayMoveBackwardPointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   True)
{
SEQAN_CHECKPOINT
	::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void 
_arrayMoveBackwardPointer(TValue * source_begin, 
						   TValue * source_end, 
						   TValue * target_begin,
						   False)
{
SEQAN_CHECKPOINT
	_arrayMoveBackwardDefault(source_begin, source_end, target_begin); 
}
template<typename TValue>
inline void 
arrayMoveBackward(TValue * source_begin, 
				  TValue * source_end, 
				  TValue * target_begin)
{
SEQAN_CHECKPOINT
	_arrayMoveBackwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

//////////////////////////////////////////////////////////////////////////////
//arrayClearSpace
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void 
_arrayClearSpacePointer(TValue * array_begin, 
						size_t array_length, 
						size_t keep_from, 
						size_t move_to,
						True)
{
	if (keep_from == move_to) return;
SEQAN_CHECKPOINT
	arrayMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
}
template <typename TValue>
inline void 
_arrayClearSpacePointer(TValue * array_begin, 
						size_t array_length, 
						size_t keep_from, 
						size_t move_to,
						False)
{
	_arrayClearSpaceDefault(array_begin, array_length, keep_from, move_to);
}
template <typename TValue>
void arrayClearSpace(TValue * array_begin, 
					 size_t array_length, 
					 size_t keep_from, 
					 size_t move_to)
{
	_arrayClearSpacePointer(array_begin, array_length, keep_from, move_to, typename IsSimple<TValue>::Type() );
}



//////////////////////////////////////////////////////////////////////////////
// IsSimple specializations
//////////////////////////////////////////////////////////////////////////////

// standard types
template <> struct IsSimple_< bool > { typedef True Type; };
template <> struct IsSimple_< char > { typedef True Type; };

template <> struct IsSimple_< unsigned char > { typedef True Type; };
template <> struct IsSimple_< unsigned short > { typedef True Type; };
template <> struct IsSimple_< unsigned int > { typedef True Type; };
template <> struct IsSimple_< unsigned long > { typedef True Type; };

template <> struct IsSimple_< signed char > { typedef True Type; };
template <> struct IsSimple_< signed short > { typedef True Type; };
template <> struct IsSimple_< signed int > { typedef True Type; };
template <> struct IsSimple_< signed long > { typedef True Type; };

template <> struct IsSimple_< float > { typedef True Type; };
template <> struct IsSimple_< double > { typedef True Type; };
template <> struct IsSimple_< long double > { typedef True Type; };

// user defined types (re-specializations are allowed here)
template <> struct IsSimple< wchar_t > { typedef True Type; };
template <> struct IsSimple< __int64 > { typedef True Type; };
template <> struct IsSimple< __uint64 > { typedef True Type; };

//////////////////////////////////////////////////////////////////////////////
// gapValue
//////////////////////////////////////////////////////////////////////////////

inline char const &
gapValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _gap = '-';
	return _gap;
}
inline char const &
gapValueImpl(char const *)
{
SEQAN_CHECKPOINT
	static char const _gap = '-';
	return _gap;
}

inline char const &
unknownValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _unknown = 'N';
	return _unknown;
}
inline char const &
unknownValueImpl(char const *)
{
SEQAN_CHECKPOINT
	static char const _unknown = 'N';
	return _unknown;
}

//////////////////////////////////////////////////////////////////////////////
// generic extreme values
//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T const &
supremumValueImpl(T *)
{
SEQAN_CHECKPOINT
	return MaxValue<T>::VALUE;
}
template <typename T>
inline T const &
infimumValueImpl(T *)
{
SEQAN_CHECKPOINT
	return MinValue<T>::VALUE;
}

//////////////////////////////////////////////////////////////////////////////
// bool 
//////////////////////////////////////////////////////////////////////////////

template <> struct BitsPerValue< bool > { enum { VALUE = 1 }; };

/*
//////////////////////////////////////////////////////////////////////////////
// char 
//////////////////////////////////////////////////////////////////////////////

inline char const &
supremumValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _value = (char) 127;
	return _value;
}
inline char const &
infimumValueImpl(char *)
{
SEQAN_CHECKPOINT
	static char const _value = (char) -128;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed char 
//////////////////////////////////////////////////////////////////////////////

inline signed char const &
supremumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
	static signed char const _value = 127;
	return _value;
}
inline signed char const &
infimumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
	static signed char const _value = -128;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned char 
//////////////////////////////////////////////////////////////////////////////

inline unsigned char const &
supremumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
	static unsigned char const _value = 255;
	return _value;
}
inline unsigned char const &
infimumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
	static unsigned char const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// wchar_t
//////////////////////////////////////////////////////////////////////////////

inline wchar_t const &
supremumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
	static wchar_t const _value = 1UL << (BitsPerValue<wchar_t>::VALUE) - 1;
	return _value;
}
inline wchar_t const &
infimumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
	static wchar_t const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed short 
//////////////////////////////////////////////////////////////////////////////

inline signed short const &
supremumValueImpl(signed short *)
{
SEQAN_CHECKPOINT
	static signed short const _value = (((1 << (BitsPerValue<signed short>::VALUE - 2)) - 1) << 1) + 1;
	return _value;
}
inline signed short const &
infimumValueImpl(signed short *dummy)
{
SEQAN_CHECKPOINT
	static signed short const _value = -supremumValueImpl(dummy) - 1;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned short 
//////////////////////////////////////////////////////////////////////////////

inline unsigned short const &
supremumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
	static unsigned short const _value = (((1 << (BitsPerValue<unsigned short>::VALUE - 1)) - 1) << 1) + 1;
	return _value;
}
inline unsigned short const &
infimumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
	static unsigned short const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed int 
//////////////////////////////////////////////////////////////////////////////

inline signed int const &
supremumValueImpl(signed int *)
{
SEQAN_CHECKPOINT
	static signed int const _value = (((1 << (BitsPerValue<signed int>::VALUE - 2)) - 1) << 1) + 1;
	return _value;
}
inline signed int const &
infimumValueImpl(signed int *dummy)
{
SEQAN_CHECKPOINT
	static signed int const _value = -supremumValueImpl(dummy) - 1;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned int 
//////////////////////////////////////////////////////////////////////////////

inline unsigned int const &
supremumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
	static unsigned int const _value = ~0ul;
	return _value;
}
inline unsigned int const &
infimumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
	static unsigned int const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed long
//////////////////////////////////////////////////////////////////////////////

inline signed long const &
supremumValueImpl(signed long *)
{
SEQAN_CHECKPOINT
	static signed long const _value = (((1 << (BitsPerValue<signed long>::VALUE - 2)) - 1) << 1) + 1;
	return _value;
}
inline signed long const &
infimumValueImpl(signed long *dummy)
{
SEQAN_CHECKPOINT
	static signed long const _value = -supremumValueImpl(dummy) - 1;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned long
//////////////////////////////////////////////////////////////////////////////

inline unsigned long const &
supremumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
	static unsigned long const _value = ~0ul;
	return _value;
}
inline unsigned long const &
infimumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
	static unsigned long const _value = 0;
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed 64bit int (cannot use long long <- no ISO C++)
//////////////////////////////////////////////////////////////////////////////

inline __int64 const &
supremumValueImpl(__int64 *)
{
SEQAN_CHECKPOINT
	static __int64 const _value = ((((__int64)1 << (BitsPerValue<__int64>::VALUE - 2)) - 1) << 1) + 1;
	return _value;
}
inline __int64 const &
infimumValueImpl(__int64 *dummy)
{
SEQAN_CHECKPOINT
	static __int64 const _value = -supremumValueImpl(dummy) - 1;
	return _value;
}
*/

//////////////////////////////////////////////////////////////////////////////
// float 
//////////////////////////////////////////////////////////////////////////////

inline float const &
supremumValueImpl(float *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static float const _value = ::std::numeric_limits<float>::infinity( );
#else
	static float const _value = 3.40282347e+38F;
#endif
	return _value;
}
inline float const &
infimumValueImpl(float *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static float const _value = -::std::numeric_limits<float>::infinity( );
#else
	static float const _value = -3.40282347e+38F;
#endif
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// double 
//////////////////////////////////////////////////////////////////////////////

inline double const &
supremumValueImpl(double *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static double const _value = ::std::numeric_limits<double>::infinity( );
#else
	static double const _value = 1.7976931348623157e+308;
#endif
	return _value;
}
inline double const &
infimumValueImpl(double *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static double const _value = -::std::numeric_limits<double>::infinity( );
#else
	static double const _value = -1.7976931348623157e+308;
#endif
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
// long double 
//////////////////////////////////////////////////////////////////////////////

inline long double const &
supremumValueImpl(long double *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static long double const _value = ::std::numeric_limits<long double>::infinity( );
#else
	static long double const _value = 1.7976931348623157e+308;
#endif
	return _value;
}
inline long double const &
infimumValueImpl(long double *)
{
SEQAN_CHECKPOINT
#ifdef PLATFORM_WINDOWS
	static long double const _value = -::std::numeric_limits<long double>::infinity( );
#else
	static long double const _value = -1.7976931348623157e+308;
#endif
	return _value;
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
