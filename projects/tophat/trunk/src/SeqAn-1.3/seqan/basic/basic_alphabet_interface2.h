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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_BASIC_ALPHABET_INTERFACE2_H
#define SEQAN_HEADER_BASIC_ALPHABET_INTERFACE2_H

#include <new>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// gapValue, gapValueImpl
//////////////////////////////////////////////////////////////////////////////
/**
.Function.gapValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.gapValue@.
..signature:gapValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A gap character.
..remarks.text:This function implements @Function.getValue@. 
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
..include:seqan/basic.h
*/

template <typename T>
inline T const &
gapValueImpl(T *)
{
SEQAN_CHECKPOINT
	static T const _gap = T();
	return _gap;
}

template <typename T>
inline T
unknownValueImpl(T *)
{
SEQAN_CHECKPOINT
	return 'N';
}

/**
.Function.gapValue:
..cat:Alphabets
..cat:Alignments
..summary:Returns reference to a value that is used as gap character.
..signature:gapValue<TValue>()
..param.TValue:Value type.
..returns:A gap character.
..remarks.text:The function is implemented in @Function.gapValueImpl@. 
Do not specialize $gapValue$, specialize @Function.gapValueImpl@ instead!
..see:Function.gapValueImpl
..include:seqan/basic.h
*/

/*
template <typename T>
inline T const &
gapValue()
{
SEQAN_CHECKPOINT
	static T * _tag = 0;
	return gapValueImpl(_tag);
}
*/
template <typename T>
inline T
gapValue()
{
SEQAN_CHECKPOINT
	static T * _tag = 0;
	return gapValueImpl(_tag);
}

template <typename T>
inline T
unknownValue()
{
SEQAN_CHECKPOINT
	static T * _tag = 0;
	return unknownValueImpl(_tag);
}


//////////////////////////////////////////////////////////////////////////////
// maxValue, supremumValueImpl
//////////////////////////////////////////////////////////////////////////////

/**
.Function.supremumValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.maxValue@.
..signature:supremumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$.
..remarks.text:This function implements @Function.maxValue@. 
It is recommended to use @Function.maxValue@ rather than $supremumValueImpl$.
..status:deprecated, will be removed in favour of @Metafunction.MaxValue@
..include:seqan/basic.h
*/

/*
template <typename T>
inline T const &
supremumValueImpl(T *)
{
	static T const _value = -1;
	return _value;
}
*/

/**
.Function.maxValue:
..cat:Alphabets
..summary:Supremum for a given type.
..signature:maxValue<T>()
..param.T:An ordered type.
..returns:A value $inf$ that holds: $inf >= i$ for all values $i$ of type $T$.
..remarks.text:The function is implemented in @Function.supremumValueImpl@. 
Do not specialize $maxValue$, specialize @Function.supremumValueImpl@ instead!
..see:Function.supremumValueImpl
..status:deprecated, will be removed in favour of @Metafunction.MaxValue@
..include:seqan/basic.h
*/

template <typename T>
inline T const &
maxValue()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

template <typename T>
inline T const &
maxValue(T)
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////
// minValue, infimumValueImpl
//////////////////////////////////////////////////////////////////////////////

/**
.Function.infimumValueImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.minValue@.
..signature:infimumValueImpl(value_pointer_tag)
..param.value_pointer_tag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$.
..remarks.text:This function implements @Function.minValue@. 
It is recommended to use @Function.minValue@ rather than $infimumValueImpl$.
..status:deprecated, will be removed in favour of @Metafunction.MinValue@
..include:seqan/basic.h
*/

/*
template <typename T>
inline T const &
infimumValueImpl(T *)
{
	static T const _value = -1;
	return _value;
}
*/

/**
.Function.minValue:
..cat:Alphabets
..summary:Infimum for a given type.
..signature:minValue<T>()
..param.T:An ordered type.
..returns:A value $inf$ that holds: $inf <= i$ for all values $i$ of type $T$.
..remarks.text:The function is implemented in @Function.infimumValueImpl@. 
Do not specialize $minValue$, specialize @Function.infimumValueImpl@ instead!
..see:Function.infimumValueImpl
..see:Function.maxValue
..status:deprecated, will be removed in favour of @Metafunction.MinValue@
..include:seqan/basic.h
*/

template <typename T>
inline T const &
minValue()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return infimumValueImpl(_tag);
}

template <typename T>
inline T const &
minValue(T)
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return infimumValueImpl(_tag);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
