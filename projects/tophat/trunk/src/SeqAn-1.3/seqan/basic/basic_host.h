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

#ifndef SEQAN_HEADER_BASIC_HOST_H
#define SEQAN_HEADER_BASIC_HOST_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Host Functions
//////////////////////////////////////////////////////////////////////////////
//these functions assume that the hosted object exports a function "_dataHost"
//that returns a reference to a holder type of Host<T>::Type & 

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline bool
emptyHost(T const & me)
{
SEQAN_CHECKPOINT
	return empty(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline bool
dependentHost(T const & me)
{
SEQAN_CHECKPOINT
	return dependent(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
clearHost(T & me)
{
SEQAN_CHECKPOINT
	clear(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
createHost(T & me)
{
SEQAN_CHECKPOINT
	create(_dataHost(me));
}

//____________________________________________________________________________

template <typename T, typename THost>
inline void
createHost(T & me,
		   THost & host_)
{
SEQAN_CHECKPOINT
	create(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
createHost(T & me,
		   THost const & host_)
{
SEQAN_CHECKPOINT
	create(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
setHost(T & me,
		THost & host_)
{
SEQAN_CHECKPOINT
	setValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
setHost(T & me,
		THost const & host_)
{
SEQAN_CHECKPOINT
	setValue(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename T>
inline typename Host<T>::Type &
host(T & me)
{
SEQAN_CHECKPOINT
	return value(_dataHost(me));
}
template <typename T>
inline typename Host<T const>::Type &
host(T const & me)
{
SEQAN_CHECKPOINT
	return value(_dataHost(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
assignHost(T & me,
		   THost & host_)
{
SEQAN_CHECKPOINT
	assignValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
assignHost(T & me,
		   THost const & host_)
{
SEQAN_CHECKPOINT
	assignValue(_dataHost(me), host_);
}
//////////////////////////////////////////////////////////////////////////////

template <typename T, typename THost>
inline void
moveHost(T & me,
		 THost & host_)
{
SEQAN_CHECKPOINT
	moveValue(_dataHost(me), host_);
}
template <typename T, typename THost>
inline void
moveHost(T & me,
		 THost const & host_)
{
SEQAN_CHECKPOINT
	moveValue(_dataHost(me), host_);
}

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


