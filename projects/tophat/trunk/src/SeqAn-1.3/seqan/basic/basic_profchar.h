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

#ifndef SEQAN_HEADER_SEQAN_BASIC_PROFCHAR_H
#define SEQAN_HEADER_SEQAN_BASIC_PROFCHAR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Profile alphabet character
//////////////////////////////////////////////////////////////////////////////


template<typename TValue, typename TCount = unsigned int, typename TSpec = Default>
class ProfileType;

template<typename TValue, typename TCount, typename TSpec>
struct ValueSize<ProfileType<TValue, TCount, TSpec> >
{
	enum { VALUE = ValueSize<TValue>::VALUE + 1};
};

template<typename TValue, typename TCount, typename TSpec>
class ProfileType {
	public:
		typedef typename Size<ProfileType>::Type TSize;

		TCount count[ValueSize<ProfileType>::VALUE];

		ProfileType() {
			memset<ValueSize<ProfileType>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
		}

		~ProfileType() {}

		ProfileType(ProfileType const& other_data) {
			for(TSize i = 0; i<ValueSize<ProfileType>::VALUE; ++i) 
				count[i] = other_data.count[i];
		}


		template <typename TOther> 
		ProfileType(TOther const& other_data) {
			memset<ValueSize<ProfileType>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
			count[ordValue(TValue(other_data))] = 1;
		}

		ProfileType const& 
		operator = (ProfileType const& other_data) 
		{
			if (this == &other_data) return *this;
			for(TSize i = 0; i<ValueSize<ProfileType>::VALUE; ++i) 
				count[i] = other_data.count[i];
			return *this;
		}

		template <typename TOther>
		ProfileType const& 
		operator = (TOther const& other_data) {
			memset<ValueSize<ProfileType>::VALUE * sizeof(TCount), (unsigned char) 0>(count);
			count[ordValue(TValue(other_data))] = 1;
			return *this;
		}

		operator char() { 
			typedef typename Size<ProfileType>::Type TSize;
			TSize maxIndex = 0;
			TCount maxCount = count[0];
			for(TSize i = 1; i<ValueSize<ProfileType>::VALUE; ++i) {
				if (count[i] > maxCount) {
					maxIndex = i;
					maxCount = count[i];
				}
			}
			return (maxIndex == ValueSize<ProfileType>::VALUE - 1) ? gapValue<char>() : (char) TValue(maxIndex);
		}

		bool
		operator==(ProfileType const& other_data) const {
			for(TSize i = 0; i<ValueSize<ProfileType>::VALUE; ++i) {
				if (count[i] != other_data.count[i]) return false;
			}
			return true;
		}

		bool 
		operator!=(ProfileType const& other_data) const {
			for(TSize i = 0; i<ValueSize<ProfileType>::VALUE; ++i) {
				if (count[i] != other_data.count[i]) return true;
			}
			return false;
		}
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct SourceValue;

template<typename TValue, typename TCount, typename TSpec>
struct SourceValue<ProfileType<TValue, TCount, TSpec> >
{
	typedef TValue Type;
};

template<typename TValue, typename TCount, typename TSpec>
struct SourceValue<ProfileType<TValue, TCount, TSpec> const>
{
	typedef TValue const Type;
};

////////////////////////////////////////////////////////////////////////////////

template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline bool
empty(ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	typedef typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type TSize;
	// Check if there are only gaps
	for(TSize i = 0; i<ValueSize<TSourceValue>::VALUE; ++i) {
		if (source.count[i]) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type
_getMaxIndex(ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	typedef ProfileType<TSourceValue, TSourceCount, TSourceSpec> TProfileType;
	typedef typename Size<TProfileType>::Type TSize;
	TSize maxIndex = 0;
	TSourceCount maxCount = source.count[0];
	for(TSize i = 1; i<ValueSize<TProfileType>::VALUE; ++i) {
		if (source.count[i] > maxCount) {
			maxIndex = i;
			maxCount = source.count[i];
		}
	}
	return maxIndex;
}


////////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	target.value = _getMaxIndex(source);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Convert<TTarget, ProfileType<TSourceValue, TSourceCount, TSourceSpec> >::Type
convertImpl(Convert<TTarget, T> const,
			ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	return (_getMaxIndex(source) == ValueSize<TSourceValue>::VALUE) ? convertImpl(Convert<TTarget, T>(), '-') : convertImpl(Convert<TTarget, T>(), TSourceValue(_getMaxIndex(source)));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStream, typename TValue, typename TCount, typename TSpec>
inline TStream& 
operator<<(TStream& os, ProfileType<TValue, TCount, TSpec> const& rhs) {
	typedef ProfileType<TValue, TCount, TSpec> TProfileType;
	typedef typename Size<TProfileType>::Type TSize;
	for(TSize i = 0; i<ValueSize<TProfileType>::VALUE; ++i) {
		os << i << ':' << rhs.count[i] << ' ' << ';';
	}
	return os;
}

}

#endif

