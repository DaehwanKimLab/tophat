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

#ifndef SEQAN_HEADER_BASIC_HOLDER_DYNAMIC_H
#define SEQAN_HEADER_BASIC_HOLDER_DYNAMIC_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

class VoidHolder
{
public:
	char * _data;
	size_t _length;
	bool _empty;

	VoidHolder()
		: _length(0)
		, _empty(true)
	{ }
	VoidHolder(VoidHolder const & other)
		: _length(other._length)
		, _empty(other._empty)
	{ 
		if (_length)
		{
			allocate(*this, _data, _length);
			arrayCopyForward(other._data, other._data + _length, _data);
		}
	}
	template <typename TValue>
	VoidHolder(TValue const & val)
		: _length(0)
		, _empty(true)
	{
		assign(*this, val);
	}

	~VoidHolder()
	{
		if (_length)
		{
			deallocate(*this, _data, _length);
			_length = 0;
		}
	}

	VoidHolder const & 
	operator = (VoidHolder const & other)
	{
		_empty = other._empty;
		if (!_empty)
		{
			if (_length < other._length)
			{
				deallocate(*this, _data, _length);
				_length = other._length;
				allocate(*this, _data, _length);
			}
			arrayCopyForward(other._data, other._data + other._length, _data);
		}
		return *this;
	}
	template <typename TValue>
	VoidHolder const & 
	operator = (TValue const & val)
	{
		assign(*this, val);
		return *this;
	}

};

//////////////////////////////////////////////////////////////////////////////

inline bool
empty(VoidHolder & me)
{
	return me._empty;
}
inline bool
empty(VoidHolder const & me)
{
	return me._empty;
}

//////////////////////////////////////////////////////////////////////////////

inline void
clear(VoidHolder & me)
{
	me._empty = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize>
inline void
resize(VoidHolder & me, 
	   TSize length)
{
	if (me._length < length)
	{
		deallocate(me, me._data, me._length);
		me._length = length;
		allocate(me, me._data, me._length);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline TValue &
create(VoidHolder & me)
{
	resize(me, sizeof(TValue));
	valueConstruct(reinterpret_cast<TValue *>(me._data));
	me._empty = false;
	return * reinterpret_cast<TValue *>(me._data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
destroy(VoidHolder & me)
{
	if (empty(me)) return;
	me._empty = true;

	SEQAN_ASSERT(me._length >= sizeof(TValue))
	valueDestruct(reinterpret_cast<TValue *>(me._data));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline TValue &
value(VoidHolder & me)
{
	if (me._empty)
	{
		return create<TValue>(me);
	}
	return * reinterpret_cast<TValue *>(me._data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline TValue &
getValue(VoidHolder & me)
{
	return * reinterpret_cast<TValue *>(me._data);
}
template <typename TValue>
inline TValue const &
getValue(VoidHolder const & me)
{
	return * reinterpret_cast<TValue const *>(me._data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
inline void
assign(VoidHolder & me, 
	   TValue const & val)
{
	assign(value<TValue>(me), val);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


