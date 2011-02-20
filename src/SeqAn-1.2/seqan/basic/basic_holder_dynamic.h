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
  $Id: $
 ==========================================================================*/

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


