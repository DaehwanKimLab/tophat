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
  $Id: align_trace.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_TRACE_H
#define SEQAN_HEADER_ALIGN_TRACE_H

namespace SEQAN_NAMESPACE_MAIN
{

//??? DIESER CODE IST NOCH UNGETESTET

//////////////////////////////////////////////////////////////////////////////
// AlignTrace

template <typename TSize = size_t, typename TSpec = void>
struct AlignTrace
{
	String<TSize> data_lengths;
	String<TSize> data_factors;
	String<unsigned int> data_dat;
	unsigned int data_bits_per_entry;
	unsigned int data_entries_per_word; //entries per word, or words per entry
	bool data_has_few_bits_per_entry;
	unsigned char data_max_sub_position;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
struct Size< AlignTrace<TSize, TSpec> >
{
	typedef TSize Type;
};

//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline unsigned int
dimension(AlignTrace<TSize, TSpec> & me)
{
	return length(me.data_lengths);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline unsigned int
setDimension(AlignTrace<TSize, TSpec> & me,
			 unsigned int _dim)
{
	return resize(me.data_lengths, _dim);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline TSize
length(AlignTrace<TSize, TSpec> const & me,
	   TSize _dim)
{
	return me.data_lengths[_dim];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline void
setLength(AlignTrace<TSize, TSpec> & me,
		  unsigned int _dim,
		  TSize _length)
{
	me.data_lengths[_dim] = _length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline unsigned int
_getBitsPerEntry(AlignTrace<TSize, TSpec> & me)
{
	return length(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec>
inline void
_initMembers(AlignTrace<TSize, TSpec> & me)
{
	SEQAN_ASSERT(dimension(me) >= 1);

	unsigned int bits_per_word = BitsPerValue<unsigned int>::VALUE;

	me.data_bits_per_entry = _computeBitsPerEntry(me);
	me.data_has_few_bits_per_entry = (me.data_bits_per_entry <= bits_per_word);

	unsigned int i;
	TSize _size;
	if (me.data_has_few_bits_per_entry)
	{
		me.data_entries_per_word = bits_per_word / me.data_bits_per_entry;
		_size = length(me, 0) / me.data_entries_per_word;
		if (_size < 1) _size = 1;
		me.data_max_sub_position = (me.data_entries_per_word - 1) * me.data_bits_per_entry;
	}
	else
	{
		me.data_entries_per_word = (me.data_bits_per_entry + bits_per_word - 1) / bits_per_word;
		_size = length(me, 0) * me.data_entries_per_word;
	}

	resize(me.data_factors, dimension(me));

	me.data_factors[0] = 1; //actually, this is not needed

	unsigned int i = 1;
	for (i = 1; i < dimension(me); ++i)
	{
		me.data_factors[i] = _size;
		_size *= length(me, i);
	}

	resize(me.data_dat, _size);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TCoordinates>
inline void
_CoordinatesToPositions(AlignTrace<TSize, TSpec> & me,
						TCoordinates const & coordinates_,
						TSize & position_,
						unsigned char & bit_position)
{
	for (unsigned int i=dimension(me)-1; i > 0; --i)
	{
		pos_ += me.data_factor[i] * coordinates[i];
	}
	pos_ += coordinates[0];

	if (me.data_has_few_bits_per_entry)
	{
		position_ = pos_ / me.data_entries_per_word;
		bit_position = me.data_bits_per_entry *(pos_ - position * me.data_entries_per_word);
	}
	else
	{
		position_ = pos_ * me.data_entries_per_word;
		bit_position = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TCoordinates>
inline void
_PositionsToCoordinates(AlignTrace<TSize, TSpec> & me,
						TCoordinates & coordinates_,
						TSize position_,
						unsigned char bit_position = 0)
{
	if (me.data_has_few_bits_per_entry)
	{
		position_ *= me.data_entries_per_word;
		position_ += bit_position / me.data_bits_per_entry;
	}
	else
	{
		position_ /= me.data_entries_per_word;
	}
	for (unsigned int i = dimension(me)-1; i > 0; --i)
	{
		TSize factor_ = me.data_factors[i];
		TSize coor_ = position_ / factor_;
		_coordinates[i] = _coor;
		position_ -= _coor * _factor;
	}

	dimension_[0] = position_;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, TSpec>
struct Value< AlignTrace<TSize, TSpec> >
{
	typedef bool Type;
};

//////////////////////////////////////////////////////////////////////////////

//???verschieben
template <typename T>
struct Navigator;

template <typename TSize, typename TSpec>
struct Navigator< AlignTrace<TSize, TSpec> >
{
	typedef Navi< AlignTrace<TSize, TSpec>, void> Type;
};


//////////////////////////////////////////////////////////////////////////////
// Navigator
//////////////////////////////////////////////////////////////////////////////

//???verschieben
template <typename TContainer, typename TSpec>
struct Navi;

template <typename TSize, typename TSpec, typename TSpec2>
struct Navi< AlignTrace<TSize, TSpec>, TSpec2>
{
	typedef AlignTrace<TSize, TSpec> TContainer;
	TContainer & data_container;
	TSize data_position;
	unsigned char data_bit_position;
};

//////////////////////////////////////////////////////////////////////////////

//???verschieben
template <typename TContainer, typename TSpec>
inline TContainer &
container(Navi<TContainer, TSpec> & me)
{
	return me.data_container;
}
template <typename TContainer, typename TSpec>
inline TContainer const &
container(Navi<TContainer, TSpec> const & me)
{
	return me.data_container;
}
template <typename TContainer, typename TSpec>
inline void
setContainer(Navi<TContainer, TSpec> & me,
			 TContainer & container_)
{
	me.data_container = container_;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline bool
getBit(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me)
{
	return (me.data_dat[me.data_position] >> me.data_bit_position) & 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
setBit(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me)
{
	me.data_dat[me.data_position] |= (1 << me.data_bit_position);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
clearBit(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me)
{
	me.data_dat[me.data_position] &= ~(1 << me.data_bit_position);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
goNext(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me,
	   unsigned int dimension_)
{
	if (dimension_)
	{
		me.data_position += me.data_factor[dimension_];
	}
	else
	{
		if (me.data_container.data_has_few_bits_per_entry)
		{
			me.data_bit_position += me.data_container.data_bits_per_entry;
			if (me.data_bit_position >= BitsPerValue<unsigned int>::VALUE)
			{
				me.data_bit_position = 0;
				++me.data_position;
			}
		}
		else
		{
			me.data_position += me.data_container.data_entries_per_word;
		}
	}	
}

//____________________________________________________________________________

template <typename TSize, typename TSpec, typename TSpec2>
inline Navi< AlignTrace<TSize, TSpec>, TSpec2> const &
operator ++ (Navi< AlignTrace<TSize, TSpec>, TSpec2> & me)
{
	goNext(me);
}

//todo??? postfix

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
goPrevious(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me,
		   unsigned int dimension_)
{
	if (dimension_)
	{
		me.data_position -= me.data_factor[dimension_];
	}
	else
	{
		if (me.data_container.data_has_few_bits_per_entry)
		{
			me.data_bit_position -= me.data_container.data_bits_per_entry;
			if (me.data_bit_position > BitsPerValue<unsigned int>::VALUE)
			{
				me.data_bit_position = me.data_container.data_max_sub_position;
				--me.data_position;
			}
		}
		else
		{
			me.data_position -= me.data_container.data_entries_per_word;
		}
	}	
}

//____________________________________________________________________________

template <typename TSize, typename TSpec, typename TSpec2>
inline Navi< AlignTrace<TSize, TSpec>, TSpec2> const &
operator -- (Navi< AlignTrace<TSize, TSpec>, TSpec2> & me)
{
	goPevious(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
moveNextBit(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me)
{
	if (me.data_bit_position >= BitsPerValue<unsigned int>::VALUE -1)
	{
		me.data_bit_position = 0;
		++me.data_position;
	}
	else
	{
		++me.data_bit_position;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSize, typename TSpec, typename TSpec2>
inline void
movePreviousBit(Navi< AlignTrace<TSize, TSpec>, TSpec2 > & me)
{
	if (me.data_bit_position > 0)
	{
		--me.data_bit_position;
	}
	else
	{
		me.data_bit_position = BitsPerValue<unsigned int>::VALUE -1;
		--me.data_position;
	}
}
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
