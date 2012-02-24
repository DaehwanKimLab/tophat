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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

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
_coordinatesToPositions(AlignTrace<TSize, TSpec> & me,
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
_positionsToCoordinates(AlignTrace<TSize, TSpec> & me,
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
