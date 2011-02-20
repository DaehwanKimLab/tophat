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
  $Id: align_cols_base.h 1769 2008-03-12 08:51:24Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_COLS_BASE_H
#define SEQAN_HEADER_ALIGN_COLS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Align Cols
//////////////////////////////////////////////////////////////////////////////
//implements Cols sequence for align implementations

/**
.Class.AlignCols:
..cat:Alignments
..summary:Pseudo columns container for row-based alignment classes.
..signature:AlignCols<TAlign>
..param.TAlign:Alignment type.
...metafunction:Metafunction.Host
..remarks:This class emulates a container of columns on alignment classes
  that store the alignment in a container of rows.
  Note that accessing a row-based alignment column-wise can be significantly
  slower than accessing the alignment row-wise.
..see:Class.Align
*/

template <typename TAlign>
struct AlignCols
{
//____________________________________________________________________________
	mutable TAlign * data_align;

	AlignCols():
		data_align(0)
	{
SEQAN_CHECKPOINT
	}

	AlignCols(AlignCols const & other):
		data_align(other.data_align)
	{
SEQAN_CHECKPOINT
	}

	AlignCols(TAlign & ali):
		data_align(& ali)
	{
SEQAN_CHECKPOINT
	}

	~AlignCols() 
	{
SEQAN_CHECKPOINT
	}

	AlignCols const & operator = (AlignCols const & other)
	{
SEQAN_CHECKPOINT
		data_align = other.data_align;
		return *this;
	}
//____________________________________________________________________________

	template <typename TPosition>
	inline typename Value<AlignCols>::Type
	operator [](TPosition _pos)
	{
SEQAN_CHECKPOINT
		return value(*this, _pos);
	}

	template <typename TPosition>
	inline typename Value<AlignCols const>::Type
	operator [](TPosition _pos) const
	{
SEQAN_CHECKPOINT
		return value(*this, _pos);
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cols.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Cols<Align<TSource, TSpec> >
{
	typedef AlignCols<Align<TSource, TSpec> > Type;
};
template <typename TSource, typename TSpec>
struct Cols<Align<TSource, TSpec> const>
{
	typedef AlignCols<Align<TSource, TSpec> const> Type;
};

//////////////////////////////////////////////////////////////////////////////
// Metafunctions for AlignCols
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.AlignCols

template <typename TAlign>
struct Host<AlignCols<TAlign> >
{
	typedef TAlign Type;
};
template <typename TAlign>
struct Host<AlignCols<TAlign> const>
{
	typedef TAlign Type;
};

//////////////////////////////////////////////////////////////////////////////
/**
.Spec.AlignColIterator:
..cat:Iterators
..summary:Iterator for @Class.AlignCols@ pseudo container.
..signature:Iter< TAlign, AlignColIterator<TSpec> >
..param.TSpec:Specialization tag.
..general:Class.Iter
..see:Class.AlignCols
*/


template <typename TSpec>
struct AlignColIterator;

///.Metafunction.Iterator.param.T.type:Class.AlignCols

template <typename TAlign, typename TIteratorSpec>
struct Iterator<AlignCols<TAlign>, TIteratorSpec>
{
	typedef Iter< TAlign, AlignColIterator<void> > Type;
};
template <typename TAlign, typename TIteratorSpec>
struct Iterator<AlignCols<TAlign> const, TIteratorSpec>
{
	typedef Iter< TAlign, AlignColIterator<void> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// Iterator is also used as Value

///.Metafunction.Value.param.T.type:Class.AlignCols

template <typename TAlign>
struct Value<AlignCols<TAlign> >:
	Iterator<AlignCols<TAlign>, Standard >
{
};
template <typename TAlign>
struct Value<AlignCols<TAlign> const>:
	Iterator<AlignCols<TAlign> const, Standard>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.AlignCols

template <typename TAlign>
struct Size<AlignCols<TAlign> >:
	Size<typename Row<TAlign>::Type>
{
};
template <typename TAlign>
struct Size<AlignCols<TAlign> const>:
	Size<typename Row<TAlign const>::Type>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Position.param.T.type:Class.AlignCols

template <typename TAlign>
struct Position<AlignCols<TAlign> >:
	Position<typename Row<TAlign>::Type>
{
};
template <typename TAlign>
struct Position<AlignCols<TAlign> const>:
	Position<typename Row<TAlign const>::Type>
{
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

///.Function.host.param.object.type:Class.AlignCols

template <typename TAlign>
inline typename Host<AlignCols<TAlign> >::Type &
host(AlignCols<TAlign> & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(me.data_align)
	return *me.data_align; 
}
template <typename TAlign>
inline typename Host<AlignCols<TAlign> const>::Type &
host(AlignCols<TAlign> const & me)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(me.data_align)
	return *me.data_align; 
}

//////////////////////////////////////////////////////////////////////////////

///.Function.iter.param.object.type:Class.AlignCols

template <typename TAlign, typename TPosition, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
iter(AlignCols<TAlign> & me,
	 TPosition pos_,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	return typename Iterator<AlignCols<TAlign>, Tag<TTag> const >::Type(host(me), pos_); 
}
template <typename TAlign, typename TPosition, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
iter(AlignCols<TAlign> const & me,
	 TPosition pos_,
	 Tag<TTag> const)
{
SEQAN_CHECKPOINT
	return typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type(host(me), pos_); 
}

//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.container.type:Class.AlignCols

template <typename TAlign, typename TPosition>
inline typename Value< AlignCols<TAlign> >::Type
value(AlignCols<TAlign> & me,
	  TPosition _pos)
{
SEQAN_CHECKPOINT
	return iter(me, _pos);
}
template <typename TAlign, typename TPosition>
inline typename Value< AlignCols<TAlign> const>::Type
value(AlignCols<TAlign> const & me,
	  TPosition _pos)
{
SEQAN_CHECKPOINT
	return iter(me, _pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.beginPosition.param.object.type:Class.AlignCols

template <typename TAlignCols>
inline typename Position<TAlignCols>::Type
beginPosition_AlignCols(TAlignCols const & me)
{
SEQAN_CHECKPOINT
	typedef typename Host<TAlignCols>::Type TAlign;
	typename Position<typename Rows<TAlign>::Type>::Type _i = length(rows(host(me)));

	if (!_i)
	{
		return 0;
	}

	--_i;
	typename Position<TAlignCols>::Type _pos = beginPosition(row(host(me), _i));

	while (_i > 0)
	{
		--_i;
		typename Position<TAlignCols>::Type _pos2 = beginPosition(row(host(me), _i));
		if (_pos2 < _pos)
		{
			_pos = _pos2;
		}
	}
	return _pos;
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
beginPosition(AlignCols<TAlign> const & me)
{
SEQAN_CHECKPOINT
	return beginPosition_AlignCols(me);
}
template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
beginPosition(AlignCols<TAlign> & me)
{
SEQAN_CHECKPOINT
	return beginPosition_AlignCols(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.begin.param.object.type:Class.AlignCols

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
begin(AlignCols<TAlign> & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, beginPosition(me), tag_); 
}
template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
begin(AlignCols<TAlign> const & me,
	  Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, beginPosition(me), tag_); 
}

//////////////////////////////////////////////////////////////////////////////

///.Function.endPosition.param.object.type:Class.AlignCols

template <typename TAlignCols>
inline typename Position<TAlignCols>::Type
endPosition_AlignCols(TAlignCols const & me)
{
SEQAN_CHECKPOINT
	typedef typename Host<TAlignCols>::Type TAlign;

	typename Position<typename Rows<TAlign>::Type>::Type _i = length(rows(host(me)));
	typename Position<TAlignCols>::Type _pos = 0;

	while (_i > 0)
	{
		--_i;
		typename Position<TAlignCols>::Type _pos2 = endPosition(row(host(me), _i));
		if (_pos2 > _pos)
		{
			_pos = _pos2;
		}
	}
	return _pos;
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
endPosition(AlignCols<TAlign> & me)
{
SEQAN_CHECKPOINT
	return endPosition_AlignCols(me);
}
template <typename TAlign>
inline typename Position<AlignCols<TAlign> const>::Type
endPosition(AlignCols<TAlign> const & me)
{
SEQAN_CHECKPOINT
	return endPosition_AlignCols(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.end.param.object.type:Class.AlignCols

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
end(AlignCols<TAlign> & me,
	 Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, endPosition(me), tag_); 
}
template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
end(AlignCols<TAlign> const & me,
	 Tag<TTag> const tag_)
{
SEQAN_CHECKPOINT
	return iter(me, endPosition(me), tag_); 
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign>
inline typename Size<AlignCols<TAlign> >::Type
length(AlignCols<TAlign> const & me)
{
SEQAN_CHECKPOINT
	return endPosition(me) - beginPosition(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign>
inline bool
operator ==(AlignCols<TAlign> const & left, 
			AlignCols<TAlign> const & right)
{
SEQAN_CHECKPOINT
	return left.data_align == right.data_align; 
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
