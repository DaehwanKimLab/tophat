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
  $Id: align_base.h 4687 2009-08-16 10:05:41Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_BASE_H
#define SEQAN_HEADER_ALIGN_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Cols:
..summary:Type of column container of an alignment.
..signature:Cols<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the columns of $T$.
*/
template <typename T> struct Cols;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Col:
..summary:Type of a column in an alignment.
..signature:Col<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The column type of $T$.
..remarks:The returned type is equivalent to $Value<Cols<T>::Type>::Type$.
..see:Metafunction.Cols
..see:Metafunction.Value
*/
template <typename T>
struct Col:
	Value<typename Cols<T>::Type>
{
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Rows:
..summary:Type of row container of an alignment.
..signature:Rows<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the rows of $T$.
..see:Metafunction.Cols
*/
template <typename T> struct Rows;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Row:
..summary:Type of a row in an alignment.
..signature:Row<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The row type of $T$.
..remarks:The returned type is equivalent to $Value<Rows<T>::Type>::Type$.
..see:Metafunction.Rows
..see:Metafunction.Value
*/
template <typename T>
struct Row:
	Value<typename Rows<T>::Type>
{
};

template <typename T>
struct Row<T const> {
	typedef typename Row<T>::Type const Type;
};



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Align
//////////////////////////////////////////////////////////////////////////////
//Default implementation: array of Gaps objects

/**
.Class.Align:
..cat:Alignments
..summary:An alignment of sequences.
..signature:Align<TSource, TSpec>
..param.TSource:Type of the ungapped sequences.
...metafunction:Metafunction.Source
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...XXdefault:@Spec.ArrayGaps@
..remarks:The default implementation of $Align$ stores the alignment in a set of @Class.Gaps.Gaps<TSource.TSpec>@ objects.
 Hence, the default implementation is row-based, so it will be faster to access the alignment row-wise than column-wise.
*/

template <typename TSource, typename TSpec = ArrayGaps>
class Align
{
public:
	typedef Gaps<TSource, TSpec> TGaps;
	typedef String<TGaps> TRows;
	typedef typename Size<TRows>::Type TRowsSize;

	TRows data_rows;

//____________________________________________________________________________
public:
	Align()
	{
SEQAN_CHECKPOINT
	}
	Align(Align const & _other):
		data_rows(_other.data_rows)
	{
SEQAN_CHECKPOINT
	}
	template <typename TString, typename TStringsetSpec>
	Align(StringSet<TString, TStringsetSpec> & stringset)
	{
SEQAN_CHECKPOINT
		setStrings(*this, stringset);
	}
	~Align()
	{
SEQAN_CHECKPOINT
	}

	Align const &
	operator = (Align const & _other)
	{
SEQAN_CHECKPOINT
		data_rows = _other.data_rows;
		return *this;
	}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
//ALIGN INTERFACE
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Value<Align<TSource, TSpec> >:
	Value<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct Value<Align<TSource, TSpec> const>:
	Value<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct GetValue<Align<TSource, TSpec> >:
	GetValue<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct GetValue<Align<TSource, TSpec> const>:
	GetValue<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Reference<Align<TSource, TSpec> >:
	Reference<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct Reference<Align<TSource, TSpec> const>:
	Reference<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

// struct Cols<Align>: see below (in align_cols_base)

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Rows.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Rows<Align<TSource, TSpec> >
{
	typedef String<Gaps<TSource, TSpec> > Type;
};
template <typename TSource, typename TSpec>
struct Rows<Align<TSource, TSpec> const>
{
	typedef String<Gaps<TSource, TSpec> > const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Source.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Source<Align<TSource, TSpec> >
{
	typedef TSource Type;
};
template <typename TSource, typename TSpec>
struct Source<Align<TSource, TSpec> const >
{
	typedef TSource Type;
};


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.StringSetType:
..summary:Return type of @Function.stringSet@ function. 
..signature:StringSetType<T>::Type
..param.T:Alignment data structure.
..param.T.Type:Spec.Alignment Graph
..param.T.Type:Class.Align
..returns.param.Type:A @Class.StringSet.string set@ type of a reference to a string set type.
..see:Function.stringSet
*/
template <typename T>
struct StringSetType;

template <typename TSource, typename TSpec>
struct StringSetType<Align<TSource, TSpec> >
{
	typedef StringSet<TSource, Dependent<> > Type;
};
template <typename TSource, typename TSpec>
struct StringSetType<Align<TSource, TSpec> const >
{
	typedef StringSet<TSource, Dependent<> > Type;
};


//Use StringSet<Owner> for SequenceGaps to store temporary source strings
template <typename TSource>
struct StringSetType<Align<TSource, SequenceGaps> >
{
	typedef StringSet<TSource, Owner<> > Type;
};
template <typename TSource>
struct StringSetType<Align<TSource, SequenceGaps> const>
{
	typedef StringSet<TSource, Owner<> > Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.rows:
..cat:Alignments
..summary:The container of rows in an alignment. 
..signature:Rows rows(align)
..param.align:An alignment.
...type:Class.Align
..returns:The container of rows in $align$. 
...metafunction:Metafunction.Rows
..see:Function.cols
..see:Metafunction.Rows
*/
template <typename TSource, typename TSpec>
inline typename Rows< Align<TSource, TSpec> >::Type &
rows(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_rows;
}
template <typename TSource, typename TSpec>
inline typename Rows< Align<TSource, TSpec> const >::Type &
rows(Align<TSource, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_rows;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.row:
..cat:Alignments
..summary:A row in an alignment. 
..signature:Row & row(align, position)
..param.align:An alignment.
...type:Class.Align
..param.position:A position in the @Function.rows@ container of $align$.
..returns:The row in @Function.rows@ container of $align$ at the given $position$. 
...metafunction:Metafunction.Row
..remarks:This function is equivalent to $value(rows(align), position)$.
..see:Function.rows
..see:Function.col
..see:Metafunction.Row
*/
template <typename TSource, typename TSpec, typename TPosition>
inline typename Row< Align<TSource, TSpec> >::Type &
row(Align<TSource, TSpec> & me, 
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(rows(me), _pos);
}
template <typename TSource, typename TSpec, typename TPosition>
inline typename Row< Align<TSource, TSpec> const>::Type &
row(Align<TSource, TSpec> const & me, 
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(rows(me), _pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.cols:
..cat:Alignments
..summary:The container of columns in an alignment. 
..signature:Cols cols(align)
..param.align:An alignment.
...type:Class.Align
..returns:The container of columns in $align$. 
...metafunction:Metafunction.Cols
..see:Metafunction.Cols
*/
template <typename TSource, typename TSpec>
inline typename Cols< Align<TSource, TSpec> >::Type
cols(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	return typename Cols< Align<TSource, TSpec> >::Type(me);
}
template <typename TSource, typename TSpec>
inline typename Cols< Align<TSource, TSpec> const>::Type
cols(Align<TSource, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return typename Cols< Align<TSource, TSpec> const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.col:
..cat:Alignments
..summary:A column in an alignment. 
..signature:Col & col(align, position)
..param.align:An alignment.
...type:Class.Align
..param.position:A position in the @Function.cols@ container of $align$.
..returns:The column in @Function.cols@ container of $align$ at the given $position$. 
...metafunction:Metafunction.Col
..remarks:This function is equivalent to $value(cols(align), position)$.
..see:Function.cols
..see:Metafunction.Col
*/
template <typename TSource, typename TSpec, typename TPosition>
inline typename Col< Align<TSource, TSpec> >::Type
col(Align<TSource, TSpec> & me,
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(cols(me), _pos);
}
template <typename TSource, typename TSpec, typename TPosition>
inline typename Col< Align<TSource, TSpec> const>::Type
col(Align<TSource, TSpec> const & me,
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(cols(me), _pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.detach.param.object.type:Class.Align

template <typename TSource, typename TSpec>
inline void
detach(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Iterator<TRows, Standard>::Type TRowsIterator;

	TRowsIterator it = begin(rows(me));
	TRowsIterator it_end = end(rows(me));

	while (it != it_end)
	{
		detach(*it);
		++it;
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSource, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Align<TSource, TSpec> const & source,
	  TIDString const &,
	  Raw)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	TRowsPosition row_count = length(rows(source));
	TPosition begin_ = beginPosition(cols(source));
	TPosition end_ = endPosition(cols(source));
	
	unsigned int baseCount=0;
	unsigned int leftSpace=6;
	while(begin_ < end_) {
		unsigned int windowSize_ = 50;
		if ((begin_ + windowSize_)>end_) windowSize_=end_ - begin_;

		// Print header line
		unsigned int offset=0;
		if (baseCount != 0) offset = (unsigned int) floor(log((double)baseCount) / log((double)10));
		for(unsigned int j = 0;j<leftSpace-offset;++j) {
			_streamPut(target, ' ');
		}
		_streamPutInt(target, baseCount);
		baseCount+=windowSize_;
		_streamPut(target, ' ');
		for(TPosition i = 1;i<=windowSize_;++i) {
			if ((i % 10)==0) _streamPut(target, ':');
			else if ((i % 5)==0) _streamPut(target, '.');
			else _streamPut(target, ' ');
		}
		_streamPut(target, ' ');
		_streamPut(target, '\n');
		
		// Print sequences
		for(TRowsPosition i=0;i<2*row_count-1;++i) {
			for(unsigned int j = 0;j<leftSpace+2;++j) _streamPut(target, ' ');
			if ((i % 2)==0) {
				TRow& row_ = row(source, i/2);
				typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
				TIter begin1_ = iter(row_, begin_);
				TIter end1_ = iter(row_, begin_ + windowSize_);
				for (; begin1_ != end1_; ++begin1_) {
					if (isGap(begin1_)) _streamPut(target, gapValue<char>());
					else _streamPut(target, *begin1_);
				}
				//_streamWriteRange(target, iter(row_, begin_), iter(row_, begin_ + windowSize_));
			} else {
				for(unsigned int j = 0;j<windowSize_;++j) {
					if ((!isGap(row(source, (i-1)/2), begin_+j)) &&
						(!isGap(row(source, (i+1)/2), begin_+j)) &&
						(row(source, (i-1)/2)[begin_+j]==row(source, (i+1)/2)[begin_+j]))
					{
						_streamPut(target, '|');
					} else {
						_streamPut(target, ' ');
					}
				} 
			}
			_streamPut(target, '\n');
		}
		_streamPut(target, '\n');
		begin_+=50;
	}
	_streamPut(target, '\n');
}

//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Align<TSource, TSpec> const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator >> (TStream & source, 
			 Align<TSource, TSpec> & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
*/
//////////////////////////////////////////////////////////////////////////////

/**
.Function.setStrings:
..cat:Alignments
..summary:Loads the sequences of a stringset into an alignment.
..signature:setStrings(align, stringset)
..param.align:An alignment.
...type:Class.Align
..param.stringset:A string set.
...type:Class.StringSet
..remarks:The function clears $align$ and creates an new global alignment between strings in $stringset$ that contains only trainling gaps.
The alignment will be dependent from the strings in the stringset; use @Function.detach@ to make $align$ the owner of its strings.
*/

template <typename TSource, typename TSpec, typename TSpec2>
inline void
setStrings(Align<TSource, TSpec> & me,
		   StringSet<TSource, TSpec2> & stringset)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> TAlign;
	typedef StringSet<TSource, TSpec2> TStringset;

	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Iterator<TRows>::Type TRowsIterator;
	typedef typename Size<TStringset>::Type TStringsetSize;

	clear(me.data_rows);
	resize(me.data_rows, length(stringset));
	TRowsIterator it = begin(rows(me));
	TStringsetSize stringset_length = length(stringset);
	for (TStringsetSize i = 0; i < stringset_length; ++i)
	{
		setSource(*it, value(stringset, i));
		++it;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec>
inline void
clearGaps(Align<TSource, TSpec> & me)
{
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Iterator<TRows>::Type TRowsIterator;
	
	for (TRowsIterator it = begin(rows(me)); it != end(rows(me)); goNext(it))
	{
		clearGaps(*it);
	}
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.align#stringSet
..param.g.type:Class.Align
*/
template <typename TSource, typename TSpec>
inline typename StringSetType<Align<TSource, TSpec> >::Type 
stringSet(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> TAlign;
	typedef typename StringSetType<TAlign>::Type TStringSet;

	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Iterator<TRows>::Type TRowsIterator;

	TStringSet ss;

	for (TRowsIterator it = begin(rows(me)); it != end(rows(me)); goNext(it))
	{
		appendValue(ss, source(*it));
	}
	return ss;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
