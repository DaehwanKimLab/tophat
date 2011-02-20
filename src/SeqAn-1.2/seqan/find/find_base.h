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
  $Id: find_base.h 4743 2009-09-01 07:34:56Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BASE_H
#define SEQAN_HEADER_FIND_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Tags

struct FindInfix; //find needle as a substring of haystack. This is the default
struct FindPrefix; //find needle as a prefix of haystack. (prefix search)

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultFinder:
..cat:Searching
..summary:Default @Class.Finder@ specialization type.
..signature:DefaultFinder<THaystack>::Type
..param.THaystack:The given haystack type.
..returns:Is $void$ by default and @Tag.Index Find Algorithm.ESA_FIND_MLR@ if $THaystack$ is an @Class.Index@.
*/
template < typename TObject >
struct DefaultFinder 
{
	typedef void Type;
};

/**
.Metafunction.DefaultPattern:
..cat:Searching
..summary:Default @Class.Pattern@ specialization type.
..signature:DefaultPattern<TNeedle>::Type
..param.TNeedle:The given needle type.
..returns:Is $void$ by default.
*/
template < typename TObject >
struct DefaultPattern 
{
	typedef void Type;
};

/**
.Metafunction.Haystack:
..summary:Returns the haystack type of a @Class.Finder@ type.
..cat:Searching
..signature:Haystack<TFinder>::Type
..param.TFinder:A @Class.Finder@ type.
...type:Class.Finder
..returns:The haystack type of $TFinder$, i.e. $THaystack$ for $Finder<THaystack, TSpec>$.
*/

template <typename TFinder>
struct Haystack 
{
	typedef typename Container<TFinder>::Type Type;
};

/**
.Metafunction.Needle:
..summary:Returns the needle type of a @Class.Pattern@ type.
..cat:Searching
..signature:Needle<TPattern>::Type
..param.TPattern:A @Class.Pattern@ type.
...type:Class.Pattern
..returns:The needle type of $TPattern$, i.e. $TNeedle$ for $Pattern<TNeedle, TSpec>$.
*/

template <typename TPattern>
struct Needle 
{
	typedef typename Host<TPattern>::Type Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> > 
{
	typedef Segment<THost, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Needle<Segment<THost, TSpec> const> 
{
	typedef Segment<THost, TSpec> const Type;
};


//////////////////////////////////////////////////////////////////////////////

/**
.Function.find:
..summary:Search for a @Class.Pattern@ in a @Class.Finder@ object.
..cat:Searching
..signature:find(finder, pattern)
..signature:find(finder, pattern, k)
..param.finder:The @Class.Finder@ object to search through.
...remarks:For online-algorithm $patterns$, finder can also be an arbitrary @Concept.Rooted Iterator@.
...type:Class.Finder
...type:Concept.Rooted Iterator
..param.pattern:The @Class.Pattern@ object to search for.
...remarks:For index $finders$, pattern can also be a Sequence.
...type:Class.Pattern
..param.k:Desired minimal score (for approximate matching).
...remarks:$k$ has to be a number <= 0.
...remarks:Differences are deletions, insertions and substitutions.
..returns:$boolean$ that indicates whether an occurence of $pattern$ was found or not.
..remarks:Repeated calls of this function iterate through all occurences of $pattern$.
*/

/**
.Class.Finder:
..summary:Holds the haystack and a current search context.
..cat:Searching
..signature:Finder<THaystack[, TSpec]>
..param.THaystack:The haystack type.
...type:Class.String
...type:Class.Index
..param.TSpec:The index-algorithm to search with (Optional).
...default:The result of @Metafunction.DefaultFinder@
...remarks:Leave empty for online pattern matching (see @Class.Pattern@).
...remarks:If $THaystack$ is an @Class.Index@, then $TSpec$ specifies the index search algorithm.
..remarks:$position(finder)$ returns the position of the current hit in the haystack.
If $THaystack$ is a set of strings or an index of a set of strings, then $position(finder)$ returns a @Class.Pair@ $(hayNo, pos)$,
in which $hayNo$ is the haystack index and $pos$ the local position of the hit.
..remarks:Use $clear(finder)$ to reset a finder object and search from the beginning.
*/

///.Function.clear.param.object.type:Class.Finder
///.Function.position.param.iterator.type:Class.Finder

template < typename THaystack, typename TSpec = typename DefaultFinder<THaystack>::Type >
class Finder
{
	typedef typename Iterator<THaystack, Rooted>::Type TIterator;
	typedef typename Position<THaystack>::Type TPosition;
	typedef typename Size<THaystack>::Type TSize;

public:
	TIterator data_iterator;
	TPosition data_endPos; //note: we need this since iterator could point to begin or end (depending on pattern type)
	TSize data_length;
	bool _needReinit;					// if true, the Pattern needs to be reinitialized
	bool _beginFind_called;					// if false, then findBegin was not yet called for this match position (see findBegin default implementation)

	Finder()
		: data_endPos(0)
		, data_length(0)
		, _needReinit(true)
		, _beginFind_called(false)
	{}

	Finder(THaystack & haystack)
		: data_iterator(begin(haystack, Rooted()))
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}

	Finder(TIterator &iter)
		: data_iterator(iter)
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}

	Finder(TIterator const &iter)
		: data_iterator(iter)
		, data_endPos(0)
		, data_length(0)
		, _needReinit(true) 
		, _beginFind_called(false)
	{}

	Finder(Finder const &orig)
		: data_iterator(orig.data_iterator)
		, data_endPos(orig.data_endPos)
		, data_length(orig.data_length)
		, _needReinit(orig._needReinit) 
		, _beginFind_called(orig._beginFind_called)
	{}

	~Finder() {}

//____________________________________________________________________________

	Finder const &
	operator = (Finder const & other)
	{
		data_iterator = other.data_iterator;
		data_endPos = other.data_endPos;
		data_length = other.data_length;
		_needReinit = other._needReinit;
		_beginFind_called = other._beginFind_called;
		return *this;
	}

//____________________________________________________________________________

	inline typename Reference<TIterator>::Type 
	operator* () 
	{
SEQAN_CHECKPOINT
		return value(hostIterator(*this));
	}

	inline typename Reference<TIterator const>::Type 
	operator* () const
	{
SEQAN_CHECKPOINT
		return value(hostIterator(*this));
	}

//____________________________________________________________________________

	operator TIterator () const
	{
SEQAN_CHECKPOINT
		return data_iterator;
	}

//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me)
{//shortcut: move end position to iterator position +1
	me._beginFind_called = false;
	me.data_endPos = position(me)+1;
}
template <typename THaystack, typename TSpec, typename TPosition>
inline void
_setFinderEnd(Finder<THaystack, TSpec> & me,
			  TPosition end_pos)
{
	me._beginFind_called = false;
	me.data_endPos = end_pos;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec, typename TSize>
inline void
_setFinderLength(Finder<THaystack, TSpec> & me,
				 TSize _length)
{
	me.data_length = _length;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.beginPosition.param.object.type:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type 
beginPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos - me.data_length;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
beginPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos - me.data_length;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.begin.param.object.type:Class.Finder

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> & me,
	  Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), beginPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
begin(Finder<THaystack, TSpec> const & me,
	  Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), beginPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.endPosition.param.object.type:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Position<THaystack>::Type 
endPosition(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<THaystack const>::Type
endPosition(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_endPos;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.end.param.object.type:Class.Finder

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> & me,
	Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), endPosition(me), tag);
}

template <typename THaystack, typename TSpec, typename TTag>
inline typename Iterator<THaystack const, Tag<TTag> const>::Type
end(Finder<THaystack, TSpec> const & me,
	Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
	return iter(haystack(me), endPosition(me), tag);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.length.param.object.type:Class.Finder

template <typename THaystack, typename TSpec>
inline typename Size<THaystack>::Type 
length(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_length;
}
template <typename THaystack, typename TSpec>
inline typename Size<THaystack const>::Type 
length(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

/**.Function.finder#infix:
..description:Segment of the last found match in haystack.
..signature:Infix infix(finder)
..param.finder:An online finder.
...type:Class.Finder
..returns:An @Metafunction.Infix.infix@ of the @Metafunction.Haystack.haystack@.
...metafunction:Metafunction.Infix
..remarks:This function works only correct if the begin position of the match was already found,
see @Function.findBegin@.
*/
template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	typedef typename Infix<THaystack>::Type TInfix;
	return TInfix(haystack(me), beginPosition(me), endPosition(me));
}

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack const>::Type
infix(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	typedef typename Infix<THaystack const>::Type TInfix;
	return TInfix(haystack(me), beginPosition(me), endPosition(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec>
inline typename _Parameter<THaystack>::Type 
host(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename _Parameter<THaystack>::Type 
host(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename _Parameter<THaystack>::Type 
container(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename _Parameter<THaystack>::Type 
container(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return container(hostIterator(me));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void
setHost(Finder<THaystack, TSpec> & me, 
		typename _Parameter<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
	setContainer(hostIterator(me), container_);
	goBegin(me);
}

template <typename THaystack, typename TSpec>
inline void
setContainer(Finder<THaystack, TSpec> & me, 
			 typename _Parameter<THaystack>::Type container_)
{
SEQAN_CHECKPOINT
	setContainer(hostIterator(me), container_);
	goBegin(me);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type &
hostIterator(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}

template <typename THaystack, typename TSpec>
inline typename Iterator<THaystack, Rooted>::Type const &
hostIterator(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_iterator;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
empty(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me._needReinit;
}

template <typename THaystack, typename TSpec>
inline void
clear(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	me._needReinit = true;
}

//____________________________________________________________________________

template <typename T>
inline void
_finderSetNonEmpty(T & me)
{
SEQAN_CHECKPOINT
	goBegin(me);
}


template <typename THaystack, typename TSpec>
inline void
_finderSetNonEmpty(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	me._needReinit = false;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline bool
atBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return (!empty(me) && atBegin(hostIterator(me)));
}

template <typename THaystack, typename TSpec>
inline bool
atEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	return (!empty(me) && atEnd(hostIterator(me)));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void
goBegin(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	//_finderSetNonEmpty(me);
	goBegin(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline void
goEnd(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	//_finderSetNonEmpty(me);
	goEnd(hostIterator(me));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	if (empty(me)) return 0;
	return position(hostIterator(me));
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, TSpec> >::Type
position(Finder<THaystack, TSpec> const & me)
{
SEQAN_CHECKPOINT
	if (empty(me)) return 0;
	return position(hostIterator(me));
}

//____________________________________________________________________________
/**
.Function.setPosition:
..cat:Searching
..summary:Sets the position of a finder.
..signature:setPosition(finder, pos)
..param.finder:A finder.
...class:Class.Finder
..param.pos:A position.
...metafunction:Metafunction.Position
..see:Function.position
*/

template <typename THaystack, typename TSpec, typename TPosition>
inline void 
setPosition(Finder<THaystack, TSpec> & me, TPosition pos_)
{
SEQAN_CHECKPOINT
	setPosition(hostIterator(me), pos_);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator--(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
	--hostIterator(me);
	return me;
}

template <typename THaystack, typename TSpec>
inline Finder<THaystack, TSpec> &
operator++(Finder<THaystack, TSpec> & me)
{
SEQAN_CHECKPOINT
/*			if (beforeBegin()) {
		goBegin(hostIterator(me));
	} else*/
		++hostIterator(me);
	return me;
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator + (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
	return Finder<THaystack, TSpec>(hostIterator(left) + right);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator += (Finder<THaystack, TSpec> & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) += right;
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> const
operator - (Finder<THaystack, TSpec> const & left, TIntegral right)
{
SEQAN_CHECKPOINT
	return Finder<THaystack, TSpec>(hostIterator(left) - right);
}

template <typename THaystack, typename TSpec, typename TIntegral>
inline typename Difference<Finder<THaystack, TSpec> const>::Type
operator - (Finder<THaystack, TSpec> const & left, Finder<THaystack, TSpec> const & right)
{
SEQAN_CHECKPOINT
	return hostIterator(left) - hostIterator(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TSpec, typename TIntegral>
inline Finder<THaystack, TSpec> &
operator -= (Finder<THaystack, TSpec> & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	hostIterator(left) -= right;
	return left;
}

//____________________________________________________________________________


/**
.Function.setHaystack:
..summary:Sets the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:setHaystack(finder, haystack)
..param.finder:The @Class.Finder@ object to search with.
...type:Class.Finder
..param.haystack:The haystack object the finder searches through.
...type:Class.String
*/

template < typename THaystack, typename TSpec >
inline void
setHaystack(Finder<THaystack, TSpec> &obj, THaystack const &hstk) {
	setHost(obj, hstk);
}

/**
.Function.haystack:
..summary:Returns the haystack of a @Class.Finder@ object.
..cat:Searching
..signature:haystack(finder)
..param.finder:The @Class.Finder@ object to search through.
...type:Class.Finder
..returns:The haystack object.
..remarks:The result type is @Metafunction.Haystack@$<TFinder>::Type$ for finder of type $TFinder$.
*/

template < typename TObject >
inline typename Haystack<TObject>::Type &
haystack(TObject &obj) {
	return container(obj);
}

template < typename TObject >
inline typename Haystack<TObject const>::Type &
haystack(TObject const &obj) {
	return container(obj);
}



//////////////////////////////////////////////////////////////////////////////


template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> > {
	typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Container< Finder<THaystack, TSpec> const> {
	typedef THaystack const Type;
};

template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> > {
	typedef THaystack Type;
};

template <typename THaystack, typename TSpec>
struct Host< Finder<THaystack, TSpec> const> {
	typedef THaystack const Type;
};


template <typename THaystack, typename TSpec>
struct Value< Finder<THaystack, TSpec> > {
	typedef typename Value<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Position< Finder<THaystack, TSpec> >:
	Position<THaystack> {};

template <typename THaystack, typename TSpec>
struct Difference< Finder<THaystack, TSpec> > {
	typedef typename Difference<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec>
struct Size< Finder<THaystack, TSpec> > {
	typedef typename Size<THaystack>::Type Type;
};

template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec>, TIteratorSpec >
{
	typedef typename Iterator<THaystack>::Type Type;
};
template <typename THaystack, typename TSpec, typename TIteratorSpec>
struct Iterator< Finder<THaystack, TSpec> const, TIteratorSpec >
{
	typedef typename Iterator<THaystack>::Type const Type;
};


template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> >:
	DefaultGetIteratorSpec< THaystack >
{
};
template <typename THaystack, typename TSpec>
struct DefaultGetIteratorSpec< Finder<THaystack, TSpec> const>:
	DefaultGetIteratorSpec< THaystack const>
{
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
