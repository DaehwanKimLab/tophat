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
  $Id: find_base.h 2377 2008-06-07 08:10:34Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_PATTERN_BASE_H
#define SEQAN_HEADER_FIND_PATTERN_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pattern:
..summary:Holds the needle and preprocessing data (depends on algorithm).
..cat:Searching
..signature:Pattern<TNeedle[, TSpec]>
..param.TNeedle:The needle type.
...type:Class.String
..param.TSpec:The online-algorithm to search with.
...remarks:Leave empty for index-based pattern matching (see @Class.Index@).
...default:The result of @Metafunction.DefaultPattern@
..remarks:If $TNeedle$ is a set of strings, then $position(pattern)$ returns the index of the currently matching needle.
*/

template < typename TNeedle, typename TSpec = typename DefaultPattern<TNeedle>::Type >
class Pattern;

//default implementation
template < typename TNeedle >
class Pattern<TNeedle, void>
{
public:
	typedef typename Position<TNeedle>::Type TNeedlePosition;

	Holder<TNeedle> data_host;
	TNeedlePosition data_begin_position;
	TNeedlePosition data_end_position;

	Pattern() {}

	template <typename _TNeedle>
	Pattern(_TNeedle & ndl):
		data_host(ndl) {}

	template <typename _TNeedle>
	Pattern(_TNeedle const & ndl):
		data_host(ndl) {}

};
//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TSpec>
struct Container< Pattern<TNeedle, TSpec> > {
	typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Container< Pattern<TNeedle, TSpec> const > {
	typedef TNeedle const Type;
};

template <typename TNeedle, typename TSpec>
struct Host< Pattern<TNeedle, TSpec> > {
	typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Host< Pattern<TNeedle, TSpec> const > {
	typedef TNeedle const Type;
};


template <typename TPattern, typename TSpec>
struct Value< Pattern<TPattern, TSpec> > {
	typedef typename Value<TPattern>::Type Type;
};

template <typename TPattern, typename TSpec>
struct Position< Pattern<TPattern, TSpec> > {
	typedef typename Position<TPattern>::Type Type;
};

template <typename TPattern, typename TSpec>
struct Difference< Pattern<TPattern, TSpec> > {
	typedef typename Difference<TPattern>::Type Type;
};

template <typename TPattern, typename TSpec>
struct Size< Pattern<TPattern, TSpec> > {
	typedef typename Size<TPattern>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.ScoringScheme:
..summary:Returns the scoring scheme of an approximate searching algorithm.
..cat:Searching
..signature:ScoringScheme<TPattern>::Type
..param.TPattern:A @Class.Pattern@ type.
...type:Class.Pattern
..returns:The scoring scheme.
...default:@Shortcut.EditDistanceScore@
*/

template <typename TNeedle>
struct ScoringScheme
{
	typedef EditDistanceScore Type;
};
template <typename TNeedle>
struct ScoringScheme<TNeedle const>:
	ScoringScheme<TNeedle>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> & 
_dataHost(Pattern<TNeedle, TSpec> & me) 
{ 
	return me.data_host;
}
template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> & 
_dataHost(Pattern<TNeedle, TSpec> const & me) 
{
	return const_cast<Holder<TNeedle> &>(me.data_host);
}

//host access: see basic_host.h


//???TODO: Diese Funktion entfernen! (sobald setHost bei anderen pattern nicht mehr eine Art "assignHost" ist)
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, TSpec> & me,
		TNeedle2 const & ndl) 
{
	 me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void 
setHost(Pattern<TNeedle, TSpec> & me,
		TNeedle2 & ndl) 
{ 
	 me.data_host = ndl; //assign => Pattern haelt eine Kopie => doof!
}
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type & 
beginPosition(Pattern<TNeedle, TSpec> & me) 
{
	return me.data_begin_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type & 
beginPosition(Pattern<TNeedle, TSpec> const & me) 
{
	return me.data_begin_position;
}


template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setBeginPosition(Pattern<TNeedle, TSpec> & me, 
				 TPosition _pos) 
{
	me.data_begin_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> >::Type & 
endPosition(Pattern<TNeedle, TSpec> & me) 
{
	return me.data_end_position;
}
template <typename TNeedle, typename TSpec>
inline typename Position<Pattern<TNeedle, TSpec> const >::Type & 
endPosition(Pattern<TNeedle, TSpec> const & me) 
{
	return me.data_end_position;
}

template <typename TNeedle, typename TSpec, typename TPosition>
inline void
setEndPosition(Pattern<TNeedle, TSpec> & me, 
			   TPosition _pos) 
{
	me.data_end_position = _pos;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type 
segment(Pattern<TNeedle, TSpec> & me) 
{
	typedef typename Infix<TNeedle>::Type TInfix;
	return TInfix(host(me), me.data_begin_position, me.data_end_position);
}
template <typename TNeedle, typename TSpec>
inline typename Infix<TNeedle>::Type 
segment(Pattern<TNeedle, TSpec> const & me) 
{
	typedef typename Infix<TNeedle>::Type TInfix;
	return TInfix(host(me), me.data_begin_position, me.data_end_position);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.host.param.object.type:Class.Pattern

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, TSpec> >::Type & 
host(Pattern<TNeedle, TSpec> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

template <typename TNeedle, typename TSpec>
inline typename Host<Pattern<TNeedle, TSpec> const>::Type & 
host(Pattern<TNeedle, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.needle:
..summary:Returns the needle of a @Class.Pattern@ object (not implemented for some online-algorithms).
..cat:Searching
..signature:needle(pattern)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..returns:The needle object to search for.
..remarks:The result type is @Metafunction.Needle@$<TPattern>::Type$ for pattern of type $TPattern$.
*/

template < typename TObject >
inline typename Needle<TObject>::Type &
needle(TObject &obj) 
{
	return obj;
}

template < typename TObject >
inline typename Needle<TObject const>::Type &
needle(TObject const &obj) 
{
	return obj;
}


///.Function.position.param.iterator.type:Class.Pattern

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> >::Type &
needle(Pattern<TNeedle, TSpec> & obj) 
{
	return host(obj);
}

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern<TNeedle, TSpec> const>::Type &
needle(Pattern<TNeedle, TSpec> const & obj) 
{
	return host(obj);
}

/**
.Function.setNeedle:
..summary:Sets the needle of a @Class.Pattern@ object and optionally induces preprocessing.
..cat:Searching
..signature:setNeedle(pattern, needle)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..param.needle:The needle object to search for.
...type:Class.String
*/

template < typename TNeedle, typename TSpec >
inline void
setNeedle(Pattern<TNeedle, TSpec> &obj, TNeedle const &ndl) {
	setHost(obj, ndl);
}


//____________________________________________________________________________

/**.Function.scoringScheme
..cat:Searching
..summary:The @glos:scoring scheme@ used for finding or aligning.
..signature:scoringScheme(obj)
..param.obj:Object that holds a @glos:scoring scheme@
...type:Class.Pattern
..returns:The @glos:scoring scheme@ used in $obj$
...default:@Shortcut.EditDistanceScore@
..see:glos:scoring scheme
..see:Metafunction.ScoringScheme
*/

template <typename TNeedle, typename TSpec>
inline typename ScoringScheme<Pattern<TNeedle, TSpec> >::Type 
scoringScheme(Pattern<TNeedle, TSpec> &)
{
SEQAN_CHECKPOINT
	return typename ScoringScheme<Pattern<TNeedle, TSpec> >::Type();
}
template <typename TNeedle, typename TSpec>
inline typename ScoringScheme<Pattern<TNeedle, TSpec> const>::Type 
scoringScheme(Pattern<TNeedle, TSpec> const &)
{
SEQAN_CHECKPOINT
	return typename ScoringScheme<Pattern<TNeedle, TSpec> const>::Type();
}

//____________________________________________________________________________

/**.Function.setScoringScheme
..cat:Searching
..summary:Sets the @glos:scoring scheme@ used for finding or aligning.
..signature:setScoringScheme(obj, score)
..param.obj:Object that holds a @glos:scoring scheme@.
...type:Class.Pattern
..param.score:The new @glos:scoring scheme@ used by $obj$.
..see:glos:scoring scheme
..see:Function.scoringScheme
*/

template <typename TNeedle, typename TSpec, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, TSpec> & me, 
				 TScore2 & score)
{
//dummy implementation for compatibility reasons
}
//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
