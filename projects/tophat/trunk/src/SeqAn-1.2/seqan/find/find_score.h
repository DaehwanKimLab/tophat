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
  $Id: find_score.h 3704 2009-03-20 12:34:55Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SCORE_H
#define SEQAN_HEADER_FIND_SCORE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// DPSearch
//////////////////////////////////////////////////////////////////////////////

//template <typename TScore, typename TSpec = FindInfix, typename TFindBeginPatternSpec = DPSearch<TScore, FindPrefix, void> >
template <typename TScore, typename TSpec = FindInfix, typename TFindBeginPatternSpec = typename DefaultFindBeginPatternSpec<TScore>::Type>
struct DPSearch;

/*
.Spec.DPSearch:
..cat:Searching
..general:Class.Pattern
..summary:A dynamic programming algorithm for approximate string-matching with a user-definable scoring function.
..signature:Pattern<TNeedle, DPSearch<TScore [, TSpec [, TFindBeginPatternSpec] ]> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TScore:The scoring function.
...type:Class.Score
..remarks.text:The algorithm is based on the Sellers/Needleman-Wunsch dynamic progamming algorithm. 
The $Pattern$ object only contains the right-most column of the dynamic programming matrix.
...note:At the moment, the algorithm only works on linear gap costs.
*/

///.Class.Pattern.param.TSpec.type:Class.Score


template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
class Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> >:
	public _FindBegin<Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > >
{
public:
	typedef typename Value<TScore>::Type TScoreValue;

	Holder<TNeedle>		data_host;
	TScore				data_score;
	TScoreValue			data_limit;
	String<TScoreValue>	data_tab;
	TScoreValue			data_maxscore;  //score of the needle matching itself (needed for banding FindPrefix)

public: 
	Pattern(): 
		data_limit(0)
	{ 
SEQAN_CHECKPOINT
	}

	Pattern(TNeedle & _needle, 
			TScore & _score_func, 
			TScoreValue _limit = 0): 
		data_score(_score_func),
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		setHost(*this, _needle);
	}

	Pattern(TNeedle & _needle,
			TScoreValue _limit = 0): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		setHost(*this, _needle);
	}

	Pattern(TScoreValue _limit): 
		data_limit(_limit)
	{ 
SEQAN_CHECKPOINT
		create(data_score);
	}

	Pattern(Pattern const & other): 
		data_host( other.data_host ),
		data_score( other.data_score ), 
		data_limit( other.data_limit ),
		data_tab( other.data_tab ),
		data_maxscore( other.data_maxscore)
	{
SEQAN_CHECKPOINT
	}

	inline Pattern & 
	operator = (Pattern const & other) 
	{ 
SEQAN_CHECKPOINT
		this->data_host = other.data_host;
		this->data_score = other.data_score;
		this->data_limit = other.data_limit;
		this->data_tab = other.data_tab;
		this->data_maxscore = other.data_maxscore;

		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
struct ScoringScheme <Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > >
{
	typedef TScore Type;
};


//DEPRECATED
//.Metafunction.ScoreValue.param.T.type:Spec.DPSearch
//template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
//struct ScoreValue <Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > >:
//	Value<TScore>
//{
//};


//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > > 
{
	typedef TFindBeginPatternSpec Type;
};
template <typename TNeedle, typename TScore, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, DPSearch<TScore, FindPrefix, TFindBeginPatternSpec> > >
{// no find begin for FindPrefix
	typedef void Type;
};

//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > >::Type & 
host(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const>::Type & 
host(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> >  const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_host);
}


//____________________________________________________________________________

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void 
setHost(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me, 
		TNeedle2 & ndl)
{
	me.data_host = ndl;
	clear(me.data_tab);
}
template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec, typename TNeedle2>
void 
setHost(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me, 
		TNeedle2 const & ndl)
{
	me.data_host = ndl;
	clear(me.data_tab);
}


//____________________________________________________________________________


template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline TScore const & 
scoringScheme(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me)
{
SEQAN_CHECKPOINT
	return me.data_score;
}

//____________________________________________________________________________


template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me, 
				 TScore2 & score)
{
SEQAN_CHECKPOINT
	me.data_score = score;
	clear(me.data_tab);
}
template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec, typename TScore2>
inline void
setScoringScheme(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me, 
				 TScore2 const & score)
{
SEQAN_CHECKPOINT
	me.data_score = score;
	clear(me.data_tab);
}

//____________________________________________________________________________


/**.Function.scoreLimit
..cat:Searching
..summary:The minimal score a match must reach in approximate searching.
..signature:scoreLimit(pattern)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..returns:The current score limit of $pattern$.
*/

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline typename Value<TScore>::Type 
scoreLimit(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > const & me)
{
SEQAN_CHECKPOINT
	return me.data_limit;
}

//____________________________________________________________________________

/**.Function.setScoreLimit
..cat:Searching
..summary:Sets the minimal score a match must reach in approximate searching.
..signature:setScoreLimit(pattern, limit)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..param.limit:The new score limit.
..see:Function.scoreLimit
*/

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void 
setScoreLimit(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
	me.data_limit = _limit;
}

//____________________________________________________________________________
// returns the score of the last hit position found (note:position = end of occurrence in haystack)

/**.Function.getScore
..cat:Searching
..summary:Score of the last found match in approximate searching.
..signature:getScore(pattern)
..param.pattern:A @Concept.Pattern|pattern@ that can be used for approximate searching.
...type:Spec.DPSearch
..returns:The score of the last match found using $pattern$.
...remarks:If no match was found, the value is undefined.
..see:Function.scoreLimit
..see:Function.setScoreLimit
..see:Function.find
*/

template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline typename Value<TScore>::Type
getScore(Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me)
{
	return front(me.data_tab);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline void _patternInit (Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me) 
{
	typedef Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > TPattern;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TPattern>::Type TSize;
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TIterator;
	typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;

	TScore const & scoring = scoringScheme(me);
	TScoreValue score_gap = scoreGapExtend(scoring);

	TTab & string_tab = me.data_tab;

	//allocate enough memory for one column of DP matrix
	TSize need_length = length(needle(me));
	SEQAN_ASSERT(need_length);

	resize(string_tab, need_length);
	SEQAN_ASSERT(length(string_tab) >= need_length);

//	if (length(_dataNeedle(me)) < got_length) throw(0); //???TODO: Throw "not enough memory" exception

	//init matrix
	//note: The column is stored in reverse order
	TIterator tab_end = begin(string_tab, Standard());
	TIterator tab = end(string_tab, Standard());
	
	TScoreValue x = score_gap;

	while (tab > tab_end)
	{
		--tab;
		*tab = x;
		x += score_gap;
	}

	if (TYPECMP<TSpec, FindPrefix>::VALUE)
	{//compute data_maxscore
		me.data_maxscore = 0;
		TNeedleIterator it = begin(needle(me), Standard());
		TNeedleIterator it_end = end(needle(me), Standard());
		for (; it != it_end; ++it)
		{
			me.data_maxscore += score(scoring, *it, *it);
		}

	}

	_findBeginInit(me);
}



//////////////////////////////////////////////////////////////////////////////
// find, findNext
//////////////////////////////////////////////////////////////////////////////

//proportional gap cost: Needleman-Wunsch 

//???TODO: Ukkonen trick?
//???TODO: Finder for affine gap costs?

template <typename TFinder, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
bool 
_find_score_simple_proportional(TFinder & finder, Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TTab;
	typedef typename Iterator<TTab, Standard>::Type TTabIterator;
	typedef typename Iterator<TNeedle const, Standard>::Type TNeedleIterator;
	typedef typename Value<typename Haystack<TFinder>::Type>::Type THaystackValue;
	typedef typename Haystack<TFinder>::Type THaystack;
	typedef typename Size<THaystack>::Type TSize;

	String<TScoreValue> & string_tab = me.data_tab;
	TScore const & scoring = scoringScheme(me);

	TScoreValue score_gap = scoreGapExtend(scoring);
//	TScoreValue score_match = scoreMatch(scoringScheme(me));
//	TScoreValue score_mismatch = scoreMismatch(scoringScheme(me));

	TSize prefix_begin_position;	

	if (empty(finder))
	{
		clear(me.data_tab);
		_patternInit(me);
		_finderSetNonEmpty(finder);
		prefix_begin_position = position(finder);
	}
	else
	{
		goNext(finder);
		prefix_begin_position = beginPosition(finder);
	}

	
	TSize haystack_length = length(container(hostIterator(finder)));

	//limit search width for prefix search
	if (TYPECMP<TSpec, FindPrefix>::VALUE && (score_gap < 0))
	{
		TSize maxlen = prefix_begin_position + length(needle(me)) + ((scoreLimit(me) - me.data_maxscore) / score_gap) + 1;
		if (haystack_length > maxlen)
		{
			haystack_length = maxlen;
		}
	}

	//start searching

	TTabIterator tab_begin = end(string_tab, Standard());

	TNeedle const & ndl = needle(me);
	TNeedleIterator it_begin = begin(ndl, Standard());
	TNeedleIterator it_end = end(ndl, Standard());

	//for each character in haystack, do...
	for (; position(finder) < haystack_length; ++finder)
	{
		//get character
		THaystackValue c = *finder;

		//init some variables
		TNeedleIterator it = it_begin;
		TScoreValue * tab = tab_begin;
		TScoreValue h = (TYPECMP<TSpec, FindPrefix>::VALUE) ? score_gap * (position(finder)-prefix_begin_position) : 0;
		TScoreValue v = (TYPECMP<TSpec, FindPrefix>::VALUE) ? h + score_gap : 0;

		//fill the column
		while (it < it_end)
		{
			--tab; //note: the column is stored in "reverse order"

//			TScoreValue m2 = (c == *it) ? h + score_match : h + score_mismatch;
//char d = *it;
			TScoreValue m2 = h + score(scoring, c, *it);
			h = *tab;
			TScoreValue m1 = (h > v) ? h + score_gap : v + score_gap;

			v = (m1 > m2) ? m1 : m2;
			*tab = v;

			++it;
		}

		if (*tab >= scoreLimit(me) )
		{//found a hit
			_setFinderEnd(finder);
			if (TYPECMP<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
			return true;
		}

	}

	//found nothing
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline bool 
find(TFinder & finder, 
	 Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me)
{
	SEQAN_ASSERT(scoreGapOpen(scoringScheme(me)) == scoreGapExtend(scoringScheme(me))) //this finder is only defined for linear gap costs
	return _find_score_simple_proportional(finder, me);
}

template <typename TFinder, typename TNeedle, typename TScore, typename TSpec, typename TFindBeginPatternSpec>
inline bool 
find(TFinder & finder, 
	 Pattern<TNeedle, DPSearch<TScore, TSpec, TFindBeginPatternSpec> > & me,
	 int const limit_)
{
	SEQAN_ASSERT(scoreGapOpen(scoringScheme(me)) == scoreGapExtend(scoringScheme(me))) //this finder is only defined for linear gap costs
	setScoreLimit(me, limit_);
	return _find_score_simple_proportional(finder, me);
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
