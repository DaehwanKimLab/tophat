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
  $Id: global_seed_chain.h 3361 2009-02-05 16:38:14Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GLOBAL_SEED_CHAIN_H
#define SEQAN_HEADER_GLOBAL_SEED_CHAIN_H

namespace SEQAN_NAMESPACE_MAIN
{

//Changed version of skiplist find so that the biggest element smaller than the searched one is found
template <typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
_findPrev(Map<TValue, Skiplist<TSpec> > & me,
	 TFind const & _find, //can be a TKey or a SkiplistElement or GoEnd
	 SkiplistPath<TValue, TSpec> & path) 
{
	typedef typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type TIterator;

	_skiplistFind(me, _find, path);
	return TIterator(path.data_elements[0]);
}

//Changed version of skiplist find so that the biggest element smaller than the searched one is found
template <typename TValue, typename TSpec, typename TFind>
inline typename Iterator< Map<TValue, Skiplist<TSpec> > >::Type
_findPrev(Map<TValue, Skiplist<TSpec> > & me,
	 TFind const & _find) //can be a TKey or a SkiplistElement or GoEnd
{
	typedef SkiplistPath<TValue, TSpec> TPath;
	TPath path;
	return _findPrev(me, _find, path);
}

/*DISABLED
.Function.globalChaining
..summary:Global chaining of a set of seeds.
..cat:Seed Handling
..cat:Chaining
..signature:globalChaining(source, result)
..signature:globalChaining(source, result, gapCost, xLength, yLength)
..param.source: The set of seeds to chain.
...type:Spec.Scored SeedSet
..param.result: Container in which the result should be stored. The chain is in reversed order.
..param.gapCost: Gap cost value.
..param.xLength: Length of the first sequence.
..param.yLength: Length of the second sequence.
..returns: The score of the chain.
*/

/**
.Function.globalChaining:
..summary:Computes the chain on a set of fragments.
..cat:Chaining
..signature:globalChaining(source, dest [, score] [, algorithm])
..param.source:The set of fragments.
...remarks:This could be either a container of seeds or a @Class.SeedSet@ object.
..param.dest:Container in which the result should be stored.
..param.score:The gap scoring scheme.
...default:@Spec.Score Zero@
...class:Spec.Score Manhattan
...class:Spec.Score Zero
...class:Spec.Score ChainSoP
..param.algorithm:A tag that identifies the algorithm which is used for chaining.
...default:$Default$
...value:$Default$: Compiler selects best algorithm.
...value:$GenericChaining$: A simple generic chaining algorithm.
...value:$RangetreeChaining$: An elaborated chaining algorithm for @Spec.Score Zero@, @Spec.Score Manhattan@, and @Spec.Score ChainSoP@ scoring schemes.
*/
template<typename TValue, typename TSeedSpec, typename TScoreSpec, typename TSpec, typename TTargetContainer>
typename ScoreType<TScoreSpec>::Type
globalChaining(SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> const &source,	//Seedset containing the seeds
			   TTargetContainer &result)									//container for the result chain
{
	SEQAN_CHECKPOINT
	typedef typename ScoreType<TScoreSpec>::Type TScore;
	typedef SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const , Standard>::Type TSeedSetIterator;
	typedef Triple<TSeedSetIterator, TScore, void*> TChainElement;
	typedef Map<Pair<TValue, TChainElement*> > TSkiplist; //sorted by y-coordinate
	typedef typename Iterator<TSkiplist, Standard>::Type TSkiplistIterator;

	typedef std::multimap<TValue, Pair<bool, TChainElement*> > TMultiMap;
	typedef typename TMultiMap::iterator TMapIterator;
	TSkiplist list;
	TChainElement* pElement =0;
	add(list,-1,pElement);

	TMultiMap pointArray, test; //sorted by x-coodinate
	
	TSeedSetIterator it_end = end(source);
	for (TSeedSetIterator it = begin(source); it != it_end; ++it)
	{
		pElement = new TChainElement(it, seedScore(it), 0);
		pointArray.insert(std::make_pair(leftDim0(*it), Pair<bool, TChainElement*>(true, pElement)));
		test.insert(std::make_pair(rightDim0(*it), Pair<bool, TChainElement*>(false, pElement)));
	}

	TMapIterator it_map_end2 = test.end();
	for (TMapIterator it = test.begin(); it != it_map_end2; ++it)
	{
		pointArray.insert(std::make_pair(it->first, it->second));
	}

	TSkiplistIterator it_test = begin(list);

	TMapIterator it_map_end = pointArray.end();
	for (TMapIterator it = pointArray.begin(); it != it_map_end; ++it)
	{
		if (it->second.i1)
		{
			TSkiplistIterator it_skip = _findPrev(list, leftDim1(*it->second.i2->i1));
			if (it_skip != it_test){
				it->second.i2->i2 += ((*it_skip).i2)->i2;   //score
				it->second.i2->i3 = (*it_skip).i2;			//predecessor
			}
		} 
		else
		{
			insertTriple(list, it->second.i2);
		}

	}
	TSkiplistIterator it_skip = _findPrev(list, supremumValue<TValue>());
	pElement = (*it_skip).i2;
	TScore best = pElement->i2;
	TChainElement* delete_pointer;

	while (pElement != 0)
	{
		appendValue(result, *pElement->i1);
		delete_pointer = pElement;
		pElement = (TChainElement*)pElement->i3;
		delete(delete_pointer);
	}
	
	::std::reverse(begin(result), end(result));
	return best;
}

template<typename TValue, typename TValue2, typename TSeedSpec, typename TScoreSpec, typename TSpec, typename TTargetContainer, typename TScore>
TScore
globalChaining(SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> const &source, //Seedset containing the seeds
				TTargetContainer &result,	//container for the result chain
				TScore gapCost,				//Value for gap costs
				TValue2 xLength,				//length of the first sequence
				TValue2 yLength)				//length of the second sequence
{
	SEQAN_CHECKPOINT
	gapCost *= -1;
	typedef SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const, Standard>::Type TSeedSetIterator;
	typedef Triple<TSeedSetIterator, TScore, void*> TChainElement;
	typedef Map<Pair<TValue, TChainElement*> > TSkiplist; //sorted by y-coordinate
	typedef typename Iterator<TSkiplist, Standard>::Type TSkiplistIterator;
	typedef std::multimap<TValue, Pair<bool, TChainElement*> > TMultiMap;
	typedef typename TMultiMap::iterator TMapIterator;

	TValue total = xLength + yLength;
	
	TSkiplist list;
	TChainElement * pElement =0;
	//insert
	add(list,-3,pElement);
	TMultiMap pointArray, test; //sorted by x-coodinate

	TSeedSetIterator it_end = end(source);
	for (TSeedSetIterator it = begin(source); it != it_end; ++it)
	{
		pElement = new TChainElement(it, seedScore(it)+ (rightDim0(*it) + rightDim1(*it) - leftDim0(*it) - leftDim1(*it)+2-total)*gapCost, 0);
		pointArray.insert(make_pair(leftDim0(*it), Pair<bool, TChainElement*>(true, pElement)));
		test.insert(make_pair(rightDim0(*it), Pair<bool, TChainElement*>(false, pElement)));
	}

	TMapIterator it_map_end2 = test.end();
	for (TMapIterator it = test.begin(); it != it_map_end2; ++it)
		pointArray.insert(make_pair(it->first, it->second));
	
	TSkiplistIterator it_test = begin(list);
	TMapIterator it_map_end = pointArray.end();
	for (TMapIterator it = pointArray.begin(); it != it_map_end; ++it)
	{
		if (it->second.i1)
		{
			TSkiplistIterator it_skip = _findPrev(list, leftDim1(*it->second.i2->i1));
			if (it_skip != it_test){
				it->second.i2->i2 += ((*it_skip).i2)->i2 + total*gapCost;	//score
				it->second.i2->i3 = (*it_skip).i2;					//predecessor
			}
		} 
		else
		{
			insertTriple(list, it->second.i2);
		}

	}

	TSkiplistIterator it_skip = _findPrev(list, total);
	pElement = (*it_skip).i2;
	TScore best = pElement->i2;
	
	TChainElement* delete_pointer;
	while (pElement != 0)
	{
		appendValue(result, *pElement->i1);
		delete_pointer = pElement;
		pElement = (TChainElement*)pElement->i3;
		delete(delete_pointer);
	}

	::std::reverse(begin(result), end(result));
	return best;
}

template<typename TValue, typename TSeedSpec, typename TScoreSpec, typename TSpec, typename TTargetContainer, typename TScore>
TScore
globalChaining(SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> const &source, //Seedset containing the seeds
				TTargetContainer &result,	//container for the result chain
				TScore gapCost)				//Value for gap costs
{
	typedef SeedSet<TValue, TSeedSpec, TScoreSpec, TSpec> TSeedSet;
	typedef typename Iterator<TSeedSet const, Standard>::Type TSeedSetIterator;
	TSeedSetIterator it_end = end(source);
	TValue x_length = 0;
	TValue y_length = 0;
	for (TSeedSetIterator it = begin(source); it != it_end; ++it)
	{
		TValue right = rightPosition(*it, 0);
		if (right > x_length) x_length = right;

		right = rightPosition(*it, 1);
		if (right > y_length) y_length = right;
	}
	return globalChaining(source, result, gapCost, x_length, y_length);
}


template<typename TValue, typename TChainElement>
void
insertTriple(Map<Pair<TValue, TChainElement*> > &list, 
			 TChainElement* pElement)
{
	typedef typename Iterator<Map<Pair<TValue, TChainElement*> > >::Type TIterator;
	TIterator it = _findPrev(list, rightDim1(*pElement->i1)+1);
	TIterator it_begin = begin(list);
	if (it != it_begin)
	{
		if (pElement->i2 > ((*it).i2)->i2)
		{
			insert(list, rightDim1(*pElement->i1), pElement);
			TIterator it_tmp = find(list, rightDim1(*pElement->i1));
			TIterator it_end = end(list);
			TIterator del;
			while (it_tmp != it_end)
			{
				if ((*it_tmp).i2->i2 < pElement->i2)
				{
					del = it_tmp;
					++it_tmp;
					erase(list, del);
				} else
					++it_tmp;
			}
		} 
		else
			delete(pElement);
	} else {
		add(list, rightDim1(*pElement->i1), pElement);
	}
}



}

#endif //#ifndef SEQAN_HEADER_...
