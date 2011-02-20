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
  $Id: graph_impl_interval_tree.h 1757 2008-02-27 16:26:20Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_H
#define SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree
//////////////////////////////////////////////////////////////////////////////


/**
.Class.IntervalTree:
..cat:Miscellaneous
..summary:A datastructure that efficiently stores intervals.
..signature:IntervalTree<TValue, TCargo>
..param.TValue:The value type.
..param.TCargo:The cargo/id type.
...default:int
...remarks:If the intervals are not associated with cargos/IDs, they will be numbered consecutively.
*/
template<typename TValue = int, typename TCargo = unsigned int>
class IntervalTree
{
public:
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalAndCargo<TValue,TCargo> TInterval;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;

	TGraph g;
	TPropertyMap pm;
	size_t interval_counter;
	
	
	IntervalTree()
	{
SEQAN_CHECKPOINT
		interval_counter = 0;
	}
	
	
	template<typename TIterator,typename TCargoIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends, 
				 TCargoIterator interval_cargos, 
				 size_t len)	
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
		{
			intervals[i].i1 = value(interval_begins);
			++interval_begins;
			intervals[i].i2 = value(interval_ends);
			++interval_ends;
			intervals[i].cargo = value(interval_cargos);
			++interval_cargos;
			++i;
		}
		interval_counter = len;
		
		createIntervalTree(g,pm,intervals);
	}

	template<typename TIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends,
				 size_t len)
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
		{
			intervals[i].i1 = value(interval_begins);
			++interval_begins;
			intervals[i].i2 = value(interval_ends);
			++interval_ends;
			intervals[i].cargo = i;
			++i;
		}
		interval_counter = len;
		createIntervalTree(g,pm,intervals);
	}
	

	IntervalTree(String<TInterval> intervals)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals);
	}

	IntervalTree(IntervalTree const & other):
		g(other.g),
		pm(other.pm),
		interval_counter(other.interval_counter)
	{
SEQAN_CHECKPOINT
	}


	IntervalTree & operator = (IntervalTree const & other)
	{
SEQAN_CHECKPOINT
		g = other.g;
		pm = other.pm;
		interval_counter = other.interval_counter;
		return *this;
	}


	~IntervalTree()
	{
SEQAN_CHECKPOINT
		
	}

};


//template<typename TValue, typename TIterator,typename TCargoIterator, typename TCargo>
//IntervalTree<TValue,TCargo>::IntervalTree<TValue,TCargo>(TIterator interval_begins,
//				 TIterator interval_ends, 
//				 TCargoIterator interval_cargos, 
//				 size_t len)
//	{
//		String<TInterval> intervals;
//		resize(intervals,len);
//		size_t i = 0;
//		while(i<len)
//		{
//			intervals[i].i1 = value(interval_begins);
//			++interval_begins;
//			intervals[i].i2 = value(interval_ends);
//			++interval_ends;
//			intervals[i].cargo = value(interval_cargos);
//			++interval_cargos;
//			++i;
//		}
//		interval_counter = len;
//		createIntervalTree(g,pm,intervals);
//	}
// 




	
///////Specs for the way interval centers are determined
//center = minbegin + (maxend-minbegin)/2
//template <typename TSpec = SpecPointAndCargo>
struct TagComputeCenter_;
typedef Tag<TagComputeCenter_> const ComputeCenter;

//center = length(sequence)/2
//template <typename TSpec = SpecPointAndCargo>
struct TagMidCenter_;
typedef Tag<TagMidCenter_> const MidCenter;

//center = center of random interval
//template <typename TSpec = SpecPointAndCargo>
struct TagRandomCenter_;
typedef Tag<TagRandomCenter_> const RandomCenter;


template<typename TIntervals, typename TIntervalPointers>
void
_makePointerInterval(TIntervals & intervals,TIntervalPointers & interval_pointers)
{
SEQAN_CHECKPOINT
	typedef typename Value<TIntervalPointers>::Type TIntervalPointer;
	typedef typename Iterator<TIntervalPointers, Rooted>::Type TIntervalPointerIterator;

	TIntervalPointer it;
	TIntervalPointerIterator iit = begin(interval_pointers);
	if(length(intervals)>0)
		for(it = &intervals[0]; it <= &intervals[length(intervals)-1]; ++it)
		{
			*iit = it;	
			++iit;
		}

}


//
//center of root node is computed by _calcIntervalTreeRootCenter
template<typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
void 
createIntervalTree(TGraph & g, 
				   TPropertyMap & pm, 
				   TIntervals & intervals,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TIntervals>::Type TInterval;
	typedef typename Value<TInterval>::Type TValue;

	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));
	
	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);
	
	TValue center =	_calcIntervalTreeRootCenter(intervals);
	
	std::sort(begin(intervals),end(intervals),_less_compI1_ITree<TInterval>);

	String<TInterval*> interval_pointers;
	resize(interval_pointers,length(intervals));

	_makePointerInterval(intervals,interval_pointers);

	_createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
		
	reserve(pm, length(pm), Exact());
	reserve(g.data_vertex, length(g.data_vertex), Exact());

}


// most user friendly interval tree construction for the moment...
// CompCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, 
				   TPropertyMap & pm, 
				   TIntervals & intervals)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,ComputeCenter());
}


//center of root is specified by user
template<typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
void 
createIntervalTree(TGraph & g,
 TPropertyMap & pm, 
				   TIntervals & intervals,
				   typename Value<typename Value<TIntervals>::Type>::Type center,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TIntervals>::Type TInterval;
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	
	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));
	
	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);
	
	TInterval a;
	typename Iterator<TIntervals, Standard>::Type begin_ = begin(intervals);
	typename Iterator<TIntervals, Standard>::Type end_ = end(intervals);
	std::sort(begin_, end_ ,_less_compI1_ITree<TInterval>);

	String<TInterval*> interval_pointers;
	resize(interval_pointers,length(intervals));

	_makePointerInterval(intervals,interval_pointers);

	if(length(intervals)==1)
		center = (rightBoundary(intervals[0])-leftBoundary(intervals[0]))/(TValue)2.0;

	_createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
		
	reserve(pm, length(pm), Exact());
	reserve(g.data_vertex, length(g.data_vertex), Exact());

}

// CompCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   TIntervals & intervals, 
				   typename Value<typename Value<TIntervals>::Type>::Type center)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,center,ComputeCenter());
}



//////////////////////////////////////////////////////////////////////////////
//remembers minimum and maximum of point values in intervals and sets the center
//of each node to min+(max-min)/2
template<typename TGraph, typename TPropertyMap, typename TIntervalPointer, typename TValue>
void 
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TIntervalPointer*> & intervals,
				   typename VertexDescriptor<TGraph>::Type & knot, 
				   TValue, 
				   TValue center, 
				   typename VertexDescriptor<TGraph>::Type len,
				   Tag<TagComputeCenter_> const tag)
{
SEQAN_CHECKPOINT
	//IntervalAndCargo<int,Fragment<>* >* checkOne = getValue(intervals,0); 
	//std::cout <<leftBoundary(*checkOne)<<"->\n";
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*intervals[0]);
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TIntervalPointer*> TIntervalPointers;

	TIntervalPointers S_left;
	TIntervalPointers S_right;

	TValue min1 = 300000000;
	TValue min2 = 300000000;
	TValue max1 = 0;
	TValue max2 = 0;

	value(pm,knot).center = center;
	
 
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	while(it != it_end)
	{
		if((**it).i2<=center)
		{
			appendValue(S_left,*it);
			if((**it).i2 > max1)
				max1 = (**it).i2;
			if((**it).i1 < min1)
				min1 = (**it).i1;
		}
		else
		{
			if((**it).i1>center)
			{
				appendValue(S_right,(*it));
				if((**it).i2 > max2)
					max2 = (**it).i2;
				if ((**it).i1 < min2)
					min2 = (**it).i1;
			}
			else
			{
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_left,vd,center,min1+(max1-min1)/2,length(S_left),tag);
	}
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_right,vd,center,min2+(max2-min2)/2,length(S_right),tag);
	}
}




//////////////////////////////////////////////////////////////////////////////
//createIntervalTree for all specs except CompCenter, the center value of each 
//node is determined by functions _calcIntervalTreeNodeCenterLeft and 
//_calcIntervalTreeNodeCenterRight
template<typename TGraph, typename TPropertyMap, typename TSpec, typename TInterval, typename TValue>
void 
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TInterval*> & intervals, 
				   typename VertexDescriptor<TGraph>::Type & knot, 
				   TValue last_center, TValue center, 
				   typename VertexDescriptor<TGraph>::Type len,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT

	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*value(intervals,0));
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TInterval*> TIntervalPointers;
	
	TIntervalPointers S_left;
	TIntervalPointers S_right;
		
	value(pm,knot).center = center;
	
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	while(it != it_end)
	{
		if((**it).i2<=center)
		{
			appendValue(S_left,*it);
		}
		else
		{
			if((**it).i1>center)
			{
				appendValue(S_right,(*it));
			}
			else
			{
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterLeft(S_left,last_center,center,tag);
		_createIntervalTree(g,pm,S_left,vd,center,next_center,length(S_left),tag);
	}
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterRight(S_right,last_center,center,tag);
		_createIntervalTree(g,pm,S_right,vd,center,next_center,length(S_right),tag);
	}
}




//the RandomCenter spec way of chosing center values:
//pick a random interval from the list and take its center as the center value 
//for the left child node (during interval tree construction)
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterLeft(TIntervals & intervals, TValue &, TValue &, Tag<TagRandomCenter_> const)
{
SEQAN_CHECKPOINT
	TValue rand_index = rand()%length(intervals);  
	return (rightBoundary(*value(intervals,rand_index))+leftBoundary(*value(intervals,rand_index)))/(TValue)2.0;
}

//the RandomCenter spec way of chosing center values:
//pick a random interval from the list and take its center as the center value 
//for the right child node (during interval tree construction)
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterRight(TIntervals & intervals, TValue &, TValue &, Tag<TagRandomCenter_> const)
{
SEQAN_CHECKPOINT
	TValue rand_index = rand()%length(intervals);  
	return (rightBoundary(*value(intervals,rand_index))+leftBoundary(*value(intervals,rand_index)))/(TValue)2.0;
}

//the MidCenter spec way of chosing center values:
//simply take the middle 
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterLeft(TIntervals &, TValue & last_center, TValue & center, Tag<TagMidCenter_> const)
{
SEQAN_CHECKPOINT
	if (center > last_center)
		return (center - (center-last_center)/(TValue)2.0);
	else
		return (center - (last_center-center)/(TValue)2.0);
}

//the MidCenter spec way of chosing center values:
//picks a random interval from the list and takes its center as the center value 
//for the right child node (during interval tree construction)
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterRight(TIntervals &, TValue & last_center, TValue & center, Tag<TagMidCenter_> const)
{
SEQAN_CHECKPOINT
	if (center > last_center)
		return (center + (center-last_center)/(TValue)2.0);
	else
		return (center + (last_center-center)/(TValue)2.0);
}




template<typename TIntervals>
typename Value<typename Value<TIntervals>::Type>::Type
_calcIntervalTreeRootCenter(TIntervals & intervals)
{
SEQAN_CHECKPOINT
	
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	typedef typename Iterator<TIntervals,Standard>::Type TIntervalIterator;

	TIntervalIterator it = begin(intervals);
	TIntervalIterator it_end = end(intervals);

	TValue min = (TValue)300000000.0;//maxTValue;
	TValue max = 0;

	while(it != it_end)
	{
		if(leftBoundary(*it)<min) min = leftBoundary(*it);
		if(rightBoundary(*it)>max) max = rightBoundary(*it);
		++it;
	}
	
	return (min+(max-min)/(TValue)2.0);

}




template<typename TGraph, typename TPropertyMap, typename TInterval>
void
addInterval(TGraph & g, TPropertyMap & pm, TInterval interval)
{
SEQAN_CHECKPOINT

	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Value<TInterval>::Type TValue;
	typedef typename ListType<TProperty>::Type TList;
	

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < leftBoundary(interval))
		{
			if(atEnd(it)){
				TVertexDescriptor vd = addVertex(g);
				resizeVertexMap(g,pm);
				addEdge(g,act_knot,vd);
				_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
				break;
			}
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)){
						TVertexDescriptor vd = addVertex(g);
						resizeVertexMap(g,pm);
						addEdge(g,act_knot,vd);
						_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
						break;
					}
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(rightBoundary(interval) <= act_prop.center)
			{
				if(atEnd(it)){
					TVertexDescriptor vd = addVertex(g);
					resizeVertexMap(g,pm);
					addEdge(g,act_knot,vd);
					_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
					break;
				}
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)){
							TVertexDescriptor vd = addVertex(g);
							resizeVertexMap(g,pm);
							addEdge(g,act_knot,vd);
							_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
							break;
						}
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				_appendIntervalTreeNodeLists(property(pm, act_knot),interval);
				std::sort(begin(property(pm,act_knot).list1),end(property(pm,act_knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
				std::sort(begin(property(pm,act_knot).list2),end(property(pm,act_knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);
				break;
			}
		}
	}

}

template<typename TValue, typename TCargo, typename TInterval>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TInterval interval)
{
SEQAN_CHECKPOINT

	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end, TCargo cargo)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = cargo;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = itree.interval_counter;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

//template<typename TValue, typename TCargo>
//void
//addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end)
//{
//
//	IntervalAndCargo<TValue,TCargo> interval;
//	interval.i1 = begin;
//	interval.i2 = end;
//	interval.cargo = itree.interval_counter;
//	++itree.interval_counter;
//	addInterval(itree.g,itree.pm,interval);
//
//}

//
template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, TPropertyMap & pm, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < query)
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)));
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query < act_prop.center)
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) <= query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)));
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)));
				break;
			}
		}
	}

}

template<typename TValue,typename TCargo>
void
findIntervals(IntervalTree<TValue,TCargo> & it, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervals(it.g,it.pm,query,result);

}




//
template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervalsExcludeTouching(TGraph & g, TPropertyMap & pm, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename Iterator<TGraph, OutEdgeIterator >::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	
	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if( (TValue) act_prop.center < query)
		{
			int i = 0;
			while(i < (int) length(act_prop.list2) && (TValue) rightBoundary(value(act_prop.list2,i)) > query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)));
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query < (TValue) act_prop.center)
			{
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)));
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)));
					++i;
				}
				break;
			}
		}
	}

}


template<typename TValue,typename TCargo>
void
findIntervalsExcludeTouching(IntervalTree<TValue,TCargo> & it, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervalsExcludeTouching(it.g,it.pm,query,result);

}


//find overlapping intervals... under construction and far from being done
//template<typename TGraph, typename TPropertyMap, typename TValue,typename TInterval>
//void
//findIntervals(TGraph & g, TPropertyMap & pm, TInterval query, String<TInterval> & result)
//{
//
//	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
//	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
//	typedef typename Value<TPropertyMap>::Type TProperty;
//	
//	// start at root
//	TVertexDescriptor act_knot = 0;
//	TProperty act_prop = property(pm,act_knot);
//	TProperty next_prop;
//	
//	while(true)
//	{
//		TOutEdgeIterator it(g, act_knot);
//		act_prop = property(pm,act_knot);
//		if(act_prop.center < leftBoundary(query))
//		{
//			int i = 0;
//			while(i < length(act_prop.list2) && rightBoundary(act_prop.list2[i]) > leftBoundary(query))
//			{
//				appendValue(result,act_prop.list2[i]);
//				++i;	
//			}
//			if(atEnd(it)) break;
//			else{
//				next_prop = property(pm,targetVertex(it));
//				if(next_prop.center <= act_prop.center)
//				{
//					goNext(it);
//					if(atEnd(it)) break;
//				}
//			}
//			act_knot = targetVertex(it);
//		}
//		else{
//			if(query < act_prop.center)
//			{
//				int i = 0;
//				while(i < length(act_prop.list1) && leftBoundary(act_prop.list1[i]) < rightBoundary(query))
//				{
//					appendValue(result,act_prop.list1[i]);
//					++i;
//				}
//				if(atEnd(it)) break;
//				else
//				{
//					next_prop = property(pm,targetVertex(it));
//					if(next_prop.center >= act_prop.center)
//					{
//						goNext(it);
//						if(atEnd(it)) break;
//					}
//				}
//				act_knot = targetVertex(it);
//			}
//			else{
//				append(result, act_prop.list1);
//				break;
//			}
//		}
//	}
//
//}






/////////////////// Metafunctions ///////////////////////
template<typename TValue, typename TCargo>
struct Value<IntervalTree<TValue,TCargo> >
{
	typedef TValue Type;
};

template<typename TValue, typename TCargo>
struct Cargo<IntervalTree<TValue,TCargo> >
{
	typedef TCargo Type;
};







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
