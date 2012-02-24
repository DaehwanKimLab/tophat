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
// Author: Anne-Katrin Emde <emde@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_MISC_INTERVAL_TREE_H
#define SEQAN_HEADER_MISC_INTERVAL_TREE_H

#include <seqan/graph_types.h>


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree Types
//////////////////////////////////////////////////////////////////////////////

///---------------------------------------------------------------///

//////////////////// Interval and ID type ///////////////////
/**
.Class.IntervalAndCargo:
..cat:Miscellaneous
..summary:A simple record type that stores an interval and a cargo value.
..signature:IntervalAndCargo<TValue, TCargo>
..param.TValue:The value type, that is the type of the interval borders.
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:The cargo type.
...default:int.
...metafunction:Metafunction.Cargo
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue = int, typename TCargo = int>
class IntervalAndCargo
{
public:
    /**
.Memvar.IntervalAndCargo#i1:
..class:Class.PointAndCargo
..summary:The first element in the interval of type i1.
     */
	TValue i1;

    /**
.Memvar.IntervalAndCargo#i2:
..class:Class.PointAndCargo
..summary:The last element in the interval of type i2.
     */
	TValue i2;

    /**
.Memvar.IntervalAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..signature:IntervalAndCargo()
     */
    IntervalAndCargo()
    {
SEQAN_CHECKPOINT
    }

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..class:Class.IntervalAndCargo
..summary:Constructor.
..signature:IntervalAndCargo(i1, i2, cargo)
..param.i1:The first element in the interval, of type TValue.
..param.i2:The last element in the interval of type TValue.
..param.cargo:The cargo value of type TCargo.
     */
	IntervalAndCargo(TValue i1, TValue i2, TCargo cargo):
		i1(i1), i2(i2), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};



/////////////////////// Point and ID type ////////////////
/**
.Class.PointAndCargo:
..cat:Miscellaneous
..summary:Simple record class storing a point (one-value interval) and a cargo.
..signature:PointAndCargo<TValue, TCargo>
..param.TValue:
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:
...default:int.
...metafunction:Metafunction.Value
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue=int, typename TCargo=int>
class PointAndCargo {
public:
    /**
.Memvar.PointAndCargo#point:
..class:Class.PointAndCargo
..summary:The stored point of type TValue.
     */
	TValue point;

    /**
.Memvar.PointAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..signature:PointAndCargo(point, cargo)
    */
	PointAndCargo() {
SEQAN_CHECKPOINT
	}

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..summary:Constructor.
..signature:PointAndCargo(point, cargo)
..param.point:
...summary:The point to store of type TValue.
..param.cargo:
...summary:The cargo to store of type TCargo.
    */
	PointAndCargo(TValue point, TCargo cargo):
		point(point), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};

///////////////////////////////////////////////////////////////////////////
/////////////////////////// IntervalTreeNode	///////////////////////////

/**
.Tag.IntervalTree Node Types
..summary:Tags to select the node type for @Class.IntervalTree@.
..cat:Miscellaneous

..tag.StorePointsOnly:The tree nodes store points.
..include:seqan/misc/misc_interval_tree.h
*/
struct StorePointsOnly {};


///..tag.StoreIntervals:The tree nodes store intervals.
struct StoreIntervals {};


/**
.Class.IntervalTreeNode:
..cat:Miscellaneous
..summary:Element of @Class.IntervalTree@.
..signature:IntervalTreeNode<TInterval, TSpec>
..param.TInterval:The type of interval to store.
..param.TSpec:The type of interval to store.
...default:StorePointsOnly.
...metafunction:Metafunction.Spec
..include: seqan/misc/misc_interval_tree.h

.Memvar.IntervalTreeNode#center:
..class:Class.IntervalTreeNode
..summary:The center of the interval of type TValue.

.Memvar.IntervalTreeNode#list1
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in ascending according to their left boundary points.

.Memvar.IntervalTreeNode#list2
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in descending according to their right boundary points.
 */
template<typename TInterval, typename TSpec=StorePointsOnly>
class IntervalTreeNode;


/**
.Spec.Interval Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:An Interval Tree Node that stores intervals explicitely in each node.
..signature:IntervalTreeNode<TInterval, StoreIntervals>
..param.TInterval:The interval type to store in the node.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StoreIntervals> {
public:
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<TInterval> list1;
	String<TInterval> list2;
};


/**
.Spec.Points Only Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:Spec for IntervalTreeNode that stores only the relevant point in each node meaning the endpoint of the interval in the list sorted by endpoints (list2) and only the beginpoint of the interval in the list sorted by beginpoints (list1).
..signature:IntervalTreeNode<TInterval, StorePointsOnly>
..param.TInterval:The interval type to store in the node.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StorePointsOnly> {
public:
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<PointAndCargo<TValue,TCargo> > list1;
	String<PointAndCargo<TValue,TCargo> > list2;

    /**
.Memfunc.IntervalTreeNode#IntervalTreeNode:
..class:Class.IntervalTreeNode
..summary:Default constructor.
..signature:IntervalTreeNode()
     */
    IntervalTreeNode()
    {
SEQAN_CHECKPOINT
    }

	IntervalTreeNode(IntervalTreeNode const & other):
		center(other.center),
		list1(other.list1),
		list2(other.list2)
	{
SEQAN_CHECKPOINT
	}
};





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
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue=int, typename TCargo=unsigned int>
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
	
	/**
.Memfunc.IntervalTree#IntervalTree
..class:Class.IntervalTree
..summary:Constructor
..signature:IntervalTree(intervalBegins, intervalEnds, intervalCargos, len)
..param.intervalBegins:Iterator pointing to begin position of first interval.
..param.intervalEnds:Iterator pointing to end position of first interval.
..param.intervalCargos:Iterator pointing to cargos/ids for intervals.
..param.len:Number of intervals to store in tree.
     */
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

	/**
..signature:IntervalTree(intervalBegins, intervalEnds, len)
     */
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
	
	/**
..signature:IntervalTree(String<TInterval> intervals)
     */
	IntervalTree(String<TInterval> intervals)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals);
	}

	/**
..signature:IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
     */
	template <typename TTagSpec>
	IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals,tag);
	}

    /**
..signature:IntervalTree(String<TInterval> intervals, TValue center)
    */
	IntervalTree(String<TInterval> intervals, TValue center)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals,center);
	}
};



///////Specs for the way interval centers are determined
/**
.Tag.IntervalTree Centers
..cat:Miscellaneous
..summary:Tag to select a specific way to compute the center of an interval tree node.
..see:Class.IntervalTree
..include:seqan/misc/misc_interval_tree.h
 */


/**
..Tag.ComputeCenter
...summary:For intervals that are more or less uniformly distributed in the value range, using the ComputeCenter tag may result in a more balanced tree compared to using the RandomCenter tag.
...signature:ComputeCenter
...remarks:center = minbegin + (maxend-minbegin)/2
 */
//template <typename TSpec = SpecPointAndCargo>
struct TagComputeCenter_;
typedef Tag<TagComputeCenter_> const ComputeCenter;


/**
..Tag.RandomCenter
...summary:The RandomCenter tag guarantees that each node contains at least one interval, therefore the size of the tree is limited by the nummer of intervals. This may lead to an unbalanced tree, but is the most space-efficient and in practice the fastest method.
...signature:RandomCenter
...remarks:center = center of random interval
 */
//template <typename TSpec = SpecPointAndCargo>
struct TagRandomCenter_;
typedef Tag<TagRandomCenter_> const RandomCenter;


///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalAndCargo functions //////////////////////////
///////////////////////////////////////////////////////////////////////////




/**
.Function.leftBoundary:
..cat:Miscellaneous
..summary:Access to the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the left boundary of the interval of type TValue&.
..see:Function.getLeftBoundary
..see:Function.rightBoundary
..see:Function.getRightBoundary
*/
template<typename TValue, typename TCargo>
TValue &
leftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}


/**
.Function.rightBoundary:
..cat:Miscellaneous
..summary:Access to the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the right boundary of the interval of type TValue&.
..see:Function.getRightBoundary
..see:Function.leftBoundary
..see:Function.getLeftBoundary
*/
template<typename TValue, typename TCargo>
TValue &
rightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}


/**
.Function.getLeftBoundary:
..cat:Miscellaneous
..summary:Get method for the left boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the left boundary of the interval of type TValue.
..see:Function.leftBoundary
..see:Function.getRightBoundary
..see:Function.rightBoundary
*/
template<typename TValue, typename TCargo>
TValue
getLeftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}


/**
.Function.getRightBoundary:
..cat:Miscellaneous
..summary:Get method for the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the right boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the right boundary of the interval of type TValue.
..see:Function.rightBoundary
..see:Function.getLeftBoundary
..see:Function.leftBoundary
*/
template<typename TValue, typename TCargo>
TValue
getRightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}


/**
.Function.cargo:
..signature:cargo(me)
..param.me:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/
template<typename TValue, typename TCargo>
TCargo &
cargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}

/**
.Function.getCargo:
..signature:getCargo(me)
..param.me:
...type:Class.IntervalAndCargo
..see:Function.cargo
*/
template<typename TValue, typename TCargo>
TCargo
getCargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}


/////////////////// Metafunctions //////////////////////
    
///.Metafunction.Value.param.T.type:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Value<IntervalAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Cargo<IntervalAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};


///////////////////////////////////////////////////////////////////////////
///////////////////// PointAndCargo functions /////////////////////////////
///////////////////////////////////////////////////////////////////////////

/**
.Function.leftBoundary:
..signature:leftBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue &
leftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.rightBoundary:
..signature:rightBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue &
rightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.getLeftBoundary:
..signature:getLeftBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue
getLeftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.getRightBoundary:
..signature:getRightBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue
getRightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.cargo:
..signature:cargo(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TCargo &
cargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}


/**
.Function.cargo:
..signature:getCargo(point)
..param.point:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/
template<typename TValue, typename TCargo>
TCargo
getCargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}



////////////////// Metafunctions //////////////////
///.Metafunction.Value.param.T.type:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Value<PointAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Cargo<PointAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};


//// Comparators
template <typename TPair>
bool _less_compI1_ITree(TPair const& p1, TPair const& p2){
SEQAN_CHECKPOINT
  return (leftBoundary(const_cast<TPair&>(p1)) < leftBoundary(const_cast<TPair&>(p2)));
}


template <typename TPair>
bool _greater_compI2_ITree(TPair const& p1, TPair const& p2){
SEQAN_CHECKPOINT
  return (rightBoundary(const_cast<TPair&>(p1)) > rightBoundary(const_cast<TPair&>(p2)));
}



///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalTreeNode functions //////////////////////////
///////////////////////////////////////////////////////////////////////////





//internal set functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval,StoreIntervals> & knot,TValue center,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	knot.center = center;
	appendValue(knot.list1,interval);
	appendValue(knot.list2,interval);

}

template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval,StoreIntervals> & knot,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	appendValue(knot.list1,interval);
	appendValue(knot.list2,interval);

}


//internal set functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval,StorePointsOnly> & knot,TValue center,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	knot.center = center;
	appendValue(knot.list1,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
	

}


template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval,StorePointsOnly> & knot,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	appendValue(knot.list1,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
	

}

/////////////////// Metafunctions ///////////////////////
///.Metafunction.Value.param.T.type:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Value<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Value<TInterval>::Type Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Cargo<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Cargo<TInterval>::Type Type;
};


/**
.Metafunction.ListType:
..cat:Miscellaneous
..signature:ListType<T>::ListType
..summary:Type of lists in tree nodes.
..param.T:The type to retrieve the list type for.
..returns:Returns the type of the the lists in @Class.IntervalTreeNode@ objects.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename T>
struct ListType;


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StorePointsOnly> >
{
	typedef String<PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StoreIntervals> >
{
	typedef String<IntervalAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};



///////////////////////////////////////////////////////////////////////////
/////////////////////// IntervalTree functions ////////////////////////////
///////////////////////////////////////////////////////////////////////////


/**
.Function.createIntervalTree
..summary:Create an interval tree.
..cat:Miscellaneous
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, Tag<TSpec> const tag)
..param.g:DirectedGraph to create interval tree in.
...type:Class.Graph
..param.pm:Property map to use for the created interval tree.
...type:Class.PropertyMap
..param.intervals:Container of intervals.
...type:Class.String
...remark:Should be a String of @Class.Interval@ or @Class.IntervalAndCargo@ objects.
..param.tag:Tag for tree construction method. @Tag.IntervalTree Centers.Tag.RandomCenter@ or @Tag.IntervalTree Centers.@Tag.ComputeCenter@
..remark:center of root node is computed by _calcIntervalTreeRootCenter
..include:seqan/misc/misc_interval_tree.h
 */
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

    if (length(intervals) > 0u) {
        TValue center =	_calcIntervalTreeRootCenter(intervals);

        std::sort(begin(intervals),end(intervals),_less_compI1_ITree<TInterval>);

        String<TInterval*> interval_pointers;
        resize(interval_pointers,length(intervals));
        _makePointerInterval(intervals,interval_pointers);

        _createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
        reserve(pm, length(pm), Exact());
        reserve(g.data_vertex, length(g.data_vertex), Exact());
    }
}


/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals)
..param.tag.default:Tag.IntervalTree Centers.tag.ComputeCenter
 */
// most user friendly interval tree construction for the moment...
// RandomCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, 
				   TPropertyMap & pm, 
				   TIntervals & intervals)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,RandomCenter());
}


/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, center, tag)
 */
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

/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, center)
 */
// RandomCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   TIntervals & intervals, 
				   typename Value<typename Value<TIntervals>::Type>::Type center)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,center,RandomCenter());
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
	//  Rekursionsanker
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*intervals[0]);
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TIntervalPointer*> TIntervalPointers;

	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;

	TValue min1 = maxValue<TValue>();
	TValue min2 = maxValue<TValue>();
	TValue max1 = minValue<TValue>();
	TValue max2 = minValue<TValue>();

	value(pm,knot).center = center;
	
 
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
			 //remember right most and left most point in left list
			if((**it).i2 > max1)
				max1 = (**it).i2;
			if((**it).i1 < min1)
				min1 = (**it).i1;
		}
		else
		{
			// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
				 //remember right most and left most point in right list
				if((**it).i2 > max2)
					max2 = (**it).i2;
				if ((**it).i1 < min2)
					min2 = (**it).i1;
			}
			else // interval belongs to this node
			{
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_left,vd,center,min1+(max1-min1)/2,length(S_left),tag);
	}
	// build subtree to the right
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
	// Rekursionsanker
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*value(intervals,0));
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TInterval*> TIntervalPointers;
	
	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;
		
	value(pm,knot).center = center;
	
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
		}
		else
		{	// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
			}
			else
			{
				// interval belongs to the current node
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterLeft(S_left,last_center,center,tag);
		_createIntervalTree(g,pm,S_left,vd,center,next_center,length(S_left),tag);
	}
	// build subtree to the right
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterRight(S_right,last_center,center,tag);
		_createIntervalTree(g,pm,S_right,vd,center,next_center,length(S_right),tag);
	}
}


// fill the container interval_pointers with pointers to the corresponding objects in intervals.
// this is done to avoid copying and passing the whole IntervalAndCargo objects during interval tree construction
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


// if the center of the root is not given, it is placed in the "ComputeCenter way": in the middle of minValue and maxValue
// where minValue is the minimum left boundary and maxValue is the maximum right boundary of all intervals
template<typename TIntervals>
typename Value<typename Value<TIntervals>::Type>::Type
_calcIntervalTreeRootCenter(TIntervals & intervals)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(intervals), 0u);
	
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	typedef typename Iterator<TIntervals,Standard>::Type TIntervalIterator;

	TIntervalIterator it = begin(intervals);
	TIntervalIterator it_end = end(intervals);

	TValue min = maxValue<TValue>();
	TValue max = minValue<TValue>();

	while(it != it_end)
	{
		if(leftBoundary(*it)<min) min = leftBoundary(*it);
		if(rightBoundary(*it)>max) max = rightBoundary(*it);
	  SEQAN_ASSERT_LEQ(min, max);
		++it;
	}

	SEQAN_ASSERT_LEQ(min, max);
	
	return (min+(max-min)/(TValue)2.0);

}



/**
.Function.addInterval
..cat:Miscellaneous
..signature:addInterval(graph, propertyMap, interval)
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree.
..param.interval:The interval to be added to the interval tree.
..summary:Adds an interval to an interval tree.
..include:seqan/misc/misc_interval_tree.h
*/
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
	

	if(empty(pm))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
		return;
		
	}
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


/**
..signature:addInterval(intervalTree, interval)
..param.intervalTree:The interval tree to add the interval to.
...type:Class.IntervalTree
 */
template<typename TValue, typename TCargo, typename TInterval>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TInterval interval)
{
SEQAN_CHECKPOINT

	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}


// TODO(holtgrewe): Is this begin/end in C++ style or is it first/last?
/**
..signature:addInterval(intervalTree, begin, end, cargo)
..param.begin:Begin position of interval of type TValue.
..param.end:End position of interval of type TValue.
..param.cargo:Cargo to attach to the interval.
 */
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


/**
..signature:addInterval(intervalTree, begin, end)
 */
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


/**
.Function.findIntervals
..summary:Find all intervals that contain the query point or overlap with the query interval.
..signature:findIntervals(graph, propertyMap, query, result)
..param.query:A query point.
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, TPropertyMap & pm, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Value<TProperty>::Type    TPropertyValue;
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
		if(act_prop.center < (TPropertyValue)query)
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > (TPropertyValue)query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
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
			if((TPropertyValue)query < act_prop.center)
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) <= (TPropertyValue)query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
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
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				break;
			}
		}
	}

}


/**
..signature:findIntervals(intervalTree, query, result)
..param.intervalTree:An interval tree
...type:Class.IntervalTree
*/
template<typename TValue,typename TCargo>
void
findIntervals(IntervalTree<TValue,TCargo> & it, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervals(it.g,it.pm,query,result);

}


/**
.Function.findIntervalsExcludeTouching
..cat:Miscellaneous
..summary::Find all intervals that contain the query point, exclude intervals that touch the query, i.e. where the query point equals the start or end point.
..signature:findIntervalsExcludeTouching(graph, propertyMap, query, result)
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree
..param.query:The TValue to query here.
..param.result:The resulting string of cargos/ids of the intervals that contain the query point.
...type:Class.String
...remark:Should be a string of TCargo.
..include:seqan/misc/misc_interval_tree.h
 */
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
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
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
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
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
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				break;
			}
		}
	}

}


/**
..signature:findIntervalsExcludeTouching(intervalTree, query, result)
..param.intervalTree:An interval tree
...type:Class.IntervalTree
*/
template<typename TValue,typename TCargo>
void
findIntervalsExcludeTouching(IntervalTree<TValue,TCargo> & tree, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervalsExcludeTouching(tree.g,tree.pm,query,result);

}




/**
.Function.findIntervals
..signature:findIntervals(intervalTree, query_begin, query_end, result)
..param.query_begin:The begin position of the query interval.
..param.query_end:The end position of the query interval.
*/
template<typename TValue,typename TCargo>
void
findIntervals(IntervalTree<TValue,TCargo> & tree, TValue query_begin, TValue query_end, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervals(tree.g,tree.pm,query_begin,query_end,result);

}



template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, TPropertyMap & pm, TValue query_begin, TValue query_end, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	findIntervals(g, pm, act_knot, query_begin, query_end, result);
}


template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, 
			  TPropertyMap & pm, 
			  typename VertexDescriptor<TGraph>::Type & act_knot, 
			  TValue query_begin, 
			  TValue query_end, 
			  String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		//
		if(act_prop.center < query_begin) // query interval is to the right of node center
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > query_begin)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
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
			if(query_end < act_prop.center) // query interval is to the left of node center
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) < query_end)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
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
			else{//node center is contained in query interval
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				
				while(!atEnd(it))
				{
					TVertexDescriptor next_knot = targetVertex(it);
					findIntervals(g,pm, next_knot, query_begin, query_end, result);
					goNext(it);
				}
				break;

				//break; //dont break! continue in both subtrees!!
			}
		}
	}

}



/////////////////// Metafunctions ///////////////////////

///.Metafunction.Value.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Value<IntervalTree<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Cargo<IntervalTree<TValue,TCargo> >
{
	typedef TCargo Type;
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  //#ifndef SEQAN_MISC_INTERVAL_TREE_H
