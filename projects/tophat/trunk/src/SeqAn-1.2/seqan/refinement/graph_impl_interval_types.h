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
  $Id: graph_impl_interval_types.h 1757 2008-02-27 16:26:20Z emde@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_TYPES_H
#define SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree Types
//////////////////////////////////////////////////////////////////////////////






	

///---------------------------------------------------------------///

//////////////////// Interval and ID type ///////////////////
template<typename TValue = int, typename TCargo = int>
class IntervalAndCargo
{
public:
	TValue i1;
	TValue i2;
	TCargo cargo;

	IntervalAndCargo()
	{
SEQAN_CHECKPOINT
	}
	
	IntervalAndCargo(TValue i1, TValue i2, TCargo cargo):
		i1(i1),
		i2(i2),
		cargo(cargo)
	{
SEQAN_CHECKPOINT
	}

	IntervalAndCargo(IntervalAndCargo const & other):
		i1(other.i1),
		i2(other.i2),
		cargo(other.cargo)
	{
SEQAN_CHECKPOINT
	}
	
	IntervalAndCargo & operator = (IntervalAndCargo const & other)
	{
SEQAN_CHECKPOINT
		i1 = other.i1;
		i2 = other.i2;
		cargo = other.cargo;
		return *this;
	}



	~IntervalAndCargo()
	{
SEQAN_CHECKPOINT
	}
};

//get by reference functions
template<typename TValue, typename TCargo>
TValue &
leftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}

template<typename TValue, typename TCargo>
TValue &
rightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}

template<typename TValue, typename TCargo>
TValue
getLeftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}

template<typename TValue, typename TCargo>
TValue
getRightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}

template<typename TValue, typename TCargo>
TCargo &
cargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}

template<typename TValue, typename TCargo>
TCargo
getCargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}


/////////////////// Metafunctions //////////////////////
template<typename TValue,typename TCargo>
struct Value<IntervalAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


template<typename TValue,typename TCargo>
struct Cargo<IntervalAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};



/////////////////////// Point and ID type ////////////////
template<typename TValue = int, typename TCargo = int>
class PointAndCargo{
public:
	TValue point;
	TCargo cargo;

	PointAndCargo()
	{
SEQAN_CHECKPOINT
	}
	
	PointAndCargo(TValue point, TCargo cargo):
		point(point),
		cargo(cargo)
	{
SEQAN_CHECKPOINT
	}

	PointAndCargo(PointAndCargo const & other):
		point(other.point),
		cargo(other.cargo)
	{
SEQAN_CHECKPOINT
	}
	PointAndCargo & operator = (PointAndCargo const & other)
	{
SEQAN_CHECKPOINT
		point = other.point;
		cargo = other.cargo;
		return *this;
	}

	

	~PointAndCargo()
	{
SEQAN_CHECKPOINT
	}
};


//get by reference
template<typename TValue, typename TCargo>
TValue &
leftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}

template<typename TValue, typename TCargo>
TValue &
rightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}

template<typename TValue, typename TCargo>
TValue
getLeftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}

template<typename TValue, typename TCargo>
TValue
getRightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}

template<typename TValue, typename TCargo>
TCargo &
cargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}

template<typename TValue, typename TCargo>
TCargo
getCargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}



////////////////// Metafunctions //////////////////
template<typename TValue,typename TCargo>
struct Value<PointAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};

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
/////////////////////////// IntervalTreeNode	///////////////////////////

struct StorePointsOnly{};
struct StoreIntervals{};


template<typename TInterval, typename TSpec = StorePointsOnly>
class IntervalTreeNode;


//Spec for IntervalTreeNode that stores the intervals explicitly in each node
template<typename TInterval>
class IntervalTreeNode<TInterval,StoreIntervals>{

public:
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<TInterval> list1;
	String<TInterval> list2;

	IntervalTreeNode()
	{
SEQAN_CHECKPOINT
	}

//	IntervalTreeNode(TValue & center, String<TInterval> & list1, String<TInterval> & list2):
//		center(center),
//		list1(list1),
//		list2(list2)
//	{
//	}

	IntervalTreeNode(IntervalTreeNode const & other):
		center(other.center),
		list1(other.list1),
		list2(other.list2)
	{
SEQAN_CHECKPOINT
	}

	IntervalTreeNode & operator = (IntervalTreeNode const & other)
	{
SEQAN_CHECKPOINT
		center = other.center;
		list1 = other.list1;
		list2 = other.list2;
		return *this;
	}


	~IntervalTreeNode()
	{
SEQAN_CHECKPOINT
	}
	
};





//Spec for IntervalTreeNode that stores only the relevant point in each node
//meaning the endpoint of the interval in the list sorted by endpoints (list2) and 
//only the beginpoint of the interval in the list sorted by beginpoints (list1)
template<typename TInterval>
class IntervalTreeNode<TInterval,StorePointsOnly>{

public:
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<PointAndCargo<TValue,TCargo> > list1;
	String<PointAndCargo<TValue,TCargo> > list2;

	IntervalTreeNode()
	{
SEQAN_CHECKPOINT
	}

//	IntervalTreeNode(TValue center, 
//		String<PointAndCargo<TValue,TCargo> > & list1, 
//		String<PointAndCargo<TValue,TCargo> > & list2):
//			center(center),
//			list1(list1),
//			list2(list2)
//	{
//	}

	IntervalTreeNode(IntervalTreeNode const & other):
		center(other.center),
		list1(other.list1),
		list2(other.list2)
	{
SEQAN_CHECKPOINT
	}

	IntervalTreeNode & operator = (IntervalTreeNode const & other)
	{
SEQAN_CHECKPOINT
		center = other.center;
		list1 = other.list1;
		list2 = other.list2;
		return *this;
	}


	~IntervalTreeNode()
	{
SEQAN_CHECKPOINT
	}
	
};




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
template<typename TInterval, typename TSpec>
struct Value<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Value<TInterval>::Type Type;
};

template<typename TInterval, typename TSpec>
struct Cargo<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Cargo<TInterval>::Type Type;
};

template<typename T>
struct ListType;


template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StorePointsOnly> >
{
	typedef String<PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};

template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StoreIntervals> >
{
	typedef String<IntervalAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};








}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
