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
  $Id: graph_impl_fragment.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H
#define SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Fragment Specs
//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = Default>
struct ExactFragment;	


template<typename TSpec = Default>
struct ExactReversableFragment;	


//////////////////////////////////////////////////////////////////////////////
// Default Fragment is the exact one
//////////////////////////////////////////////////////////////////////////////

template<typename TSize = typename Size<String<char> >::Type, typename TSpec = ExactFragment<> >
class Fragment;


//////////////////////////////////////////////////////////////////////////////
// Size Metafunction
//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> > {
	typedef TSize Type;
};


template<typename TSize, typename TSpec>
struct Size<Fragment<TSize, TSpec> const> {
	typedef TSize Type;
};


//////////////////////////////////////////////////////////////////////////////
// Exact Fragment
//////////////////////////////////////////////////////////////////////////////
	

template<typename TSize, typename TSpec>
class Fragment<TSize, ExactFragment<TSpec> > {
 public:
	 typedef typename Id<Fragment>::Type TId;
	 TId seqId1;
	 TSize begin1;
	 TId seqId2;
	 TSize begin2;
	 TSize len;
  
	 Fragment()	 {}

	 Fragment(TId sqId1, TSize beg1, TId sqId2, TSize beg2, TSize l) :
	 seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l) 
	 {
		 SEQAN_CHECKPOINT
	 }

};


//////////////////////////////////////////////////////////////////////////////
// Exact Fragment that is a forward or reverse match
//////////////////////////////////////////////////////////////////////////////
	

template<typename TSize, typename TSpec>
class Fragment<TSize, ExactReversableFragment<TSpec> > {
 public:
	 typedef typename Id<Fragment>::Type TId;
	 TId seqId1;
	 TSize begin1;
	 TId seqId2;
	 TSize begin2;
	 TSize len;
	 bool reversed;
  
	 Fragment()	 {}

	 Fragment(TId sqId1, TSize beg1, TId sqId2, TSize beg2, TSize l) :
	 seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(false) 
	 {
		 SEQAN_CHECKPOINT
	 }

 	 Fragment(TId sqId1, TSize beg1, TId sqId2, TSize beg2, TSize l, bool rev) :
	 seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l), reversed(rev) 
	 {
		 SEQAN_CHECKPOINT
	 }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TStringSet, typename TVal>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Fragment<TSize, TSpec> const& f,
      TStringSet& str,
      TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == (f.seqId1)) ? infix(getValueById(str, (TId) seqId), f.begin1, f.begin1 + f.len) : infix(getValueById(str, (TId) seqId), f.begin2, f.begin2 + f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TVal>
inline typename Id<Fragment<TSize, TSpec> >::Type
sequenceId(Fragment<TSize, TSpec> const& f,
		   TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == 0) ? f.seqId1 : f.seqId2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentBegin(Fragment<TSize, TSpec> const& f,
			  TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	return ((TId) seqId == f.seqId1) ? const_cast<TSize&>(f.begin1) : const_cast<TSize&>(f.begin2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f,
			   TVal const)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec>
inline TSize&
fragmentLength(Fragment<TSize, TSpec> const& f)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactFragment<TSpec> > const& f,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	
	if ((TId) seqId == f.seqId1) {
		SEQAN_TASSERT((TPosition1)f.begin1<=pos)
		SEQAN_TASSERT(pos - f.begin1 < f.len)	
		pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_TASSERT((TPosition1)f.begin2<=pos)
		SEQAN_TASSERT(pos - f.begin2 < f.len)
		pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	typedef typename Id<Fragment<TSize, TSpec> >::Type TId;
	
	if ((TId) seqId == f.seqId1) {
		SEQAN_TASSERT((TPosition1)f.begin1<=pos)
		SEQAN_TASSERT(pos - f.begin1 < f.len)	
		if (f.reversed) pos2 = (f.begin2 + f.len - 1) - (pos - f.begin1);
		else pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_TASSERT((TPosition1)f.begin2<=pos)
		SEQAN_TASSERT(pos - f.begin2 < f.len)
		if (f.reversed) pos2 = (f.begin1 + f.len - 1) - (pos - f.begin2);
		else pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize, typename TSpec>
inline bool
isReversed(Fragment<TSize, ExactReversableFragment<TSpec> > const& f)
{
	SEQAN_CHECKPOINT
	return f.reversed;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
