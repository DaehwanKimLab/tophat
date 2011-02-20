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
  $Id: adapt_std_list.h 3354 2009-02-04 16:11:00Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/


//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_ADAPT_STD_LIST_H
#define SEQAN_HEADER_ADAPT_STD_LSIT_H

#include <list>


namespace SEQAN_NAMESPACE_MAIN
{

template <typename TValue>
struct Value<std::list<TValue> >
{
	typedef TValue Type;
};

template <typename TValue>
struct Iterator<std::list<TValue> >
{
	typedef typename std::list<TValue>::iterator Type;
};

template <typename TValue>
struct Iterator<const std::list<TValue> >
{
	typedef typename std::list<TValue>::const_iterator Type;
};

template<typename TValue>
inline typename Iterator<std::list<TValue>, Standard >::Type
begin(std::list<TValue> &list_){
	return list_.begin();
}

template<typename TValue>
inline typename Iterator<std::list<TValue>, Standard >::Type
end(std::list<TValue> &list_){
	return list_.end();
}

template<typename TValue>
inline typename Iterator<const std::list<TValue>, Standard >::Type
begin(const std::list<TValue> &list_){
	return list_.begin();
}

template<typename TValue>
inline typename Iterator<const std::list<TValue>, Standard >::Type
end(const std::list<TValue> &list_){
	return list_.end();
}

template <typename T, typename T2>
inline void
appendValue(std::list<T> & list, 
			T2 value)
{
	list.push_back(value);
}

}

#endif //#ifndef SEQAN_HEADER_
