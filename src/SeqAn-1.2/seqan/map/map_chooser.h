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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MAP_CHOOSER_H
#define SEQAN_HEADER_MAP_CHOOSER_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Set meta-function to choose an efficient implementation

template <typename TValue, size_t SIZE>
struct _ChooseMap2
{
	typedef Map<TValue, Skiplist< > > Type;
};
template <typename TValue>
struct _ChooseMap2<TValue, 1>
{
	typedef Map<TValue, VectorSet< > > Type;
};


template <typename TValue>
struct ChooseMap:
	_ChooseMap2<TValue, sizeof(typename Key<TValue>::Type)>
{};


//////////////////////////////////////////////////////////////////////////////

}

#endif
