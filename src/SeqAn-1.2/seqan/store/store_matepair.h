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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_STORE_MATEPAIR_H
#define SEQAN_HEADER_STORE_MATEPAIR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Mate Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct MatePairStoreElement
{
	typedef typename Id<MatePairStoreElement>::Type TId;

	static const TId INVALID_ID;
	
	TId		readId[2];	// refers to the two reads of a mate-pair, INVALID_ID if this is a singleton fragment (e.g. in afg: reads refer to fragments (mate pairs) and these refer to libraries, singletons refer to an empty fragment)
	TId		libId;

	MatePairStoreElement() : libId(INVALID_ID) 
	{
		readId[0] = INVALID_ID;
		readId[1] = INVALID_ID;
	}
};


//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
const typename Id<MatePairStoreElement<TSpec> >::Type
MatePairStoreElement<TSpec>::INVALID_ID = SupremumValue<typename Id<MatePairStoreElement<TSpec> >::Type>::VALUE;


//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
