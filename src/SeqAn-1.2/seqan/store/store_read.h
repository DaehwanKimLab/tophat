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

#ifndef SEQAN_HEADER_STORE_READ_H
#define SEQAN_HEADER_STORE_READ_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Read Store
//////////////////////////////////////////////////////////////////////////////

template <typename TReadSeq, typename TPos, typename TSpec = void>
struct ReadStoreElement
{
	typedef typename Id<ReadStoreElement>::Type TId;
	
	static const TId INVALID_ID;
	
	
	TId matePairId;				// refers to the mate-pair, INVALID_ID if not part of a mate-pair

	ReadStoreElement() : matePairId(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TReadSeq, typename TPos, typename TSpec>
const typename Id<ReadStoreElement<TReadSeq, TPos, TSpec> >::Type
ReadStoreElement<TReadSeq, TPos, TSpec>::INVALID_ID = SupremumValue<typename Id<ReadStoreElement<TReadSeq, TPos, TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
