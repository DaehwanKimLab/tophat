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
  $Id: rmt_skip_base_element.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_RMT_SKIP_BASE_ELEMENT_H
#define SEQAN_RMT_SKIP_BASE_ELEMENT_H


namespace seqan
{

		// get the score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	weight( SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return weight( *getObject( me ) );
	}

		// get the chain score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	priority( SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return priority( *getObject( me ) );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //SEQAN_RMT_SKIP_BASE_ELEMENT_H


