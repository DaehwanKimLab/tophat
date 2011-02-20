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
  $Id: rt_skip_base_element.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef _RT_MAX_SKIP_BASE_ELEMENT_H
#define _RT_MAX_SKIP_BASE_ELEMENT_H


namespace seqan
{


//___________________________ struct SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > > _______________________
// Adaption of the struct SkipBaseElement for use in a RangeTree
// Instead saving it's own theKey and value, it has a pointer to the related 
// RTEntry object, which stores the theKey and value
// To get the correct theKey, SkipBaseElement needs information about the dimension
// of the RangeTree it is part of.


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	dump( SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		typename Size< TObject >::Type dim = dimension( *getObject( &me ) );
		typename Size< TObject >::Type counter = 0;
		while( counter < dim )
		{
			if( key( &getObject( &me ), counter ) == infimumValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
			else
				std::cout << key( getObject( &me ), counter );
			++counter;
			std::cout << " ";
		}
		std::cout << std::endl;
	}


	
} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //_RT_SKIP_BASE_ELEMENT_H
