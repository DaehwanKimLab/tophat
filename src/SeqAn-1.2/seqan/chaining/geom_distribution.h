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
  $Id: geom_distribution.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

/*	
*	Random bit generator
*	
*	from Numerical Recipes in C
*	
*/


#ifndef SEQAN_HEADER_RAND_GEOM
#define SEQAN_HEADER_RAND_GEOM


#define SEQAN_RG_IB1 1
#define SEQAN_RG_IB2 2
#define SEQAN_RG_IB5 16
#define SEQAN_RG_IB18 131072
#define SEQAN_RG_MASK ( SEQAN_RG_IB1 + SEQAN_RG_IB2 + SEQAN_RG_IB5 )

namespace seqan {

	template< typename T > inline
	T
	_geomRand( )
	{
		static unsigned long seed = rand();
		T value = 0;
		while( true )
		{
			if( ( seed & SEQAN_RG_IB18 ) ){
				seed = ( ( seed ^ SEQAN_RG_MASK ) << 1 ) | SEQAN_RG_IB1;
				++value;
			}
			else {
				seed <<= 1;
				break;
			}
		}
		return value;
	}

}

#endif
