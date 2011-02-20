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
  $Id: rmt_def_algos.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_RMT_SL_DEF_ALGOS_H
#define SEQAN_HEADER_RMT_SL_DEF_ALGOS_H

namespace seqan{


/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

		// construct the towers for the complete skip list
		// -> the deferred RMT comletely builds the list of the lowest layer
	template< typename TObject, typename TSpec > inline
	void
	_completeBuild( SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > & list,
					typename Size< TObject >::Type dim )
	{
		if( dim == 0 )
		{
			_sortRecursive( list, _getBaseStore( list ), _getBaseStore( list ) + length( list ) + 1, dim );
			_setHeight( *_getBaseStore( list ), 1 );	
			_buildMaxTowers( list );
		}
	}
		
/////////////////////////////////////////////////////////////////////////////////////////
//
//	range maximum tree algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

		// _performRMQ adaption for the deferred RMT
		// The defferd RMT needs 2 searches per layer:
		// The first one assures that the element with the appropriate key is in an ordered state,
		// the second one is equivalent to _performRMQ of the common case
	template< typename TObject, typename TSpec, typename TBorder, typename TSize, typename TKey >
	void
	_performRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * list,
					TBorder & borderObj,
					TSize dim,
					TKey searchKey,
					TObject *& maxObject )
	{
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * layer_element = _getRoot( *list );
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > >::Type height = _getCurrentLayer( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > ** search_path = _getSearchPath( *_getMainTree( *list ), 0 ) + height - 1;
				// first search operation
		_searchFrom( *list, layer_element, searchKey, search_path, dim );
				// search again
		height = _getCurrentLayer( *list );
		layer_element = _getRoot( *list );
		while( height > 0 )
		{
			while( key( *_getRight( *layer_element ) ) < searchKey )
			{ 
				if( !_hasAssocStruct( layer_element ) )
					_createAssocStruct( layer_element, list, dim );
				if( _hasSmallAssocStruct( layer_element ) )
					_performSmallRMQ( layer_element, dim - 1, borderObj, maxObject );
				else
					_processRMQ( _getAssocStruct( layer_element ), dim - 1, borderObj, maxObject );
				layer_element = _getRight( *layer_element );
			}
			--layer_element;
			--height;
		}
		++layer_element;
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * base_element = _getDown( *layer_element );
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Deferred > * base_buffer = base_element;
		while( key( *base_buffer, dim ) < searchKey )
		{
			base_buffer = _getRight( *base_buffer );
		}
			// in the lowest layer, searching to the right
		while( base_element != base_buffer )
		{
			_testRangeMax( getObject( base_element), &borderObj, maxObject, dim - 1 );
			goNext( base_element );
		}
	}

	
}

#endif // SEQAN_HEADER_RMT_SL_DEF_ALGOS_H
