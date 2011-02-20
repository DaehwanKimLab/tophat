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
  $Id: rt_sl_compl_algos.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_RT_SL_COMPL_ALGOS_H
#define SEQAN_HEADER_RT_SL_COMPL_ALGOS_H

namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	
		// quick sort recursive step
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_sortRecursive( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
					SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * elem,
					SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		if( _getCount( *elem ) != 0 )
		{
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * pivot = _sort( list, elem, right, dim );
			_sortRecursive( list, elem, pivot, dim );
			_sortRecursive( list, pivot, right, dim );
		}	
	}
		// construct the towers for the complete skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize >
	void
	_buildTowers( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * buffer = _getBaseStore( list );
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * end = buffer + length( list ) + 1;
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > ** search_path = new SkipElement< TObject, TModus, RT< TSpec >, TStructuring >*[ _getMaximalSLTowerHeight( length( list ) ) ];
		*search_path = &_getUp( *buffer ); 
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring  > >::Type height;
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring  > >::Type max_height = _getMaximalSLTowerHeight( list );
		typename Key< TObject >::Type buffer_theKey = key( *_getBaseStore( list ), dim );
		goNext( buffer );

		while( buffer != end )
		{
			SEQAN_CHECK( buffer_theKey <= key( *buffer, dim ) )
			if( key( *buffer, dim ) != buffer_theKey )
			{
				height = _throwCoin< TObject, TModus, RT< TSpec >, TStructuring  >( list, max_height );
				if( height > 0 ){
					_add( list, buffer, height, search_path );
					_connect_actualize( list, buffer, height, search_path, dim );
				}
				buffer_theKey = key( *buffer, dim );
			}
			goNext( buffer  );
		}
		
		delete[] search_path;
	}

	

/////////////////////////////////////////////////////////////////////////////////////////
//
//	search algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

	
	
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet, typename TSize >
	void
	_fingerSearch(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
					TObject * left_border,
					TObject * right_border,
					TSize dim,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			
			// search for the left base element
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, dim );

		base = _checkBaseElementsLeft( base, left_border, right_border, dim, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;
	
				//	1 ) searching for highest layer,
				//		on-line search in associated structures
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );
				//	2 ) on-line search in associated structures of the higher layers
		base = _collectAssocStructs( tower_buffer, left_border, right_border, dim, right_theKey, list, results );
				//	3 ) check the remaining base elements
		_checkBaseElementsRight( base, left_border, right_border, dim, right_theKey, results );
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet >
	void
	_bottomSearch(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
					TObject * left_border,
					TObject * right_border,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type dim = 0;
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, dim );
		
		base = _checkBaseElementsLeftBottom( base, left_border, right_border, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;

		while( key( *base ) <= right_theKey )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
	}

}

#endif // SEQAN_HEADER_RT_SL_COMPL_ALGOS_H
