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
  $Id: rmt_compl_algos.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_RMT_SL_COMPL_ALGOS_H
#define SEQAN_HEADER_RMT_SL_COMPL_ALGOS_H

namespace seqan{

/*
 *  rmt_compl_algos.h
 *  rmt
 *
 *  Created by Hendrik Woehrle
 *
 *	Specializations for the complete RMT
 *
 */


/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

		// construct the towers for the complete and semi deferred skip list layer in a rmt

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_completeBuild( SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & list,
					typename Size< TObject >::Type dim )
	{
		_sortRecursive( list, _getBaseStore( list ), _getBaseStore( list ) + length( list ) + 1, dim );
		_setHeight( *_getBaseStore( list ), 1 );
		if( dim != 0 )
			_buildTowers( list, dim );
		else _buildMaxTowers( list );
	}
	
		// activation of the max pointers during initialization
	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_activateScoreBuild( SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * base,
						SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
						SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** search_path )
	{
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = 0;
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type max_height = _getMaximalSLTowerHeight( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > ** layer_element = search_path;
		TObject * max_obj = getObject( base );
		typename Weight< TObject >::Type max_score = priority( *max_obj );
		while( height < max_height )
		{
			if( priority( *layer_element ) <= max_score )
				_setMaxObject( *layer_element, max_obj );
			else
				break;
			++height;
			++layer_element;
		}
	}

	

		// reconnect the pointers after building
		// special for the complete case
	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect_actualize(	SkipList< TObject, TModus, RT< MaxTree< TSpec > >, Complete > & list,
						SkipBaseElement< TObject, TModus, RT< MaxTree< TSpec > >, Complete > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< MaxTree< TSpec > >, Complete > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< MaxTree< TSpec > >, Complete > * buffer = &_getUp( *base );
		typename Key< TObject >::Type searchKey = key( *base, dim );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > * tower_top = buffer + height - 1;
		while( buffer != tower_top ){
			new( buffer ) SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete >( _getRight( **search_path ), base, searchKey );
			_setRight( **search_path, buffer );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
		new( buffer ) SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete >( _getRight( **search_path ), base, searchKey );
		_setRight( **search_path, buffer );	
		_createAssocStruct( *search_path, &list, dim );
		*search_path = buffer;
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	range maximum tree algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

			// perform a rmq in higher layers, spec for semi deferred rmt
	template< typename TObject, typename TSpec, typename TStructuring, typename TBorder, typename TSize, typename TKey >
	void
	_performRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
					TBorder & borderObj,
					TSize dim,
					TKey searchKey,
					TObject *& maxObject )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( searchKey < supremumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds supremum" ) 
		SEQAN_CHECK2( searchKey > infimumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds infimum" ) 
	
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * layer_element = _getRoot( *list );
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = _getCurrentLayer( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * right_buffer = _getRight( *layer_element );
		// search in higher layers		
		while( height > 0 )
		{
			while( key( *right_buffer ) < searchKey )
			{ 
				if( !_hasAssocStruct( layer_element ) )
					_createAssocStruct( layer_element, list, dim );
				if( _hasSmallAssocStruct( layer_element ) )
					_performSmallRMQ( layer_element, dim - 1, borderObj, maxObject );
				else
					_processRMQ( _getAssocStruct( layer_element ), dim - 1, borderObj, maxObject );
				layer_element = right_buffer;
				right_buffer = _getRight( *right_buffer );
			}
			--layer_element;
			right_buffer = _getRight( *layer_element );
			--height;
		}
		++layer_element;
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * base_element = _getDown( *layer_element );
//		typename Weight< TObject >::Type max_score = priority( *maxObject );
			// in the lowest layer, searching to the right
		while( key( *base_element, dim ) < searchKey )
		{
			_testRangeMax( getObject( base_element ), &borderObj, maxObject, dim - 1 );
			goNext( base_element );
		}
	}

		// perform a rmq in higher layers, spec for complete rmt
	template< typename TObject, typename TSpec, typename TBorder, typename TSize, typename TKey >
	void
	_performRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > * list,
					TBorder & borderObj,
					TSize dim,
					TKey searchKey,
					TObject *& maxObject )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( searchKey < supremumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds supremum" ) 
		SEQAN_CHECK2( searchKey > infimumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds infimum" ) 
	
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > * layer_element = _getRoot( *list );
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > >::Type height = _getCurrentLayer( *list );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > * right_buffer = _getRight( *layer_element );
		
			// search in higher layers		
		while( height > 0 )
		{
			while( key( *right_buffer ) < searchKey )
			{ 
				if( _hasSmallAssocStruct( layer_element ) )
					_performSmallRMQ( layer_element, dim - 1, borderObj, maxObject );
				else
					_processRMQ( _getAssocStruct( layer_element ), dim - 1, borderObj, maxObject );
				layer_element = right_buffer;
				right_buffer = _getRight( *right_buffer );
			}
			--layer_element;
			right_buffer = _getRight( *layer_element );
			--height;
		}
		++layer_element;
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, Complete > * base_element = _getDown( *layer_element );
		typename Weight< TObject >::Type max_score = priority( *maxObject );
			// in the lowest layer, searching to the right
		while( key( *base_element, dim ) < searchKey )
		{
			_testRangeMax( getObject( base_element ), &borderObj, maxObject, dim - 1 );
			goNext( base_element );
		}
	}

	template< typename TObject, typename TSpec, typename TStructuring, typename TBorder >
	void
	_performRMQ(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
					TBorder & /*obj*/,
					typename Key< TObject >::Type searchKey,
					TObject *& maxObject )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( searchKey < supremumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds supremum" ) 
		SEQAN_CHECK2( searchKey > infimumValue< typename Key< TObject >::Type >( ), "search theKeyexceeds infimum" ) 
	
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * layer_element = _getRoot( *list );
		typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type height = _getCurrentLayer( *list );
		TObject * temp_obj;
		typename Weight< TObject >::Type max_score = priority( *maxObject );
		SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * right_buffer = _getRight( *layer_element );
			// search in higher layers		
		while( height > 0 )
		{
			while( key( *right_buffer ) < searchKey )
			{	
				temp_obj = _getMaxObject( layer_element );
				if( priority( *temp_obj ) >= max_score ){
					maxObject = temp_obj;
					max_score = priority( *maxObject );
				}
				layer_element = right_buffer;
				right_buffer = _getRight( *right_buffer );
			}
			--layer_element;
			--height;
			right_buffer = _getRight( *layer_element );
		}
		++layer_element;
			// lowest layer
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * base_element = _getDown( *layer_element );
		while( key( *base_element ) < searchKey )
		{
			if( priority( *maxObject ) < priority( base_element ) )
				maxObject = getObject( base_element );
			goNext( base_element );
		}
	}

		// perform a rmq for elements in a range which is below threshold
	template< typename TObject, typename TSpec, typename TStructuring, typename TBorder > inline
	void
	_performSmallRMQ(	SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * elem,
						typename Size< TObject >::Type dim,
						TBorder & border,
						TObject *& max_object )
	{
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * end = _getDown( *_getRight( *elem ) );
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * buffer = _getDown( *elem );
		while( buffer != end )
		{
			_testRangeMax( getObject( buffer ), &border, max_object, dim );
			goNext( buffer );
		}
	}

	

}

#endif // SEQAN_HEADER_RMT_SL_COMPL_ALGOS_H

