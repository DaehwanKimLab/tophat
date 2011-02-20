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
  $Id: rt_sl_def_algos.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_RT_SL_DEF_ALGOS_H
#define SEQAN_HEADER_RT_SL_DEF_ALGOS_H

namespace seqan{


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > ** 
	_getSearchPath( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
					TParam dim )
	{
		SEQAN_CHECKPOINT
		return _getSearchPath( *_getMainTree( me ), dim );
	}



/////////////////////////////////////////////////////////////////////////////////////////
//
//	search algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////


	template< typename TObject, typename TModus, typename TSpec, typename TResultSet, typename TSize >
	void
	_fingerSearch(	SkipList< TObject, TModus, RT< TSpec >, Deferred > * list,
					TObject * left_border,
					TObject * right_border,
					TSize dim,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path = _getSearchPath( *list, dim );
		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, search_path, dim );
		
		base = _checkBaseElementsLeft( base, left_border, right_border, dim, left_theKey, right_theKey, results );
		if( key( *base, dim ) > right_theKey )
			return;

			// search for the right border
			//	the cooresponding element might be sorted or unsorted
			//	=> search for the corresponding element

				// I. search for the smallest element, which is outside the search intervall
				//	1 ) search for the element, of which the right pointer of the highest tower
				//		points to an element outside of the search intervall
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );					
		
		base = _getDown( *tower_buffer );

				//  2 ) search downwards from highest element, sort all equal keys behind found object
		_sort_equals( *list, _searchFrom( *list, tower_buffer, right_theKey, search_path, dim ), right_theKey );
		
				// II. the range is now limited by a sorted element
				//	1 ) searching for highest layer again( may have changed ),
				//		on-line search in associated structures
		tower_buffer = _findTowerTop( base, list, left_border, right_border, dim, right_theKey, results );

				//	2 ) on-line search in associated structures of the higher layers
		base = _collectAssocStructs( tower_buffer, left_border, right_border, dim, right_theKey, list, results );

				//	3 ) check the remaining base elements
		_checkBaseElementsRight( base, left_border, right_border, dim, right_theKey, results );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TResultSet >
	void
	_bottomSearch(	SkipList< TObject, TModus, RT< TSpec >, Deferred > * list,
					TObject * left_border,
					TObject * right_border,
					TResultSet & results )
	{
		SEQAN_CHECKPOINT
		typename Size< TObject >::Type dim = 0;
		typename Key< TObject >::Type left_theKey = key( *left_border, dim );
		typename Key< TObject >::Type right_theKey = key( *right_border, dim );
			// search for the left base element
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path = _getSearchPath( *list, 0 );
		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base = _searchFrom( *list, _getRoot( *list ), left_theKey, search_path, dim );

		while( key( *base ) < left_theKey )
			goNext( base );
		if( key( *base ) > right_theKey )
			return;

		SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * right_base  = _searchFrom( *list, _getRoot( *list ), right_theKey, search_path, dim );
		_sort_equals( *list, right_base, right_theKey );

			//  2 ) search downwards from highest element
		while( base != right_base )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
		
		while( key( *base ) <= right_theKey )
		{
			_pushBack( results, getObject( base ) );
			goNext( base );
		}
	}

}

#endif // SEQAN_HEADER_RT_SL_DEF_ALGOS_H
