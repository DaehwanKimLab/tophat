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
  $Id: rt_sl_base.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_RT_SKIP_LIST_BASE_H
#define SEQAN_RT_SKIP_LIST_BASE_H

#include "skip_list.h"
#include "rt_skip_element.h"
#include "rt_skip_base_element.h"

namespace seqan{


/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic accessor functions
//
/////////////////////////////////////////////////////////////////////////////////////////

		// get thr right border element
		// e.g. the skip base element with supremum key
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > *
	_getRightBorder( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me );
		
	
		// specialization for static case
	template< typename TObject, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > * 
	_getRightBorder( SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &_getUp( me._baseStore[ me._numOfElements + 1 ] );
	}

		// specialization for dynamic case
	//template< typename TObject, typename TSpec, typename TStructuring > inline
	//SkipElement< TObject, SkipListDynamic, RT< TSpec >, TStructuring > * 
	//_getRightBorder( SkipList< TObject, SkipListDynamic, RT< TSpec >, TStructuring > & me )
	//{
	//	SEQAN_CHECKPOINT
	//	SkipElement< TObject, SkipListDynamic, RT< TSpec >, TStructuring > * border = _getRoot( me );
	//	while( key( *border ) < supremumValue< typename Key< TObject >::Type >() )
	//		border = _getRight( *border );
	//	return border;
	//}

		// get the element allocator of a skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited, SimpleAllocator > > &
	_getElementAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getElementAlloc( *_getMainTree( me ) );
	}

		// get the list allocator of a skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited, SimpleAllocator > > &
	_getListAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getListAlloc( *_getMainTree( me ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getBaseAlloc( *_getMainTree( me ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	RangeTree< TObject, TModus, RT< TSpec >, TStructuring > *
	_getMainTree( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me );

	
///////////////////////////////////////////////////////
//
//	initialization
//
///////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TContainer, typename TSize >
	void
	_create( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
				TContainer & data,
				RangeTree< TObject,	TModus, RT< TSpec >, TStructuring > & Tree,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		me._mainTree = &Tree;
		me._dim = dim;
		me._numOfElements = length( data );

		typename Iterator< TContainer >::Type first = begin( data );
		typename Iterator< TContainer >::Type last = end( data );

		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base_right;

		_initBases( me, me._baseStore, base_right, first, last, me._numOfElements, dim ); 

		_initSL( me, me._baseStore, base_right, me._numOfElements, me._dim );

		_completeBuild( me, dim );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize >
	void
	_create( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * first,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * last,
				RangeTree< TObject,	TModus, RT< TSpec >, TStructuring > & Tree,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		me._mainTree = &Tree;
		me._dim = dim;
		me._numOfElements = last - first;

		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * base_right;

		_initBases( me, me._baseStore, base_right, first, last, me._numOfElements, dim ); 

		_initSL( me, me._baseStore, base_right, me._numOfElements, dim );

		_completeBuild( me, dim );
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TIter, typename TSize > inline
	void
	_initBases( SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > & list, 
				SkipBaseElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > * firstBase,
				SkipBaseElement< TObject, SkipListStatic, RT< TSpec >, TStructuring > *& lastBase,
				TIter & firstData,
				TIter & lastData,
				TSize numEntries,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		allocate( _getBaseAlloc( list ), firstBase, numEntries + 2 );
		list._baseStore = firstBase;

		valueConstruct( firstBase, _getLBorderObj( *_getMainTree( list ) ), infimumValue< typename Key< TObject >::Type >() );
		
		++firstBase;

		while( firstData != lastData )
		{
			valueConstruct( firstBase, &value( firstData ), key( value( firstData ), dim ) );
			++firstBase;
			++firstData;
		}
		lastBase = firstBase;
		firstBase = list._baseStore;
		valueConstruct( lastBase, _getRBorderObj( *_getMainTree( list ) ), supremumValue< typename Key< TObject >::Type >() );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_initSL(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * first_base,
				SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * /*last_base*/,
				TSize numEntries,
				TSize /*dim*/ )
	{
		SEQAN_CHECKPOINT
			// allocate space for bording elements
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring  > * _rightBorder;
		allocate( _getElementAlloc( list ), _rightBorder, 2 );
		arrayConstruct( _rightBorder, _rightBorder + 2 );
		allocate( _getElementAlloc( list ), list._leftSideStore, _getMaximalSLTowerHeight( numEntries ) );
				
			// set values
			// left side ...
		_setHeight( *first_base, 1 );
		_setUp( *first_base, *list._leftSideStore );

		typename Key< TObject >::Type left_border_key = infimumValue< typename Key< TObject >::Type >();
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * buffer = list._leftSideStore;
		for( typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type i = 0; i < _getMaximalSLTowerHeight( numEntries ); ++i )
		{
			valueConstruct( buffer, _rightBorder, first_base, left_border_key );
			++buffer;
		}
			// ... right side
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * r_base = &list._baseStore[ numEntries + 1 ];
		_setDown( *_rightBorder, r_base);
		_setHeight( *r_base, 1 );
		_setRight( * _rightBorder, _rightBorder );
		setKey( *_rightBorder, supremumValue< typename Key< TObject >::Type >() );
		_setUp( *r_base, *_rightBorder );

		_setCount( *first_base, numEntries );
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	structure building algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect_actualize(	SkipList< TObject, TModus, RT< TSpec >, Deferred > & list,
						SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Deferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			_deleteAssocStruct( *search_path );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
	}

	// reconnect the pointers after building a tower
	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect(	SkipList< TObject, TModus, RT< TSpec >, Deferred > & list,
				SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > * base,
				TSize height,
				SkipElement< TObject, TModus, RT< TSpec >, Deferred > ** search_path,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, Deferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Deferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			_deleteAssocStruct( *search_path );
			++buffer;
			++search_path;
		}
	}

	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect_actualize(	SkipList< TObject, TModus, RT< TSpec >, SemiDeferred > & /*list*/,
						SkipBaseElement< TObject, TModus, RT< TSpec >, SemiDeferred > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred > * top = buffer + height;
		typename Key< TObject >::Type theKey = key( *base, dim );
		while( buffer != top )
		{
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, SemiDeferred >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			*search_path = buffer;
			++buffer;
			++search_path;
		}
	}


		// reconnect the pointers after building special for
		// the complete case
	template< typename TObject, typename TModus, typename TSpec, typename TSize > inline
	void 
	_connect_actualize(	SkipList< TObject, TModus, RT< TSpec >, Complete > & list,
						SkipBaseElement< TObject, TModus, RT< TSpec >, Complete > * base,
						TSize height,
						SkipElement< TObject, TModus, RT< TSpec >, Complete > ** search_path,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, TModus, RT< TSpec >, Complete > * buffer = &_getUp( *base );
		typename Key< TObject >::Type theKey = key( *base, dim );
		typename Size< SkipList< TObject, TModus, RT< TSpec >, Complete > >::Type current_height = 0;
		while( current_height < height - 1 ){
			new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Complete >( _getRight( **search_path ), base, theKey );
			_setRight( **search_path, buffer );
			if( _getHeight( **search_path ) == current_height + 1)
				_createAssocStruct( *search_path, &list, dim );
			*search_path = buffer;
			++current_height;
			++buffer;
			++search_path;
		}
		new( buffer ) SkipElement< TObject, TModus, RT< TSpec >, Complete >( _getRight( **search_path ), base, theKey );
		_setRight( **search_path, buffer );	
		_createAssocStruct( *search_path, &list, dim );
		*search_path = buffer;
	}


	// linker rand muss getestet werden -> darf nicht mit eigefügt werden
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_buildAssocStruct_left(  SkipList< TObject, TModus, RT< TSpec >, TStructuring >  * me,
							 SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * left,
							 SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
							 TSize dim )
	{
		SEQAN_CHECKPOINT
		if( _getDown( *left ) != _getBaseStore( *me ) )
			_createAssocStruct( left, me, _getDown( *left ), right, dim );
		else
			_createAssocStruct( left, me, _getSucc( *_getDown( *left ) ), _getDown( *_getRight( *_getRight( *left ) ) ), dim );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_buildAssocStruct_right(  SkipList< TObject, TModus, RT< TSpec >, TStructuring >  * me,
							  SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * left,
							  SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right,
							  TSize dim )
	{
		SEQAN_CHECKPOINT
		_createAssocStruct( left, me, _getDown( *left ), right, dim );
	}


}

#endif
