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
  $Id: skip_list_dynamic.h 3420 2009-02-12 12:10:09Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

/*	2006 Hendrik Woehrle
*
*	SkipListDynamic Specialization of the Skip List
*
*	Contains the functions for the SkipListDynamic specialization of the Skip List
*
*/

//SEQAN_NO_DDDOC: do not generate documentation for this file


#ifndef SEQAN_HEADER_SKIP_LIST_DYNAMIC_H
#define SEQAN_HEADER_SKIP_LIST_DYNAMIC_H


namespace seqan
{

	
/////////////////////////////////////////////////////////////////////////////////
//
//								 struct SkipList
//
/////////////////////////////////////////////////////////////////////////////////


///////////////////////////// SkipListDynamic< True, TSpec >  specialization ////////////////////////


		// Function to renew the left bording tower.
		// If the number of elements in the list increases, the _getMaximalSLTowerHeight 
		// value might increase as well.
		// Therefore, the height of the left bording tower must be increased.
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > 
	void
	_renewLeftTower( SkipList< TObject, TModus, TSpec, TStructuring > & list )
	{
		SEQAN_CHECKPOINT
			// allocate the new tower
		SkipElement< TObject, TModus, TSpec, TStructuring > * new_left;
		typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type height = _getMaximalSLTowerHeight( length( list ) );
		allocate( _getElementAlloc( list ), new_left, height + 1 );
		
			// construct tower elements, adjust connections and delete old tower elements 
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base = list._baseStore;
		SkipElement< TObject, TModus, TSpec, TStructuring > * buffer = &_getUp( *base );
		SkipElement< TObject, TModus, TSpec, TStructuring > * tower_top = buffer + height;
		typename Key< TObject >::Type theKey= key( *base );
		while( buffer != tower_top )
		{
			new( new_left ) SkipElement< TObject, TModus, TSpec, TStructuring >( _getRight( *buffer ), base, theKey );
			valueDestruct( buffer );
			++buffer;
			++new_left;
		}
		new( new_left ) SkipElement< TObject, TModus, TSpec, TStructuring >( &_getUp( *list._rightBorder ), base, theKey );
		deallocate( _getElementAlloc( list ), &_getUp( *base ), height );
			
			//	adjust values of the skip list and base layer
		new_left -= height;
		list._leftSideStore = new_left;
		_setUp( *base, *( new_left ) );
		_renewSearchPath( list, height, height+1 );
	}

		// calculate the count values of the deferred skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setDeferredCounts( SkipList< TObject, TModus, TSpec, TStructuring > & list,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * preceding_elem,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * new_elem )
	{
		// do nothing (non-deferred sl's don't have count values)
	}

	template< typename TObject, typename TModus, typename TSpec > inline
	void
	_setDeferredCounts( SkipList< TObject, TModus, TSpec, Deferred > & list,
						SkipBaseElement< TObject, TModus, TSpec, Deferred > * preceding_elem,
						SkipBaseElement< TObject, TModus, TSpec, Deferred > * new_elem )
	{
		SEQAN_CHECKPOINT
			_setCount( *new_elem, _getCount( *preceding_elem ) );
			_setCount( *preceding_elem, 0 );
	}

	
/**
.Function.insert:
..summary:Inserts an object into a SkipList.
..cat:SkipList
..signature:insert(list, object)
..param.list:The list.
...type:Class.SkipList
..param.object:The object.
...remarks:The TObject parameter of the SkipList.
..remarks:If the list is in a initial state, i.e. no search operation has been performed since the time of construction, the objects in the list are in an unsorted state.
..The object is therefore inserted at the end of the list. After the first search operation, the objects are sorted. Pbjects from following insert operations will then be inserted at the correct place. 
*/
		// two cases:
			// 1. list is empty apart from bording element or is not sorted 
			//	-> insert new element at the end 
			//	-> _insertBack is used
			// 2. list contains elements 
			//	-> insert at correct place
			//	-> _insertInPlace is used

		// static spec -> do nothing
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	insert( SkipList< TObject, TModus, TSpec, TStructuring > & list, 
			TObject & obj )
	{
		SEQAN_ASSERT2( false, "SkipListStatic Skip Lists don't provide insertion operations" )
	}


		// dynamic case (wrapper for additional parameters)
	template< typename TObject, typename TSpec, typename TStructuring > inline
	void 
	insert( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
			TObject & obj )
	{
		SEQAN_CHECKPOINT
		insert( list, obj, list );
	}


		// insertion at the end of the skip list
	template< typename TObject, typename TSpec, typename TStructuring, typename TParam >
	void 
	_insertBack( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
				TObject & obj,
				TParam & param )
	{
		SEQAN_CHECKPOINT
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * last_elem = _getPred( *list._rightBorder );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * new_elem;
		allocate( _getBaseAlloc( list ), new_elem, 1 );
		valueConstruct( new_elem, &obj, key( obj ) );
			// inserting values of member variables
		_setPred( *new_elem, last_elem );
		_setSucc( *last_elem, new_elem );
		_setPred( *list._rightBorder, new_elem );
		_setSucc( *new_elem, list._rightBorder );
		_setCount( *_getBaseStore( list ), _getCount( *_getBaseStore( list ) ) +1 );
		typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type list_size = length( list );
		if( _getMaximalSLTowerHeight( list_size ) < _getMaximalSLTowerHeight( list_size + 1 ) )
			_renewLeftTower( list );
		_setLength( list, list_size + 1 );
	}

		// insertion at the correct place of the skip list
	template< typename TObject, typename TSpec, typename TStructuring, typename TParam >
	void 
	_insertInPlace( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
				TObject & obj,
				TParam & param )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** search_path = _getSearchPath( list );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * preceding_elem = _searchFrom( list, _getRoot( list ), key( obj ), search_path, list );
		
		_sort_equals( list, preceding_elem );
		while( key( *_getSucc( *preceding_elem ), param ) < key( obj ) )
			goNext( preceding_elem );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * new_elem;
		allocate( _getBaseAlloc( list ), new_elem, 1 );
		valueConstruct( new_elem, &obj, key( obj ) );

			// inserting values of member variables
		_setPred( *new_elem, preceding_elem );
		_setSucc( *new_elem, _getSucc( *preceding_elem ) );
		_setSucc( *preceding_elem, new_elem );
		_setPred( *_getSucc( *new_elem ), new_elem );
		
		_renewDynConnects( *preceding_elem, *new_elem, _getCount( *preceding_elem ) + 1, static_cast< typename Size< SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > >::Type >( 0 ) );

		typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type height = _throwCoin( list, _getMaximalSLTowerHeight( list ) );
		if( key( *preceding_elem, param ) != key( *new_elem, param ) && height > 0 )
		{
			_add( list, new_elem, height, search_path );
			_connect( list, new_elem, height, search_path, param );
		}
		if( _getMaximalSLTowerHeight( length( list ) ) < _getMaximalSLTowerHeight( length( list ) + 1 ) )
			_renewLeftTower( list );
		_setLength( list, length( list ) + 1 );
	}

		// the insertion function itself
	template< typename TObject, typename TSpec, typename TStructuring, typename TParam > inline
	void 
	insert( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
			TObject & obj,
			TParam & param )
	{
		SEQAN_CHECKPOINT
			// case 1
		if( _getInitialState( list ) )
		{		
			_insertBack( list, obj, param );
		}
		else 
		{	// case 2
			_insertInPlace( list, obj, param );
		}
	}

	

/**
.Function.erase:
..summary:Deletes object from a SkipList.
..cat:SkipList
..signature:erase(list, iterator)
..signature:erase(list, object)
..signature:erase(list, theKey)
..param.list:The list containing the element.
...type:Class.SkipList
..param.iterator:Iterator to an entry, that should be deleted.
...concept:Concept.Iterator
..param.object:An object in the list, that should be deleted.
..param.key:A key.
...type:A key. All entries with the given theKeywill be deleted.
..remarks.text:Returns true, if at least one element is deleted, false otherwise.
*/

	template< typename TObject, typename TSpec, typename TStructuring, typename TParam >
	SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *
	_deleteSearchFrom(	SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
						typename Key< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type theKey,
						SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** search_path,
						TParam & param )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( theKey != supremumValue< typename Key< TObject >::Type >( ), "search key is supremum" ) 
		SEQAN_CHECK2( theKey != infimumValue< typename Key< TObject >::Type >( ), "search key is infimum" ) 
		
		if( _getInitialState( list ) )
		{
			_completeBuild( list );
			_noInitialState( list );
		}
		
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > * layer_element = _getRoot( list );
		typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type height = layer_element - &_getUp( *_getDown( *layer_element ) ) + 1;
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > * temp_right = _getRight( *layer_element );
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** sp = &search_path[ height - 1];
			// search in higher layers		
		while( height > 0 ){
			while( key( *temp_right ) < theKey ){ 
				layer_element = temp_right;
				right( temp_right );
			}
			--height;
			*sp = layer_element;
			--sp;
			--layer_element;
			temp_right = _getRight( *layer_element );
		}
		++layer_element;
			// in the lowest layer, searching to the right
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * base_element = _getDown( *layer_element );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * right_base = _getRight( *base_element );
		while( key( *right_base, param ) < theKey ){
			base_element = right_base;
			right( right_base );
		}
		if( key( *right_base, param ) == theKey && key( *base_element, param ) != theKey )
		{
			base_element = right_base;
		}

			// if the element found is in correct place and no unsorted elements follow,
			// it is the correct element. Otherwise search following intervall
		if( _getCount( *base_element ) > 0 )
		{
			base_element = _splitAction( list, base_element, theKey, search_path, param );
		}
			// handling for multiple elements
		right_base = _getRight( *base_element );
		if( key( *right_base, param ) <= theKey && key( *base_element, param ) < theKey )
			return right_base;
		return  base_element;
	}

		// deletes the tower of skip elements that is related to a skip base element and 
		// adjusts the connections of the elements in the search path
	template< typename TObject, typename TSpec, typename TStructuring >
	void
	_deleteTower( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
					SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** search_path, 
					SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * base, 
					typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type height )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK( &_getUp( *base ) != NULL )
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > * buffer = &_getUp( *base );
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > * top = buffer + height;
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** sp = search_path;
			// adjusting pointers of the elements in the search path
		while( buffer != top )
		{
			_setRight( **sp, _getRight( *buffer ) );
			valueDestruct( buffer );
			++sp;
			++buffer;
		}
			// free memory
		deallocate( _getElementAlloc( list ), &_getUp( *base ), height );
	}

		// deletes the base element in a skip list and adjusts the
		// connections between surrounding elements
	template< typename TObject, typename TSpec, typename TStructuring >
	void
	_deleteBase( SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list,
					SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > *& base )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( base != list._rightBorder, "Can't delete right border" )
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * base_buffer1 = _getPred( *base );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * base_buffer2 = _getSucc( *base );
		_setSucc( *base_buffer1, base_buffer2 );
		_setPred( *base_buffer2, base_buffer1 );
		_setDefConnects( _getLeft( *base ), _getRight( *base ) );
		valueDestruct( base );
		deallocate( _getBaseAlloc( list ), base, 1 );
		_setLength( list, length( list ) - 1 );
		base = base_buffer1;
	}

	

////////////////////////////////////////////////////////////////////////////////////////////

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey >
	bool 
	erase(	SkipList< TObject, TModus, TSpec, TStructuring > & list, 
			TKey theKey)
	{
		SEQAN_ASSERT2( false, "No deletion in static lists")
		return false;
	}

	template< typename TObject, typename TSpec, typename TStructuring, typename TKey >
	bool 
	erase(	SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list, 
			TKey theKey )
	{
		SEQAN_CHECKPOINT
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** search_path = _getSearchPath( list );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * buffer = _deleteSearchFrom( list, theKey, search_path, list );
		if( key( *buffer ) != theKey )
			return false;
		if( _getHeight( *buffer ) )
			_deleteTower( list, search_path, buffer, _getHeight( *buffer ) );
		_sort_equals( list, buffer, theKey );
		while( key( *buffer ) == theKey )
		{
			_deleteBase( list, buffer );
			goNext( buffer );
		}
		return true;
	}
	
	// deleting an element with a given object
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	bool 
	erase(	SkipList< TObject, TModus, TSpec, TStructuring > & list, 
			TObject & obj )
	{
		SEQAN_ASSERT2( false, "No deletion in static lists")
		return false;
	}


		// deleting an element with a given object
	template< typename TObject, typename TSpec, typename TStructuring >
	bool 
	erase(	SkipList< TObject, SkipListDynamic, TSpec, TStructuring > & list, 
			TObject & obj )
	{
		SEQAN_CHECKPOINT
		if( key( obj ) == infimumValue< typename Key< TObject >::Type >() || key( obj ) == supremumValue< typename Key< TObject >::Type >() )
			return false;
		SkipElement< TObject, SkipListDynamic, TSpec, TStructuring > ** search_path = _getSearchPath( list );
		SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * buffer = _deleteSearchFrom( list, key( obj ), search_path, list );
		if( key( *buffer, list ) != key( obj ) )
			return false;
		_sort_equals( list, buffer );

		while( getObject( buffer ) != &obj && key( *buffer ) == key( obj ) )
			goNext( buffer );

		if( &obj == getObject( buffer ) && buffer != list._baseStore )
		{
			typename Size< SkipList< TObject, SkipListDynamic, TSpec, TStructuring > >::Type height = _getHeight( *buffer );
			if( height > 0 )
			{ 
				if( key( *_getSucc( *buffer ) ) == key( *buffer ) ){
					SkipBaseElement< TObject, SkipListDynamic, TSpec, TStructuring > * succ_base = _getSucc( *buffer );
					_swapBases( *buffer, *succ_base );
					_deleteBase( list, succ_base );
				}
				else{
					_deleteTower( list, search_path, buffer, height );
					_deleteBase( list, buffer );
				}
			}
			else{
				_deleteBase( list, buffer );
			}
			return true;
		}
		
			// searching rightmost element that may have the given object
		while( key( *buffer, list ) == key( obj ) )
		{
			goNext( buffer );
		}
		if( getObject( buffer ) == &obj ){
			_deleteBase( list, buffer );
			return true;
		}
		return false;
	}

	
} // namespace seqan


#endif // _SKIP_LIST_DYNAMIC_H
