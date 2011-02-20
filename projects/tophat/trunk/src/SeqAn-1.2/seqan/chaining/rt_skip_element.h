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
  $Id: rt_skip_element.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef _RT_SKIP_ELEMENT_H
#define _RT_SKIP_ELEMENT_H

namespace seqan
{

//_________________________________struct SkipElement< TObject, TModus, RT< TSpec >, Deferred > >_______________________
// Adaption of the struct SkipBaseElement for use in a Deferred, RT< TSpec >, TStructuring

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct _RangeCargo
	{

		SkipList< TObject, TModus, TSpec, TStructuring > * _assocStruct;

		_RangeCargo()
			: _assocStruct( NULL )
		{}

	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, TSpec, TStructuring > *
	_getAssoc( _RangeCargo< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._assocStruct;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setAssoc( _RangeCargo< TObject, TModus, TSpec, TStructuring > & me,
				SkipList< TObject, TModus, TSpec, TStructuring > * list)
	{
		SEQAN_CHECKPOINT
		me._assocStruct = list;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > 
	struct Cargo< SkipElement< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef _RangeCargo< TObject, TModus, RT< TSpec >, TStructuring > Type;
	};
	
		//  Documentation in SkipElement.h
		//  Adaption for Deferred, RT< TSpec >, TStructuring

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return _getAssoc( *cargo( *me ) );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me,
						SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list )
	{
		SEQAN_CHECKPOINT
		_setAssoc( *cargo(* me ), list );
	}
	

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_buildAssocStruct(	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me,
						SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( _getAssocStruct( me ) == NULL, "List overwritten" )
		if( _getDown( *me ) != _getBaseStore( *list ) )
		{
			SkipList< TObject, TModus, RT< TSpec >, TStructuring > * new_list;
			allocate( _getListAlloc( *list ), new_list, 1 );
			valueConstruct( new_list );
			_create( *new_list, _getDown( *me ), _getDown( *_getRight( *me ) ), *_getMainTree( *list ), dim - 1 );
			_setAssocStruct( me, new_list );
		}
		else{
			if( _checkAssocThresh( _getDown( *me ), _getDown( *_getRight( *me ) ) ) )
			{
				SkipList< TObject, TModus, RT< TSpec >, TStructuring > * new_list;
				allocate( _getListAlloc( *list ), new_list, 1 );
				valueConstruct( new_list );
				_create( *new_list, _getSucc( *_getDown( *me ) ), _getDown( *_getRight( *me ) ), *_getMainTree( *list ), dim - 1 );
				_setAssocStruct( me, new_list );
			}
			else _setAssocStruct( me, (SkipList< TObject, TModus, RT< TSpec >, TStructuring > * ) 2 );
		}
	}
		
		//	_deleteNextLayer:
		// delete the associated structure <=> the skip list over the elements in the subtree of a lower dimension
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	void 
	_deleteAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		if( _hasBigAssocStruct( me ) )
		{
			SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > * list = _getAssocStruct( me );
			Allocator< ClassPool< SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring >, Unlimited > > * _listAlloc = &_getListAlloc( *list );
			valueDestruct( list );
			deallocate( *_listAlloc, list, 1 );
			_setAssocStruct< TObject, TModus, TSpec, TStructuring >( me, NULL );
		}
		else
			_setAssocStruct< TObject, TModus, TSpec, TStructuring >( me, NULL );
	}

		//	_hasNextLayer:
		// delete the associated structure <=> the skip list over the elements in the subtree of a lower dimension
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	bool 
	_hasAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return _getAssoc( *cargo( *me ) ) != NULL;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	bool 
	_hasSmallAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return ( _getAssoc( *cargo( *me ) ) == ( SkipList< TObject, TModus, RT< TSpec >, TStructuring > * ) 2 );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	bool 
	_hasBigAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me )
	{
		SEQAN_CHECKPOINT
		return ( _getAssoc( *cargo( *me ) ) > ( SkipList< TObject, TModus, RT< TSpec >, TStructuring > * ) 2 );
	}


			// build an associated structure in element elem over the base elements left and right
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void
	_createAssocStruct( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me,
						SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
						TSize dim )
	{
		SEQAN_CHECKPOINT
		if( dim == 0 )
			return;
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * left = _getDown( *me );
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * right = _getDown( *_getRight( *me ) );
		if( _checkAssocThresh( left, right ) )
			_buildAssocStruct( me, list, dim );
		else
			_setAssocStruct( me, ( SkipList< TObject, TModus, RT< TSpec >, TStructuring > *) 2 );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultsSet, typename TSize > inline
	void
	_scanAssocStruct(	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * me,
						SkipList< TObject, TModus, RT< TSpec >, TStructuring > * list,
						TObject * left_border,
						TObject * right_border,
						TSize dim, 
						TResultsSet & results )
	{
		SEQAN_CHECKPOINT
		if( ! _hasAssocStruct( me ) )
		{
			_createAssocStruct( me, list, dim );				
		}
		
		if( _hasSmallAssocStruct( me ) )
		{
			typename Size< TObject >::Type l_dim = dim - 1;
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * buffer = _getDown( * me );
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * end = _getDown( * _getRight( * me ) );
			while( buffer != end )
			{
				if( _testRange( *getObject( buffer ), *left_border, *right_border, l_dim ) ){
					_pushBack( results, getObject( buffer ) );
				}
				goNext( buffer );
			}
			return;
		}
		if( dim > 1 )
			_fingerSearch( _getAssocStruct( me ), left_border, right_border, dim - 1, results );
		else
			_bottomSearch( _getAssocStruct( me ), left_border, right_border, results );
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_performDestructorAction( SkipElement< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		if( _hasBigAssocStruct( &me ) ) 
			_deleteAssocStruct( &me );
		else
			_setAssocStruct< TObject, TModus, TSpec, TStructuring > ( &me, NULL );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //_RT_SKIP_ELEMENT_H
