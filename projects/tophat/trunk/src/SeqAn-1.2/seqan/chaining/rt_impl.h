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
  $Id: rt_impl.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_COMPLETE_RANGE_TREE_H
#define SEQAN_HEADER_COMPLETE_RANGE_TREE_H


namespace seqan
{


/*DISABLED
.Class.RangeTree:
..cat:Range Tree
..summary:The RangeTree is a data structure to solve the orthogonal range searching problem.
..signature:RangeTree< TObject, [ TModus, TSpec, TStructuring] >
..param.TObject:Type of stored objects.
..param.TModus:Modus of operation of a RangeTree. A RangeTree is static.
..param.TSpec:Specialization of the RangeTree.
..param.TStructuring:Parameter to specify whether the RangeTree uses Deferred Data Structuring or not.
..remarks:The object given to the RangeTree should offer the following functions:
..remarks:$key( obj, dim )$: returns the key of the object for dimension $dim$.
..remarks:$setKey( obj, dim, k )$: set the key of the object to $k$ for dimension $dim$.
..remarks:In contrast to STL-like containers, the objects are not cloned by the RangeTree. It only supports searching operations on a set of objects. This set must be handled by the user.
..remarks:The $MaxTree$ specialization offers the abbility to perform Range Maximum Queries.
*/

		
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
class RangeTree< TObject, TModus, RT< TSpec >, TStructuring >
{
public:
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > * _list;

	TObject _RBorderObj;
	TObject _LBorderObj;

	typename Size< TObject >::Type _dim;
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type _numOfElems;

	_SearchPath< TObject, TModus, RT< TSpec >, TStructuring > _sp;
	_RTreeAllocators< TObject, TModus, RT< TSpec >, TStructuring > _allocs;
	
	friend inline
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type 
	length( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._numOfElems;
	}
/*
	friend
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited > > &
	_getElementAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited > > & 
	_getListAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	TObject * 
	_getRBorderObj<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	TObject * 
	_getLBorderObj<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );

	friend
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getList<>( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me );
*/
	template< typename TSize > friend inline
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > **
	_getSearchPath( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me,
					TSize dim )
	{
		SEQAN_CHECKPOINT
		return me._sp._searchPath + _getMaximalSLTowerHeight( me ) * dim;
	}
	
	friend inline
	typename Size< TObject >::Type
	dimension( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._dim;
	}

	RangeTree()
	{}

public:
	
	template< typename TContainer, typename TSize >
	RangeTree(	TContainer & data,
				TSize dim )
	: _dim( dim )
	, _numOfElems( length( data ) )
	, _sp( _getMaximalSLTowerHeight( length( data ) ) * _dim )
	, _allocs( 2 * length( data ) * dim * dim * dim, length( data ) * dim * dim )
	{
		SEQAN_CHECKPOINT
		_setMinInfty( _LBorderObj, _dim );
		_setMaxInfty( _RBorderObj, _dim );

		allocate( _getListAlloc( *this ), _list, 1 );
		valueConstruct( _list );
		_create( *_list, data, *this, _dim - 1 );
	}


	~RangeTree()
	{	
		SEQAN_CHECKPOINT
		_clearSearchPath( this->_sp, _getMaximalSLTowerHeight( length( *this ) ) * _dim );
	}

}; // struct RangeTree

}

#endif // SEQAN_HEADER_COMPLETE_RANGE_TREE_H




