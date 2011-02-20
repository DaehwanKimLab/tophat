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
  $Id: rt_base.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file


#ifndef SEQAN_RT_BASE_H
#define SEQAN_RT_BASE_H



namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	declarations
//
/////////////////////////////////////////////////////////////////////////////////////////
		
		// standard tag struct for range tree
	template< typename TSpec = Default >
	struct RT
	{};

		// tag structs for grade of deferredness of the range tree
	struct SemiDeferred
	{};

	struct Complete
	{};

	struct Deferred
	{};

		// main classes
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct SkipList;
	
	template< typename TObject, typename TModus, typename TSpec = Default, typename TStructuring = Complete >
	class RangeTree;


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef _Empty Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef typename Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type Type;
	};

/////////////////////////////////////////////////////////////////////////////////////////
//
//	dependent members
//
/////////////////////////////////////////////////////////////////////////////////////////

		// the memory allocators of the range tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	_RTreeAllocators;

	template< typename TObject, typename TSpec, typename TStructuring >
	struct
	_RTreeAllocators< TObject, SkipListStatic, TSpec, TStructuring >
	{
		
		Allocator< ClassPool< SkipElement< TObject, SkipListStatic, TSpec, TStructuring >, Limited > > _elementAlloc;
		Allocator< ClassPool< SkipList< TObject, SkipListStatic, TSpec, TStructuring >, Unlimited > > _listAlloc;
		Allocator< SimpleAlloc<> > _baseAlloc;

		_RTreeAllocators()
			: _elementAlloc( NULL )
		{}

		template< typename TSize >
		_RTreeAllocators( TSize size1, TSize size2 )
			:_elementAlloc( size1 )
			,_listAlloc( size2 )
		{
			SEQAN_CHECKPOINT
		}

		~_RTreeAllocators( )
		{
			SEQAN_CHECKPOINT
		}
	};



/////////////////////////////////////////////////////////////////////////////////////////
//
//	utilities
//
/////////////////////////////////////////////////////////////////////////////////////////
		

	const size_t _rt_thresh = 16;

		
	template< typename TObject, typename TSpec, typename TStructuring > inline
	bool
	_checkAssocThresh( SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * first,
						SkipBaseElement< TObject, SkipListStatic, TSpec, TStructuring > * second )
	{
		SEQAN_CHECKPOINT
		return ( ( second - first ) > _rt_thresh );
	}

	template< typename TTarget, typename TSource > inline
	void
	_pushBack( TTarget & target, TSource const & source )
	{
		SEQAN_CHECKPOINT
#ifdef RTTIMETEST
		volatile TSource temp = source;
#else
		appendValue( target, source );
#endif
	}


			// setting the element to - infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMinInfty(	TObject & me,  
					TSize dim )
	{
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type infValue = infimumValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, infValue );
		}
	}

		// setting the element to + infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMaxInfty(	TObject & me,  
					TSize dim )
	{
		
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type supValue = supremumValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, supValue );
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic accessor functions
//
/////////////////////////////////////////////////////////////////////////////////////////
	

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited, SimpleAllocator > > & 
	_getElementAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._elementAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited, SimpleAllocator > > & 
	_getListAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._listAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._baseAlloc;
	}

	
		// accessor für grenzobjekte
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	TObject * 
	_getLBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._LBorderObj;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	TObject * 
	_getRBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._RBorderObj;
	}

		// skip list of the main tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getList( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._list;
	}



	
/////////////////////////////////////////////////////////////////////////////////////////
//
//	algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	
	
/*DISABLED
.Function.rangeQuery:
..summary:Get the object with maximum priority in the RMT in a given intervall
..cat:Range Tree
..signature:rangeMaxQuery(tree, lower_border, upper_border, dest)
..param.tree:A Range Tree.
...type:RangeTree
..param.lower_border:The object that stores the lower borders for all dimensions, i.e. $key( lower_border ) <= key( point in range )$
..param.lower_border:The object that stores the upper borders for all dimensions, i.e. $key( point in range ) <= key( upper_border )$
..param.dest:A container to save the objects.
..remarks:The size of $dest$ should be sufficient.
*/

		// perform a range query
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet >
	void
	rangeQuery( RangeTree< TObject, TModus, TSpec, TStructuring > & me, 
				TObject & first, 
				TObject & second,
				TResultSet & results )
	{
		SEQAN_CHECKPOINT
		if( dimension( me ) > 1 )
			_fingerSearch( _getList( me ), &first, &second, dimension( me ) - 1, results );
		else
			_bottomSearch( _getList( me ), &first, &second, results );
	}

	template< typename TObject, typename TSize > inline
	bool 
	_testBruteForce( TObject & elem,
					 TObject & first,
					 TObject & second,
					 TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Size< TObject >::Type _dim = 0;
		while( in_range && _dim <= dim )
		{
			in_range = ( ( key( first, _dim ) <= key( elem, _dim ) ) && ( key( elem, _dim ) <= key( second, _dim ) ) );
			++_dim;
		}
		return in_range;
	}

		// test if an element is in range
		// from dim to dim - 1 to 0
	template< typename TObject, typename TSize > inline
	bool 
	_testRange(	TObject & elem,
				TObject & first,
				TObject & second,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Key< TObject >::Type theKey;
		while ( in_range && dim > 0 )
		{
			theKey = key( elem, dim );
			in_range = ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
			--dim;
		}
		theKey = key( elem, dim );
		in_range &= ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
		return in_range;
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type
	_getMaximalSLTowerHeight( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & rt )
	{
		SEQAN_CHECKPOINT
		return log2( length( rt ) ) + 2;
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TSize >
	void 
	printLayerScores(	SkipList< TObject, SkipListStatic, RT< TSpec >, TStructuring > * list,
						TSize layer,
						TSize _dim )
	{}

}

#endif // SEQAN_RT_BASE
