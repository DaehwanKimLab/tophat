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
  $Id: rmt_base.h 1351 2007-11-30 15:00:38Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_RMT_MAX_BASE_H
#define SEQAN_RMT_MAX_BASE_H

/*
 *  rmt_base.h
 *  rmt
 *
 *  Created by Hendrik Woehrle
 *
 *	Contains basic decalrations, definitions, types
 *
 */

namespace seqan{

/*DISABLED
.Spec.RangeMaximumTree:
..cat:Range Tree
..summary:The RangeMaximumTree is a data structure to find the object with maximal score in a hypercorner.
..signature:RangeTree< TObject, [ TModus, TSpec, TStructuring] >
..param.TObject:Type of stored objects.
..param.TModus:Modus of operation of a RangeTree. A RangeTree is static.
..param.TSpec:Specialization of the RangeTree. Must be $MaxTree<>$.
..param.TStructuring:Parameter to specify whether the RangeTree uses Deferred Data Structuring or not.
..remarks:The object given to the RangeTree should offer the following functions:
..remarks:$key( obj, dim )$: returns the key of the object for dimension $dim$.
..remarks:$setKey( obj, dim, k )$: set the key of the object to $k$ for dimension $dim$.
..remarks:In contrast to STL-like containers, the objects are not cloned by the RangeTree. It only supports searching operations on a set of objects. This set must be handled by the user.
*/

/////////////////////////////////////////////////////////////////////////////////////////
//
//	declarations, types, tags
//
/////////////////////////////////////////////////////////////////////////////////////////

		// spec for the RMT
	template< typename TSpec = Default >
	struct MaxTree
	{};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	class RangeTree;

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipList< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > >
	{
		typedef _Empty Type;
	};

	
/////////////////////////////////////////////////////////////////////////////////////////
//
//	metafunctions
//
/////////////////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Weight:
..summary:Type of the weight of an item. 
..signature:Weight<T>::Type
..param.T:Type for which the weight type is determined.
..returns.param.Type:Weight type of $T$.
..remarks.text:The weight type of an item $T$ is the type which is used by $T$ for scores or priorities.
*/

	template< typename T >
	struct Weight
	{
		typedef int Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< RangeTree< TObject, TModus, TSpec, TStructuring > >
	{
		typedef typename Weight< TObject >::Type Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< RangeTree< TObject, TModus, TSpec, TStructuring > const >
	{
		typedef typename Weight< TObject >::Type const Type;
	};
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipList< TObject, TModus, TSpec, TStructuring > >
	{
		typedef typename Weight< TObject >::Type Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipList< TObject, TModus, TSpec, TStructuring > const >
	{
		typedef typename Weight< TObject >::Type const Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipElement< TObject, TModus, TSpec, TStructuring > >
	{
		typedef typename Weight< TObject >::Type Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipElement< TObject, TModus, TSpec, TStructuring > const >
	{
		typedef typename Weight< TObject >::Type const Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
	{
		typedef typename Weight< TObject >::Type Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Weight< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
	{
		typedef typename Weight< TObject >::Type const Type;
	};

	
/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic data utility functions
//
/////////////////////////////////////////////////////////////////////////////////////////


		// tests, if an object lies in the area of interest and has maximum priority
		// by brute force
	template< typename TObject, typename TSize > inline
	void 
	_testBruteForceMax(	TObject * elem,
						TObject * border,
						TObject *& max_object,
						TSize _dim )
	{
		bool in_range = true;
		TSize dim = 0;
		while( in_range && dim < _dim ){
			in_range = ( key( *elem, dim ) < key( *border, dim ) );
			++dim;
		}
		if( in_range && priority( *elem ) > priority( *max_object ) )
			max_object = elem;
	}

		// test if an element is in range and has a larger priority
	template< typename TObject, typename TBorder, typename TSize > inline
	void 
	_testRangeMax(	TObject * elem,
					TBorder * border,
					TObject *& max_object,
					TSize dim )
	{
		bool in_range = true;
		while( in_range && dim > 0 ){
			in_range = ( key( *elem, dim ) < key( *border, dim ) );
			--dim;
		}
		in_range = in_range && ( key( *elem, dim ) < key( *border, dim ) );
		if( in_range && priority( *elem ) > priority( *max_object ) )
			max_object = elem;
	}


	template< typename TObject, typename TSpec, typename TStructuring > inline
	bool
	_checkAssocThresh( SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * first,
						SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * second )
	{
		SEQAN_CHECKPOINT
		return ( ( second - first ) > 16 );
	}

		// dump of the lowest layer of a RMT
	template< typename TObject, typename TSpec, typename TStructuring, typename TSize1, typename TSize2 >
	void 
	printLayerScores(	SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * list,
						TSize1 layer,
						TSize2 _dim )
	{
		if( layer == 0 ){
			for( typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type j = 0; j < length( *list ) + 1; ++j ){
				std::cout.width(3);
				std::cout<< "___";
			}
			std::cout<<std::endl;
			typename Iterator< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type temp;
			for( typename Size< TObject >::Type dim = 0; dim < dimension( *_getMainTree( *list ) ); ++dim )
			{
				temp = begin( *list );
				--temp;
				while( _pointsTo( temp ) != _pointsTo( end( *list ) ) ){
					std::cout.width(3);
					if( key( *temp, dim ) == infimumValue< typename Key< TObject >::Type >( ) )
						std::cout << std::left << "L";
					else
						std::cout << std::left << key( *temp, dim );
					++temp;
				}
				if( dim == _dim )
					std::cout << "<-";
				std::cout << std::endl;
			}
			std::cout<<std::endl;
			printScores( *list );
		}
		else {
			SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * temp = _getBaseStore( *list );
			while( temp != NULL && temp != end( *list ) ){
				if( _getHeight( *temp ) >= layer ){
					if( _getRight( *( &_getUp( *temp ) + layer - 1) ) ){
						std::stringstream s;
						if( priority( &_getUp( *temp ) + layer - 1 ) == infimumValue< typename Weight< TObject >::Type >() )
							s << std::left << "L";
						else
							s << priority( &_getUp( *temp ) + layer - 1 );
						std::cout.width(3);
						std::cout << std::left << s.str();
					}
					else{
						std::cout.width(3);
						std::cout << priority( &_getUp( *temp ) + layer - 1 );
					}
				}
				else if( temp == _getBaseStore( *list ) )
				{
					std::stringstream s;
					if( priority( &_getUp( *temp ) + layer - 1 ) == infimumValue< typename Weight< TObject >::Type >() )
						s << std::left << "L";
					else
						s << priority( &_getUp( *temp ) + layer - 1 );
					std::cout.width(3);
					std::cout << std::left << s.str();
				}
				else {
					std::cout.width(3);
					std::cout << " ";
				}
				++temp;
			} 
		}
		std::cout << std::endl;
	}

	template< typename TObject, typename TSpec, typename TStructuring >
	void
	printScores( SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > & list )
	{
		SkipBaseElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring >* temp = _getBaseStore( list );
		for( typename seqan::Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type i = 0; i < length( list ) + 1; ++i ){
			std::cout.width(3);
			if( temp == _getBaseStore( list ) )
				std::cout << std::left << "L";
			else
				std::cout << weight( temp );
			++temp;
		}
		std::cout << std::endl;
		temp = _getBaseStore( list );
		for( typename Size< SkipList< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > >::Type i = 0; i < length( list ) + 1; ++i ){
			std::cout.width(3);
			if( priority( temp ) == infimumValue< typename Weight< TObject >::Type >() )
				std::cout << std::left << "L";
			else
				std::cout << priority( temp );
			++temp;
		}
		std::cout << std::endl;
	}

}


#endif // SEQAN_RMT_MAX_BASE_H
