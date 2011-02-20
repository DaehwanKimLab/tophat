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
  $Id: rt_sl_impl.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_CRT_SKIPLIST_H
#define SEQAN_HEADER_CRT_SKIPLIST_H



namespace seqan
{
	

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct SkipList< TObject, TModus, RT< TSpec >, TStructuring >
{	
	// private members:
		// container for elements on the right side
	SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * _leftSideStore;
		// container for elements in the base(lowest) layer
	union{
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * _baseStore;	
		SkipList< TObject, TModus, RT< TSpec >, TStructuring > * _next;	
	};
		// the main container class containing meta-information and interface functions
	RangeTree< TObject, TModus, RT< TSpec >, TStructuring > * _mainTree;
		// stores number of elements used
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type _numOfElements;
		// dimension of the list
	typename Size< TObject >::Type _dim;


//*************************************** internal help functions: ************************************************


	// **********************  Special SkipList< TObject, TModus, RT< TSpec >, TStructuring > accessor functions ***************************

	friend inline
	typename Size< TObject >::Type
	dimension( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._dim;
	}

	friend inline
	RangeTree< TObject, TModus, RT< TSpec >, TStructuring > *
	_getMainTree( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._mainTree;
	}

	SkipList & operator=( const SkipList & old )
	{}

	SkipList( const SkipList & old )
	{}
		
public:

		// **********************  Constructors ***************************
	SkipList()
	{
		mtRandInit();
	};


	~SkipList()
	{
		Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited > > & elem_pool_buffer = _getElementAlloc( *_getMainTree( *this ) );
		
			// destroy towers
		_setHeight( *_baseStore, _getMaximalSLTowerHeight( _numOfElements ) );
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * up_elem = &_getUp( *_baseStore );
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * up_buffer = up_elem;
		SkipElement< TObject, TModus, RT< TSpec >, TStructuring > * right_elem = &_getUp( _baseStore[ _numOfElements + 1 ] );
		typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type height;
		while( up_elem != right_elem ){
			height = _getHeight( *_getDown( *up_elem ) );
			up_buffer = _getRight( *up_elem );
			arrayDestruct( up_elem, up_elem + height );
			deallocate( elem_pool_buffer, up_elem, height );
			up_elem = up_buffer;
		}
		arrayDestruct( up_elem, up_elem + 2 );
		deallocate( elem_pool_buffer, up_elem, 2 );
			
			// destroy base layer
		arrayDestruct( _baseStore, _baseStore + _numOfElements + 2 );
		deallocate( _getBaseAlloc( *this ), _baseStore, _numOfElements + 2 );
	}
	
	
	friend
	void 
	printLayer(	SkipList< TObject, TModus, RT< TSpec >, TStructuring > & me,
				typename seqan::Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type layer,
				int column )
	{
		if( layer == 0 )
		{
			for( typename seqan::Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type j = 0; j < 11; ++j )
			{
				std::cout<< "______";
			}
			std::cout<<std::endl;
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * temp = me._baseStore;
			temp += column;
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * tempEnd = temp;
			tempEnd += me._numOfElements;
			while( temp != tempEnd )
			{
				std::cout.width(7);
				if( key( *temp ) == infimumValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
				else
					std::cout << std::left << key( *temp );
				goNext( temp );
			}
			//std::cout<<std::endl;
			//printCounts( me );
			//std::cout<<std::endl;
			//printSorting( me );
			//std::cout<<std::endl;
			//printHeights( me );
			//std::cout<<std::endl;
			//printValues( me );
		}
		else if( layer <= _getCurrentLayer( me ) )
		{
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * temp = me._baseStore;
			temp += column;
			SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > * tempEnd = temp;
			tempEnd += me._numOfElements;
			while( temp != tempEnd )
			{
				if( _getHeight( *temp ) >= layer )
				{
					if( _getRight( *( &_getUp( *temp ) + layer - 1) ) )
					{
						std::stringstream s;
						if( key( *temp ) == infimumValue< typename Key< TObject >::Type >( ) )
							s << std::left << "L";
						else
							s << key( *temp );
						s << ">";
						std::cout.width(7);
						std::cout << std::left << s.str();
					}
					//else
					//	dump( *( &_getUp( *hostIterator( temp ) ) + layer - 1 ) );
				}
				else 
				{
					std::cout.width(7);
					std::cout << " ";
				}
				goNext( temp );
			} 
		}
		std::cout<<std::endl;
	}

	void 
	dump(  )
	{
		int column = 0;

		while( column < length( *this ) )
		{
			int j = _getMaximalSLTowerHeight( *this );
			while( j > 0 ){
				printLayer( *this, --j, column );
			}
			std::cout<< std::endl;
			column += 14;
		}
	}

	friend
	void 
	printCounts( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list )
	{
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring  >* temp = &list._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type i = 0; i < list._numOfElements + 1; ++i ){
			std::cout.width(3);
			std::cout << std::left << _getCountImpl( *temp );
			++temp;
		}
	}


	friend
	void
	printHeights( SkipList< TObject, TModus, RT< TSpec >, TStructuring > & list )
	{
		SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring  >* temp = &list._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type i = 0; i < list._numOfElements + 1; ++i ){
			std::cout.width(3);
			std::cout << std::left << _getHeight( *temp );
			++temp;
		}
	}


}; // struct SkipList


}

#endif	// SEQAN_HEADER_CRT_SKIPLIST_H
