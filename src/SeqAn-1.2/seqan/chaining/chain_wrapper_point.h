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
  $Id: chain_wrapper_point.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_H_CHAIN_WRAPPER_POINT
#define SEQAN_H_CHAIN_WRAPPER_POINT

namespace seqan{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class _WrapperPoint
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
.Class._WrapperPoint:
..summary:Data structure to represent a begin or end point of a fragment
..cat:Chaining
..signature:_WrapperPoint<TBorder>
..param.TBorder:Type of the class that represents the limits (multidimensional point) of an fragment.
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			general functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

			// get the Key of the bording points of the related fragment depends on _end 
	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( _WrapperPoint< TFragType > & me )
	{
		return me._key;
	}

	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( const _WrapperPoint< TFragType > & me )
	{
		return me._key;
	}


		// get the related fragment
	
	template< typename TFragType > inline
	TFragType & 
	_getFrag( _WrapperPoint< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		SEQAN_CHECK( &_getFrag( *me._meta ) == me._frag )
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const _WrapperPoint< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		SEQAN_CHECK( _getFrag( *me._meta ) == me._frag )
		return *me._frag;
	}

	template< typename TFragType > inline
	_MetaFragment< TFragType > & 
	_meta( _WrapperPoint< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		return *me._meta;
	}

	template< typename TFragType > inline
	void
	_setMeta( _WrapperPoint< TFragType > & me,
				_MetaFragment< TFragType > & meta )
	{
		me._meta = &meta;
		me._frag = _getFrag( meta );
	}

	
	template< typename TFragType > inline
	bool 
	_isEnd( _WrapperPoint< TFragType > & me )
	{
		return me._end;
	}

	template< typename TFragType > inline 
	bool 
	_isEnd( const _WrapperPoint< TFragType > & me )
	{
		return me._end;
	}
	
	template< typename TFragType > inline 
	bool
	_isBegin( _WrapperPoint< TFragType > & me )
	{
		return !me._end;
	}

	template< typename TFragType > inline 
	bool
	_isBegin( const _WrapperPoint< TFragType > & me )
	{
		return !me._end;
	}

	template< typename TFragType > inline 
	typename Size< TFragType >::Type
	dimension( const _WrapperPoint< TFragType > & me )
	{
		return dimension( *me._frag );
	}

	
	template< typename TFragType >
	struct _WrapperPoint
	{
			// related triple, which stores the preceding frgament, chain value...
		TFragType * _frag;
			// the related key
		typename Key< TFragType >::Type _key;
			// meta information struct for that point
		_MetaFragment< TFragType > * _meta;
			// the point is either the end ( _end == true ) or the beginning of a fragment
		bool _end;
			
		
	#ifdef _SEQAN_CHAIN_DEBUG
		friend inline
		void 
		dump( _WrapperPoint< TFragType > & me )
		{
			std::cout << "[ ";
			typename Size< TFragType >::Type dim = 0;
			if( me._end )
				std::cout << rightPosition( *me._frag, dim );
			else
				std::cout << leftPosition( *me._frag, dim );
			++dim;
			while( dim != dimension( *me._frag ) )
			{
				if( me._end )
					std::cout << " , " << rightPosition( *me._frag, dim );
				else
					std::cout << " , " << leftPosition( *me._frag, dim );
				++dim;
			}
			std::cout << " ]" << std::endl;
		}

	#endif // _SEQAN_CHAIN_DEBUG

		_WrapperPoint( )
			: _frag( NULL )
			, _key( 0 )
			, _meta( NULL )
			, _end( false )
		{
		}


		_WrapperPoint( const _WrapperPoint & old )
			: _frag( old._frag )
			, _key( old._key )
			, _meta( old._meta )
			, _end( old._end )
		{
		}


		_WrapperPoint( _MetaFragment< TFragType > & meta, 
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( end ? rightPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) : leftPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) )
			, _meta( & meta )
			, _end( end )
		{}

		template< typename TKey >
		_WrapperPoint( _MetaFragment< TFragType > & meta,
						TKey key,
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( key )
			, _meta( & meta )
			, _end( end )
		{}


		_WrapperPoint & operator=( const _WrapperPoint & old )
		{
			if ( this == &old ) 
				return *this;
			_key = old._key;
			_end = old._end;
			_frag = old._frag;
			_meta = old._meta;
			return *this;
		}

	};

}
#endif // SEQAN_H_...
