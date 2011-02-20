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
  $Id: chain_point.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CHAINPOINT_H
#define SEQAN_HEADER_CHAINPOINT_H

namespace seqan
{


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class _ChainPoint
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Class._ChainPoint:
..summary:Basic data structure to represent the end-point of a fragment for use in a RMT
..cat:Chaining
..signature:_ChainPoint< TFragType, TSpec >
..param.TFragType:Type of the class which represents the limits (multidimensional point) of an fragment.
..param.TSpec:Spec of the _ChainPoint.
*/
	template< typename TFragType, typename TSpec >
	struct Spec< _ChainPoint< TFragType, TSpec > >
	{
		typedef TSpec Type;
	};

	template< typename TFragType, typename TSpec, typename TSize > inline
	typename Key< TFragType >::Type 
	key( _ChainPoint< TFragType, TSpec > & me,
			TSize dim )
	{
SEQAN_CHECK2( me._coords != NULL, "point not initialized" )
SEQAN_CHECK2( dim <= me._dim, "dim corrupted" )
		return me._coords[ dim ];
	}

	template< typename TFragType, typename TSpec, typename TSize > inline 
	typename Key< TFragType >::Type 
	key( const _ChainPoint< TFragType, TSpec > & me,
			TSize dim )
	{
SEQAN_CHECK2( me._coords != NULL, "point not initialized" )
SEQAN_ASSERT( dim <= me._dim )
		return me._coords[ dim ];
	}

	template< typename TFragType, typename TSpec, typename TSize, typename TKey > inline 
	void 
	setKey( _ChainPoint< TFragType, TSpec > & me,
			TSize dim,
			TKey val )
	{
SEQAN_ASSERT( dim <= me._dim )
		me._coords[ dim ] = val;
	}

	template< typename TFragType, typename TSpec, typename TSize, typename TKey > inline 
	void 
	setKey( const _ChainPoint< TFragType, TSpec > & me,
				TSize dim,
				TKey val )
	{
SEQAN_ASSERT( dim <= me._dim )
		me._coords[ dim ] = val;
	}

	template< typename TFragType, typename TSpec, typename TSize > inline 
	void 
	_incKey(  _ChainPoint< TFragType, TSpec > & me,
				TSize dim )
	{
SEQAN_ASSERT( dim <= me._dim )
		++me._coords[ dim ];
	}

		// get the related fragment
	template< typename TFragType, typename TSpec > inline
	TFragType & 
	_getFrag( _ChainPoint< TFragType, TSpec > & me )
	{
SEQAN_ASSERT( me._meta != NULL )
		return _getFrag( *me._meta );
	}

	template< typename TFragType, typename TSpec > inline
	TFragType & 
	_getFrag( const _ChainPoint< TFragType, TSpec > & me )
	{
SEQAN_ASSERT( me._meta != NULL )
		return _getFrag( *me._meta );
	}

		// get the related metainformation struct
	template< typename TFragType, typename TSpec > inline
	_MetaFragment< TFragType > & 
	_meta( _ChainPoint< TFragType, TSpec > & me )
	{
		return *me._meta;
	}

	template< typename TFragType, typename TSpec > inline
	_MetaFragment< TFragType > & 
	_meta( const _ChainPoint< TFragType, TSpec > & me )
	{
		return *me._meta;
	}

		// set the related metainformation struct
	template< typename TFragType, typename TSpec > inline
	void 
	_setMeta( _ChainPoint< TFragType, TSpec > & me,
				_MetaFragment< TFragType > & meta )
	{
		me._meta = meta;
	}

	template< typename TFragType, typename TSpec > inline
	void 
	_setMeta( const _ChainPoint< TFragType, TSpec > & me,
				_MetaFragment< TFragType > & meta )
	{
		me._meta = meta;
	}

	template< typename TFragType, typename TSpec > inline
	typename Size< TFragType >::Type
	dimension( _ChainPoint< TFragType, TSpec > & me )
	{
		return me._dim;
	}

	template< typename TFragType, unsigned int ISize > inline
	typename Size< TFragType >::Type
	dimension( _ChainPoint< TFragType, Array< ISize > > &)
	{
		return ISize;
	}

	template< typename TFragType, typename TSpec > inline
	typename Weight< TFragType >::Type
	priority( _ChainPoint< TFragType, TSpec > & me )
	{
		return me._prio;
	}

	template< typename TFragType, typename TSpec, typename TPrio > inline
	void
	setPriority( _ChainPoint< TFragType, TSpec > & me,
					TPrio prio )
	{
		me._prio = prio;
	}


	template< typename TFragType, typename TSpec >
	struct _ChainPoint
	{
		
		_MetaFragment< TFragType > * _meta;
		typename Key< TFragType >::Type * _coords;
		typename Size< TFragType >::Type _dim;
		typename Weight< TFragType >::Type _prio;

	public:

			// standard constructor for use in skiplist
		_ChainPoint(  )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( 1 )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
		}

		_ChainPoint( typename Size< TFragType >::Type dim )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( dim )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
		}

		_ChainPoint( _MetaFragment< TFragType > & meta, 
						bool begin = false )
		: _meta( &meta )
		, _coords( NULL )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			_dim = dimension( *frag );
			allocate( *this, _coords, _dim );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		template< typename TSize >
		_ChainPoint( _MetaFragment< TFragType > & meta,
						TSize dim, 
						bool begin = false )
		: _meta( &meta )
		, _coords( NULL )
		, _dim( dim )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			allocate( *this, _coords, _dim );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		_ChainPoint( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim,
						_MetaFragment< TFragType > * meta, 
						bool begin = false )
		: _meta( meta )
		, _coords( NULL )
		, _dim( dim )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}

		_ChainPoint( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim )
		: _meta( NULL )
		, _coords( NULL )
		, _dim( dim )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}


		~_ChainPoint()
		{
			deallocate( *this, _coords, _dim );
			_coords = NULL;
			_meta = NULL;
			_dim = 0;
		}

		_ChainPoint( const _ChainPoint & old )
		{
			_meta = old._meta;
			_dim = old._dim;
			_prio = old._prio;
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
		}


		_ChainPoint & operator=( const _ChainPoint & old )
		{
			if ( this == &old ) 
				return *this;
			_meta = old._meta;
			if( _coords )
				deallocate( *this, _coords, _dim );
			_dim = old._dim;
			_prio = old._prio;
			allocate( *this, _coords, _dim );
			for( typename Size< TFragType >::Type i = 0; i < _dim; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
			return *this;
		}

		friend inline
		void
		dump( _ChainPoint & me )
		{
			std::cout << "[ ";
			typename Size< Seed< TFragType > >::Type dim = 0;
			std::cout << key( me, dim );
			++dim;
			while( dim != me._dim )
			{
				std::cout << " , " << key( me, dim );
				++dim;
			}
			std::cout << " ] "<< me._prio << std::endl;
		}


	};


	template< typename TFragType, unsigned int ISize>
	struct _ChainPoint< TFragType, Array< ISize > >
	{
		
		typename Key< TFragType >::Type _coords[ISize];
		_MetaFragment< TFragType > * _meta;
		typename Weight< TFragType >::Type _prio;

	public:

			// standard constructor for use in skiplist
		_ChainPoint(  )
		: _meta( NULL )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
		}

		_ChainPoint( typename Size< TFragType >::Type )
		: _meta( NULL )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
		}

		_ChainPoint( _MetaFragment< TFragType > & meta, 
						bool begin = false )
		: _meta( &meta )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		template< typename TSize >
		_ChainPoint( _MetaFragment< TFragType > & meta,
						TSize, 
						bool begin = false )
		: _meta( &meta )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			TFragType * frag = &_getFrag( meta );
			if( begin ){
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = leftPosition( *frag, i );
				}
			}
			else{
				for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
				{
					_coords[ i ] = rightPosition( *frag, i );
				}
			}
		}

		_ChainPoint( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim,
						_MetaFragment< TFragType > * meta, 
						bool begin = false )
		: _meta( meta )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}

		_ChainPoint( typename Key< TFragType >::Type * coords,
						typename Size< TFragType >::Type dim )
		: _meta( NULL )
		, _prio( infimumValue< typename Weight< TFragType >::Type >() )
		{
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = coords[ i ];
			}
		}


		~_ChainPoint()
		{
			_meta = NULL;
		}

		_ChainPoint( const _ChainPoint & old )
		{
			_meta = old._meta;
			_prio = old._prio;
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
		}


		_ChainPoint & operator=( const _ChainPoint & old )
		{
			if ( this == &old ) 
				return *this;
			_meta = old._meta;
			_prio = old._prio;
			for( typename Size< TFragType >::Type i = 0; i < ISize; ++i )
			{
				_coords[ i ] = old._coords[ i ];
			}
			return *this;
		}

		friend inline
		void
		dump( _ChainPoint & me )
		{
			std::cout << "[ ";
			typename Size< Seed< TFragType > >::Type dim = 0;
			std::cout << key( me, dim );
			++dim;
			while( dim != ISize )
			{
				std::cout << " , " << key( me, dim );
				++dim;
			}
			std::cout << " ] "<< me._prio << std::endl;
		}


	};


}	// namespace

#endif // SEQAN_HEADER_CHAINPOINT_H
