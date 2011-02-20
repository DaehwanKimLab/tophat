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
  $Id: tree_chain_utils.h 1045 2007-08-21 16:41:41Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_TREECHAIN_UTILS_H
#define SEQAN_HEADER_TREECHAIN_UTILS_H



namespace seqan{

//////////////////////////////////////////////////////////////////////////////////////////
//	helper functions
//////////////////////////////////////////////////////////////////////////////////////////

		// calculate faculty of a number
	template< typename TSimpleType > inline
	TSimpleType
	fac( TSimpleType data )
	{
		TSimpleType value = 1;
		for( TSimpleType i = 2; i <= data; ++i )
		{
			value *= i;
		}
		return value;
	}

		// initailize a permutation to 0, 1, ..., k
	template< typename TPerm, typename TSize > inline
	void
	_initPerm( TPerm & perm, 
				TSize length )
	{
		for( TSize i = 0; i < length; ++i )
		{
			appendValue( perm, i );
		}
	}
		
		// reset a permutation to 0, 1, ..., k
	template< typename TPerm > inline
	void
	_resetPerm( TPerm & perm )
	{
		typename Iterator< TPerm, Rooted >::Type permIt = begin( perm );
		int i = 0;
		while( permIt != end( perm ) )
		{
			assignValue( permIt, i );
			++i;
			goNext( permIt );
		}
	}

		// get the permutation e.g. order of to sequences of numbers
	template< typename TData, typename TItPerm, typename TSize > inline
	void
	_getPerm( TData * values,
				TItPerm perm,
				TSize length )
	{
		TData max = infimumValue< TData>();
		TSize maxIndex = 0;
		for( TSize i = 0; i < length; ++i )
		{
			for( TSize j = 0; j < length; ++j )
			{
				if( values[ j ] > max )
				{
					max = values[ j ];
					maxIndex = j;
					values[ j ] = infimumValue< TData>();
				}
			}
			assignValue( perm, maxIndex );
			maxIndex = 0;
			max = infimumValue< TData>();
			goNext( perm );
		}		
	}

		// get the differences of coordinates of two fragments
	template< typename TData, typename TSize, typename FragType >
	void
	_getPermDifference( TData * values,
						TSize dim,
						FragType & upper,
						FragType & lower )
	{
		for( TSize i = 0; i < dim; ++i )
		{
			values[ i ] = leftPosition( upper, i ) - rightPosition( lower, i );
		}		
	}

//////////////////////////////////////////////////////////////////////////////////////////
//	transformations
//////////////////////////////////////////////////////////////////////////////////////////

		// transform the coordinates of a chain point
	template< typename FragType, typename SpecType, typename TPerm > inline
	void
	_chainTransformCoords( _ChainPoint< FragType, SpecType > & point_src,
							_ChainPoint< FragType, SpecType > & point_dst,
							TPerm & perm )
	{
		typename Iterator< TPerm >::Type permIt = begin( perm );
		typename Iterator< TPerm >::Type permEnd = end( perm );
		goPrevious( permEnd );
		while( permIt != permEnd )
		{
			setKey( point_dst, *permIt, key( point_src, *permIt ) - key( point_src, *( permIt + 1 ) ) );
			++permIt;
		}
		setKey( point_dst, *permIt, key( point_src, *permIt ) );
	}

			// transform the coordinates of a chain point for searching
	template< typename FragType, typename SpecType, typename TPerm > inline
	void
	_chainTransformCoordsSearch( _ChainPoint< FragType, SpecType > & point_src,
									_ChainPoint< FragType, SpecType > & point_dst,
									TPerm & perm )
	{
		typename Iterator< TPerm >::Type permIt = begin( perm );
		typename Iterator< TPerm >::Type permEnd = end( perm );
		goPrevious( permEnd );
		while( permIt != permEnd )
		{
			setKey( point_dst, *permIt, key( point_src, *permIt ) - key( point_src, *( permIt + 1 ) ) + 1 );
			++permIt;
		}
		setKey( point_dst, *permIt, key( point_src, *permIt ) );
	}


//////////////////////////////////////////////////////////////////////////////////////////
// cost functions
//////////////////////////////////////////////////////////////////////////////////////////
	
		// calulate the score
	template< typename FragType, typename SpecType, typename TCostModell, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( _MetaFragment< FragType > & last_meta,
						_MetaFragment< FragType > & current_meta,
						_ChainPoint< FragType, SpecType > & point,
						TCostModell cost, 
						TScore const & score_,
						TSize dim );

	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( _MetaFragment< FragType > & ,
					_MetaFragment< FragType > & current_meta,
					_ChainPoint< FragType, SpecType > & point,
					G_0_Cost, 
					TScore const &,
					TSize )
	{
		return weight( current_meta ) + priority( point );
	}


		// calculate the maximum priority
	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_maxPriority( _MetaFragment< FragType > &,
					_MetaFragment< FragType > & current_meta,
					_ChainPoint< FragType, SpecType > & point,
					G_1_Cost, 
					TScore const & score_,
					TSize dim )
	{
		typename Weight< FragType >::Type prio = weight( current_meta );
		prio += score( _meta( point ) );
		prio -= _costG1( current_meta, _meta( point ), score_, dim );
		return prio;
	}

	
		// calculate the priority for activation
	template< typename FragType, typename SpecType, typename TCostModell, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( _MetaFragment< FragType > & last_meta,
						_ChainPoint< FragType, SpecType > & point,
						TCostModell cost, 
						TScore const & score_,
						TSize dim );


	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( _MetaFragment< FragType > &,
						_ChainPoint< FragType, SpecType > &,
						G_0_Cost, 
						TScore const &,
						TSize )
	{
		return 0;
	}


	template< typename FragType, typename SpecType, typename TScore, typename TSize > inline
	typename Weight< FragType >::Type
	_activatePriority( _MetaFragment< FragType > & last_meta,
						_ChainPoint< FragType, SpecType > & point,
						G_1_Cost, 
						TScore const & score_,
						TSize dim )
	{
		return _costG1( last_meta, _meta( point ), score_, dim );
	}

	
		// the cost function for manhattan metric
	template< typename FragType, typename TScore, typename TSize >
	typename Weight< FragType >::Type
	_costG1( _MetaFragment< FragType > & upper,
				_MetaFragment< FragType > & lower,
				TScore const & score,
				TSize dim )
	{
		typename Weight< FragType >::Type weight = 0;
		for( typename Size< FragType >::Type i = 0; i < dim; ++i )
		{
//!!!Change dist semantics
//			weight += ( scoreGapExtend( score ) * ( leftPosition( _getFrag( upper ), i ) - rightPosition( _getFrag( lower ), i ) ) );
			weight += ( scoreGap( score ) * ( leftPosition( _getFrag( upper ), i ) - rightPosition( _getFrag( lower ), i ) - 1) );
		}
		return weight;
	}

		// the cost function for SoP metric
	template< typename FragType, typename TItPerm, typename TScoreValue, typename TScoreType, typename TSize >
	typename Weight< FragType >::Type
	_costGSoP( _MetaFragment< FragType > & upper,
				_MetaFragment< FragType > & lower,
				Score< TScoreValue, TScoreType > const & score,
				TItPerm permBeg,
				TItPerm,
				TSize dim )
	{
		typename Weight< FragType >::Type weight = 0;
		typename Weight< FragType >::Type weight_buffer = 0;
		typename Weight< FragType >::Type delta;
		typename Size< FragType >::Type dim_factor = dim - 1;
		for( typename Size< FragType >::Type i = 0; i < dim; ++i )
		{
//!!!Change dist semantics
//			delta = static_cast< typename Weight< FragType >::Type >( leftPosition( _getFrag( upper ), value( permBeg ) ) - rightPosition( _getFrag( lower ), value( permBeg ) ) );
			delta = static_cast< typename Weight< FragType >::Type >( leftPosition( _getFrag( upper ), value( permBeg ) ) - rightPosition( _getFrag( lower ), value( permBeg ) ) -1);
			weight_buffer = ( scoreGap( score ) * delta );
			weight_buffer *= static_cast< typename Weight< FragType >::Type >( ( dim_factor - i ) );
			weight += weight_buffer;
			weight_buffer = ( delta * ( scoreMismatch( score ) - scoreGap( score ) ) );
			weight_buffer *= static_cast< typename Weight< FragType >::Type >( i );
			weight += weight_buffer;
			goNext( permBeg );
		}
		return weight;
	}

//////////////////////////////////////////////////////////////////////////////////////////
//	sorting
//////////////////////////////////////////////////////////////////////////////////////////


		// struct for std::sort of wrapper points
	template< typename T >
	struct
	_ChainSorter
	{
		inline bool
		operator()( T & first, T & second  )
		{
			if ( key( first ) < key( second ) )
				return true;
			else if ( key( first ) == key( second ) && _isBegin( first ) && _isEnd( second ) )
				return true;
			return false;
		}
		inline bool 
		operator()( const T & first, const T & second  )
		{
			if ( key( first ) < key( second ) )
				return true;
			else if ( key( first ) == key( second ) && _isBegin( first ) && _isEnd( second ) )
				return true;
			return false;
		}

		_ChainSorter()
		{}	

	};
	

//////////////////////////////////////////////////////////////////////////////////////////
//	dynamic programming helper functions
//////////////////////////////////////////////////////////////////////////////////////////


		// backtracking
	template< typename TDest, typename TMetas >
	typename Weight< typename Value< TDest >::Type >::Type
	_chain_trace( TDest & dest,
					TMetas & metas )
	{
		typedef typename Value< TDest >::Type FragType;
		typename Iterator< TMetas >::Type meta = end( metas );
		goPrevious( meta );
		_MetaFragment< FragType > * pMeta = & value( meta );
		typename Weight< FragType >::Type chain_score = score( *pMeta );
		//pMeta = &_getPred( *pMeta );
		while( pMeta != &value( begin( metas ) ) )
		{
			SEQAN_CHECK( &_getFrag( *pMeta ) != 0 )
			appendValue( dest, _getFrag( *pMeta ) );
			pMeta = &_getPred( *pMeta );
		}
		appendValue( dest, _getFrag( *pMeta ) );
		std::reverse( begin( dest ), end( dest ) );
//		typename Iterator< TDest >::Type destIt = begin( dest );
		return chain_score;
	}


//////////////////////////////////////////////////////////////////////////////////////////
//	initialization helper functions
//////////////////////////////////////////////////////////////////////////////////////////

		// init of starting frag (the origin)
	template< typename FragType, typename TSize > inline
	void
	_init_starting_frag( FragType & frag,
							TSize dim )
	{
		TSize dim_counter = 0;
		while( dim_counter != dim )
		{
			_setLeftPosition( frag, dim_counter, 0 );
			_setRightPosition( frag, dim_counter, 0 );
			++dim_counter;
		}
		setWeight( frag, 0 );
	}

		// init of needed variables
		// * the wrapper points
		// * the chain points
		// * the metainformation structures
		// spec for G0 and G1 metric
	template< typename FragType, typename TSource, typename TMetas, typename TWPoints, typename TCPoints, typename TSpec > inline
	void
	_build_chain_environment( TSource & source,  
								TMetas & metas, 
								TWPoints & wPoints,
								TCPoints & cPoints,
								FragType & startingFrag,
								FragType & endFrag,
								TSpec )
	{
		typedef typename Key< FragType >::Type KeyType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;
		SizeType dim = dimension( value( begin( source ) ) );
		SizeType lower_dim  = dim - 1;
				
			// initialize starting frag
		_init_starting_frag( startingFrag, dim );
		appendValue( metas, _MetaFragment< FragType >( startingFrag ) );
		wPoints.push_back( _WrapperPoint< FragType >( value( begin( metas ) ), true ) );
		appendValue( cPoints, _ChainPoint< FragType, SpecType >( value( begin( metas ) ) ) );

			// buffers to find the maximal coordinates
		KeyType * maxCoords;
		allocate( maxCoords, maxCoords, dim );
		for( SizeType i = 0; i < dim; ++i )
		{
			maxCoords[ i ] = 0;
		}

			// traverse the set of fragments
		typename Iterator< TSource, Rooted >::Type sourceIt = begin( source );
		typename Iterator< TMetas, Rooted >::Type metaIt = begin( metas );
		goNext( metaIt );
		while( sourceIt != end( source ) )
		{
			appendValue( metas, _MetaFragment< FragType >( value( sourceIt ) ) );
			wPoints.push_back( _WrapperPoint< FragType >( value( metaIt ), leftPosition( value( sourceIt ), lower_dim ), false ) );
			wPoints.push_back( _WrapperPoint< FragType >( value( metaIt ), rightPosition( value( sourceIt ), lower_dim ), true ) );
			appendValue( cPoints, _ChainPoint< FragType, SpecType >( value( metaIt ) ) );

			for( SizeType i = 0; i < dim; ++i )
			{
				if( maxCoords[ i ] < rightPosition( value( sourceIt ), i ) )
					maxCoords[ i ] = rightPosition( value( sourceIt ), i );
			}
			
			goNext( metaIt );
			goNext( sourceIt );
		}

			// set the coordinates of the terminus
		for( SizeType dim_counter = 0; dim_counter < dim; ++dim_counter )
		{
			_setLeftPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
			_setRightPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
		}
		appendValue( metas, _MetaFragment< FragType >( endFrag ) );
		wPoints.push_back( _WrapperPoint< FragType >( *metaIt, false ) );
		deallocate( maxCoords, maxCoords, dim );
		appendValue( cPoints, _ChainPoint< FragType, SpecType >( *metaIt ) );
	}


		// spec for G_SoP metric
	template< typename FragType, typename TSource, typename TMetas, typename TWPoints, typename TPoints, typename TPerm, typename TSpec > inline
	void
	_build_chain_environment( TSource & source, 
								TMetas & metas, 
								TWPoints & wPoints,
								TPoints  & tPoints,
								FragType & startingFrag,
								FragType & endFrag,
								TPerm & perm,
								typename Size< FragType >::Type fac,
								TSpec &)
	{
		typedef typename Key< FragType >::Type KeyType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;
		SizeType dim = dimension( value( begin( source ) ) );
		SizeType lower_dim = dim - 1;
		SizeType dim_counter = 0;
		
		_init_starting_frag( startingFrag, dim );

		appendValue( metas, _MetaFragment< FragType >( startingFrag ) );
		wPoints.push_back( _WrapperPoint< FragType >( value( begin( metas ) ), true ) );
		
		KeyType * maxCoords;
		allocate( maxCoords, maxCoords, dim );
		for( SizeType i = 0; i < dim; ++ i )
			maxCoords[ i ] = 0;

		typename Iterator< TSource, Rooted >::Type sourceIt = begin( source );
		typename Iterator< TMetas, Rooted >::Type metaIt = begin( metas );
		goNext( metaIt );

		while( sourceIt != end( source ) )
		{
			appendValue( metas, _MetaFragment< FragType >( value( sourceIt ) ) );
			wPoints.push_back( _WrapperPoint< FragType >( value( metaIt ), leftPosition( value( sourceIt ), lower_dim ), false ) );
			wPoints.push_back( _WrapperPoint< FragType >( value( metaIt ), rightPosition( value( sourceIt ), lower_dim ), true ) );

			for( dim_counter = 0; dim_counter < dim; ++dim_counter )
			{
				if( maxCoords[ dim_counter ] < rightPosition( value( sourceIt ), dim_counter ) )
					 maxCoords[ dim_counter ] = rightPosition( value( sourceIt ), dim_counter );
			}

			goNext( metaIt );
			goNext( sourceIt );
		}
		dim_counter = 0;
		while( dim_counter != dim )
		{
			_setLeftPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1);
			_setRightPosition( endFrag, dim_counter, maxCoords[ dim_counter ] + 1 );
			++dim_counter;
		}
		appendValue( metas, _MetaFragment< FragType >( endFrag ) );
		deallocate( maxCoords, maxCoords, dim );
		wPoints.push_back( _WrapperPoint< FragType >( *metaIt, false ) );
		metaIt = begin( metas );
		
			// transformate the point coordinates
		typename Iterator< TPoints >::Type transIt =  begin( tPoints );
			// traverse all dim! permutations
		for( SizeType i = 0; i < fac; ++i )
		{
			assignValue( transIt, typename Value< TPoints >::Type() );
			reserve( value( transIt ), length( source ) + 2 );
			
				// transform points
			typename Iterator< typename Value< TPoints >::Type >::Type cPointIt = begin( value( transIt ) );
			for( SizeType pointCount = 0; pointCount < ( length( source ) + 2 ); ++pointCount )
			{
				_ChainPoint< FragType, SpecType > buffer( value( metaIt ), dim );
				appendValue( value( transIt ), buffer );
				_chainTransformCoords( buffer, value( cPointIt ), perm );
				goNext( metaIt );
				goNext( cPointIt );
			}
			std::next_permutation( begin( perm ), end( perm ) );
			metaIt = begin( metas );
			goNext( transIt );
		}
	}

		// construct dim! range trees for gSoP metric
	template< typename TTrees, typename TTPoints, typename TSize >
	void
	_build_chain_trees( TTrees & trees, 
						TTPoints & tPoints, 
						TSize dim,
						TSize facValue )
	{
		typename Iterator< TTPoints >::Type tPointIt = begin( tPoints );
		typename Iterator< TTrees >::Type treeIt = begin( trees );
		for( TSize i = 0; i < facValue; ++i )
		{
			allocate( value( treeIt ), value( treeIt ), 1 );
			new( value( treeIt ) ) typename Value< typename Value< TTrees >::Type >::Type ( value( tPointIt ), dim );
			goNext( tPointIt );
			goNext( treeIt );
		}
	}

		// delete all RMT's
	template< typename TTrees, typename TSize >
	void
	_delete_chain_trees( TTrees & trees, 
							TSize facValue )
	{
		typename Iterator< TTrees >::Type treeIt = begin( trees );
		for( TSize i = 0; i < facValue; ++i )
		{
			delete ( value( treeIt ) );
			goNext( treeIt );
		}
	}

}

#endif
