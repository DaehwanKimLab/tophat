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
  $Id: fragment.h 3038 2008-11-12 21:07:25Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FRAGMENT_H
#define SEQAN_HEADER_FRAGMENT_H


namespace seqan
{



/**
.Spec.MultiSeed:
..summary:Data structure which represents a seed of multiple sequences.
..cat:Chaining
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, MultiSeed>
..param.TPosition:Type of the class which represents the limits (multidimensional point) of a seed.
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class Seed implementation
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct MultiSeed;

template <typename TBorder>
struct Seed <TBorder, MultiSeed>
{
	typedef MultiSeed TSpec;
	typename Key< Seed< TBorder, TSpec > >::Type * _left;
	typename Key< Seed< TBorder, TSpec > >::Type * _right;
	typename Size< Seed< TBorder, TSpec > >::Type _dim;
	typename Weight< Seed< TBorder, TSpec > >::Type _weight;

/*
#ifdef _SEQAN_CHAIN_DEBUG // some debugging variables to identify fragments while debugging
	char _id;

	friend inline 
	char 
	_getID( Seed & me )
	{
		return me._id;
	}

	friend inline 
	char 
	_getID( const Seed & me )
	{
		return me._id;
	}

#endif // _SEQAN_CHAIN_DEBUG


#ifdef _SEQAN_CHAIN_DEBUG
	static int _frag_counter;
#endif // _SEQAN_CHAIN_DEBUG
*/

	Seed( )
		: _left( NULL )
		, _right( NULL )
		, _dim( 0 )
		, _weight( 0 )
	{}

	Seed( typename Size< Seed< TBorder, TSpec > >::Type dim )
		: _left( NULL )
		, _right( NULL )
		, _dim( dim )
		, _weight( 0 )
	{
		allocate( *this, _left, _dim );
		allocate( *this, _right, _dim );
	}

	Seed( TBorder * left,
				TBorder * right,
				typename Size< Seed< TBorder, TSpec > >::Type dim )
		: _left( NULL )
		, _right( NULL )
		, _dim( dim )
		, _weight( 0 )
	{
		allocate( *this, _left, _dim );
		allocate( *this, _right, _dim );
		_weight = 0;
		for( typename Size< Seed< TBorder, TSpec > >::Type i = 0; i < dim; ++i )
		{
			_left[ i ] = left[ i ];
			_right[ i ] = right[ i ];
		}
	}

	Seed( TBorder * left,
				TBorder * right,
				typename Size< Seed< TBorder, TSpec > >::Type dim,
				typename Weight< Seed< TBorder, TSpec > >::Type weight )
		: _left( NULL )
		, _right( NULL )
		, _dim( dim )
		, _weight( weight )
	{
		allocate( *this, _left, _dim );
		allocate( *this, _right, _dim );
		for( typename Size< Seed< TBorder, TSpec > >::Type i = 0; i < dim; ++i )
		{
			_left[ i ] = left[ i ];
			_right[ i ] = right[ i ];
		}
	}

	
	Seed & operator=( const Seed & old )
	{
		if( this == &old) 
			return *this;
		if( _left )
			deallocate( *this, _left, _dim );
		if( _right )
			deallocate( *this, _right, _dim );
		_dim = old._dim;
		_weight =  old._weight;
		allocate( *this, _left, _dim );
		allocate( *this, _right, _dim );
		for( typename Size< Seed< TBorder, TSpec > >::Type i = 0; i < _dim; ++i )
		{
			_left[ i ] = old._left[ i ];
			_right[ i ] = old._right[ i ];
		}			
	#ifdef _SEQAN_CHAIN_DEBUG
		_id = old._id;
	#endif
		return *this;
	}

	Seed( const Seed & old )
	{
		_dim = old._dim;
		_weight =  old._weight;
		allocate( *this, _left, _dim );
		allocate( *this, _right, _dim );
		for( typename Size< Seed< TBorder, TSpec > >::Type i = 0; i < _dim; ++i )
		{
			_left[ i ] = old._left[ i ];
			_right[ i ] = old._right[ i ];
		}			
	#ifdef _SEQAN_CHAIN_DEBUG
		_id = old._id;
	#endif
	}

	~Seed()
	{
		deallocate( *this, _left, _dim );
		deallocate( *this, _right, _dim );
		_dim = 0;
		_weight = 0;
	}

	friend inline
	void 
	dump( Seed & me )
	{
		if( me._left != NULL )
		{
			std::cout << "[ ";
			typename Size< Seed< TBorder, TSpec > >::Type dim = 0;
			std::cout << leftPosition( me, dim );
			++dim;
			while( dim != me._dim )
			{
				std::cout << " , " << leftPosition( me, dim );
				++dim;
			}
			std::cout << " ] * [ ";
			dim = 0;
			std::cout << rightPosition( me, dim );
			++dim;
			while( dim != me._dim )
			{
				std::cout << " , "<< rightPosition( me, dim );
				++dim;
			}
			std::cout << " ] "  << weight( me ) << std::endl;
		}
	}

	#ifdef _SEQAN_CHAIN_DEBUG
	friend inline
	void 
	dump( Seed & me, std::ostream & os )
	{
		if( me._left != NULL )
		{
			typename Size< Seed< TBorder, TSpec > >::Type dim = 0;
			os << "# " << weight( me )/10 << std::endl;
			while( dim != me._dim )
			{
				os << "[" << leftPosition( me, dim ) << "," << rightPosition( me, dim ) << "]";
				++dim;
			}
			os << std::endl;
		}
	}
	#endif // _SEQAN_CHAIN_DEBUG

};

template< typename TBorder > inline
typename Size< Seed< TBorder, MultiSeed > >::Type
dimension( Seed< TBorder, MultiSeed > & me )
{
	return me._dim;
}

template< typename TBorder, typename TSize > inline 
TBorder 
leftPosition( Seed< TBorder, MultiSeed >  & me, 
				TSize dim )
{
	return me._left[ dim ];
}

template< typename TBorder, typename TSize > inline 
TBorder 
leftPosition( const Seed< TBorder, MultiSeed >  & me, 
				TSize dim )
{
	return me._left[ dim ];
}

template< typename TBorder, typename TSize > inline 
TBorder 
rightPosition( Seed< TBorder, MultiSeed >  & me, 
				TSize dim )
{
	return me._right[dim];
}

template< typename TBorder, typename TSize > inline 
TBorder 
rightPosition( const Seed< TBorder, MultiSeed >  & me,
				TSize dim )
{
SEQAN_ASSERT( me._left->size() == me._right->size() )
	return me._right[dim];
}



template< typename TBorder, typename TSpec > inline 
typename Weight< Seed< TBorder, TSpec > >::Type 
weight( Seed< TBorder, TSpec > & me )
{
	return me._weight;
}

template< typename TBorder, typename TSpec > inline 
typename Weight< Seed< TBorder, TSpec > >::Type 
weight( const Seed< TBorder, TSpec > & me )
{
	return me._weight;
}


template< typename TBorder, typename TSpec, typename TWeight > 
inline void 
setWeight( Seed< TBorder, TSpec > & me,
			TWeight weight )
{
	me._weight = weight;
}

template< typename TBorder, typename TSpec, typename TSize, typename TPosition > 
inline void 
_setLeftPosition( Seed< TBorder, TSpec > & me,
					TSize dim,
					TPosition value )
{
	SEQAN_CHECK2( dim >= 0 && dim < me._dim, "Dimension index out of bounds" );
	me._left[ dim ] = value;
}
template< typename TPosition, typename TSize, typename TPosition2> 
inline void
setLeftPosition(Seed<TPosition, MultiSeed>  & me, 
				TSize dim,
				TPosition2 new_pos)
{
	_setLeftPosition(me, dim, new_pos);
}

template< typename TBorder, typename TSpec, typename TSize, typename TPosition > 
inline void 
_setRightPosition( Seed< TBorder, TSpec > & me,
					TSize dim, 
					TPosition value )
{
	SEQAN_CHECK2( dim >= 0 && dim < me._dim, "Dimension index out of bounds" );
	me._right[ dim ] = value;
}
template< typename TPosition, typename TSize, typename TPosition2> 
inline void 
setRightPosition(Seed<TPosition, MultiSeed>  & me, 
				 TSize dim,
				 TPosition2 new_pos)
{
	_setRightPosition(me, dim, new_pos);
}




//set a fragment to "top" position, that is the starting fragment of a global chain
template <typename TFragment> 
inline void 
makeBeginFragment(TFragment & me)
{
	unsigned int dim = dimension(me);
	for (unsigned int i = 0; i < dim; ++i)
	{
		_setLeftPosition(me, i, 0);		//probably dont need left positions for top fragment
		_setRightPosition(me, i, 0);
	}
}

//set a fragment to "bottom" position in respect to a fragment set
template <typename TFragment, typename TFragments> 
inline void 
makeEndFragment(TFragment & me,
				TFragments & fragments)
{
	typedef typename Size<TFragment>::Type TSize;

	unsigned int dim = dimension(me);
	String<TSize> maxes;
	resize(maxes, dim);
	arrayFill(begin(maxes, Standard()), end(maxes, Standard()), 0);

	typedef typename Iterator<TFragments, Standard>::Type TIterator;
	for (TIterator it = begin(fragments, Standard()); it < end(fragments, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			TSize rpos = rightPosition(*it, i);
			if (rpos > maxes[i])
			{
				maxes[i] = rpos;
			}
		}
	}
	for (unsigned int i = 0; i < dim; ++i)
	{
		_setLeftPosition(me, i, maxes[i] );	
		_setRightPosition(me, i, maxes[i] );		//probably dont need right positions for bottom fragments
	}
}




} // namespace seqan

#endif // SEQAN_HEADER_FRAGMENT_H

