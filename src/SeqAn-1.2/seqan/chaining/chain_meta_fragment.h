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
  $Id: chain_meta_fragment.h 1448 2007-12-20 15:56:43Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_META_FRAGMENT
#define SEQAN_HEADER_META_FRAGMENT

namespace seqan{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class _MetaFragment
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Class._MetaFragment:
..summary:Basic data which associates fragments with prededing fragments and stores chain score informations
..cat:Chaining
..signature:_MetaFragment< TFragType >
..param.TFragType:Type of the fragment
*/

		// get/set the weight of the related fragment
	template< typename TFragType > inline
	typename Weight< TFragType >::Type
	weight( _MetaFragment< TFragType > & me )
	{
		return weight( *me._frag );
	}

	template< typename TFragType, typename TWeight> inline
	void
	setWeight( _MetaFragment< TFragType > & me,
				TWeight weight )
	{
		setWeight( *me._frag, weight );
	}

		// get/set the score of the chain
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	score( _MetaFragment< TFragType > & me )
	{
		return me._score;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setScore( _MetaFragment< TFragType > & me,
						TWeight score )
	{
		me._score = score;
	}

		// get/set the priority
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	priority( _MetaFragment< TFragType > & me )
	{
		return me._priority;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setPriority( _MetaFragment< TFragType > & me,
					TWeight prio )
	{
		me._priority = prio;
	}

		// get the associated fragment
	template< typename TFragType > inline
	TFragType & 
	_getFrag( _MetaFragment< TFragType > & me )
	{
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const _MetaFragment< TFragType > & me )
	{
		return *me._frag;
	}

		// get preceding fragment
	template< typename TFragType > inline
	_MetaFragment< TFragType > & 
	_getPred( _MetaFragment< TFragType > & me )
	{
		return *me._pred;
	}

	template< typename TFragType > inline 
	_MetaFragment< TFragType > & 
	_getPred( const _MetaFragment< TFragType > & me )
	{
		return *me._pred;
	}

		// set preceding fragment
	template< typename TFragType > inline
	void
	_setPred( _MetaFragment< TFragType > & me, 
				_MetaFragment< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void
	_setPred( const _MetaFragment< TFragType > & me, 
				_MetaFragment< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void 
	dump( _MetaFragment< TFragType > & me )
	{
		if( me._frag )
			dump( *me._frag );
		std::cout << me._priority << " " << me._score << std::endl;
	}


	template< typename TFragType >
	struct _MetaFragment
	{
		TFragType * _frag;
			// preceding element in a chain
		typename Weight< TFragType >::Type _priority;
		typename Weight< TFragType >::Type _score;
		_MetaFragment< TFragType > * _pred;

		_MetaFragment()		
			: _frag( NULL )
			, _priority( infimumValue< typename Weight< TFragType >::Type >() )
			, _score( infimumValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		_MetaFragment( TFragType & frag )
			: _frag( &frag )
			, _priority( infimumValue< typename Weight< TFragType >::Type >() )
			, _score( infimumValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		_MetaFragment( const _MetaFragment & old )
			: _frag( old._frag)
			, _priority( old._priority )
			, _score( old._score )
			, _pred( old._pred )
		{}

		_MetaFragment & operator=( const _MetaFragment & old )
		{
			if ( this == &old ) 
				return *this;
			_frag = old._frag;
			_priority = old._priority;
			_score = old._score;
			_pred = old._pred;
			return *this;
		}

	};


}

#endif // SEQAN_HEADER_META_FRAGMENT
