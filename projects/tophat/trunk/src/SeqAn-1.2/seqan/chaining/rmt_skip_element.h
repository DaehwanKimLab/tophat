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
  $Id: rmt_skip_element.h 954 2007-07-27 11:48:23Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_RMT_SKIP_ELEMENT_H
#define SEQAN_RMT_SKIP_ELEMENT_H

namespace seqan
{

// Modifications of the struct SkipElement for use in a SkipListStatic< True, RT< MaxTag > >

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct _RangeMaxCargo
	{

		SkipList< TObject, TModus, TSpec, TStructuring > * _assocStruct;
		TObject * _maxObj;

		_RangeMaxCargo()
			: _assocStruct( NULL )
			, _maxObj( NULL )
		{}

	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, TSpec, TStructuring > *
	_getAssoc( _RangeMaxCargo< TObject, TModus, TSpec, TStructuring > & me )
	{
		return me._assocStruct;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setAssoc( _RangeMaxCargo< TObject, TModus, TSpec, TStructuring > & me,
				SkipList< TObject, TModus, TSpec, TStructuring > * list )
	{
		me._assocStruct = list;
	}

	
		// handling for max objects

	template< typename TObject, typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( _RangeMaxCargo< TObject, SkipListStatic, TSpec, TStructuring > & me )
	{
		return me._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( _RangeMaxCargo< TObject, SkipListStatic, TSpec, TStructuring > & me,
				TObject * obj )
	{
		me._maxObj = obj;
	}

	template< typename TObject,typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return me->_cargo._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me,
					TObject * maxObj )
	{
		me->_cargo._maxObj = maxObj;
	}
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipElement< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > >
	{
		typedef _RangeMaxCargo< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > Type;
	};

		// get the score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	weight( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return weight( *me->_cargo._maxObj );
	}

		// get the chain score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	priority( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return priority( *me->_cargo._maxObj );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif // SEQAN_RMT_SKIP_ELEMENT_H
