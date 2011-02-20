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
  $Id: skip_list_type.h 3348 2009-02-03 17:14:40Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SKIPLIST_TYPE_H
#define SEQAN_HEADER_SKIPLIST_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{

//SEQAN_NO_DDDOC: do not generate documentation for this file


/**
.Metafunction.Key:
..summary:Type of the theKeyattribute of an object. 
..signature:Key<T>::Type
..param.T:Type for which the key type is determined.
..returns.param.Type:Key type of $T$.
..remarks.text:The theKeytype of an object is used for sorting and searching.
*/

/* moved to basic_h, see #6
	template< typename T >
	struct Key
	{
		typedef T Type;
	};
*/
		// default for Pair
	template < typename _T1, typename _T2, typename Compressed >
    struct Pair;

	// auxiliary functions for Objects that are plugged into the 
		// SkipList
		// Specialized for std::pair
	template< typename TKey, typename TVal > inline
	TKey & key( Pair<TKey, TVal> & p )
	{
		return p.i1;
	}

	template< typename TKey, typename TVal > inline
	TVal & getValue( Pair<TKey, TVal> & p )
	{
		return p.i2;
	}

	template< typename TKey2, typename TVal >
	void setKey( Pair<TKey2, TVal> & p, TKey2 theKey)
	{
		p.i1 = theKey;
	}

	template< typename TKey, typename TVal >
	struct Value< Pair< TKey, TVal> >
	{
		typedef TVal Type;
	};

	template< typename TKey, typename TVal >
	struct Key< Pair< TKey, TVal > >
	{
		typedef TKey Type;
	};

		// specialization for std::pair
	template< typename TKey, typename TVal > inline
	TKey key( std::pair< TKey, TVal > & p )
	{
		return p.first;
	}

	template< typename TKey, typename TVal >
	void setKey( std::pair<TKey, TVal> & p, TKey theKey )
	{
		p.first = theKey;
	}

	template< typename TKey, typename TVal >
	struct Value< std::pair< TKey, TVal> >
	{
		typedef TVal Type;
	};

	template< typename TKey, typename TVal >
	struct Key< std::pair< TKey, TVal > >
	{
		typedef TKey Type;
	};

	
	/**
	.Metafunction.Cargo:
	..summary:An additional cargo member of an object. 
	..signature:Cargo<T>::Type
	..param.T:Type for which the cargo type is determined.
	..returns.param.Type:Cargo type of $T$.
	..remarks.text:The cargo type is used for additional cargo information.
	*/

	struct _Empty
	{};
/*
	template< typename T >
	struct Cargo
	{
		typedef _Empty Type;
	};
*/

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

struct SkipListDynamic;

struct SkipListStatic;

struct Complete;

struct Deferred;

//////////////////////////////////////////////////////////////////////////////
// Forward declarations
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring = Complete >
struct SkipElement;

template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring  = Complete >
struct SkipBaseElement;

template< typename TObject, typename TModus = SkipListDynamic, typename TSpec = Default, typename TStructuring  = Complete >
struct SkipList;


//////////////////////////////////////////////////////////////////////////////
//
// METAFUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef size_t Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Position Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Position< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > * Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		GetValue Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > * >
{
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Value Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef TObject Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Key Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< TObject >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//
//		Cargo Type
//
///////////////////////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef _Empty Type;	// default: no cargo
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};


template< typename TTag, typename TCargo > inline
void
_initCargo( TTag * /*tag*/, 
		   TCargo & /*_cargo*/ )
{}


}

#endif // SEQAN_HEADER_SKIPLIST_TYPE_H
