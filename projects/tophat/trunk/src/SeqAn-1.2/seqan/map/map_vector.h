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
  $Id: $
 ==========================================================================*/


#ifndef SEQAN_HEADER_MAP_VECTOR_H
#define SEQAN_HEADER_MAP_VECTOR_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.VectorSet:
..cat:Map
..general:Class.Map
..summary:A fast map for small key types.
..signature:Map<TValue, VectorSet<TSpec> >
..param.TValue:The type of value stored in the map.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...default:@Spec.Alloc String.Alloc@
...remarks:$TSpec$ is used as specialization for the @Class.String@ type that is used to store the values.
..remarks:The memory needed is linear to the number different values the @Metafunction.Key.key type@ of $TValue$ can get.
So do not use this map for large key types.
*/
//////////////////////////////////////////////////////////////////////////////
// 
//	VectorSet
// 
//////////////////////////////////////////////////////////////////////////////


template <typename TSpec = Alloc<> >
struct VectorSet;


//////////////////////////////////////////////////////////////////////////////
// internal set meta-functions

template <typename TCargo>
struct _VectorSetElement
{
	bool data_valid;
	TCargo data_cargo;
};
template <>
struct _VectorSetElement<Nothing>
{
	bool data_valid;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct _VectorSetElements
{
	typedef char * Type; //dummy implementation for VC++
};
template <typename TValue, typename TSpec>
struct _VectorSetElements<Map<TValue, VectorSet<TSpec> > >
{
	typedef Map<TValue, VectorSet<TSpec> > TMap;
	typedef typename Cargo<TMap>::Type TCargo;
	typedef _VectorSetElement<TCargo> TElement;
	typedef String<TElement, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////
// VectorSet class

template <typename TValue, typename TSpec>
class Map<TValue, VectorSet<TSpec> >
{
public:
	typedef typename Key<Map>::Type TKey;
	typedef typename Size<Map>::Type TSize;
	typedef typename _VectorSetElements<Map>::Type TElements;

	TElements	data_elements;
	TSize		data_length;

	Map()
	{
		resize(data_elements, (TSize) ValueSize<TKey>::VALUE); 
		clear(*this);
	}

	Map(TSize _vectorSize)
	{
		resize(data_elements, (TSize) _vectorSize);
		clear(*this);
	}
	
	template <typename TKey>
	inline typename MapValue<Map>::Type
	operator [] (TKey const & key)
	{
		return mapValue(*this, key);
	}

};


//////////////////////////////////////////////////////////////////////////////
// set meta-functions

template <typename TValue, typename TSpec>
struct Value< Map<TValue, VectorSet<TSpec> > > 
{
	typedef TValue Type;
};
template <typename TValue, typename TSpec>
struct Reference< Map<TValue, VectorSet<TSpec> > > 
{
	typedef Map<TValue, VectorSet<TSpec> > TMap;
	typedef typename Iterator<TMap>::Type TIterator;
	typedef Proxy<IteratorProxy<TIterator> > Type;
};


template <typename TValue, typename TSpec>
struct Size< Map<TValue, VectorSet<TSpec> > >:
	Size<typename _VectorSetElements< Map<TValue, VectorSet<TSpec> > >::Type> {};

template <typename TValue, typename TSpec>
struct Key< Map<TValue, VectorSet<TSpec> > > :
	Key<TValue> {};

template <typename TValue, typename TSpec>
struct Cargo< Map<TValue, VectorSet<TSpec> > > :
	Cargo<TValue> {};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size< Map<TValue, VectorSet<TSpec> > >::Type 
length(Map<TValue, VectorSet<TSpec>  > const &set) 
{
	return set.data_length;
}

//////////////////////////////////////////////////////////////////////////////
//

struct VectorSetIterator;

template <typename TMap>
class Iter<TMap,  VectorSetIterator> 
{
public:
	typedef typename Size<TMap>::Type TSize;
	typedef typename _VectorSetElements<TMap>::Type TElements;
	typedef typename Iterator<TElements, Rooted>::Type TElementsIterator;
	typedef typename Value<TElementsIterator>::Type TValue;

	TElementsIterator data_iterator;

	Iter() 
	{	}
	Iter(TElementsIterator it) 
		: data_iterator(it)
	{	
		if (!atEnd(*this) && !(*data_iterator).data_valid) 
		{
			goNext(*this);
		}
	}
	Iter(Iter const & other)
		: data_iterator(other.data_iterator)
	{	}
	~Iter()
	{	}

	Iter const & operator = (Iter const & other)
	{
		data_iterator = other.data_iterator;
		return *this;
	}

	operator bool()
	{
		return data_iterator;
	}

	TValue *
	operator -> () const
	{
		return & value(data_iterator);
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TObject, typename TSpec, typename TIteratorSpec>
struct Iterator< Map<TObject, VectorSet<TSpec> >, TIteratorSpec > 
{
	typedef Map<TObject, VectorSet<TSpec> > TMap;
	typedef Iter<TMap, VectorSetIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void 
clear(Map<TValue, VectorSet<TSpec> > & me) 
{
	typedef Map<TValue, VectorSet<TSpec> > TMap;
	typedef typename _VectorSetElements<TMap>::Type TElements;
	typedef typename Iterator<TElements>::Type TIterator;
	TIterator it_end = end(me.data_elements);
	for (TIterator it = begin(me.data_elements); it != it_end; ++it)
	{
		(*it).data_valid = false;
	}
	me.data_length = 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
struct _VectorSet_Insert
{
	template <typename TValue, typename TSpec, typename TElement>
	static inline void 
	insert_(Map<TValue, VectorSet<TSpec> > & set,
		TElement const &element) 
	{
		if (!set.data_elements[ordValue(key(element))].data_valid) 
		{
			++set.data_length;
			set.data_elements[ordValue(key(element))].data_valid = true;
		}
		set.data_elements[ordValue(key(element))].data_cargo = cargo(element);
	}
};
template <>
struct _VectorSet_Insert<Nothing>
{
	template <typename TValue, typename TSpec, typename TElement>
	static inline void 
	insert_(Map<TValue, VectorSet<TSpec> > & set,
		TElement const &element) 
	{
		if (!set.data_elements[ordValue(key(element))].data_valid) 
		{
			++set.data_length;
			set.data_elements[ordValue(key(element))].data_valid = true;
		}
	}
};

template <typename TValue, typename TSpec, typename TElement>
inline void 
insert(Map<TValue, VectorSet<TSpec> > & set,
	   TElement const &element) 
{
	typedef Map<TValue, VectorSet<TSpec> > TMap;
	typedef typename Cargo<TMap>::Type TCargo;
	_VectorSet_Insert<TCargo>::insert_(set, element);
}
template <typename TValue, typename TSpec, typename TKey, typename TCargo>
inline void 
insert(Map<TValue, VectorSet<TSpec> > &set,
	   TKey const & _key,
	   TCargo const & _cargo) 
{
	typedef Map<TValue, VectorSet<TSpec> > TMap;
	typedef typename Value<TMap>::Type TValue2;
	TValue2 new_val;
	key(new_val) = _key;
	cargo(new_val) = _cargo;
	_VectorSet_Insert<TCargo>::insert_(set, new_val);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename T>
inline void 
erase(Map<TValue, VectorSet<TSpec> > &set, 
	  T const & tokill)  //can be an element, a key, or an iterator
{
	if (set.data_elements[ordValue(key(tokill))].data_valid) 
	{
		--set.data_length;
		set.data_elements[ordValue(key(tokill))].data_valid = false;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline bool 
hasKey(Map<TValue, VectorSet<TSpec> > const &set, TKey const &key) 
{
	return set.data_elements[ordValue(key)].data_valid;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TKey2, typename TSpec>
inline typename Iterator< Map<TKey2, VectorSet<TSpec> > >::Type 
find(Map<TKey2, VectorSet<TSpec> > & set, 
	 TKey const &key) 
{
	typedef Map<TKey2, VectorSet<TSpec> > TMap;
	typedef typename Iterator<TMap>::Type TIterator;
	return TIterator(begin(set.data_elements) + ordValue(key));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TSpec, typename TKey2>
inline typename Cargo< Map<TKey, VectorSet<TSpec> > >::Type &
cargo(Map<TKey, VectorSet<TSpec> > & set, 
	  TKey2 const & _key) 
{
	if (!hasKey(set, _key))
	{
		typedef Map<TKey, VectorSet<TSpec> > TMap;
		typedef typename Value<TMap>::Type TValue;
		TValue new_value;
		key(new_value) = _key;
		insert(set, new_value);
	}
	return set.data_elements[ordValue(_key)].data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSpec>
inline typename Iterator< Map<TElement, VectorSet<TSpec> > >::Type 
begin(Map<TElement, VectorSet<TSpec> > & set) 
{
	typedef Map<TElement, VectorSet<TSpec> > TMap;
	typedef typename Iterator<TMap>::Type TIterator;
	return TIterator(begin(set.data_elements));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSpec>
inline typename Iterator< Map<TElement, VectorSet<TSpec> > >::Type 
end(Map<TElement, VectorSet<TSpec> > &set) 
{
	typedef Map<TElement, VectorSet<TSpec> > TMap;
	typedef typename Iterator<TMap>::Type TIterator;
	return TIterator(end(set.data_elements));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline bool 
operator==(Iter<TSet, VectorSetIterator> const &a, 
		   Iter<TSet, VectorSetIterator> const &b) 
{
	return a.data_iterator == b.data_iterator;
}
template <typename TSet>
inline bool 
operator!=(Iter<TSet, VectorSetIterator> const &a, 
		   Iter<TSet, VectorSetIterator> const &b) 
{
	return a.data_iterator != b.data_iterator;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline bool 
atEnd(Iter<TSet, VectorSetIterator> & a) 
{
	return atEnd(a.data_iterator);
}
template <typename TSet>
inline bool 
atEnd(Iter<TSet, VectorSetIterator> const & a) 
{
	return atEnd(a.data_iterator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Key<TSet>::Type
key(Iter<TSet, VectorSetIterator> &it) 
{
	return position(it.data_iterator);
}
template <typename TSet>
inline typename Key<TSet>::Type
key(Iter<TSet, VectorSetIterator> const &it) 
{
	return position(it.data_iterator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Iter<TSet, VectorSetIterator> &it) {
	return (*it.data_iterator).data_cargo;
}
template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Iter<TSet, VectorSetIterator> const &it) {
	return (*it.data_iterator).data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Reference<TSet>::Type 
value(Iter<TSet, VectorSetIterator> &it) 
{
	typedef typename Reference<TSet>::Type TReference;
	return TReference(it);
}
template <typename TSet>
inline typename Reference<TSet>::Type 
value(Iter<TSet, VectorSetIterator> const &it) 
{
	typedef typename Reference<TSet>::Type TReference;
	return TReference(it);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSet>
inline void
goNext(Iter<TSet, VectorSetIterator> & it) 
{
	++it.data_iterator;
	while (!atEnd(it.data_iterator))
	{
		if ((*it.data_iterator).data_valid) 
			break;
		++it.data_iterator;
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Key<TSet>::Type
key(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > & pr) 
{
	return key(iter(pr));
}
template <typename TSet>
inline typename Key<TSet>::Type
key(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > const & pr) 
{
	return key(iter(pr));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > >& pr) 
{
	return cargo(iter(pr));
}
template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > const & pr) 
{
	return cargo(iter(pr));
}


//////////////////////////////////////////////////////////////////////////////

}

#endif
