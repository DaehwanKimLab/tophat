 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2008

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: map_adapter_stl.h 1988 2008-05-14 12:02:05Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MAP_ADAPTER_STL_H
#define SEQAN_HEADER_MAP_ADAPTER_STL_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Value< ::std::map<TKey, TCargo, TCompare, TAlloc> > 
{
    typedef Pair<TKey,TCargo> Type;
};

template <typename TKey, typename TCompare, typename TAlloc>
struct Value< ::std::set<TKey, TCompare, TAlloc> > 
{
    typedef bool Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Key< ::std::map<TKey,TCargo, TCompare, TAlloc>  >
{
    typedef TKey Type;
};

template <typename TKey, typename TCompare, typename TAlloc>
struct Key< ::std::set<TKey, TCompare, TAlloc>  >
{
    typedef TKey Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Cargo< ::std::map<TKey,TCargo, TCompare, TAlloc> >
{
    typedef TCargo Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCompare, typename TAlloc>
struct Size< ::std::set<TKey, TCompare, TAlloc> >
{
    //typedef ::std::set<TKey>::size_type Type;
    typedef unsigned Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Size< ::std::map<TKey,TCargo, TCompare, TAlloc> >
{
    //typedef ::std::map<TKey,TCargo>::size_type Type;
    typedef unsigned Type;
};

//////////////////////////////////////////////////////////////////////////////

// traits for stl usage
/*
template <typename TKey, typename TCompare, typename TAlloc>
struct AllocatorType< ::std::set<TKey, TCompare, TAlloc> >
{
    typedef TAlloc Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct AllocatorType< ::std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef TAlloc Type;
};

template <typename T>
struct _STLComparator;

template <typename TKey, typename TCompare, typename TAlloc>
struct _STLComparator< ::std::set<TKey, TCompare, TAlloc> >
{
    typedef TCompare Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct _STLComparator< ::std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef TCompare Type;
};

*/


template <typename T>
struct _STLIterator
{
	typedef int Type;
};
template <typename TKey, typename TCompare, typename TAlloc>
struct _STLIterator< ::std::set<TKey, TCompare, TAlloc> >
{
	typedef typename ::std::set<TKey, TCompare, TAlloc>::iterator Type;
};
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct _STLIterator< ::std::map<TKey, TCargo, TCompare, TAlloc> >
{
	typedef typename ::std::map<TKey, TCargo, TCompare, TAlloc>::iterator Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCompare, typename TAlloc>
inline void
assign(::std::set<TKey, TCompare, TAlloc> & target,
	   ::std::set<TKey, TCompare, TAlloc> const & source)
{
    target = source;
}

template <typename TKey,typename TCargo, typename TCompare, typename TAlloc>
inline void
assign(::std::map<TKey,TCargo, TCompare, TAlloc> & target,
	   ::std::map<TKey,TCargo, TCompare, TAlloc> const & source)
{
    target = source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(::std::set<TValue,TCompare,TAlloc> & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

template <typename TValue, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(::std::set<TValue,TCompare,TAlloc> const & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}


template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(::std::map<TKey, TCargo, TCompare, TAlloc> & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(::std::map<TKey, TCargo, TCompare, TAlloc> const & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc>
inline typename Size< ::std::set<TValue, TCompare, TAlloc> >::Type
length(::std::set<TValue, TCompare, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return me.size();
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Size< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type
length(::std::map<TKey,TCargo, TCompare, TAlloc> const & me)
{
SEQAN_CHECKPOINT
	return me.size();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TValue2>
inline void
insert(::std::set<TValue, TCompare, TAlloc> & me,TValue2 const & _value)
{
SEQAN_CHECKPOINT
    me.insert(_value);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2, typename TCargo2>
inline void
insert(::std::map<TKey,TCargo, TCompare, TAlloc> & me, 
	   TKey2 const & _key, 
	   TCargo2 const & _cargo)
{
SEQAN_CHECKPOINT

	me[_key] = _cargo;
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline void
insert(::std::map<TKey,TCargo, TCompare, TAlloc> & me,TKey2 const & _key)
{
SEQAN_CHECKPOINT

	insert(me, _key, TCargo());
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2, typename TCargo2, typename TSpec>
inline void
insert(::std::map<TKey,TCargo, TCompare, TAlloc> & me, Pair<TKey2,TCargo2,TSpec> const & _value)
{
SEQAN_CHECKPOINT
	insert(me, _value.i1, _value.i2);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc>
inline void
clear(::std::set<TValue, TCompare, TAlloc> & me)
{
SEQAN_CHECKPOINT
    me.clear();
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline void
clear(::std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
SEQAN_CHECKPOINT
    me.clear();
}

//////////////////////////////////////////////////////////////////////////////

//template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
//inline typename Value< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
//value(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
//	  TKey2 const & _key)
//{    
//    return me[_key];
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline typename Cargo< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
cargo(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
	  TKey2 const & _key)
{
SEQAN_CHECKPOINT
	return me[_key];
}

//////////////////////////////////////////////////////////////////////////////

struct STLSetIterator;
struct STLMapIterator;

//////////////////////////////////////////////////////////////////////////////
//                            MapIterator                                   //
//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
struct Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc> , TIteratorSpec >
{
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TSTLMap;
	typedef Iter<TSTLMap, STLMapIterator> Type;
};


template <typename TSTLMap>
class Iter< TSTLMap, STLMapIterator>
{
public:    
    //typedef typename ::std::map<typename Key<TSTLMap>::Type,
    //                            typename Cargo<TSTLMap>::Type,
    //                            typename _STLComparator<TSTLMap>::Type,
    //                            typename AllocatorType<TSTLMap>::Type >::iterator THostIter;

	typedef typename _STLIterator<TSTLMap>::Type THostIter;
	THostIter _iter;
    Holder<TSTLMap> host_map_holder;

	Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter(Iter const & other)
		: _iter(other._iter)
	{
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
	}

	Iter(TSTLMap & map)
		: _iter(map.begin())
    {
SEQAN_CHECKPOINT
        setValue(host_map_holder,map);
	}

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & 
	operator = (Iter const & other)
	{
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
	}
	operator bool () const
	{
SEQAN_CHECKPOINT
        return (_iter != value(host_map_holder).end());
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
operator == (Iter<TSTLMap, STLMapIterator> const & left,
			 Iter<TSTLMap, STLMapIterator> const & right)
{
SEQAN_CHECKPOINT
	return left._iter == right._iter;
}

template <typename TSTLMap>
inline bool
operator != (Iter<TSTLMap, STLMapIterator> const & left,
			 Iter<TSTLMap, STLMapIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc>, TIteratorSpec>::Type
begin(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
	  TIteratorSpec)
{
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TSTLMap;
    typedef typename Iterator<TSTLMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type
begin(::std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TSTLMap;
    typedef typename Iterator<TSTLMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type
end(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
	  TIteratorSpec)
{
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TSTLMap;
	typedef typename Iterator<TSTLMap, TIteratorSpec>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type
end(::std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TSTLMap;
	typedef typename Iterator<TSTLMap>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
atEnd(Iter<TSTLMap, STLMapIterator> & it)
{
SEQAN_CHECKPOINT
	return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline void
goNext(Iter<TSTLMap, STLMapIterator> & it)
{
SEQAN_CHECKPOINT
	it._iter++;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Value<TSTLMap>::Type &
value(Iter<TSTLMap, STLMapIterator> & it)
{
SEQAN_CHECKPOINT
    return it._iter->second;
}
template <typename TSTLMap>
inline typename Value<TSTLMap>::Type &
value(Iter<TSTLMap, STLMapIterator> const & it)
{
SEQAN_CHECKPOINT
	return it._iter->second;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLMapIterator> & it)
{
SEQAN_CHECKPOINT
	return it._iter->first;
}
template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLMapIterator> const & it)
{
SEQAN_CHECKPOINT
	return it._iter->first;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Cargo<TSTLMap>::Type &
cargo(Iter<TSTLMap, STLMapIterator> & it)
{
SEQAN_CHECKPOINT
	return it._iter->second;
}
template <typename TSTLMap>
inline typename Cargo<TSTLMap>::Type &
cargo(Iter<TSTLMap, STLMapIterator> const & it)
{
SEQAN_CHECKPOINT
	return it._iter->second;
}

////////////////////////////////////////////////////////////////////////////////
//                             SetIterator                                    //
////////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
struct Iterator< ::std::set<TValue, TCompare, TAlloc> , TIteratorSpec >
{
    typedef ::std::set<TValue, TCompare, TAlloc> TSTLMap;
	typedef Iter<TSTLMap, STLSetIterator> Type;
};

template <typename TSTLMap>
class Iter< TSTLMap, STLSetIterator>
{
public:    
    //typedef typename ::std::set<typename Key<TSTLMap>::Type,
    //                            typename _STLComparator<TSTLMap>::Type,
    //                            typename AllocatorType<TSTLMap>::Type >::iterator THostIter;
	typedef typename _STLIterator<TSTLMap>::Type THostIter;
	THostIter _iter;
    Holder<TSTLMap> host_map_holder;

	Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter(Iter const & other)
		: _iter(other._iter)
	{
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
	}

	Iter(TSTLMap & map)
		: _iter(map.begin())
    {
SEQAN_CHECKPOINT
        setValue(host_map_holder,map);
	}

	~Iter()
	{
SEQAN_CHECKPOINT
	}

	Iter const & 
	operator = (Iter const & other)
	{
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
	}
	operator bool () const
	{
SEQAN_CHECKPOINT
        return (_iter != value(host_map_holder).end());
	}

};

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
operator == (Iter<TSTLMap, STLSetIterator> const & left,
			 Iter<TSTLMap, STLSetIterator> const & right)
{
SEQAN_CHECKPOINT
	return left._iter == right._iter;
}

template <typename TSTLMap>
inline bool
operator != (Iter<TSTLMap, STLSetIterator> const & left,
			 Iter<TSTLMap, STLSetIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< ::std::set<TValue,TCompare,TAlloc>, TIteratorSpec>::Type
begin(::std::set<TValue,TCompare,TAlloc> & me,
	  TIteratorSpec)
{
    typedef ::std::set<TValue,TCompare,TAlloc> TSTLMap;
    typedef typename Iterator<TSTLMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TValue, typename TCompare, typename TAlloc>
inline typename Iterator< ::std::set<TValue,TCompare,TAlloc> >::Type
begin(::std::set<TValue,TCompare,TAlloc> & me)
{
    typedef ::std::set<TValue,TCompare,TAlloc> TSTLMap;
    typedef typename Iterator<TSTLMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< ::std::set<TValue, TCompare, TAlloc>, TIteratorSpec>::Type
end(::std::set<TValue, TCompare, TAlloc> & me,
	  TIteratorSpec)
{
    typedef ::std::set<TValue,TCompare,TAlloc> TSTLMap;
	typedef typename Iterator<TSTLMap, TIteratorSpec>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TValue, typename TCompare, typename TAlloc>
inline typename Iterator< ::std::set<TValue,TCompare,TAlloc> >::Type
end(::std::set<TValue,TCompare,TAlloc> & me)
{
    typedef ::std::set<TValue,TCompare,TAlloc> TSTLMap;
	typedef typename Iterator<TSTLMap>::Type TIterator;
	TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline bool
atEnd(Iter<TSTLMap, STLSetIterator> & it)
{
SEQAN_CHECKPOINT
	return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline void
goNext(Iter<TSTLMap, STLSetIterator> & it)
{
SEQAN_CHECKPOINT
    it._iter++;
}

//////////////////////////////////////////////////////////////////////////////
//
//template <typename TSTLMap>
//inline typename Value<TSTLMap>::Type &
//value(Iter<TSTLMap, STLSetIterator> & it)
//{
//    return hasKey(*it);
//}
//template <typename TSTLMap>
//inline typename Value<TSTLMap>::Type &
//value(Iter<TSTLMap, STLSetIterator> const & it)
//{
//	return hasKey(*it);
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLSetIterator> & it)
{
SEQAN_CHECKPOINT
    return (*it._iter);
}
template <typename TSTLMap>
inline typename Key<TSTLMap>::Type const &
key(Iter<TSTLMap, STLSetIterator> const & it)
{
SEQAN_CHECKPOINT
	return (*it._iter);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc,typename TFind>
inline typename Iterator< ::std::set<TValue, TCompare, TAlloc> >::Type
find(::std::set<TValue, TCompare, TAlloc> & me,
	 TFind const & _find)
{
SEQAN_CHECKPOINT
    typedef ::std::set<TValue, TCompare, TAlloc> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;  

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find); 
    }
    return _iter;
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TFind>
inline typename Iterator< ::std::map<TKey,TCargo, TCompare, TAlloc> >::Type
find(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
	 TFind const & _find)
{
SEQAN_CHECKPOINT
    typedef ::std::map<TKey,TCargo, TCompare, TAlloc> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;  

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find); 
    }
    return _iter;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TMap2>
inline void
erase(::std::set<TValue, TCompare, TAlloc> & me,
	  Iter<TMap2, STLSetIterator> const & it)
{
SEQAN_CHECKPOINT
    me.erase(it._iter);
}

template <typename TKey, typename TCargo ,typename TMap2>
inline void
erase(::std::map<TKey,TCargo> & me,
	  Iter<TMap2, STLMapIterator> const & it)
{
SEQAN_CHECKPOINT
    me.erase(it._iter);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TToRemove>
inline void
erase(::std::set<TValue, TCompare, TAlloc> & me,
	  TToRemove const & to_remove)
{
SEQAN_CHECKPOINT
    me.erase(to_remove);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TToRemove>
inline void
erase(::std::map<TKey,TCargo, TCompare, TAlloc> & me,
	  TToRemove const & to_remove)
{
SEQAN_CHECKPOINT
    me.erase(to_remove);
}

}
#endif // #ifndef SEQAN_HEADER
