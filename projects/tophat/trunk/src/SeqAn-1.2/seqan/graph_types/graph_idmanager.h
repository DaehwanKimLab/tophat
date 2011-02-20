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
  $Id: graph_idmanager.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IDMANAGER_H
#define SEQAN_HEADER_GRAPH_IDMANAGER_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// IdManager
//////////////////////////////////////////////////////////////////////////////

/**
.Class.IdManager:
..cat:Graph
..summary:Id manager that provides unique ids for vertices and edges.
..signature:IdManager<TIdType,TSpec>
..param.TIdType:The id type of the managed ids.
...metafunction:Metafunction.Value
...remarks:Use the Value Metafunction to get the id type managed by a given id manager.
...default:$unsigned int$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
*/
template<typename TIdType, typename TSpec>
class IdManager 
{
	public:
		String<TIdType> data_freeIds;  
		String<bool> data_in_use;   //1 = in use, 0 = not in use

//____________________________________________________________________________	
	public:
		IdManager()
		{
			SEQAN_CHECKPOINT
			clear(data_in_use);
			clear(data_freeIds);
		}

		~IdManager() 
		{
			SEQAN_CHECKPOINT
		}

		IdManager(IdManager const & _other)
		{
			SEQAN_CHECKPOINT
			data_freeIds = _other.data_freeIds;
			data_in_use = _other.data_in_use;
		}

		IdManager const& 
		operator = (IdManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_freeIds = _other.data_freeIds;
			data_in_use = _other.data_in_use;
			return *this;
		}

//____________________________________________________________________________
};
	

//////////////////////////////////////////////////////////////////////////////
// IdManager - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.IdManager

template<typename TIdType, typename TSpec> 
struct Value<IdManager<TIdType, TSpec> > 
{
	typedef TIdType Type;
};

template<typename TIdType, typename TSpec> 
struct Value<IdManager<TIdType, TSpec> const> 
{
	typedef TIdType Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.IdManager

template<typename TIdType, typename TSpec> 
struct Spec<IdManager<TIdType, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TIdType, typename TSpec> 
struct Spec<IdManager<TIdType, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.obtainId:
..cat:Graph
..summary:Obtains a new id from the id manager.
..signature:obtainId(idm)
..param.idm:The IdManager.
...type:Class.IdManager
..returns:Returns a new unique id.
..remarks:If it is a dummy id manager, i.e., IdManager<void>, the return type is (void*) 0.
..see:Function.releaseId
*/

template<typename TIdType, typename TSpec>
inline typename Value<IdManager<TIdType, TSpec> >::Type 
obtainId(IdManager<TIdType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT

	TIdType id;
	if (!empty(idm.data_freeIds)) {
		id = getValue(idm.data_freeIds, length(idm.data_freeIds) - 1);
		resize(idm.data_freeIds, length(idm.data_freeIds) - 1, Generous());
		assignValue(idm.data_in_use, id, true);
	} else {
		if (empty(idm.data_in_use)) id = 0;
		else id = (TIdType) length(idm.data_in_use);
		resize(idm.data_in_use, id + 1, Generous());
		assignValue(idm.data_in_use, id, true);
	}
	return id;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.releaseId:
..cat:Graph
..summary:Releases a given id so it can be redistributed later on.
..signature:releaseId(idm, id)
..param.idm:The IdManager.
...type:Class.IdManager
..param.id:The id that is to be released.
..returns:void
..see:Function.obtainId
*/

template<typename TIdType, typename TSpec, typename TId>
inline void 
releaseId(IdManager<TIdType, TSpec>& idm, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(idm,id) == true)
	if (id == (TId) length(idm.data_in_use) - 1) {
		resize(idm.data_in_use, length(idm.data_in_use) - 1, Generous());
	} else {
		assignValue(idm.data_in_use, id, false);
		appendValue(idm.data_freeIds, id, Generous());
	}
	if (idCount(idm)==0) {
		releaseAll(idm);
	}
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.releaseAll:
..cat:Graph
..summary:Releases all ids handled by this id manager at once.
..signature:releaseAll(idm)
..param.idm:The IdManager.
...type:Class.IdManager
..returns:void
..see:Function.releaseId
*/


template<typename TIdType, typename TSpec>
inline void 
releaseAll(IdManager<TIdType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	clear(idm.data_freeIds);
	clear(idm.data_in_use);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getIdUpperBound:
..cat:Graph
..summary:Returns the largest distributed id plus 1. That is, the return value is guaranteed to be an upper bound on all distributed ids.
..signature:getIdUpperBound(idm)
..param.idm:The IdManager.
...type:Class.IdManager
..returns:An upper bound on all distributed ids.
..see:Function.getIdLowerBound
*/

template<typename TIdType, typename TSpec>
inline typename Value<IdManager<TIdType, TSpec> >::Type 
getIdUpperBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (empty(idm.data_in_use)) ? 0 : (typename Value<IdManager<TIdType, TSpec> >::Type) length(idm.data_in_use);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getIdLowerBound:
..cat:Graph
..summary:Returns the smallest distributed id. That is, the return value is guaranteed to be the smallest id obtained so far.
..signature:getIdLowerBound(idm)
..param.idm:The IdManager.
...type:Class.IdManager
..returns:The smallest obtained id.
..see:Function.getIdUpperBound
*/

template<typename TIdType, typename TSpec>
inline typename Value<IdManager<TIdType, TSpec> >::Type 
getIdLowerBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	for(TIdType it = 0; it < length(idm.data_in_use); ++it) {
		if (getValue(idm.data_in_use, it)) return it;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.idCount:
..cat:Graph
..summary:Determines the number of ids that were obtained.
..signature:idCount(idm)
..param.idm:The IdManager.
...type:Class.IdManager
..returns:Number of ids in use.
*/

template<typename TIdType, typename TSpec>
inline typename Value<IdManager<TIdType, TSpec> >::Type 
idCount(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (length(idm.data_in_use) - length(idm.data_freeIds));
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.idInUse:
..cat:Graph
..summary:Checks whether the given id is in use or not.
..signature:idInUse(idm, id)
..param.idm:The IdManager.
...type:Class.IdManager
..param.id:The given id.
..returns:True if the id was distributed, false otherwise.
*/

template<typename TIdType, typename TSpec, typename TId>
inline bool 
idInUse(IdManager<TIdType, TSpec> const& idm, 
		TId const id)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(id < (TId) length(idm.data_in_use))
	return idm.data_in_use[id];
}


//////////////////////////////////////////////////////////////////////////////
// Dummy IdManager
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Counting IdManager:
..cat:Graph
..general:Class.IdManager
..summary:Id Manager that just counts the number of ids in use.
..signature:IdManager<void, TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
*/
template<typename TSpec>
class IdManager<void, TSpec> 
{
	public:
		typedef typename Id<IdManager>::Type TIdType;
		TIdType data_idCount;

//____________________________________________________________________________	
	public:
		IdManager() : data_idCount(0) 
		{
			SEQAN_CHECKPOINT
		}

		~IdManager() 
		{
			SEQAN_CHECKPOINT
		}

		IdManager(IdManager const & _other) : data_idCount(_other.data_idCount) 
		{
			SEQAN_CHECKPOINT
		}

		IdManager const& 
		operator = (IdManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_idCount = _other.data_idCount;
			return *this;
		}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Dummy IdManager - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.IdManager
template<typename TSpec> 
struct Value<IdManager<void, TSpec> > {
	typedef typename Size<IdManager<void, TSpec> >::Type Type;
};

template<typename TSpec> 
struct Value<IdManager<void, TSpec> const> {
	typedef typename Size<IdManager<void, TSpec> const>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IdManager<void, TSpec> >::Type 
obtainId(IdManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	++idm.data_idCount;
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TId>
inline void 
releaseId(IdManager<void, TSpec>& idm, 
		  TId const) 
{
	SEQAN_CHECKPOINT
	if (idm.data_idCount > 0) --idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline void 
releaseAll(IdManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	idm.data_idCount = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IdManager<void, TSpec> >::Type 
getIdUpperBound(IdManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	// Must be data_idCount in order to resize property maps!!!
	// Don't change to 0
	return idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Value<IdManager<void, TSpec> >::Type 
getIdLowerBound(IdManager<void, TSpec> const&)
{
	SEQAN_CHECKPOINT
	return 0;
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSpec>
inline typename Value<IdManager<void, TSpec> >::Type 
idCount(IdManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return idm.data_idCount;
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSpec, typename TId>
inline bool 
idInUse(IdManager<void, TSpec> const&, 
		TId const) 
{
	SEQAN_CHECKPOINT
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
