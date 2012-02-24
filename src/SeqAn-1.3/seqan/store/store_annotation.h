// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ANNOTATION_H
#define SEQAN_HEADER_STORE_ANNOTATION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Annotation Store
//////////////////////////////////////////////////////////////////////////////

template <typename TPos_, typename TSpec = void>
struct AnnotationStoreElement
{
	typedef typename Id<AnnotationStoreElement>::Type		TId;
	typedef TPos_											TPos;
	typedef StringSet<CharString, Owner< ConcatDirect<> > >	TValues;

	static const TId  INVALID_ID;
	static const TPos INVALID_POS;
	
	TId					parentId;
	TId					contigId;
	TId					countId;
	TId					typeId;			// gene, intron, ...

	TPos				beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos				endPos;			// end position of ..., for reverse aligned reads holds end < begin
	
	TId					lastChildId;	// generated back links to child
	TId					nextSiblingId;	// and sibling
	
	TValues				values;			// stores values for each keyId of (key,value) pairs

	AnnotationStoreElement() : 
		parentId(INVALID_ID), contigId(INVALID_ID), countId(INVALID_ID), typeId(INVALID_ID), 
		beginPos(INVALID_POS), endPos(INVALID_POS),
		lastChildId(INVALID_ID), nextSiblingId(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TSpec>
const typename Id<AnnotationStoreElement<TPos, TSpec> >::Type
AnnotationStoreElement<TPos, TSpec>::INVALID_ID = MaxValue<typename Id<AnnotationStoreElement<TPos, TSpec> >::Type>::VALUE;

template <typename TPos, typename TSpec>
const TPos
AnnotationStoreElement<TPos, TSpec>::INVALID_POS = MaxValue<TPos>::VALUE;

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct AnnotationTree {};

template <typename TFragmentStore, typename TSpec>
class Iter<TFragmentStore, AnnotationTree<TSpec> >
{
public:
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TFragmentStore *store;
	TId id;
	
	Iter():
		store(NULL),
		id(TAnnotation::INVALID_ID) {}

	Iter(TFragmentStore &_store):
		store(&_store),
		id(0) {}

	Iter(TFragmentStore &_store, MinimalCtor):
		store(&_store),
		id(TAnnotation::INVALID_ID) {}

	inline Iter const &
	operator = (Iter const &_origin)
	{
		store = &container(_origin);
		id = _origin.id;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Iterator< TFragmentStore, AnnotationTree<TSpec> > {
	typedef Iter< TFragmentStore, AnnotationTree<TSpec> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Value< Iter< TFragmentStore, AnnotationTree<TSpec> > >:
	VertexDescriptor<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Size< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Size<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Position< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Position<TFragmentStore> {};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return it.id;
}

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type const & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return it.id;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return *it.store;
}

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return *it.store;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename Reference<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return it.store->annotationStore[it.id];
}

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoName(*it.store, it.id);
}

template <typename TFragmentStore, typename TSpec, typename TName>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
setName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TName & name)
{
	if (length(it.store->annotationNameStore) <= it.id)
		resize(it.store->annotationNameStore, it.id + 1);
	it.store->annotationNameStore[it.id] = name;
}

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getParentName(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TId parentId = it.store->annotationStore[it.id].parentId;
	if (parentId == TAnnotation::INVALID_ID) parentId = it.id;
	return getAnnoName(*it.store, parentId);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationTypeStore>::Type
getType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoType(*it.store, it.id);
}

template <typename TFragmentStore, typename TSpec, typename TTypeName>
inline void
setType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TTypeName & typeName)
{
	_storeAppendType(*it.store, getAnnotation(it).typeId, typeName);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline CharString
getUniqueName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoUniqueName(*it.store, it.id);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void 
clearValues(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	clear(getAnnotation(it).values);
}

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline void 
assignValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > & it,
	TKey const & key,
	TValue const & value)
{
	annotationAssignValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline bool 
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key,
	TValue & value)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey>
inline CharString
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it.id = 0;
}

template <typename TFragmentStore, typename TSpec>
inline void
goEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it.id = TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline void
clear(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it.id = TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline bool
atBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	return it.id == 0;
}

template <typename TFragmentStore, typename TSpec>
inline bool
atEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	return it.id == TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goNext(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goDown(it) && !goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void
goNextRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void
goNextUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it.id = 0;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId parentId = getAnnotation(it).parentId;
	if (parentId != TAnnotation::INVALID_ID)
	{
		it.id = parentId;
		return true;
	}
	return false;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goDown(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId lastChildId = getAnnotation(it).lastChildId;
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		it.id = it.store->annotationStore[lastChildId].nextSiblingId;
		return true;
	}
	return false;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TAnnotation &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		if (it.id != lastChildId)
		{
			it.id = nextSiblingId;
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeUp(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goUp(tmp);
	return tmp;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeDown(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goDown(tmp);
	return tmp;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeRight(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goRight(tmp);
	return tmp;
}

//////////////////////////////////////////////////////////////////////////////

// insert a new id into a cyclic list and returns new last child id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListFrontInsert(TAnnotationStore & annotationStore, TId newId, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TId nextId, newLastId;
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		// get last node in the cycle
		TAnnotation &lastChild = annotationStore[lastChildId];
		// last child points to first child
		nextId = lastChild.nextSiblingId;
		// insert new node between last and first
		lastChild.nextSiblingId = newId;
		// last child remains the same
		newLastId = lastChildId;
	} else
		// cyclic list was empty
		newLastId = nextId = newId;
	
	// link new node to former first node
	annotationStore[newId].nextSiblingId = nextId;
	
	return newLastId;
}

// delete an id from a cyclic list and returns new last child id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListSearchPrev(TAnnotationStore & annotationStore, TId id, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	if (lastChildId == TAnnotation::INVALID_ID)
		return TAnnotation::INVALID_ID;
	
	TId prevId, i = lastChildId;
	do {
		prevId = i;
		i = annotationStore[i].nextSiblingId;
		if (i == id) break;
	} while (i != lastChildId);

	if (i == id)
		return prevId;
	else
		return TAnnotation::INVALID_ID;
}

// delete an id from a cyclic list and returns new last child id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListRemove(TAnnotationStore & annotationStore, TId id, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TId prevId = _cyclicListSearchPrev(annotationStore, id, lastChildId);
	
	if (prevId != TAnnotation::INVALID_ID)
	{
		annotationStore[prevId].nextSiblingId = annotationStore[id].nextSiblingId;
		
		if (id == lastChildId)
		{
			if (prevId != id)
				return prevId;
			else
				return TAnnotation::INVALID_ID;
		} else
			return lastChildId;
	}
	return lastChildId;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createLeftChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	parentAnno.lastChildId = _cyclicListFrontInsert(it.store->annotationStore, childId, parentAnno.lastChildId);
	childAnno.parentId = it.id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter.id = childId;
	return childIter;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createRightChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	_cyclicListFrontInsert(it.store->annotationStore, childId, parentAnno.lastChildId);
	parentAnno.lastChildId = childId;
	childAnno.parentId = it.id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter.id = childId;
	return childIter;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createSibling(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &anno = getAnnotation(it);

	TId siblingId = length(it.store->annotationStore) - 1;

	TAnnotation &parentAnno = it.store->annotationStore[anno.parentId];
	if (parentAnno.lastChildId == it.id)
		parentAnno.lastChildId = siblingId;

	TAnnotation &siblingAnno = it.store->annotationStore[siblingId];
	siblingAnno.nextSiblingId = anno.nextSiblingId;
	siblingAnno.parentId = anno.parentId;
	siblingAnno.lastChildId = TAnnotation::INVALID_ID;
	anno.nextSiblingId = siblingId;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > siblingIter(it);
	siblingIter.id = siblingId;
	return siblingIter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline bool
isRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it.id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it.id].parentId == TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline bool
isLeaf(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it.id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it.id].lastChildId == TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline bool
isLastChild(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TAnnotation &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		return it.id == lastChildId;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAnnotationStore>
inline void
_storeClearAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;

	TAnnoIter it = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	
	for (; it != itEnd; ++it)
	{
		(*it).lastChildId = TAnnotation::INVALID_ID;
		(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

template <typename TAnnotationStore>
inline void
_storeCreateAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename TAnnotation::TId							TId;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;
	
	TAnnoIter itBegin = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	TId id = (itEnd - itBegin) - 1;
	TAnnoIter it = itBegin + id;
	
	for (; itBegin <= it; --it, --id)
	{
		if ((*it).parentId != TAnnotation::INVALID_ID)
		{
			TAnnoIter parent = itBegin + (*it).parentId;
			if ((*parent).lastChildId == TAnnotation::INVALID_ID)
			{
				(*parent).lastChildId = id;
				(*it).nextSiblingId = id;
			}

			if ((*it).nextSiblingId == TAnnotation::INVALID_ID)
			{
				TAnnoIter lastChild = itBegin + (*parent).lastChildId;
				(*it).nextSiblingId = (*lastChild).nextSiblingId;
				(*lastChild).nextSiblingId = id;
			}
		}
		else 
			(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
