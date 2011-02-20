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
  $Id: memoryManager_int.h 2334 2008-06-06 13:13:23Z kemena@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MEMORYMANAGER_INT_H
#define SEQAN_HEADER_MEMORYMANAGER_INT_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


/**
.Spec.MemoryManagerIndex:
..general:Class.MemoryManager
..cat:Seed Handling
..summary: Faster than the pointer version but needs data size of at least size_t.
..signature:MemoryManager<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the manager.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size has to be a power of 2, e.g., 1024.
*/

//////////////////////////////////////////////////////////////////////////////

struct _Free_Int;
typedef Tag<_Free_Int> const FreeMemoryInt;	

template<typename TValue, unsigned int SPACE>
class MemoryManager<TValue, Block<SPACE>, FreeMemoryInt> 
{
	typedef String<TValue, Array<SPACE> >				TBlock;
	typedef TBlock*										PBlock;
	typedef Allocator< SinglePool<sizeof(TBlock)> >		TAllocator;

public:

	typedef typename Iterator<TBlock, Standard>::Type	TBlockIter;
	typedef String<PBlock>								TBlockTable;

	TBlockTable		blocks;
	TBlockIter		blockFirst, blockLast;	// current block boundaries
	TBlockIter		lastValue;				// pointer to top value
	TAllocator		alloc;
	bool			change;

    size_t pEmptySpace;
	unsigned int exponent;
	
	//____________________________________________________________________________
	      
	public:
	MemoryManager():
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter()) 
	{
		SEQAN_CHECKPOINT
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		pEmptySpace = -1;
		change = false;
	}

	template<typename TSource>
	MemoryManager(TSource const& source):
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter())
	{
		SEQAN_CHECKPOINT
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		assign(*this, source);
		change = false;
	} 

	MemoryManager(MemoryManager const & source):
		blockFirst(TBlockIter()),
		blockLast(TBlockIter()),
		lastValue(TBlockIter())
	{
		SEQAN_CHECKPOINT 
		unsigned int x = SPACE;
		exponent = 0;
		while (x != 1){
			x >>=1;
			++exponent;
		}
		assign(*this, source);
	}

	~MemoryManager() 
	{
		clear(*this);
	}

	//____________________________________________________________________________



	public:
		template<typename TPos>
		inline typename Reference<MemoryManager>::Type 
			operator[] (TPos pos) 
		{
		SEQAN_CHECKPOINT	
			return value(*this, pos);
		}

		template<typename TPos>
		inline typename Reference<MemoryManager const>::Type 
			operator[] (TPos pos) const 
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


template<typename TValue, unsigned int SPACE>
struct DefaultOverflowImplicit< MemoryManager<TValue, Block<SPACE> , FreeMemoryInt> >
{
	typedef Generous Type;
};


///.Metafunction.Value.param.T.type:Class.MemoryManager
template <typename TValue, unsigned int SPACE>
struct Value<MemoryManager<TValue,Block<SPACE>, FreeMemoryInt > >
{
	typedef TValue Type;
};


template<typename TValue, unsigned int SPACE,  typename TSource>
inline void 
assign(
	MemoryManager<TValue, Block<SPACE>, FreeMemoryInt>& target, 
	TSource const& source) 
{
	SEQAN_CHECKPOINT
	clear(target);
	_raiseMemory(target,length(source));
	for (unsigned int i = 0; i < length(source);++i){
		target[i] = source[i];
	}

	target.lastValue = &target[length(source)-1];

	target.pEmptySpace = source.pEmptySpace;
	size_t tmpSource = source.pEmptySpace;
	size_t tmpTarget = source.pEmptySpace;
	//hmmmmmmmmmmm
	size_t x = supremumValue<size_t>();
	while (tmpSource != x){
		tmpSource = source[tmpSource];
		target[tmpTarget] = tmpSource;
		tmpTarget = target[tmpTarget];
	}
	target.change = source.change;
}

template<typename TValue, unsigned int SPACE,  typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(
	MemoryManager<TValue, Block<SPACE>, FreeMemoryInt >& stack, 
	TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos & (SPACE-1));
}


template<typename TValue, unsigned int SPACE,  typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type 
value(
	MemoryManager<TValue, Block<SPACE>, FreeMemoryInt > const& stack, 
	TPos const pos) 
{
	SEQAN_CHECKPOINT
	return value(*(stack.blocks[pos>>stack.exponent]), pos & (SPACE-1));
}

template<typename TValue, unsigned int SPACE>
inline void 
clear(MemoryManager<TValue, Block<SPACE>, FreeMemoryInt >& me)
{
SEQAN_CHECKPOINT
	typedef String<TValue, Block<SPACE>	>			TBlockString;
	typedef typename TBlockString::TBlockTable		TBlockTable;
	typedef typename Iterator<TBlockTable, Standard>::Type	TIter;
	
	TIter it = begin(me.blocks), itEnd = end(me.blocks);
	while (it != itEnd) {
		deallocate(me.alloc, *it, 1);
		++it;
	}
	me.pEmptySpace = -1;
	clear(me.blocks);
	me.lastValue = me.blockLast = typename TBlockString::TBlockIter();
}


template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
length(MemoryManager<TValue, Block<SPACE>, FreeMemoryInt > const & me) 
{
	SEQAN_CHECKPOINT
	if (length(me.blocks))
		return (length(me.blocks) - 1) * SPACE + (me.lastValue - me.blockFirst) + 1;
	else
		return 0;
}


template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
capacity(MemoryManager<TValue, Block<SPACE>, FreeMemoryInt > const & me) 
{
SEQAN_CHECKPOINT
	if (length(me.blocks))
		return length(me.blocks) * SPACE;
	else
		return 0;
}


template<typename TValue, unsigned int SPACE> 
void
_raiseMemory(MemoryManager<TValue,Block<SPACE>, FreeMemoryInt > &manager, 
			 typename Size<String<TValue, Block<SPACE> > >::Type leng)
{
	SEQAN_CHECKPOINT
	while (capacity(manager) < leng){
		typename Size< String<TValue, Block<SPACE> > >::Type last = length(manager.blocks);
		resize(manager.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
		allocate(manager.alloc, manager.blocks[last], 1);
		manager.lastValue = manager.blockFirst = begin(*manager.blocks[last]);
		manager.blockLast = (manager.blockFirst + (SPACE - 1));
	}
}

template<typename TValue, unsigned int SPACE,  typename TPosition> 
void
releaseID(MemoryManager<TValue,Block<SPACE>, FreeMemoryInt > &manager, 
		  TPosition pos)
{
	SEQAN_CHECKPOINT
	void* pTmp = &manager[pos];
	*(size_t*)pTmp = manager.pEmptySpace;
	manager.pEmptySpace = pos;
}


template<typename TValue, unsigned int SPACE> 
inline typename Size<String<TValue, Block<SPACE> > >::Type
obtainID(MemoryManager<TValue,Block<SPACE>,FreeMemoryInt > &manager)
{
	SEQAN_CHECKPOINT
	manager.change = false;
	if (manager.pEmptySpace == supremumValue<size_t>()){//-1){
		if (length(manager) == capacity(manager)){
			typename Size< String<TValue, Block<SPACE> > >::Type last = length(manager.blocks);
			resize(manager.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
			allocate(manager.alloc, manager.blocks[last], 1);
			manager.lastValue = manager.blockFirst = begin(*manager.blocks[last]);
			manager.blockLast = (manager.blockFirst + (SPACE - 1));
			manager.change = true;
		} else
			++manager.lastValue;
		return length(manager)-1;
		
	}else{
		size_t tmp = manager.pEmptySpace;
		void* pTmp = &manager[tmp];
		
		manager.pEmptySpace = *(size_t*) pTmp;
		return tmp;
	}
}


template<typename TValue, unsigned int SPACE,  typename TPosition> 
void
assignValue(MemoryManager<TValue,Block<SPACE>, FreeMemoryInt > &manager, 
			TPosition pos, 
			TValue value)
{
	SEQAN_CHECKPOINT
	manager[pos] = value;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
