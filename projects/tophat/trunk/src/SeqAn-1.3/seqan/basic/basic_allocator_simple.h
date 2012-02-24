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

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_SIMPLE_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// SimpleAlloc Allocator
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Simple Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:General purpose allocator.
..signature:Allocator< SimpleAlloc<ParentAllocator> >
..param.ParentAllocator:An allocator that is by the simple allocator used to allocate memory.
...default:@Tag.Default@
...remarks:@Tag.Default@ used as allocator means that the default implementations
of @Function.allocate@ and @Function.deallocate@ are used.
..include:seqan/basic.h
*/

template <typename TParentAllocator = Default>
struct SimpleAlloc;

//////////////////////////////////////////////////////////////////////////////


typedef Allocator<SimpleAlloc<Default> > SimpleAllocator;

template <typename TParentAllocator>
struct Allocator<SimpleAlloc<TParentAllocator> >
{
	struct Header
	{
		Header * left;
		Header * right;
		size_t size;
	};

	Header * data_storages;
	Holder<TParentAllocator> data_parent_allocator;

	Allocator():
		data_storages(0)
	{
SEQAN_CHECKPOINT
	}

	Allocator(TParentAllocator & parent_alloc):
		data_storages(0)
	{
SEQAN_CHECKPOINT
		setValue(data_parent_allocator, parent_alloc);
	}

	//Dummy copy
	Allocator(Allocator const &):
		data_storages(0)
	{
	}
	inline Allocator &
	operator = (Allocator const &)
	{
		clear(*this);
		return *this;
	}

	~Allocator()
	{
SEQAN_CHECKPOINT
		clear(*this);
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SimpleAlloc<TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	return value(me.data_parent_allocator);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.Allocator#clear:
..cat:Memory
..summary:Deallocates all memory blocks.
..signature:clear(allocator)
..param.allocator:Allocator object.
...type:Class.Allocator
...concept:Concept.Allocator
..remarks:This function deallocates all memory blocks 
that was allocated using @Function.allocate@ for $allocator$.
The memory is not pooled but directly passed back to the heap manager.
..see:Function.allocate
..see:Function.deallocate
..include:seqan/basic.h
*/
template <typename TParentAllocator>
void
clear(Allocator<SimpleAlloc<TParentAllocator> > & me)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;

	while (me.data_storages)
	{
		typename TAllocator::Header * next_storage = me.data_storages->right;
		deallocate(parentAllocator(me), reinterpret_cast<char *>(me.data_storages), me.data_storages->size);
		me.data_storages = next_storage;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<SimpleAlloc<TParentAllocator> > & me, 
		 TValue * & data,
		 TSize count,
		 Tag<TUsage> const)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;
	typedef typename TAllocator::Header THeader;

	//compute needed bytes
	size_t bytes_needed = count * sizeof(TValue) + sizeof(THeader);

	//allocate storage from parent
	char * ptr;
	allocate(parentAllocator(me), ptr, bytes_needed, TagAllocateStorage());

	THeader * new_block = reinterpret_cast<THeader *>(ptr);
	new_block->left = 0;
	new_block->right = me.data_storages;
	new_block->size = bytes_needed;

	if (me.data_storages)
	{
		me.data_storages->left = new_block;
	}
	me.data_storages = new_block;

	//return data
	data = reinterpret_cast<TValue *>(ptr + sizeof(THeader));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SimpleAlloc<TParentAllocator> > & me,
		   TValue * data, 
		   TSize,
		   Tag<TUsage> const)
{
SEQAN_CHECKPOINT
	typedef Allocator<SimpleAlloc<TParentAllocator> > TAllocator;
	typedef typename TAllocator::Header THeader;

	//update links
	THeader & header = *(reinterpret_cast<THeader *>(data) - 1);
	if (header.left)
	{
		header.left->right = header.right;
	}
	else
	{
		me.data_storages = header.right;
	}
	if (header.right)
	{
		header.right->left = header.left;
	}

	//deallocate storage using parent
	char * ptr = reinterpret_cast<char *>(& header);
	deallocate(parentAllocator(me), ptr, header.size);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
