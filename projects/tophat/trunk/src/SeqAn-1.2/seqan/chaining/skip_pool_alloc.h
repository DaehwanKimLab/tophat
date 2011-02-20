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
  $Id: skip_pool_alloc.h 3420 2009-02-12 12:10:09Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

/*	Hendrik Woehrle
*
*	Deferred Skip List Datastructure
*
*	ClassPool -
*
*	Memory pool allocator for multiple objects at once
*
*/

#ifndef SEQAN_HEADER_SKIP_POOL_ALLOC_H
#define SEQAN_HEADER_SKIP_POOL_ALLOC_H

namespace seqan
{

struct Unlimited
{};

struct Limited
{};

/*DISABLED
.Spec.Class Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks for a specific class.
..signature:Allocator< ClassPool<Class, Type, ParentAllocator> >
..param.Class:The class.
..param.Type:A specialization. The Class Pool Allocator can manage either slices of a fixed size or of 
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap allocator is reduced and that speeds up memory management. The Class Pool Allocator
pools only memory blocks of a fixed size, namely $sizeof(Class)$. The Class Pool Allocator only pools the memory, constructor and destructor have to be called manually.
*/

template< typename TClass, typename TType, typename TParentAllocator = SimpleAllocator >
struct ClassPool;

template< typename TClass, typename TSpec, typename TParentAlloc >
void clear( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );

template< typename TClass, typename TSpec, typename TParentAlloc, typename TSize >
void allocate( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me, 
				TClass *& dest, 
				TSize number );

template< typename TClass, typename TSpec, typename TParentAlloc, typename TSize >
void deallocate( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me, 
				 TClass * location );

template< typename TClass, typename TSpec, typename TParentAlloc >
TParentAlloc & 
parentAlloc( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );

template< typename TClass, typename TSpec, typename TParentAlloc >
void
setParentAlloc( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me );


template< typename TClass, typename TSpec, typename TParentAlloc > inline
TClass * 
_getNextBlock( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & me,
				TClass & block )
{
	return _getNext( block );
}

template< typename TClass, typename TSpec, typename TParentAlloc > inline
void 
_setNextBlock( Allocator< ClassPool< TClass, TSpec, TParentAlloc > > & /*me*/,
				TClass & dest,
				TClass * block )
{
	_setNext( dest, block );
}


//********************************** Partially specialised Template Class of the memory allocator for SkipBaseElement
template< typename TClass, typename TParentAllocator  >
struct Allocator< ClassPool< TClass, Unlimited, TParentAllocator > >
{
		// Last free block in queue, which can be used again
	TClass * _freeBlock;
		// Pointer to end of "used" memory Block
	TClass * _end;
		// Pointer to single end of reserved memory
	TClass * _terminal_end;
		// size of blocks
	typename Size< TClass >::Type _blockSize;
		// parent Allocator
	Holder<TParentAllocator> data_parent_allocator;
		

	Allocator(Allocator const &):
		_freeBlock(0),
		_end(0),
		_terminal_end(0),
		_blockSize(0)
	{
	}

	inline Allocator &
	operator = (Allocator const &)
	{
		clear(*this);
		return *this;
	}


	Allocator( typename Size< TClass >::Type numElements = 1000 )
	{
		TClass * block;
		_blockSize = numElements > 20 ? numElements : 20;
		allocate( parentAllocator( *this ), block, _blockSize );
		_end = block;
		_terminal_end = _end + _blockSize;
		_freeBlock = NULL;
	}

	~Allocator(void)
	{
		clear( *this );
	}

};

	template< typename TClass, typename TParentAlloc >
	inline 
	TParentAlloc &
	parentAllocator( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me )
	{
		return value( me.data_parent_allocator );
	}

	template< typename TClass, typename TParentAlloc >
	inline 
	void
	setParentAllocator( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me,
					   TParentAlloc & alloc_)
	{
		setValue( me.data_parent_allocator, alloc_ );
	}


	template< typename TClass, typename TParentAlloc, typename TSize >
	void allocate( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me, 
					TClass *& dest, 
					TSize number )
	{
		SEQAN_CHECK2( number <= me._blockSize, "tried to allocate more elements than available in block");
		SEQAN_CHECK2( number > 0, "tried to allocate 0 elements");
		if( number == 1 )
		{
			if( me._end < me._terminal_end )
			{		// enough place available in current memory pool
				dest = me._end;
				++me._end;
				return;
			}
			else{
				allocate( parentAllocator( me ), me._end, me._blockSize );
				me._terminal_end = me._end + me._blockSize;
				dest = me._end;
				++me._end;
				return;
			}
		}
		else
		{
			if( me._end + number < me._terminal_end )
			{		// enough place available in current memory pool
				dest = me._end;
				me._end += number;
				return;
			}
			else{		
					// not enough memory in current pool avaiable: allocate new, adjust pointers
				if( me._end < me._terminal_end ){
					_setNext( * me._end, me._freeBlock );
					++me._end;
					while( me._end < me._terminal_end){
						_setNext( *me._end, me._end );
						++me._end;
					}
					me._freeBlock = me._end-1;
				}
				allocate( parentAllocator( me ), me._end, me._blockSize );
				me._terminal_end = me._end + me._blockSize;
				dest = me._end;
				me._end+=number;
				return;
			}
		}	
	}

	
	template< typename TClass, typename TParentAlloc, typename TSize > inline
	void deallocate( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me, 
					 TClass * location,
					 TSize count )
	{
		SEQAN_CHECK2( location != NULL, "Tried to free NULL-pointer" )
		_setNext( *location, me._freeBlock );
		me._freeBlock = location;
	}


	template< typename TClass, typename TParentAlloc > inline
	void
	clear( Allocator< ClassPool< TClass, Unlimited, TParentAlloc > > & me )
	{
	SEQAN_CHECKPOINT

		me._freeBlock = NULL;
		me._end = NULL;
		me._terminal_end = NULL;
		me._blockSize = 0;

		clear( parentAllocator( me ) );
	}


//***************************** Partially specialised Template Class of the memory allocator for SkipElement
template< typename TClass, typename TParentAllocator  >
struct Allocator< ClassPool< TClass, Limited, TParentAllocator > >
{
		// Map for freed Blocks, which can be used again
	TClass * _freeBlocks[ sizeof( typename Size< TClass >::Type ) * 8];
		// Pointer to end of "used" memory Block
	TClass * _end;
		// Pointer to single end of reserved memory
	TClass * _terminal_end;
		// size of Elements
	typename Size< TClass >::Type _blockSize;
		// parent Allocator
	Holder<TParentAllocator> data_parent_allocator;


public:

	Allocator( typename Size< TClass >::Type numElements = 1000 )
	{
		_blockSize = numElements > 32 ? numElements : 32;
		allocate( parentAllocator( *this ), _end, _blockSize );
		_terminal_end = _end + _blockSize;
		for( typename Size< TClass >::Type i = 0; i < sizeof( typename Size< TClass >::Type ) * 8; i++ )
			_freeBlocks[i] = NULL;
	}

	~Allocator(void)
	{
		clear( *this );
	}

};

	template< typename TClass, typename TParentAlloc >
	inline 
	TParentAlloc &
	parentAllocator( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me )
	{
		return value(me.data_parent_allocator);
	}

	template< typename TClass, typename TParentAlloc >
	inline 
	void
	setParentAllocator( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me,
					   TParentAlloc & alloc_)
	{
		setValue( me.data_parent_allocator, alloc_ );
	}

	template< typename TClass, typename TParentAlloc, typename TSize >
	void allocate( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me, 
					TClass *& dest, 
					TSize number )
	{
		SEQAN_CHECK2( number <= me._blockSize, "tried to allocate more elements than available in block")
		SEQAN_CHECK( number != 0 )		
			// recycle old memory block
		if( me._freeBlocks[ number-1 ] != NULL )
		{
			dest = me._freeBlocks[number-1];
			me._freeBlocks[number-1] = _getNext( *dest );
			return;
		}
			// use in memory from current block
		else if( me._end + number < me._terminal_end )
		{
			dest = me._end;
			me._end += number;
			return;
		}
		else {	// allocate new memory block
			TClass * block;
			allocate( parentAllocator( me ), block, me._blockSize );
			typename Size< TClass >::Type rest_space = me._terminal_end - me._end;
			if( rest_space != 0 )
			{
				_setNextBlock( me, *me._end, me._freeBlocks[ rest_space - 1 ] );
				me._freeBlocks[ rest_space - 1 ] = me._end;
			}
			me._end = block;
			me._terminal_end = me._end + me._blockSize;
			dest = me._end;
			me._end += number;
		}
	}

	template< typename TClass, typename TParentAlloc, typename TSize >
	void 
	deallocate( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me,
				TClass * location,
				TSize count )
	{
		SEQAN_CHECK( count != 0 )
		_setNextBlock( me, *location, me._freeBlocks[ count - 1 ] );
		me._freeBlocks[ count - 1 ] = location;
	}

	template< typename TClass, typename TParentAlloc >
	void
	clear( Allocator< ClassPool< TClass, Limited, TParentAlloc > > & me )
	{
	SEQAN_CHECKPOINT

		for( size_t i = 0; i < sizeof( typename Size< TClass >::Type ) * 8; ++i )
			me._freeBlocks[i] = NULL;
		me._end = NULL;
		me._terminal_end = NULL;
		me._blockSize = 0;

		clear( parentAllocator( me ) );
	}

} // namespace seqan

#endif // SKIPPOOLALLOC_H
