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

#ifndef SEQAN_HEADER_BASIC_ALLOCATOR_TO_STD_H
#define SEQAN_HEADER_BASIC_ALLOCATOR_TO_STD_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Filter that adapts seqan allocator zu std allocator
/**
.Class.ToStdAllocator:
..summary:Emulates standard conform allocator.
..signature:ToStdAllocator<THost, TValue>
..param.THost:Type of the host allocator object.
...text:This object is used to call @Function.allocate@ and @Function.deallocate@.
..param.TValue:Type of allocated items.
..remarks:The member functions $allocate$ and $deallocate$ of $ToStdAllocator$ call
the (globale) functions @Function.allocate@ and @Function.deallocate@, respectively. The globale functions
get an allocator object as their first arguments. This allocator object is not the $ToStdAllocator$ object itself,
but the host object that was given to the constructor. 
..cat:Basic
..remarks:
..see:Function.allocate
..see:Function.deallocate
..include:seqan/basic.h
*/
template <typename THost, typename TValue>
struct ToStdAllocator
{
	typedef TValue value_type;
	typedef value_type * pointer;
	typedef value_type & reference;
	typedef value_type const * const_pointer;
	typedef value_type const & const_reference;

//	typedef typename THost::Size size_type;
//	typedef typename THost::Difference difference_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

/**
.Memfunc.ToStdAllocator:
..summary:Constructor
..signature:ToStdAllocator(host)
..class:Class.ToStdAllocator
..param.host:The host object that is used as allocator for @Function.allocate@ and @Function.deallocate@.
*/
	ToStdAllocator(THost & host): m_host(& host)
	{
	}
	template <typename TValue2>
	ToStdAllocator(ToStdAllocator<THost, TValue2> const & alloc): m_host(alloc.m_host)
	{
	}
	ToStdAllocator & operator= (ToStdAllocator const & alloc)
	{
		m_host = alloc.m_host;
		return *this;
	}
	~ToStdAllocator()
	{
	}

// TODO(holtgrew): Move to basic_host.h?
/**
.Function.host:
..summary:The object a given object depends on.
..cat:Dependent Objects
..signature:host(object)
..param.object:An object.
...type:Class.ToStdAllocator
..returns:The host object.
..include:seqan/basic.h
*/
	pointer allocate(size_type count)
	{
		value_type * ptr;
		seqan::allocate(*m_host, ptr, count);
		return pointer(ptr);
	}
	pointer allocate(size_type count, const void *)
	{
		value_type * ptr;
		seqan::allocate(*m_host, ptr, count);
		return pointer(ptr);
	}

	void deallocate(pointer data, size_type count)
	{
		seqan::deallocate(*m_host, data, count);
	}

	void construct(pointer ptr, const_reference data)
	{
		new(ptr) TValue(data);
	}

	void destroy(pointer ptr)
	{
		ptr->~TValue();
	}

	pointer address(reference value) const
	{
		return (&value);
	}
	const_pointer address(const_reference value) const
	{
		return (&value);
	}

	size_type max_size() const
	{
		return ~0UL / sizeof(value_type);
	}

	template<class TValue2>
	struct rebind
	{
		typedef ToStdAllocator<THost, TValue2> other;
	};

	template <typename THost2, typename TValue2>
	friend
	struct ToStdAllocator;

	private:
		THost * m_host;
};

template <typename THost, typename TValue>
THost & 
host(ToStdAllocator<THost, TValue> & me)
{
   return *me.m_host;
}


//////////////////////////////////////////////////////////////////////////////



//returns std-allocator type (for allocators)
template <typename T, typename TData>
struct StdAllocator
{
	typedef ToStdAllocator<T, TData> Type;
};


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN


//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
