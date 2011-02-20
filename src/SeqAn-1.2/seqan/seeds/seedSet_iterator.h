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
  $Id: seedSet_iterator.h 2334 2008-06-06 13:13:23Z kemena@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEEDSET_ITERATOR_H
#define SEQAN_HEADER_SEEDSET_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Iter
//////////////////////////////////////////////////////////////////////////////



struct SeedIterator;

template <typename TSeedSet>
class Iter<TSeedSet, SeedIterator>
{
public:
	typedef typename Size<TSeedSet>::Type TSize;
	typename std::set<TSize>::iterator data_ptr;
	TSeedSet *set;

	Iter()
	{
		this->set = 0;
	}

	Iter(TSeedSet &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet const &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet &set, Iter & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
		
	Iter(TSeedSet &set, typename std::set<TSize>::iterator other_)//:data_ptr(other_)
	{
		data_ptr = other_;
		this->set = &set;
	}
        
        
        
	Iter(TSeedSet &set, TSize &other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TSeedSet2>
	Iter(TSeedSet &set, Iter<TSeedSet2, SeedIterator> & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
	~Iter()
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
                this->set = other_.set;
		return *this;
	}
	
        /*
        Iter const &
	operator = (TValue * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}*/
        
	template <typename TContainer2>
	Iter const &
	operator = (Iter<TContainer2, SeedIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
                this->set = other_.set;
		return *this;
	}

	operator typename Value<TSeedSet>::Type * ()
	{
		return set->manager[*data_ptr];
	}

};


template <typename TSeedSet>
class Iter<TSeedSet const, SeedIterator>
{
public:
	typedef typename Size<TSeedSet>::Type TSize;
	typename std::set<TSize>::const_iterator data_ptr;
	TSeedSet const *set;

	Iter()
	{
		this->set = 0;
	}

	Iter(TSeedSet const &set)
	{
		this->set = &set;
	}

	Iter(TSeedSet &set, Iter const & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}


	Iter(TSeedSet const &set, typename std::set<TSize>::const_iterator const & other_):
		data_ptr(other_)
	{
		//data_ptr = other_;
		this->set = &set;
	}

	Iter(TSeedSet const  &set, TSize * other_data_ptr):
		data_ptr(other_data_ptr)
	{
	}
	template <typename TSeedSet2>
	Iter(TSeedSet const &set, Iter<TSeedSet2 const, SeedIterator> const & other_):
		data_ptr(other_.data_ptr)
	{
		this->set = &set;
	}
		
		
	~Iter()
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		this->data_ptr = other_.data_ptr;
		this->set = other_.set;
		return *this;
	}
	Iter const &
	operator = (TSize * other_data_ptr)
	{
		data_ptr = other_data_ptr;
		return *this;
	}
	template <typename TSeedSet2>
	Iter const &
	operator = (Iter<TSeedSet2 const, SeedIterator> const & other_)
	{
		this->data_ptr = other_.data_ptr;
		this->set = other_.set;
		return *this;
	}

};


template <typename TSeedSet>
struct Value<Iter<TSeedSet const, SeedIterator> >
{
	typedef typename Value<TSeedSet>::Type Type;
};



template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
operator ++(Iter<TSeedSet, SeedIterator > &me)
{
SEQAN_CHECKPOINT
	++me.data_ptr;
	return me;
}

template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
goNext(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	++(it.data_ptr);
	return it;
}

template<typename TSeedSet>
inline Iter<TSeedSet, SeedIterator >&
operator --(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	--(it.data_ptr);
	return it;
}

template<typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> >::Type
value(Iter<TSeedSet, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	return it.set->manager[*it.data_ptr];
}

template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> const>::Type 
 value(Iter<TSeedSet, SeedIterator> const &me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}



template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet, SeedIterator> const>::Type 
operator * (Iter<TSeedSet, SeedIterator> & me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}

template <typename TSeedSet>
inline typename Reference<Iter<TSeedSet const, SeedIterator> const>::Type 
operator * (Iter<TSeedSet const, SeedIterator> & me)
{
SEQAN_CHECKPOINT
	return me.set->manager[*me.data_ptr];
}

template<typename TSeedSet>
inline bool
operator !=(Iter<TSeedSet, SeedIterator > it1, Iter<TSeedSet, SeedIterator > it2)
{
	return it1.data_ptr != it2.data_ptr;
}

/**
.Function.seedScore:
..summary:Returns the score of a seed.
..cat:Seed Handling
..signature:seedScore(it);
..param.it: The seedSet iterator.
..return:Score of the seed.
*/
template<typename TSeedSet>
inline typename ScoreType<typename ScoringScheme<TSeedSet>::Type>::Type &
seedScore(Iter<TSeedSet, SeedIterator > it)
{
SEQAN_CHECKPOINT
    return it.set->scoreMap[*it.data_ptr];
}



template<typename TSeedSet>
inline typename ScoreType<typename ScoringScheme<TSeedSet>::Type>::Type const&
seedScore(Iter<TSeedSet const, SeedIterator > &it)
{
SEQAN_CHECKPOINT
	return it.set->scoreMap[*it.data_ptr];
}

template<typename TSeedSet, typename TSize>
void
setScore(Iter<TSeedSet, SeedIterator > &it, TSize score)
{
SEQAN_CHECKPOINT
	it.set->scoreMap[*it.data_ptr] = score;
}
//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
