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
 $Id: blast_hsp_iterator.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_HSP_ITERATOR_H
#define SEQAN_HEADER_BLAST_HSP_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Hsp Iterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



/**
.Spec.HspIterator:
..cat:Blast
..summary:Hsp iterator for @Class.BlastHit@.
..signature:Iterator<TBlastHit, HspIterator>
..param.TBlastHit:A Blast hit.
...type:Class.BlastHit
..general:Class.Iter
*/
template<typename TBlastHsp, typename TStoreSpec>
class Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > 
{
public:

	BlastHit<TBlastHsp, TStoreSpec> * data_host;
	unsigned int data_pos;

	Iter()	
	{
	}
	
	Iter(BlastHit<TBlastHsp, TStoreSpec> & blast) : 
		data_host(&blast), 
		data_pos(0) 
	{
	SEQAN_CHECKPOINT
	}

	Iter(Iter const& it) : 
		data_host(it.data_host), 
		data_pos(it.data_pos) 
	{
	SEQAN_CHECKPOINT
	}

	~Iter() 
	{
	SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & other) 
	{
	SEQAN_CHECKPOINT
		if (this == &other) return *this;
		data_host = other.data_host;
		data_pos = other.data_pos;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<HspIterator> - Metafunctions
//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastHit<TBlastHsp,StoreReport<TInfoSpec> >, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StoreReport<TInfoSpec> >, SimpleBlastIterator<HspIterator> > Type;
};

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastHit<TBlastHsp,StoreReport<TInfoSpec> > const, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StoreReport<TInfoSpec> > const, SimpleBlastIterator<HspIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TInfoSpec>
struct Value<Iter<BlastHit<TBlastHsp,StoreReport<TInfoSpec> >, SimpleBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TInfoSpec>
struct Value<Iter<BlastHit<TBlastHsp,StoreReport<TInfoSpec> > const, SimpleBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};





//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<HspIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline typename GetValue<Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > >::Type
getValue(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hsps[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHsp, typename TStoreSpec>
inline typename Reference<Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > >::Type
value(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hsps[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline typename Reference<Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > >::Type
operator * (Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

////**
//.Function.hostHit:
//..cat:Blast
//..summary:The BlastHit this iterator is working on.
//..signature:hostHit(it)
//..param.it:A Blast hsp iterator.
//...type:Spec.HspIterator
//..returns:A pointer to the host Blast hit.
//*/
//template<typename TBlastHsp, typename TStoreSpec>
//inline typename Host<Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > >::Type &
//hostHit(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
//{
//	return *it.data_host;
//} 

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline bool
atBegin(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == 0);	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline void
goBegin(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = 0;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHsp, typename TStoreSpec>
inline bool
atEnd(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == length(it.data_host->hsps));	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline void
goEnd(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	it.data_pos = length(it.data_host->hsps);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline void
goNext(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) ++it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >&
operator ++(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline void
goPrevious(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) --it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >&
operator --(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastHsp, typename TStoreSpec>
//inline Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >
//operator --(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it, int)
//{
//	Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> > ret = it;
//	goPrevious(it);
//	return ret;
//}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline bool
operator ==(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it1,
			Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
inline bool
operator !=(Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it1,
			Iter<BlastHit<TBlastHsp, TStoreSpec> , SimpleBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
