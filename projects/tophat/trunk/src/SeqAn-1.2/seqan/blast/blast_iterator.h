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
 $Id: blast_iterator.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_ITERATOR_H
#define SEQAN_HEADER_BLAST_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Iterators
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
struct SimpleBlastIterator;	// for BlastReport<TBlastHsp,StoreBlast,TSpec>


// Hit Iterator
struct HitIterator;

// Hsp Iterator
struct HspIterator;


//////////////////////////////////////////////////////////////////////////////



template <typename TSpec>
struct StreamBlastIterator;	// for BlastReport<TBlastHsp,StreamBlast,TSpec>



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
