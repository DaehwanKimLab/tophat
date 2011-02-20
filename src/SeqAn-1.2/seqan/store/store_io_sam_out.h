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
  $Id$
 ==========================================================================*/

#include <iostream>

#ifndef SEQAN_HEADER_STORE_IO_SAM_OUT_H
#define SEQAN_HEADER_STORE_IO_SAM_OUT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TSpec, typename TConfig>
	inline void _writeHeader(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 SAM)
    {
		typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
		typedef typename TFragmentStore::TLibraryStore					TLibraryStore;
		typedef typename TFragmentStore::TContigStore					TContigStore;
		typedef typename TFragmentStore::TNameStore						TNameStore;

        typedef typename Value<TContigStore>::Type						TContig;
        typedef typename Iterator<TLibraryStore, Standard>::Type		TLibraryIter;
        typedef typename Iterator<TContigStore, Standard>::Type			TContigIter;
        typedef typename Iterator<TNameStore, Standard>::Type			TContigNameIter;
        typedef typename Id<TContig>::Type								TId;

        TContigIter it = begin(store.contigStore, Standard());
        TContigIter itEnd = end(store.contigStore, Standard());
		TContigNameIter nit = begin(store.contigNameStore, Standard());
		TContigNameIter nitEnd = end(store.contigNameStore, Standard());
		
		_streamWrite(target, "@HD\tVN:1.0\n");
        for(; it != itEnd; ++it)
		{
			_streamWrite(target, "@SQ\tSN:");
			if (nit != nitEnd)
			{
				_streamWrite(target, *nit);
				++nit;
			}
			_streamWrite(target, "\tLN:");
			_streamPutInt(target, length((*it).seq));
            _streamPut(target, '\n');
		}

		TLibraryIter lit = begin(store.libraryStore, Standard());
		TLibraryIter litEnd = end(store.libraryStore, Standard());
        for(TId id = 0; lit != litEnd; ++lit, ++id)
		{
			_streamWrite(target, "@RG\tID:");
			_streamPutInt(target, id + 1);
			_streamWrite(target, "\tLB:");
			_streamWrite(target, store.libraryNameStore[id]);
			_streamWrite(target, "\tPI:");
			_streamPutInt(target, (int)/*std::round*/(store.libraryStore[id].mean));
			_streamWrite(target, "\tSM:none");	// sample name needs to be included into fragment store
            _streamPut(target, '\n');
		}
		_streamWrite(target, "@PG\tID:SeqAn\n");
	}
	
	
//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TSpec, typename TConfig>
	inline void _writeAlignments(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 SAM)
    {
		typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

		typedef typename TFragmentStore::TReadStore						TReadStore;
		typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
		typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
		typedef typename TFragmentStore::TContigStore					TContigStore;
		typedef typename TFragmentStore::TReadSeq						TReadSeq;

        typedef typename Value<TReadStore>::Type						TRead;
        typedef typename Value<TReadSeqStore>::Type						TReadSeqStored;
        typedef typename Value<TContigStore>::Type						TContig;
        typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

		typedef typename TContig::TContigSeq							TContigSeq;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignIter;
        typedef typename Iterator<TReadSeqStored, Standard>::Type		TReadSeqIter;
        typedef typename Id<TAlignedRead>::Type							TId;

		typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
		typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;

		String<int> mateIndex;	// store outer library size for each pair match (indexed by pairMatchId)
		calculateMateIndices(mateIndex, store);
		
        TAlignIter it = begin(store.alignedReadStore, Standard());
        TAlignIter itEnd = end(store.alignedReadStore, Standard());
		TAlignIter mit = it;
		CharString cigar;
		TReadSeq readSeq;
		
        for(; it != itEnd; ++it)
		{
            TId alignedId = (*it).id;
			TId readId = (*it).readId;
			TId mateIdx = TRead::INVALID_ID;
			if ((*it).pairMatchId != TRead::INVALID_ID)
				mateIdx = mateIndex[2*(*it).pairMatchId + getMateNo(store, (*it).readId)];

			TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
			__int64 pos = positionGapToSeq(contigGaps, _min((*it).beginPos, (*it).endPos)) + 1;
			__int64 mpos = 0;
			int isize = 0;
			unsigned short flag = 0;

			if ((*it).beginPos > (*it).endPos)
				flag |= 0x0010;			

			// calculate flags, mpos, isize
			if (mateIdx < length(store.alignedReadStore))
			{
				mit = begin(store.alignedReadStore, Standard()) + mateIdx;
				if ((*it).contigId == (*mit).contigId)
				{
					mpos = positionGapToSeq(contigGaps, _min((*mit).beginPos, (*mit).endPos)) + 1;
					if ((*it).beginPos < (*mit).beginPos)
						isize = positionGapToSeq(contigGaps, _max((*mit).beginPos, (*mit).endPos) - 1) + 2 - pos;
					else
						isize = mpos - positionGapToSeq(contigGaps, _max((*it).beginPos, (*it).endPos) - 1) - 2;
				}
				flag |= 0x0002;
				if ((*mit).beginPos > (*mit).endPos)
					flag |= 0x0020;				
			}
			unsigned char mateNo = getMateNo(store, readId);
			if (mateNo == 0) flag |= 0x0040;
			if (mateNo == 1) flag |= 0x0080;

			if (readId < length(store.readStore))
			{
				TRead &read = store.readStore[readId];
				if (read.matePairId != TRead::INVALID_ID)
					flag |= 0x0001;
			}
			
			// <qname>
			if (readId < length(store.readNameStore))
				_streamWrite(target, store.readNameStore[readId]);
            _streamPut(target, '\t');
            
            // <flag>
            _streamPutInt(target, flag);
            _streamPut(target, '\t');
            
			// <rname>
			if ((*it).contigId < length(store.contigNameStore))
				_streamWrite(target, store.contigNameStore[(*it).contigId]);
            _streamPut(target, '\t');
            
			// <pos>
            _streamPutInt(target, pos);
            _streamPut(target, '\t');
            
			// <mapq>
			if (alignedId < length(store.alignQualityStore))
				_streamPutInt(target, store.alignQualityStore[alignedId].score);
            _streamPut(target, '\t');
            
			// get read sequence
			if (readId < length(store.readSeqStore))
			{
				readSeq = store.readSeqStore[readId];
				if ((*it).beginPos <= (*it).endPos) 
				{
					setBeginPosition(contigGaps, (*it).beginPos);
					setEndPosition(contigGaps, (*it).endPos);
				} else
				{
					setBeginPosition(contigGaps, (*it).endPos);
					setEndPosition(contigGaps, (*it).beginPos);
					reverseComplementInPlace(readSeq);
				}
			} else
				clear(readSeq);
			
            // <cigar>
			TReadGaps readGaps(readSeq, (*it).gaps);
			getCigarString(cigar, contigGaps, readGaps);
			
			_streamWrite(target, cigar);
            _streamPut(target, '\t');
            
            // <mrnm>
			if ((mateIdx < length(store.alignedReadStore)))
			{
				if ((*it).contigId == (*mit).contigId)
					_streamWrite(target, '=');
				else
					if ((*mit).contigId < length(store.contigNameStore))
						_streamWrite(target, store.contigNameStore[(*mit).contigId]);
			} else
				_streamWrite(target, '*');
				
            _streamPut(target, '\t');
            
            // <mpos>
			_streamPutInt(target, (int)mpos);
            _streamPut(target, '\t');
            
            // <isize>
			_streamPutInt(target, isize);
            _streamPut(target, '\t');

            // <seq>
			_streamWrite(target, readSeq);
            _streamPut(target, '\t');
            
            // <qual>
			TReadSeqIter it = begin(store.readSeqStore[readId], Standard());
			TReadSeqIter itEnd = end(store.readSeqStore[readId], Standard());
			for (; it != itEnd; ++it)
				_streamPut(target, (char)(getQualityValue(*it) + 33));
            
			// <tags>
			if (alignedId < length(store.alignedReadTagStore) && !empty(store.alignedReadTagStore[alignedId]))
			{
				_streamPut(target, '\t');
				_streamWrite(target, store.alignedReadTagStore[alignedId]);
			}
            
            _streamPut(target, '\n');
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// write
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void write(TFile & target,
                      FragmentStore<TSpec, TConfig> & store,
                      SAM)
    {
        // write header
		_writeHeader(target, store, SAM());
        
        // write aligments
        _writeAlignments(target, store, SAM());
    }
    
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
