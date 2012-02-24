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

#include <iostream>

#ifndef SEQAN_HEADER_STORE_IO_SAM_H
#define SEQAN_HEADER_STORE_IO_SAM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Sam:
	Sam mapping file.
..include:seqan/store.h
*/
struct Sam_;
typedef Tag<Sam_> const Sam;


//////////////////////////////////////////////////////////////////////////////
// CIGAR struct
//////////////////////////////////////////////////////////////////////////////

#if SEQAN_HAS_SAMTOOLS
    struct FromBam_;
    typedef Tag<FromBam_> FromBam;
#endif  // #if SEQAN_HAS_SAMTOOLS
    
    template <typename TOperation_ = char, typename TCount_ = unsigned>
	struct CigarElement
	{
		typedef TOperation_ TOperation;
		typedef TCount_		TCount;

		TOperation			operation;
		TCount				count;

		CigarElement() : operation(0), count(0) {}
		
		CigarElement(TOperation o, TCount c):
			operation(o),
			count(c) {}

#if SEQAN_HAS_SAMTOOLS
        CigarElement(__uint32 bamCigarElement, FromBam const &)
        {
            SEQAN_ASSERT_LEQ(bamCigarElement & BAM_CIGAR_MASK, 8u);
            operation = "MIDNSHP=X"[bamCigarElement & BAM_CIGAR_MASK];
            count = bamCigarElement >> 4;
        }
#endif  // #if SEQAN_HAS_SAMTOOLS
	};

template <typename TOperation, typename TCount>
__uint32 toBamCigarElement(CigarElement<TOperation, TCount> const & cigarElement)
{
    char operation = 0;
    switch (cigarElement.operation) {
        case 'X': operation += 1;
        case '=': operation += 1;
        case 'P': operation += 1;
        case 'H': operation += 1;
        case 'S': operation += 1;
        case 'N': operation += 1;
        case 'D': operation += 1;
        case 'I': operation += 1;
        case 'M': break;
    }
    return (cigarElement.count << 4) | operation;
}

//____________________________________________________________________________

template <
	typename TCigar,
	typename TGaps1,
	typename TGaps2,
	typename TThresh>
inline void
getCigarString(
	TCigar &cigar,
	TGaps1 &gaps1,
	TGaps2 &gaps2,
	TThresh splicedGapThresh)
{
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	clear(cigar);
	char op, lastOp = ' ';
	unsigned numOps = 0;

//	std::cout << gaps1 << std::endl;
//	std::cout << gaps2 << std::endl;
	for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
	{
		if (isGap(it1))
		{
			if (isGap(it2))
				op = 'P';
			else if (isClipped(it2))
				op = '?';
			else
				op = 'I';
		} 
		else if (isClipped(it1))
		{
			op = '?';
		}
		else 
		{
			if (isGap(it2))
				op = 'D';
			else if (isClipped(it2))
				op = 'S';
			else
				op = 'M';
		}
		if (lastOp != op)
		{
			if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
				lastOp = 'N';
			if (numOps > 0)
			{
				std::stringstream num;
				num << numOps;
				append(cigar, num.str());
				appendValue(cigar, lastOp);
			}
			numOps = 0;
			lastOp = op;
		}
		++numOps;
	}
	SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
	if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
		lastOp = 'N';
	if (numOps > 0)
	{
		std::stringstream num;
		num << numOps;
		append(cigar, num.str());
		appendValue(cigar, lastOp);
	}
}

template <
	typename TCigar,
	typename TGaps1,
	typename TGaps2>
inline void
getCigarString(
	TCigar &cigar,
	TGaps1 &gaps1,
	TGaps2 &gaps2)
{
	return getCigarString(cigar, gaps1, gaps2, 20);
}

template <
    typename TOperation,
    typename TCount,
	typename TSpec,
	typename TGaps1,
	typename TGaps2,
	typename TThresh>
inline void
getCigarString(
        String<CigarElement<TOperation, TCount>, TSpec> &cigar,
        TGaps1 &gaps1,
        TGaps2 &gaps2,
        TThresh splicedGapThresh)
{
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	clear(cigar);
	char op, lastOp = ' ';
	unsigned numOps = 0;

//	std::cout << gaps1 << std::endl;
//	std::cout << gaps2 << std::endl;
	for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
	{
		if (isGap(it1))
		{
			if (isGap(it2))
				op = 'P';
			else if (isClipped(it2))
				op = '?';
			else
				op = 'I';
		} 
		else if (isClipped(it1))
		{
			op = '?';
		}
		else 
		{
			if (isGap(it2))
				op = 'D';
			else if (isClipped(it2))
				op = 'S';
			else
				op = 'M';
		}
		if (lastOp != op)
		{
			if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
				lastOp = 'N';
			if (numOps > 0)
				appendValue(cigar, CigarElement<>(op, numOps));
			numOps = 0;
			lastOp = op;
		}
		++numOps;
	}
	SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
	if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
		lastOp = 'N';
	if (numOps > 0)
        appendValue(cigar, CigarElement<>(op, numOps));
}

//////////////////////////////////////////////////////////////////////////////
// _getClippedLength
    
    template <typename TCigarString, typename TNum>
    inline void _getClippedLength(TCigarString const & cigar, TNum & sum)
    {
        typedef typename Iterator<TCigarString, Standard>::Type TCigarIter;
        
        TCigarIter it = begin(cigar, Standard());
        TCigarIter itEnd = end(cigar, Standard());
        
        sum = 0;        
        for (; it != itEnd; ++it)
            if (getValue(it).operation != 'S' && getValue(it).operation != 'H')
                sum += getValue(it).count;
    }

//////////////////////////////////////////////////////////////////////////////
// convert CIGAR to gaps

    template<typename TCigarString, typename TGaps>
    inline void cigarToGapAnchorRead(TCigarString const & cigar, TGaps & gaps)
    {
		typename Iterator<TGaps>::Type it = begin(gaps);
		for (unsigned i = 0; i < length(cigar); ++i)
		{
			switch (cigar[i].operation)
			{
				case 'D':
				case 'N':
				case 'P':
					insertGaps(it, cigar[i].count);
				case 'I':
				case 'M':
        case 'S':
					it += cigar[i].count;
			}
		}
	}

    template<typename TCigarString, typename TGaps>
    inline void cigarToGapAnchorContig(TCigarString const & cigar, TGaps & gaps)
    {
		typename Iterator<TGaps>::Type it = begin(gaps);
		for (unsigned i = 0; i < length(cigar); ++i)
		{
			switch (cigar[i].operation)
			{
				case 'I':
				case 'P':
					insertGaps(it, cigar[i].count);
				case 'D':
				case 'M':
				case 'N':
        case 'S':
					it += cigar[i].count;
			}
		}
	}



//////////////////////////////////////////////////////////////////////////////
// Parsing functions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// _parseReadCigar
    
    template <typename TFile, typename TCigarString, typename TChar>
    inline void
    _parseReadCigar(TFile & file, TCigarString & cigar, TChar & c)
    {
		typedef typename Value<TCigarString>::Type	TCigarElement;
		typedef typename TCigarElement::TOperation	TOperation;
		typedef typename TCigarElement::TCount		TCount;

		clear(cigar);
		
        // if the CIGAR is not set and '*'
        if (c == '*')
		{
            c = _streamGet(file);
            return;
        }
        
        while (!_streamEOF(file)) 
		{
            TCount count = _parseReadNumber(file, c);
			if (c >= 'a' && c <= 'z')
				c = c + 'A' - 'a';
            appendValue(cigar, TCigarElement(c, count));
            
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') break;
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadSamIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadSamIdentifier(TFile & file, TString & str, TChar& c)
    {
        if (c == ' ' || c == '\t' || c == '\n') return;
        appendValue(str, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') return;
            appendValue(str, c);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseIsDna
    
    template<typename TChar>
    inline bool
    _parseIsDna(TChar const & c)
    {
        char x = TChar(Dna5(c));
        return (c == x) || (c + 'A' - 'a' == x);
    }
    
//////////////////////////////////////////////////////////////////////////////
//_parseReadDnaSeq
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parseReadDnaSeq(TFile & file, TString & str, TChar & c)
    {
		TChar first = c;
		if (!_streamEOF(file)) 
			c = _streamGet(file);

        if (!_parseIsDna(first))
			return;
        appendValue(str, first, Generous());
        
        for (; !_streamEOF(file) && _parseIsDna(c); c = _streamGet(file))
            appendValue(str, c, Generous());
    }
        
//////////////////////////////////////////////////////////////////////////////
// _parseIsPhredQual
    
    template <typename TChar>
    inline bool
    _parseIsPhredQual(TChar c)
    {
        return c >= '!' && c <= '~';
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadSeqQual
//
    
    template<typename TFile, typename TQualString, typename TChar>
    inline void
    _parseReadSeqQual(TFile & file, TQualString & str, TChar & c)
    {
        typedef typename Size<TQualString>::Type				TSize;
        typedef typename Iterator<TQualString, Standard>::Type	TIter;
        
        TIter itBegin = begin(str, Standard());
        TSize rest = length(str);
        
        int q = 0;
        for (TIter it = itBegin; rest != 0 && _parseIsPhredQual(c); --rest, ++it)
        {
            q = c - '!';
			if (!_streamEOF(file)) 
				c = _streamGet(file);
			else
				if (rest > 1)
					rest = 1;
			
			if (q == '*' - '!' && !_parseIsPhredQual(c) && it == itBegin)
				return;
			
            assignQualityValue(*it, q);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parseReadCharsUntilEndOfLine
//
// Reads all symbols till the next '\n' and writes them in the CharString str
// the c is the first character after the '\n'.
    
    template<typename TFile, typename TChar>
    inline void
    _parseReadCharsUntilEndOfLine(TFile & file, String<char> & str, TChar& c)
    {
        // read all chars till '\n'
        while (c != '\n')
        {
            appendValue(str, c, Generous());
			if (_streamEOF(file)) return;
	        c = _streamGet(file);
        }
        
        // read the first char after the '\n'
		if (!_streamEOF(file))
	        c = _streamGet(file);
    }

//////////////////////////////////////////////////////////////////////////////
// appendAlignment
    
    template<typename TSpec, typename TConfig, typename TId, typename TPos, typename TGaps>
    inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
	appendAlignment(
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId readId, 
		TId contigId, 
		TPos beginPos, 
		TPos endPos, 
		TGaps const & gaps)
	{
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        
        TId id = length(fragStore.alignedReadStore);
        TAlignedElement alignedElem = TAlignedElement(id, readId, contigId, beginPos, endPos, gaps);
        appendValue(fragStore.alignedReadStore, alignedElem);
		
		return id;
    }
    


//////////////////////////////////////////////////////////////////////////////
// read functions for Sam
//////////////////////////////////////////////////////////////////////////////

    
//////////////////////////////////////////////////////////////////////////////
// _generatePairMatchIds
//
	template <typename TPos, typename TId>
	struct MatchMateInfo_
	{
		TId		readId;
		TId		contigId;
		TId		pairMatchId;
		TPos	beginPos;
	};
    
    struct MatchMateInfoLess_
	{
        template <typename TMInfo>
        inline bool 
        operator() (TMInfo const &a, TMInfo const &b) const 
		{
            return (a.contigId < b.contigId) || (a.contigId == b.contigId && a.beginPos < b.beginPos);
        }
    };

    template<typename TSpec, typename TConfig, typename TMatchMateInfos>
    inline void 
	_generatePairMatchIds (
		FragmentStore<TSpec, TConfig> & fragStore,
		TMatchMateInfos & matchMateInfos)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        
        typedef typename TFragmentStore::TReadStore			TReadStore;
        typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
        typedef typename TFragmentStore::TContigPos			TContigPos;

        typedef typename Value<TReadStore>::Type						TRead;
        typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIter;    
		typedef typename Iterator<TMatchMateInfos, Standard>::Type		TMIter;    
        typedef typename Id<TFragmentStore>::Type						TId;
                
        TIter it = begin(fragStore.alignedReadStore, Standard());
		TIter itEnd = end(fragStore.alignedReadStore, Standard());
		TMIter mit = begin(matchMateInfos, Standard());
		TMIter mitEnd = end(matchMateInfos, Standard());

		if (it == itEnd || mit == mitEnd) return;

        // sort the aligned read store by: begin position, contig name
        sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());
        sortAlignedReads(fragStore.alignedReadStore, SortContigId());
		std::sort(mit, mitEnd, MatchMateInfoLess_());

		TMIter mitNext = mit;
		while (mitNext != mitEnd && (*mit).beginPos == (*mitNext).beginPos)
			++mitNext;
		
		TContigPos pos = _min((*it).beginPos, (*it).endPos);
		typename Size<TReadStore>::Type readStoreLength = length(fragStore.readStore);

		while (mit != mitEnd)
		{
			int cmp = 0;
			if ((*it).contigId < (*mit).contigId) cmp = -1;
			else if ((*it).contigId > (*mit).contigId) cmp = 1;
			else if (pos < (*mit).beginPos) cmp = -1;
			else if (pos > (*mit).beginPos) cmp = 1;

			if (cmp == 1)
			{
				mit = mitNext;
				while (mitNext != mitEnd && (*mit).beginPos == (*mitNext).beginPos)
					++mitNext;
				continue;
			}
			if (cmp == 0)
			{
				for (TMIter m = mit; m != mitNext; ++m)
					if ((*it).readId != (*m).readId && (*it).readId < readStoreLength && (*m).readId < readStoreLength)		// correct position found
						if (fragStore.readStore[(*it).readId].matePairId == fragStore.readStore[(*m).readId].matePairId)	// correct mate found
							(*it).pairMatchId = (*m).pairMatchId;															// link mate
			}
			if (++it == itEnd) break;
			pos = _min((*it).beginPos, (*it).endPos);
		}
    }    

//////////////////////////////////////////////////////////////////////////////
// read

///.Function.read.param.tag.type:Tag.File Format.tag.Sam
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void 
    read(
		TFile & file,
		FragmentStore<TSpec, TConfig> & fragStore,
		Sam)
    {
        typedef Value<FILE>::Type TValue;
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
		typedef typename TFragmentStore::TContigPos TContigPos;
        typedef typename Id<TFragmentStore>::Type TId;
        
        // data structure to temporarily store the gaps that need to be inserted in the contig sequences
        typedef MatchMateInfo_<TContigPos, TId> TMatchMateInfo;
        typedef String<TMatchMateInfo> TMatchMateInfos;
        typedef StringSet<String<typename TFragmentStore::TContigGapAnchor> > TContigAnchorGaps;

        // data structure to temporarily store information about match mates
        TMatchMateInfos matchMateInfos;
		TContigAnchorGaps contigAnchorGaps;
        
        if (_streamEOF(file)) return;

		// get first character from the stream
        char c = _streamGet(file);
        
        // Read in header section
        _readHeader(file, fragStore, c, Sam());
        
        // Read in alignments section
        _readAlignments(file, fragStore, contigAnchorGaps, matchMateInfos, c, Sam());
        
        // set the match mate IDs using the information stored in matchMateInfos
        _generatePairMatchIds(fragStore, matchMateInfos);
        
		convertPairWiseToGlobalAlignment(fragStore, contigAnchorGaps);
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readHeader

    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readHeader (
		TFile & file,
		FragmentStore<TSpec, TConfig> &,
		TChar & c,
		Sam)
    {
        // skip header for now
        while (c == '@')
            _parseSkipLine(file, c);
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readAlignments
//
// reads in alignement sections from a Sam file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TChar>
    inline void 
    _readAlignments (
		TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        TContigAnchorGaps & contigAnchorGaps,	
        TMatchMateInfos & matchMateInfos,
        TChar & c,
        Sam)
    {
        // create dummy entries in Sam specific aligned read quality store and aligned read tag store
        // is needed so the ID in the aligned store can be use to access the other stores
        // even if there exists previous entries without
		typedef FragmentStore<TSpec, TConfig> TFragmentStore;
		typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
		typedef typename TFragmentStore::TNameStore TNameStore;
		typedef typename Value<TAlignQualityStore>::Type TAlignQuality;
		
		TAlignQuality q;
		q.score = maxValue(q.score);
        int diff = length(fragStore.alignedReadStore) - length(fragStore.alignQualityStore);
        for(int i = 0; i < diff; ++i)
            appendValue(fragStore.alignQualityStore, q, Generous());
		
        diff = length(fragStore.alignedReadStore) - length(fragStore.alignedReadTagStore);
        for(int i = 0; i < diff; ++i)
            appendValue(fragStore.alignedReadTagStore, "", Generous());
        
        // read in alignments
		Nothing contextSAM;
		refresh(fragStore.contigNameStoreCache);
		refresh(fragStore.readNameStoreCache);

        while (!_streamEOF(file))
            _readOneAlignment(file, fragStore, contigAnchorGaps, matchMateInfos, c, Sam(), contextSAM);
    }
    
    
//////////////////////////////////////////////////////////////////////////////
// _readOneAlignment
//
// reads in one alignement section from a Sam file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TChar, typename TContextSAM>
    inline void 
    _readOneAlignment (
		TFile & file,
		FragmentStore<TSpec, TConfig> & fragStore,
		TContigAnchorGaps & contigAnchorGaps,
		TMatchMateInfos & matchMateInfos,
		TChar & c,
		Sam,
		TContextSAM & contextSAM)
    {
        // Basic types
        typedef FragmentStore<TSpec, TConfig>										TFragmentStore;
        typedef typename Id<TFragmentStore>::Type									TId;
        typedef typename Size<TFragmentStore>::Type									TSize;
        
        // All fragment store element types
        typedef typename Value<typename TFragmentStore::TContigStore>::Type			TContigElement;
        typedef typename Value<typename TFragmentStore::TLibraryStore>::Type		TLibraryStoreElement;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type		TMatePairElement;
        typedef typename Value<typename TFragmentStore::TReadStore>::Type			TReadStoreElement;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
        typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type	TAlignQualityElement;
        typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;
		
        // Type for sequence in readstore
        typedef typename TFragmentStore::TReadSeq TReadSeq2;
        
        // Type for gap anchor
        typedef typename TFragmentStore::TContigPos									TContigPos;
		typedef Gaps<TReadSeq2, AnchorGaps<TReadGapAnchors> >						TReadGaps;
		typedef Gaps<Nothing, AnchorGaps<typename Value<TContigAnchorGaps>::Type> >	TContigGapsPW;
        
        // Type to temporarily store information about match mates
        typedef typename Value<TMatchMateInfos>::Type								TMatchMateInfo;
        
        // read fields of alignments line        
        _parseSkipWhitespace(file, c);

		// Read the query name.  The letters until the first
		// whitespace will be read into qname.  Then, we skip until we
		// hit the first tab character.
        String<char> qname;
        _parseReadSamIdentifier(file, qname, c);
        _parseSkipUntilChar(file, '\t', c);

		// read the flag
        int flag;
        flag = _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);
		bool reverse = (flag & (1 << 4)) == (1 << 4);

		// Read reference name.  Same behaviour as for query name:  Read up to
        // the first whitespace character and skip to next tab char.
        String<char> rname;
        _parseReadSamIdentifier(file, rname, c);
        _parseSkipUntilChar(file, '\t', c);

		// read begin position
        TContigPos beginPos;
        beginPos = _parseReadNumber(file, c);
        --beginPos; // Sam stores positions starting at 1 the fragment store starting at 0
        _parseSkipWhitespace(file, c);

        // read map quality
        TAlignQualityElement mapQ;
        mapQ.score = _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);

		// read CIGAR
        String<CigarElement<> > cigar;
        _parseReadCigar(file, cigar, c);
        _parseSkipWhitespace(file, c);
        
        // calculate the end position
        TContigPos endPos;
        _getClippedLength(cigar, endPos);
        endPos = beginPos + endPos;
        // if the read is on the antisense strand switch begin and end position
        if (reverse)
		{
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }

		// generate gap anchor string for the read
        TReadGapAnchors readGapAnchors;
        
        // read mate reference name
        String<char> mrnm;
        _parseReadSamIdentifier(file, mrnm, c);
        _parseSkipWhitespace(file, c);

		// read mate position
        TContigPos mPos;
        mPos = _parseReadNumber(file, c);
        --mPos; // Sam stores positions starting at 1 the fragment store starting at 0
        _parseSkipWhitespace(file, c);

		// read iSize
        _parseReadNumber(file, c);
        _parseSkipWhitespace(file, c);

		// read in sequence
        TReadSeq2 readSeq;
        _parseReadDnaSeq(file, readSeq, c);
        SEQAN_ASSERT_GT(length(readSeq), 0u);
		if (reverse)  // TODO(holtgrew): Why not after parsing quality? ib. for BAM...
			reverseComplement(readSeq);
        _parseSkipWhitespace(file, c);

		// and associated qualities
        _parseReadSeqQual(file, readSeq, c);

		// insert alignment gaps
		TReadGaps readGaps(readSeq, readGapAnchors);
        cigarToGapAnchorRead(cigar, readGaps);
        
        // read in Sam tags
        String<char> tags;
        _parseSkipSpace(file, c);
        _parseReadCharsUntilEndOfLine(file, tags, c);
		
		if (empty(qname) || empty(rname))
			return;
        
        // check if read sequence is already in the store.
        // if so get the ID, otherwise create new entries in the
        // read, read name and mate pair store
        
        TId readId = 0;
        _storeAppendRead(fragStore, readId, qname, readSeq, flag, contextSAM);
        
        // check if the contig is already in the store
        // get its ID or create a new one otherwise
        TId contigId = 0;
        _storeAppendContig(fragStore, contigId, rname);

		if (empty(cigar)) return;
		
        // create a new entry in the aligned read store
        TId pairMatchId = appendAlignment(fragStore, readId, contigId, beginPos, endPos, readGapAnchors);
		resize(contigAnchorGaps, length(fragStore.alignedReadStore), Generous());
		TContigGapsPW contigGaps(back(contigAnchorGaps));
        cigarToGapAnchorContig(cigar, contigGaps);
		
		// create entries in Sam specific stores
        appendValue(fragStore.alignQualityStore, mapQ, Generous());
        appendValue(fragStore.alignedReadTagStore, tags, Generous());
        
        // store additional data about match mate temporarily
        // used in the end of the read function to generate match mate IDs
		TId mcontigId = contigId;
        if (mrnm != "*")
		{
			if (mrnm != "=")
				_storeAppendContig(fragStore, mcontigId, mrnm);

			if (flag & 0x40)	// store mate info only for the first read in the pair
			{
				TMatchMateInfo matchMateInfo = {readId, mcontigId, pairMatchId, mPos};
				appendValue(matchMateInfos, matchMateInfo);
				back(fragStore.alignedReadStore).pairMatchId = pairMatchId;
			}
		}
    }


//////////////////////////////////////////////////////////////////////////////
// write functions for Sam
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TSpec, typename TConfig>
	inline void _writeHeader(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 Sam)
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
                                 Sam)
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
		typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;

		String<int> mateIndex;	// store outer library size for each pair match (indexed by pairMatchId)
		calculateMateIndices(mateIndex, store);
		
        TAlignIter it = begin(store.alignedReadStore, Standard());
        TAlignIter itEnd = end(store.alignedReadStore, Standard());
		TAlignIter mit = it;
		CharString cigar;
		TReadSeq readSeq;

//        int i = 0;
        for(; it != itEnd; ++it)
		{
//            std::cerr << "alignment " << i++ << std::endl;
            TId alignedId = (*it).id;
			TId readId = (*it).readId;
			TId mateIdx = TRead::INVALID_ID;

			if ((*it).pairMatchId != TRead::INVALID_ID)
				mateIdx = mateIndex[2*(*it).pairMatchId + getMateNo(store, (*it).readId)];

			TContigGaps	contigGaps(/*store.contigStore[(*it).contigId].seq, */store.contigStore[(*it).contigId].gaps);
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
			else
				flag |= 0x0008;					// mate is unmapped (actually we should check if the mate has no match at all)
			
			signed char mateNo = getMateNo(store, readId);
			if (mateNo == 0) flag |= 0x0040;	// this read is the first in the pair
			if (mateNo == 1) flag |= 0x0080;	// this read is the second in the pair

			if (readId < length(store.readStore))
			{
				TRead &read = store.readStore[readId];
				if (read.matePairId != TRead::INVALID_ID)
					flag |= 0x0001;
			}
			
			// <qname>
			if (readId < length(store.readNameStore)) {
                typedef typename Iterator<CharString, Standard>::Type TCharStringIterator;
                for (TCharStringIterator it = begin(store.readNameStore[readId]); it != end(store.readNameStore[readId]); ++it) {
                    if (*it == ' ' || *it == '\t' || *it == '\n' || *it == '\r')
                        break;
                    _streamPut(target, *it);
                }
            }
            _streamPut(target, '\t');
            
            // <flag>
            _streamPutInt(target, flag);
            _streamPut(target, '\t');
            
			// <rname>
			if ((*it).contigId < length(store.contigNameStore))
				_streamWrite(target, store.contigNameStore[(*it).contigId]);
            else
                _streamWrite(target, '.');  // No reference name given.  Standard says field must not be empty but gives no "NULL" value.
            _streamPut(target, '\t');
            
			// <pos>
            _streamPutInt(target, pos);
            _streamPut(target, '\t');
            
			// <mapq>
			if (alignedId < length(store.alignQualityStore))
				_streamPutInt(target, store.alignQualityStore[alignedId].score);
            else
                _streamPutInt(target, 255);
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
					reverseComplement(readSeq);
				}
			} else
				clear(readSeq);
			
            // <cigar>
			TReadGaps readGaps(readSeq, (*it).gaps);
            typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps2;
			// TContigGaps	contigGaps2(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
            // if (i == 4)
            //     printf("It's it!\n");
            // std::cerr << "read gaps:  " << readGaps << std::endl;
            // std::cerr << "contig gaps:" << contigGaps2 << std::endl;
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

///.Function.write.param.tag.type:Tag.File Format.tag.Sam
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void write(TFile & target,
                      FragmentStore<TSpec, TConfig> & store,
                      Sam)
    {
        // write header
		_writeHeader(target, store, Sam());
        
        // write aligments
        _writeAlignments(target, store, Sam());
    }
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
