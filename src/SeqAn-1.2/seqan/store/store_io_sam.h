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

#ifndef SEQAN_HEADER_STORE_IO_SAM_H
#define SEQAN_HEADER_STORE_IO_SAM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.SAM:
	SAM mapping file.
*/
//struct TagSAM_;
//typedef Tag<TagSAM_> const SAM;

    
//////////////////////////////////////////////////////////////////////////////
// Cigar struct
//////////////////////////////////////////////////////////////////////////////
    
    template<typename TChar = char, typename TNum = unsigned int>
    struct Cigar
    {
        String<TChar> operationType;
        String<TNum> operationCount;
        
        Cigar(){}
    };

//////////////////////////////////////////////////////////////////////////////
// append
    
    template<typename TChar, typename TNum>
    inline void append(Cigar<TChar, TNum> & cigar, TChar const c, TNum const i)
    {
        append(cigar.operationType, c);
        append(cigar.operationCount, i);
    }

//////////////////////////////////////////////////////////////////////////////
// _getClippedLength
    
    template<typename TChar, typename TNum, typename TNum2>
    inline void _getClippedLength(Cigar<TChar, TNum> const & cigar, TNum2 & sum)
    {
        typedef typename Iterator< String<TNum> >::Type TNumIter;
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        
        TNumIter it_c = begin(cigar.operationCount);
        TCharIter it_t = begin(cigar.operationType);
        
        sum = 0;
        
        while(it_c != end(cigar.operationCount)){
            if(value(it_t) != 'S' && value(it_t) != 's' && value(it_t) != 'H' && value(it_t) != 'h'){
                sum += value(it_c);
            }
            
            ++it_c;
            ++it_t;
        }

    }

//////////////////////////////////////////////////////////////////////////////
// cigarToGapAnchorRead
    
    template<typename TChar, typename TNum, typename TGapAnchor>
    inline void cigarToGapAnchorRead(Cigar<TChar, TNum> const & cigar, String<TGapAnchor> & gaps)
    {
        typedef typename Iterator< String<TNum> >::Type TNumIter;
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        
        // Iterators on cigar
        TNumIter it_c = begin(cigar.operationCount);
        TCharIter it_t = begin(cigar.operationType);
        
        // delete information in gap anchor string
        gaps = String<TGapAnchor>();
        
        // boolean that keeps track on type of the last cigar operation
        // is true if match/mismatch or insertion
        // if false if deletion, padding or clipped
        bool was_mi = true;
        
        // positions in the un-gapped and gapped sequence
        int pos_seq = 0;
        int pos_gapped = 0;
        
        // if the CIGAR is empty
        if(length(cigar.operationCount) == 0){
            return;
        }
        
        // check if there is a (soft) clipped sequence in the begining
        if(value(it_t) == 'S' || value(it_t) == 's'){
            
            // skip first 'count' characters
            pos_seq += value(it_c);
            
            // insert in GapAnchor
            TGapAnchor anchor = TGapAnchor(pos_seq, pos_gapped);
            append(gaps, anchor);
        }
        
        while(it_c != end(cigar.operationCount)){
            
            // If the operation type is match/mismatch or a Deletion in the reference there is no gap in the read sequence.
            if(value(it_t) == 'M' || value(it_t) == 'D' || value(it_t) == 'm' || value(it_t) == 'd'){
                
                // If this is the first m/i operation after a d/p operation (first one after a gap)
                if(!was_mi){
                    
                    // insert in gap anchor
                    TGapAnchor anchor = TGapAnchor(pos_seq, pos_gapped);
                    append(gaps, anchor);
                    
                    // switch operation type to match/mismatch or insertion
                    was_mi = true;
                }
                
                // Increment positions in the gapped and un-gapped sequnece
                pos_gapped += value(it_c);
                pos_seq += value(it_c);
            }
            
            // If the operation type is insertion or skipped in the reference or padding there is a gap in the read sequence.
            if(value(it_t) == 'I' || value(it_t) == 'P' || value(it_t) == 'N' || value(it_t) == 'i' || value(it_t) == 'p' || value(it_t) == 'n'){

                // switch operation type to deletion or padding
                was_mi = false;
                
                // Increment only the position in the gapped sequence
                pos_gapped += value(it_c);
            }
            
            // Iterate.
            ++it_c;
            ++it_t;
        }
        
        // following (soft) klipped characters are encode by the exceeding length of the sequence
    }
    //////////////////////////////////////////////////////////////////////////////
    // cigarToContigGaps
    
    template<typename TChar, typename TNum, typename TPos, typename TFlag, typename TGap>
    inline void cigarToContigGaps(Cigar<TChar, TNum> const & cigar, TPos beginPos, TFlag flag, String<Pair<TGap, TGap> > & gaps)
    {
        typedef typename Iterator< String<TNum> >::Type TNumIter;
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        typedef Pair<TGap, TGap> TPair;
        
        // Iterators on cigar
        TNumIter it_c = begin(cigar.operationCount);
        TCharIter it_t = begin(cigar.operationType);
        
        TNumIter stop = end(cigar.operationCount);
        
        int direction = 1;
        
        if (flag & (1 << 4) == (1 << 4)){
            it_c = end(cigar.operationCount);
            it_t = end(cigar.operationType);
            stop = begin(cigar.operationCount);
            direction = -1;
        }
        
        // boolean that keeps track on type of the last cigar operation
        // is true if it was one introducing a gap
        // if false otherwise
        bool was_gap = false;
        
        // positions and length of the current gap
        int gapLength = 0;
        int pos = beginPos;
        
        // if the CIGAR is empty
        if(length(cigar.operationCount) == 0){
            return;
        }
        
        // check if there is a (soft) clipped sequence in the begining
        if(value(it_t) == 'S' || value(it_t) == 's'){
            
            // skip first 'count' characters
            pos += value(it_c);
        }
        
        while(it_c != stop){
            
            // If the operation type is insertion in the reference, match/mismatch or skipped there is no gap in the read sequence.
            if(value(it_t) == 'I' || value(it_t) == 'M' || value(it_t) == 'N' || value(it_t) == 'i' || value(it_t) == 'm' || value(it_t) == 'n'){
                
                // If this is the first i/p/n operation after a m/d operation (first one after a gap)
                if(was_gap){
                    
                    // insert in gap pair
                    TPair pair;
                    pair.i1 = pos; pair.i2 = gapLength;
                    append(gaps, pair);
                    
                    // set gap length back to 0
                    gapLength = 0;
                    
                    // switch operation type to does not introduce gap
                    was_gap = false;
                }
                
                // Increment positions
                pos += value(it_c);
            }
            
            // If the operation type is deletion or padding there is a gap in the reference sequence.
            if(value(it_t) == 'P' || value(it_t) == 'D' || value(it_t) == 'p' || value(it_t) == 'd'){
                
                // switch operation type to introduc gap
                was_gap = true;
                
                // Increment only the position and gap length
                pos += value(it_c);
                gapLength  += value(it_c);
            }
            
            // Iterate.
            it_c += direction;
            it_t += direction;
        }
        
        // following (soft) klipped characters are encode by the exceeding length of the sequence
    }
    
//////////////////////////////////////////////////////////////////////////////
// some helping functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _print_gapAnchor
    
    template<typename TChar, typename TGapAnchor>
    inline void _print_gapAnchor(String<TChar> seq, String<TGapAnchor> gaps)
    {
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        typedef typename Iterator< String<TGapAnchor> >::Type TGapIter;
        
        // create iterators and set them on the begining of the strings
        TCharIter it_c = begin(seq);
        TGapIter  it_g = begin(gaps);
        
        //
        int gapped = value(it_g).gapPos, last_gapped = gapped;
        int ungapped = value(it_g).seqPos, last_ungapped = ungapped;
        
        for(int i = 0; i < ungapped; ++i)
            ++it_c;
        
        while(it_g != end(gaps)){
            gapped = value(it_g).gapPos;
            ungapped = value(it_g).seqPos;
            //std::cout << "(" << ungapped << ", " << gapped << ")";
            
            int count_chars = ungapped - last_ungapped;
            int count_gaps = (gapped - last_gapped) - count_chars;
            
            for(int i = 0; i < count_chars && it_c != end(seq); ++i, ++it_c){
                _streamPut(std::cout, value(it_c));
            }
            
            for(int i = 0; i < count_gaps; ++i){
                _streamPut(std::cout, '-');
            }
            
            // iterate
            last_gapped = gapped;
            last_ungapped = ungapped;
            ++it_g;
        }
        
        while(it_c != end(seq)){
            _streamPut(std::cout, value(it_c));
            ++it_c;
        }
    }

//////////////////////////////////////////////////////////////////////////////
// _getMate
//
// checks if a read name is already in read name store and if corresponding 
// begin position in the aligned read strore fits the mate position
// returns the ID or -1 otherwise
    
    template<typename TSpec, typename TConfig, typename TPos, typename TID>
    inline void 
    _getMate(FragmentStore<TSpec, TConfig> & fragStore, String<char> & readName, TPos & mPos, TID & mateID)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
        
        // set mate ID = -1. Will be replaced if mate is found
        mateID = -1;
        
        // Iterator over read names
        TNameStoreIter it_rm = begin(fragStore.readNameStore);
        
        while (it_rm != end(fragStore.readNameStore)){
            // if the read name was found
            if(readName == value(it_rm)){
                // check in the aligned read store if the begin position is correct
                int suggestedID = position(it_rm);
                
                if(value(fragStore.alignedReadStore, suggestedID).beginPos == mPos){
                    mateID = suggestedID;
                    return;
                }                
            }
            
            ++it_rm;
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// _getID
    
    template<typename TType, typename TID>
    inline bool 
    _getID(StringSet<TType> & store, TType & elem, TID & elem_id)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
        
        // Iterator over read names
        TNameStoreIter iter = begin(store);
        
        while (iter != end(store)){
            // if the element was found
            if(elem == value(iter)){
                // set the ID
                elem_id = position(iter);
                // and break the loop
                return true;               
            }
            
            ++iter;
        }
        return false;
    }
    
//////////////////////////////////////////////////////////////////////////////
// appendAlignment
    
    template<typename TSpec, typename TConfig, typename TId, typename TPos, typename TGaps>
    inline void
    appendAlignment(
                    FragmentStore<TSpec, TConfig> & fragStore, 
                    TId & readId, 
                    TId & contigId, 
                    TPos & beginPos, 
                    TPos & endPos, 
                    TGaps & gaps){
        
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        
        TId alignId = length(fragStore.alignedReadStore);
        
        TAlignedElement alignedElem = TAlignedElement();
        alignedElem.id = alignId;
        alignedElem.readId = readId;
        alignedElem.contigId = contigId;
        alignedElem.beginPos = beginPos;
        alignedElem.endPos = endPos;
        alignedElem.gaps = gaps;
        
        append(fragStore.alignedReadStore, alignedElem);
    }
    
//////////////////////////////////////////////////////////////////////////////
// _appendRead
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
    template<typename TSpec, typename TConfig, typename TId, typename TName, typename TString, typename TFlag>
    inline void
    _appendRead(FragmentStore<TSpec, TConfig> & fragStore, 
                TId & readId, 
                TName & qname,
                TString & readSeq,
                TFlag & flag){
        
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

        if(_getID(fragStore.readNameStore, qname, readId)){
            // if the read is paired
            if((flag & 1) == 1){
                // check the mate pair store if it is the same mate of the pair
                // assuming that only one flag 0x040 or 0x0080 is 1
                int inPair = ((flag & (1 << 7)) >> 7);
                
                TId matePairId = value(fragStore.readStore, readId).matePairId;
                
                readId = value(fragStore.matePairStore, matePairId).readId[inPair];
                
                if(readId == TMatePairElement::INVALID_ID){
                    // create new entry in read and read name store
                    readId = length(fragStore.readStore);
                    // set sequence and mate pair ID in new read store element
                    appendRead(fragStore, readSeq, matePairId);
                    
                    // add the identifyer to the read name store
                    appendValue(fragStore.readNameStore, qname);
                    
                    // set the ID in the mate pair store
                    value(fragStore.matePairStore, matePairId).readId[inPair] = readId;
                }
            }
        } else { // if the read name is not in the store
            // create new entry in read and read name store
            
            
            // if the read is paired
            if((flag & 1) == 1){
                TMatePairElement mateElem = TMatePairElement();
                // set the first or second read ID in the mate pair element
                readId = length(fragStore.readStore);
                mateElem.readId[((flag & (1 << 7)) >> 7)] = readId;
                // get a new mate pair ID and add the new mate pair element
                TId matePairId = length(fragStore.matePairStore);
                append(fragStore.matePairStore, mateElem);
                // set the new mate pair ID in the read element
                appendRead(fragStore, readSeq, matePairId);
            } 
            // if read is not paired
            else {
                appendRead(fragStore, readSeq);
            }
            
            appendValue(fragStore.readNameStore, qname);
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// _appendContig
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
    template<typename TSpec, typename TConfig, typename TId, typename TName, typename TContigGaps>
    inline void
    _appendContig(FragmentStore<TSpec, TConfig> & fragStore, 
                TId & contigId, 
                TName & rName,
                TContigGaps & contigGap)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        typedef typename Value<TContigGaps>::Type TContigGap;
        
        contigId = 0;
        
        // if the contig already exsists
        if(_getID(fragStore.contigNameStore, rName, contigId)){
            
        } 
        // if the contig is not in the store yet
        else {
            // set the ID on the last entry after appending
            contigId = length(fragStore.contigStore);
            // append contig store
            appendValue(fragStore.contigNameStore, rName);
            append(fragStore.contigStore, TContigElement());
            // append temp contig gaps
            appendValue(contigGap, TContigGap());
        }
    }

//////////////////////////////////////////////////////////////////////////////
// _generatePairMatchIds
//
// shift is the first ID for which mateInfo contains information about a aligned read
// the function generates new pair match IDs for these entries. It uses the ID of the first found mate for this.
    
    template<typename TSpec, typename TConfig, typename TMateInfo, typename TSize>
    inline void
    _generatePairMatchIds(FragmentStore<TSpec, TConfig> & fragStore,
                          TMateInfo mateInfos,
                          TSize shift)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        
        typedef typename TFragmentStore::TAlignedReadStore TAligned;
        typedef typename Iterator<TAligned>::Type TIter;
        
        typedef typename Id<TFragmentStore>::Type TId;
        typedef typename TFragmentStore::TContigPos TContigPos;
        
        typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        
        // sort the aligned read store by: begin position, contig name
        sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());
        
        // Iterate over all entries in the aligned read store
        // for each new entry (ID >= shift) check if it is paired
        // if search for the mate and get its ID
        TIter alignIt = begin(fragStore.alignedReadStore);
        
        for(; alignIt != end(fragStore.alignedReadStore); ++alignIt){
            if(value(alignIt).id >= shift){
                TId shiftedId = value(alignIt).id - shift;
                
                // if the aligned read is paired and the pair match ID is not set yet
                if(value(fragStore.readStore, shiftedId).matePairId != TReadStoreElement::INVALID_ID & value(alignIt).pairMatchId != TAlignedElement::INVALID_ID){
                    
                    // get mates position
                    TContigPos & mpos = value(mateInfos, shiftedId).i1;
                    TId & mrnm = value(mateInfos, shiftedId).i2;
                    
                    // search for the mate
                    TIter mateIt = lowerBoundAlignedReads(fragStore.alignedReadStore, mpos, SortBeginPos());
                    while(value(mateIt).beginPos == mpos & value(mateIt).contigId != mrnm) ++mateIt;
                    
                    if(value(mateIt).beginPos == mpos){
                        // set the pair match IDs
                        value(alignIt).pairMatchId = value(alignIt).id;
                        value(mateIt).pairMatchId = value(alignIt).id;
                    }
                    
                }
            }
        }
    }    

//////////////////////////////////////////////////////////////////////////////
// comparePosLengthPair
    
    struct comparePosLengthPair{
        template<typename TPos>
        inline bool 
        operator() (TPos const & first, TPos const & second) const {
            return (first.i1 < second.i1);
        }
    };
    
//////////////////////////////////////////////////////////////////////////////
// _writeContigGapInStore
    
    template<typename TFragmentStore, typename TContigPos>
    inline void
    _writeContigsGapsInStore(TFragmentStore & fragStore, StringSet<String<Pair<TContigPos, TContigPos> > > const & contigsGaps)
    {
        
        typedef Pair<TContigPos, TContigPos> TPosLengthPair;
        typedef String<TPosLengthPair> TContigGaps;
        typedef typename Iterator<TContigGaps>::Type TGapsIter;
        typedef StringSet<TContigGaps> TContigsGaps;
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigStoreElem;
        typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
        typedef typename Iterator<String<TContigGapAnchor> >::Type TContigGapAnchorIter;
        
        // for all contigs
        for(int i = 0; i < length(contigsGaps); ++i){
            TContigGaps gaps = value(contigsGaps, i);
            
            // Sort the contig gaps according to their start position
            sort(begin(gaps), end(gaps), comparePosLengthPair());
            
            // create iterator over contigGaps
            TGapsIter gapIt = begin(gaps);
            
            // get corresponding gapAnchor and iterator on it
            TContigStoreElem contigStoreElem = value(fragStore.contigStore, i);
            TContigGapAnchorIter gapAnchorIt = begin(contigStoreElem.gaps);
            int anchorCount = 0;
            
            TContigPos shift = 0;
            TContigPos pos = 0;
            TContigPos length = 0;
            
            // for all pos-length-pairs
            for(; gapIt != end(gaps); ++gapIt){
                pos = value(gapIt).i1;
                length = value(gapIt).i2;
                
                // while position is not reached yet and not last position
                for(; value(gapAnchorIt).gapPos < pos & gapAnchorIt != end(contigStoreElem.gaps); ++gapAnchorIt, ++anchorCount){
                    // gapped position += shift
                    value(gapAnchorIt).gapPos += shift;                    
                }
                
                // difference between position and gapped postion in last gap anchors
                TContigPos diff = pos - value(gapAnchorIt).gapPos;
                TContigGapAnchor anchor((value(gapAnchorIt).seqPos + diff), (pos + length));
                insertValue(contigStoreElem.gaps, anchorCount, anchor, Generous());
                ++gapAnchorIt; // because the insertion effects the iterator
                
                shift += length;
            }
            
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// parsing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _parse_readCigar
    
    template<typename TFile, typename TChar, typename TNum>
    inline void
    _parse_readCigar(TFile & file, Cigar<TChar, TNum> & cigar, TChar& c)
    {
        TChar type;
        TNum count;
        
        // if the CIGAR is not set and '*'
        if(c == '*'){
            c = _streamGet(file);
            return;
        }
        
        while (!_streamEOF(file)) {
            count = _parse_readNumber(file, c);
            type = c;
            append(cigar, type, count);
            
            c = _streamGet(file);
            if (c== ' ' || c== '\t' || c == '\n') break;
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSamIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readSamIdentifier(TFile & file, TString & str, TChar& c)
    {
        append(str, c);
        while (!_streamEOF(file)) {
            c = _streamGet(file);
            if (c== ' ' || c== '\t' || c == '\n') break;
            append(str, c);
        }

    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_is_dna
    
    template<typename TChar>
    inline bool
    _parse_is_dna(TChar const & c)
    {
        return ((c == 'a') || (c == 'c') || (c == 'g') || (c == 't') || (c == 'n') || (c == 'A') || (c == 'C') || (c == 'G') || (c == 'T') || (c == 'N'));
    }
    
//////////////////////////////////////////////////////////////////////////////
//_parse_readDnaSeq
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readDnaSeq(TFile & file, TString & str, TChar& c)
    {
        if (!_parse_is_dna(c)) return;
        append(str, c, Generous());
        
        while (!_streamEOF(file)) {
            c = _streamGet(file);
            if (!_parse_is_dna(c)) break;
            append(str, c, Generous());
        }
        
    }
        
//////////////////////////////////////////////////////////////////////////////
// _parse_is_PhredQual
    
    template<typename TChar>
    inline bool
    _parse_is_PhredQual(TChar const & c)
    {
        return ( ((unsigned) c > 32) && ((unsigned) c < 127) );
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSeqQual
//
// As the DNA5Q data structure can only store quality values till 40, all
// qualities over this threshold are changed to 40
    
    template<typename TFile, typename TChar>
    inline void
    _parse_readSeqQual(TFile & file, String<Dna5Q> & str, TChar& c)
    {
        typedef typename Iterator< String<Dna5Q> >::Type TIter;
        
        TIter it = begin(str);
        int q = 0;
        
        do
        {
            if (!_parse_is_PhredQual(c)) break;
            
            q = c - 32;
            if (q > 39) q = 40;
            assignQualityValue(value(it), q);
            
            c = _streamGet(file);
        }
        while (!_streamEOF(file) && it != end(str));
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readCharsTillEndOfLine
//
// Reads all symbols till the next '\n' and writes them in the CharString str
// the c is the first character after the '\n'.
    
    template<typename TFile, typename TChar>
    inline void
    _parse_readCharsTillEndOfLine(TFile & file, String<char> & str, TChar& c)
    {
        // read all chars till '\n'
        do
        {
            if (c == '\n') break;
            append(str, c, Generous());
            c = _streamGet(file);
        }
        while (!_streamEOF(file));
        
        // read the first char after the '\n'
        c = _streamGet(file);
    }
    
//////////////////////////////////////////////////////////////////////////////
// read functions for SAM
//////////////////////////////////////////////////////////////////////////////
    
//////////////////////////////////////////////////////////////////////////////
// read
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void 
    read(TFile& file,
         FragmentStore<TSpec, TConfig>& fragStore,
         SAM)
    {
        typedef Value<FILE>::Type TValue;
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        
        // data structure to temporarily store the gaps that need to be inserted in the contig sequences
        typedef Pair<typename TFragmentStore::TContigPos, typename TFragmentStore::TContigPos> TPosLengthPair;
        typedef String<TPosLengthPair> TContigGaps;
        typedef StringSet<TContigGaps> TContigsGaps;
        TContigsGaps contigsGaps;
        resize(contigsGaps, length(fragStore.contigStore));
        
        // data structure to temporarily store information about match mates
        typedef typename Id<TFragmentStore>::Type TId;
        typedef Pair<typename TFragmentStore::TContigPos, TId> TMatchMateInfo;
        typedef String<TMatchMateInfo> TMatchMateInfos;
        TMatchMateInfos matchMateInfos;
        
        // if the aligned store already contains data sets
        // save their number so the temporarily saved information
        // can be assiciated with the correct entry
        typename TFragmentStore::TContigPos shift = length(fragStore.alignedReadStore);
                
        if (!_streamEOF(file)){
            // get first character from the stream
            char c = _streamGet(file);
            
            // Read in header section
            _readHeader(file, fragStore, c, SAM());
            
            // Read in alignments section
            _readAlignments(file, fragStore, contigsGaps, matchMateInfos, c, SAM());
            
            // set the match mate IDs using the information stored in matchMateInfos
            _generatePairMatchIds(fragStore, matchMateInfos, shift);
            
            // insert gaps in the contigs using the information stored in contigsGaps
            _writeContigsGapsInStore(fragStore, contigsGaps);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readHeader

    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readHeader(TFile& file,
                FragmentStore<TSpec, TConfig>& fragStore,
                TChar & c,
                SAM)
    {
        // skip header for now
        while(c == '@'){
            _parse_skipLine(file, c);
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readAlignments
//
// reads in alignement sections from a SAM file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigGaps, typename TMateInfo, typename TChar>
    inline void 
    _readAlignments(TFile& file,
        FragmentStore<TSpec, TConfig>& fragStore,
        TContigGaps & contigGaps,
        TMateInfo & mateInfos,
        TChar & c,
        SAM)
    {
        // create dummy entries in SAM specific aligned read quality store and aligned read tag store
        // is needed so the ID in the aligned store can be use to access the other stores
        // even if there exists previous entries without
		typedef FragmentStore<TSpec, TConfig> TFragmentStore;
		typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
		typedef typename Value<TAlignQualityStore>::Type TAlignQuality;
		
        int diff = length(fragStore.alignedReadStore) - length(fragStore.alignQualityStore);
        for(int i = 0; i < diff; ++i){
			TAlignQuality q;
			q.score = 255;
            append(fragStore.alignQualityStore, q, Generous());
        }
        diff = length(fragStore.alignedReadStore) - length(fragStore.alignedReadTagStore);
        for(int i = 0; i < diff; ++i){
            appendValue(fragStore.alignedReadTagStore, "", Generous());
        }
        
        // read in alignments
        int k = 0;
        while(!_streamEOF(file)){
            std::cout << k << std::endl;
            _readOneAlignment(file, fragStore, contigGaps, mateInfos, c, SAM());
            ++k;
        }
    }
    
    
//////////////////////////////////////////////////////////////////////////////
// _readOneAlignment
//
// reads in one alignement section from a SAM file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigGaps, typename TMateInfo, typename TChar>
    inline void 
    _readOneAlignment(TFile& file,
                    FragmentStore<TSpec, TConfig>& fragStore,
                    TContigGaps & contigsGaps,
                    TMateInfo & mateInfos,
                    TChar & c,
                    SAM)
    {
        SEQAN_CHECKPOINT
        // Basic types
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Id<TFragmentStore>::Type TId;
        typedef typename Size<TFragmentStore>::Type TSize;
        
        // All fragment store element types
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
        typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignQualityElement;
        
        // Type for sequence in readstore
        typedef typename TFragmentStore::TReadSeq TReadSeq2;
        
        // Type for gap anchor
        typedef typename TFragmentStore::TReadGapAnchor TReadGapAnchor2;
        typedef typename TFragmentStore::TContigPos TContigPos;
        
        // Types to temporarily store the gaps that need to be inserted in the contig sequences
        typedef Pair<TContigPos, TContigPos> TPosLengthPair;
        
        // Type to temporarily store information about match mates
        typedef Pair<TContigPos, TId> TMatchMateInfo;
        
        // read fields of alignments line
        
        // read teh query name
        String<char> qname;
        _parse_readSamIdentifier(file, qname, c);
        _parse_skipWhitespace(file, c);
        std::cout << "qname: \t" << qname << std::endl;
        // read the flag
        int flag;
        flag = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        std::cout << "flag: \t" << flag << ":"<< (flag & 1) << (flag & (1 << 4))  << (flag & (1 << 7)) << std::endl;
        // read reference name
        String<char> rname;
        _parse_readSamIdentifier(file, rname, c);
        _parse_skipWhitespace(file, c);
        std::cout << "rname: \t" << rname << std::endl;
        // read begin position
        TContigPos beginPos;
        beginPos = _parse_readNumber(file, c);
        --beginPos; // SAM stores positions starting at 1 the fragment store starting at 0
        _parse_skipWhitespace(file, c);
        std::cout << "pos: \t" << beginPos << std::endl;
                
        // read map quality
        TAlignQualityElement mapQ;
        mapQ.score = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        std::cout << "mapQ: \t" << mapQ.score << std::endl;
        // read CIGAR
        Cigar<> cigar = Cigar<>();
        _parse_readCigar(file, cigar, c);
        _parse_skipWhitespace(file, c);
        
        // calculate the end position
        TContigPos endPos;
        _getClippedLength(cigar, endPos);
        endPos = beginPos + endPos;
        // if the read is on the antisense strand switch begin and end position
        if((flag & (1 << 4)) == (1 << 4)){
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }
        std::cout << "end pos:" << endPos << std::endl;
        // generate gap anchor string for the read
        String<TReadGapAnchor2 > readGaps;
        cigarToGapAnchorRead(cigar, readGaps);
        
        // read mate reference name
        String<char> mrnm;
        _parse_readSamIdentifier(file, mrnm, c);
        _parse_skipWhitespace(file, c);
        std::cout << "mrnm: \t" << mrnm << std::endl;
        // read mate position
        TContigPos mPos;
        mPos = _parse_readNumber(file, c);
        --mPos; // SAM stores positions starting at 1 the fragment store starting at 0
        _parse_skipWhitespace(file, c);
        std::cout << "mPos: \t" << mPos << std::endl;
        // read iSizs
        int iSize;
        iSize = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        std::cout << "iSize: \t" << iSize << std::endl;
        // read in sequence
        TReadSeq2 readSeq;
        _parse_readDnaSeq(file, readSeq, c);
        _parse_skipWhitespace(file, c);
        std::cout << "seq: \t" << readSeq << std::endl;
        // and associated qualities
        _parse_readSeqQual(file, readSeq, c);
        
        // read in SAM tags
        String<char> tags;
        _parse_readCharsTillEndOfLine(file, tags, c);
        
        
        // check if read sequence is already in the store.
        // if so get the ID, otherwise create new entries in the
        // read, read name and mate pair store
        
        TId readId = 0;
        _appendRead(fragStore, readId, qname, readSeq, flag);
        
        // check if the contig is already in the store
        // get its ID or create a new one otherwise
        TId contigId = 0;
        _appendContig(fragStore, contigId, rname, contigsGaps);

        // insert gaps in the contigGaps
        cigarToContigGaps(cigar, beginPos, flag, value(contigsGaps, contigId));
        
        // create a new entry in the aligned read store
        appendAlignment(fragStore, readId, contigId, beginPos, endPos, readGaps);
        
        // create entries in SAM specific stores
        append(fragStore.alignQualityStore, mapQ, Generous());
        appendValue(fragStore.alignedReadTagStore, tags, Generous());
        
        // store additional data about match mate temporarily
        // used in the end of the read function to generate match mate IDs
        TMatchMateInfo mateInfo;
        mateInfo.i1 = mPos;
        if(mrnm != "="){
            _appendContig(fragStore, contigId, mrnm, contigsGaps);
        }
        mateInfo.i2 = contigId;
        append(mateInfos, mateInfo);
        
        std::cout << "======================================" << std::endl;
        
    }
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
