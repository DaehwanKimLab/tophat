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

#ifndef SEQAN_HEADER_STORE_IO_H
#define SEQAN_HEADER_STORE_IO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Amos message file:
	Amos message file.
*/
struct TagAmos_;
typedef Tag<TagAmos_> const Amos;


//////////////////////////////////////////////////////////////////////////////
// Auxillary functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



template <typename TSpec, typename TConfig, typename TPos, typename TGapAnchor, typename TSpecAlign, typename TBeginClr, typename TEndClr>
inline void
getClrRange(FragmentStore<TSpec, TConfig> const& fragStore,
			AlignedReadStoreElement<TPos, TGapAnchor, TSpecAlign> const& alignEl,
			TBeginClr& begClr,		// Out-parameter: left / begin position of the clear range
			TEndClr& endClr)		// Out-parameter: right / end position of the clear range
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Iterator<String<TGapAnchor>, Standard>::Type TGapIter;
	
	TSize lenRead = length(fragStore.readSeqStore[alignEl.readId]);
	TGapIter itGap = begin(alignEl.gaps, Standard() );
	TGapIter itGapEnd = end(alignEl.gaps, Standard() );
	
	// Any gaps or clipped characters?
	if (itGap == itGapEnd) {
		begClr = 0;
		endClr = lenRead;
	} else {
		// Begin clear range
		begClr = (itGap->gapPos == 0) ? itGap->seqPos : 0;
		// End clear range
		--itGapEnd;
		if (itGapEnd->seqPos != lenRead) endClr = lenRead;
		else {
			int diff = (itGap != itGapEnd) ? (*(itGapEnd - 1)).gapPos - (*(itGapEnd-1)).seqPos : 0;
			int newDiff = itGapEnd->gapPos - itGapEnd->seqPos;
			endClr = (newDiff < diff) ? lenRead - (diff - newDiff) : lenRead;	
		}
	}

	// For reverse reads adapt clear ranges
	if (alignEl.beginPos > alignEl.endPos) {
		TBeginClr tmp = begClr;
		begClr = lenRead - endClr;
		endClr = lenRead - tmp;
	}
}



//////////////////////////////////////////////////////////////////////////////
// Read / Write of AMOS message files (*.afg)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig>
inline void 
read(TFile & file,
	 FragmentStore<TSpec, TConfig>& fragStore,
	 Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TReadSeq TReadSeq;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// All maps to mirror file ids to our ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New block?
		if (c == '{') {
			c = _streamGet(file);
			String<char> blockIdentifier;
			_parse_readIdentifier(file, blockIdentifier, c);
			_parse_skipLine(file, c);

			// Library block
			if (blockIdentifier == "LIB") {
				TLibraryStoreElement libEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous());
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "mea") {
						c = _streamGet(file);
						libEl.mean = _parse_readDouble(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "std") {
						c = _streamGet(file);
						libEl.std = _parse_readDouble(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl, Generous() );
				appendValue(fragStore.libraryNameStore, eid, Generous() );
			} else if (blockIdentifier == "FRG") {  // Fragment block
				TMatePairElement matePairEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				bool foundRds = false;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "lib") {
						c = _streamGet(file);
						matePairEl.libId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "rds") {
						foundRds = true;
						c = _streamGet(file);
						matePairEl.readId[0] = _parse_readNumber(file, c);
						c = _streamGet(file);
						matePairEl.readId[1] = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Only insert valid mate pairs
				if (foundRds) {
					frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl, Generous() );
					appendValue(fragStore.matePairNameStore, eid, Generous() );
				}
			} else if (blockIdentifier == "RED") {   // Read block
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> qual;
				TId matePairId = 0;
				TReadSeq seq;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "frg") {
						c = _streamGet(file);
						matePairId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "seq") {
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							_parse_readSequenceData(file,c,seq);
							_parse_skipWhitespace(file, c);
						}
					} else if (fieldIdentifier == "qlt") {
						clear(qual);
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) appendValue(qual, c, Generous() );
							c = _streamGet(file);
						}
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Set quality
				typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
				typedef typename Iterator<String<char> >::Type TQualIter;
				TReadIter begIt = begin(seq, Standard() );
				TQualIter qualIt = begin(qual);
				TQualIter qualItEnd = end(qual);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));

				// Insert the read
				readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
				appendRead(fragStore, seq, matePairId);
				appendValue(fragStore.readNameStore, eid, Generous() );
			} else if (blockIdentifier == "CTG") {   // Contig block
				TContigElement contigEl;
				TSize fromAligned = length(fragStore.alignedReadStore);
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> contigSeq;
				String<char> contigQual;
				while (c != '}') {
					// Are we entering a TLE block
					if (c == '{') {
						TAlignedElement alignEl;
						String<char> fdIdentifier;
						typedef typename TFragmentStore::TContigPos TContigPos;
						TContigPos offsetPos = 0;
						TContigPos clr1 = 0;
						TContigPos clr2 = 0;
						String<TContigPos> gaps;
						while (c != '}') {
							clear(fdIdentifier);
							_parse_readIdentifier(file, fdIdentifier, c);
							if (fdIdentifier == "src") {
								c = _streamGet(file);
								alignEl.readId = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "off") {
								c = _streamGet(file);
								if (c != '-') offsetPos = _parse_readNumber(file, c);
								else offsetPos = 0;
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "clr") {
								c = _streamGet(file);
								clr1 = _parse_readNumber(file, c);
								c = _streamGet(file);
								clr2 = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "gap") {
								c = _streamGet(file);
								_parse_skipWhitespace(file, c);
								while (c != '.') {
									if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
										TSize nextGap = _parse_readNumber(file, c);
										appendValue(gaps, nextGap, Generous() );
									}
									c = _streamGet(file);
								}
							} else {
								_parse_skipLine(file, c);
							}
						}
						_parse_skipLine(file, c);

						// Get the length of the read
						TId readId = (readIdMap.find(alignEl.readId))->second;
						TSize lenRead = length(value(fragStore.readSeqStore, readId));

						// Create the gap anchors
						typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
						int offset = 0;
						if ((clr1 < clr2) && (clr1>0)) offset = clr1;
						else if ((clr1 > clr2) && (clr1 < lenRead)) offset = lenRead - clr1;
						int diff = -1 * (int) (offset);
						// Clipped begin
						if (offset != 0) appendValue(alignEl.gaps, TContigGapAnchor(offset, 0), Generous() );
						// Internal gaps
						typedef typename Iterator<String<TContigPos>, Standard>::Type TPosIter;
						TPosIter posIt = begin(gaps, Standard() ); 
						TPosIter posItEnd = end(gaps, Standard() );
						TContigPos lastGap = 0;
						TSize gapLen = 0;
						TSize totalGapLen = 0;
						for(;posIt!=posItEnd; goNext(posIt)) {
							if (gapLen == 0) {
								++gapLen; ++totalGapLen;
								++diff;
								lastGap = value(posIt);
							} 
							else if (lastGap == value(posIt)) {
								++gapLen; ++totalGapLen;
								++diff;
							}
							else {
								appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
								gapLen = 1; ++totalGapLen;
								lastGap = value(posIt);
								++diff;
							}
						}
						if (gapLen > 0) appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
						// Clipped end
						if ((clr1 < clr2) && (clr2 < lenRead)) {
							diff -= (lenRead - clr2);				
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						} else if ((clr1 > clr2) && (clr2 > 0)) {
							diff -= clr2;
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						}
						
						// Set begin and end position
						if (clr1 < clr2) {
							alignEl.beginPos = offsetPos;
							alignEl.endPos = offsetPos + totalGapLen + (clr2 - clr1);
						} else {
							alignEl.beginPos = offsetPos + totalGapLen + (clr1 - clr2);
							alignEl.endPos = offsetPos;
						}
		
						// Append new align fragment, note: contigId must still be set
						alignEl.id = length(fragStore.alignedReadStore);
						appendValue(fragStore.alignedReadStore, alignEl, Generous() );
					} else {
						clear(fieldIdentifier);
						_parse_readIdentifier(file, fieldIdentifier, c);
						if (fieldIdentifier == "iid") {
							c = _streamGet(file);
							id = _parse_readNumber(file, c);
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "eid") {
							c = _streamGet(file);
							while ((c != '\n') && (c != '\r')) {
								appendValue(eid, c, Generous() );
								c = _streamGet(file);
							}
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "seq") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								do {
									_parse_readSequenceData(file,c,contigSeq);
								} while (c == '-');
								_parse_skipWhitespace(file, c);
							}
						} else if (fieldIdentifier == "qlt") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
									appendValue(contigQual, c, Generous() );
								}
								c = _streamGet(file);
							}
						} else {
							_parse_skipLine(file, c);
						}
					}
				}

				// Create the gap anchors
				char gapChar = gapValue<char>();
				typedef typename Iterator<String<char> >::Type TStringIter;
				TStringIter seqIt = begin(contigSeq);
				TStringIter seqItEnd = end(contigSeq);
				TStringIter qualIt = begin(contigQual);
				typedef typename TFragmentStore::TReadPos TPos;
				typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
				TPos ungappedPos = 0;
				TPos gappedPos = 0;
				bool gapOpen = false;
				for(;seqIt != seqItEnd; goNext(seqIt), goNext(qualIt), ++gappedPos) {
					if (value(seqIt) == gapChar) gapOpen = true;				
					else {
						if (gapOpen) {
							appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );
							gapOpen = false;
						}
						Dna5Q letter = value(seqIt);
						assignQualityValue(letter, value(qualIt));
						appendValue(contigEl.seq, letter, Generous() );
						++ungappedPos;
					}
				}
				if (gapOpen) appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );

				// Set the contigId in all aligned reads
				TSize toAligned = length(fragStore.alignedReadStore);
				TId newContigId = length(fragStore.contigStore);
				for(; fromAligned < toAligned; ++fromAligned) {
					(value(fragStore.alignedReadStore, fromAligned)).contigId = newContigId;
				}

				// Insert the contig
				appendValue(fragStore.contigStore, contigEl, Generous() );
				appendValue(fragStore.contigNameStore, eid, Generous() );
			} else {
				_parse_skipLine(file, c);
			}	
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end()) mateIt->libId = libIdPos->second;
			else mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end()) mateIt->readId[0] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end()) mateIt->readId[1] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for(;readIt != readItEnd; goNext(readIt)) {
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end()) readIt->matePairId = mateIdPos->second;
			else readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}
	TId myPairMatchId = 0;  // Dummy variable to count the matches
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		alignIt->pairMatchId = myPairMatchId++;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & target,
	  FragmentStore<TSpec, TConfig>& fragStore,
	  Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Write Header
	_streamWrite(target,"{UNV\niid:1\neid:seqan\ncom:\nafg file created with SeqAn\n.\n}\n");
	
	// Write Libraries
	typedef typename Iterator<typename TFragmentStore::TLibraryStore, Standard>::Type TLibIter;
	TLibIter libIt = begin(fragStore.libraryStore, Standard() );
	TLibIter libItEnd = end(fragStore.libraryStore, Standard() );
	bool noNamesPresent = (length(fragStore.libraryNameStore) == 0);
	for(TSize idCount = 0;libIt != libItEnd; goNext(libIt), ++idCount) {
		_streamWrite(target,"{LIB\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.libraryNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"{DST\n");
		_streamWrite(target,"mea:");
		_streamPutFloat(target, libIt->mean);
		_streamPut(target, '\n');
		_streamWrite(target,"std:");
		_streamPutFloat(target, libIt->std);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");	
		_streamWrite(target,"}\n");
	}

	// Write Fragments / mate pairs
	typedef typename Iterator<typename TFragmentStore::TMatePairStore, Standard>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore, Standard() );
	TMateIter mateItEnd = end(fragStore.matePairStore, Standard() );
	noNamesPresent = (length(fragStore.matePairNameStore) == 0);
	for(TSize idCount = 0;mateIt != mateItEnd; goNext(mateIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.matePairNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"lib:");
		_streamPutInt(target, mateIt->libId + 1);
		_streamPut(target, '\n');
		if ((mateIt->readId[0] != TMatePairElement::INVALID_ID) && (mateIt->readId[1] != TMatePairElement::INVALID_ID)) {
			_streamWrite(target,"rds:");
			_streamPutInt(target, mateIt->readId[0] + 1);
			_streamPut(target, ',');
			_streamPutInt(target, mateIt->readId[1] + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Get clear ranges
	typedef Pair<typename TFragmentStore::TReadPos, typename TFragmentStore::TReadPos> TClrRange;
	String<TClrRange> clrRange;
	fill(clrRange, length(fragStore.readStore), TClrRange(0,0));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore, Standard() );
	TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		typename TFragmentStore::TReadPos begClr = 0;
		typename TFragmentStore::TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clrRange, alignIt->readId) = TClrRange(begClr, endClr);
	}

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TReadStore, Standard>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore, Standard() );
	TReadIter readItEnd = end(fragStore.readStore, Standard() );
	noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; ++readIt, ++idCount) {
		_streamWrite(target,"{RED\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(value(fragStore.readSeqStore, idCount));
		TSeqIter seqItEnd = end(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			_streamWrite(target,"frg:");
			_streamPutInt(target, readIt->matePairId + 1);
			_streamPut(target, '\n');
		}
		if ((value(clrRange, idCount)).i1 != (value(clrRange, idCount)).i2) {
			_streamWrite(target,"clr:");
			_streamPutInt(target, (value(clrRange, idCount)).i1);
			_streamPut(target, ',');
			_streamPutInt(target, (value(clrRange, idCount)).i2);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Sort aligned reads according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());

	// Write Contigs
	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	alignIt = begin(fragStore.alignedReadStore);
	alignItEnd = end(fragStore.alignedReadStore);
	noNamesPresent = (length(fragStore.contigNameStore) == 0);
	for(TSize idCount = 0;contigIt != contigItEnd; goNext(contigIt), ++idCount) {
		_streamWrite(target,"{CTG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.contigNameStore, idCount));
			_streamPut(target, '\n');
		}
		String<char> qlt;
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TContigSeq>::Type TContigIter;
		TContigIter seqContigIt = begin(contigIt->seq);
		TContigIter seqContigItEnd = end(contigIt->seq);
		typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor> >::Type TGapsIter;
		TGapsIter itGaps = begin(contigIt->gaps);
		TGapsIter itGapsEnd = end(contigIt->gaps);
		int diff = 0;
		char gapChar = gapValue<char>();
		typename TFragmentStore::TContigPos mySeqPos = 0;
		TSize k = 0;
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			while (mySeqPos < itGaps->seqPos) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, value(seqContigIt));
				Ascii c = ' ';
				convertQuality(c, getQualityValue(value(seqContigIt)));
				appendValue(qlt, c, Generous() );
				goNext(seqContigIt);++mySeqPos;
			}
			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, gapChar);
				appendValue(qlt, '0', Generous() );
			}
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;seqContigIt != seqContigItEnd; goNext(seqContigIt)) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			++k;
			_streamPut(target, value(seqContigIt));
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqContigIt)));
			appendValue(qlt, c, Generous() );
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		for(TSize k = 0;k<length(qlt); k+=60) {
			TSize endK = k + 60;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		
		while ((alignIt != alignItEnd) && (idCount < alignIt->contigId)) goNext(alignIt);
		for(;(alignIt != alignItEnd) && (idCount == alignIt->contigId); goNext(alignIt)) {
			_streamWrite(target,"{TLE\n");
			_streamWrite(target,"src:");
			_streamPutInt(target, alignIt->readId + 1);
			_streamPut(target, '\n');
			typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor> >::Type TReadGapsIter;
			TReadGapsIter itGaps = begin(alignIt->gaps);
			TReadGapsIter itGapsEnd = end(alignIt->gaps);

			// Create the gaps string and the clear ranges
			typename TFragmentStore::TReadPos lenRead = length(value(fragStore.readSeqStore, alignIt->readId));
			TSize clr1 = 0;
			TSize clr2 = lenRead;
			// Create first clear range
			if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) clr1 = itGaps->seqPos;
			int diff = clr1;
			String<unsigned int> gaps;
			for(;itGaps != itGapsEnd; goNext(itGaps)) {
				for(int i = 0; i< diff - ((int) itGaps->seqPos - (int) itGaps->gapPos); ++i) {
					appendValue(gaps, itGaps->seqPos - clr1, Generous() );
				}
				// Clipped sequence
				if (diff - ((int) itGaps->seqPos - (int) itGaps->gapPos) < 0) {
					clr2 = lenRead + diff - ((int) itGaps->seqPos - (int) itGaps->gapPos);
				}
				diff = ((int) itGaps->seqPos - (int) itGaps->gapPos);
			}
			if (alignIt->beginPos > alignIt->endPos) {
				clr1 = lenRead - clr1;
				clr2 = lenRead - clr2;
			}
			_streamWrite(target,"off:");
			if (alignIt->beginPos < alignIt->endPos) _streamPutInt(target, alignIt->beginPos);
			else _streamPutInt(target, alignIt->endPos);
			_streamPut(target, '\n');
			_streamWrite(target,"clr:");
			_streamPutInt(target, clr1);
			_streamPut(target, ',');
			_streamPutInt(target, clr2);
			_streamPut(target, '\n');
			if (length(gaps)) {
				_streamWrite(target,"gap:\n");
				for(TSize z = 0;z<length(gaps); ++z) {
					_streamPutInt(target, value(gaps, z));
					_streamPut(target, '\n');
				}
				_streamWrite(target, ".\n");
			}
			_streamWrite(target,"}\n");
		}
		_streamWrite(target,"}\n");
	}
}

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences from multiple files
template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, StringSet<CharString> const &fileNameList)
{
	unsigned seqOfs = length(store.contigStore);
	for (unsigned filecount = 0; filecount < length(fileNameList); ++filecount)
	{
		MultiFasta multiFasta;
		if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY))
			return false;

		AutoSeqFormat format;
		guessFormat(multiFasta.concat, format);					// guess file format
		split(multiFasta, format);								// divide into single sequences

		unsigned seqCount = length(multiFasta);
		resize(store.contigStore, seqOfs + seqCount, Generous());
		resize(store.contigNameStore, seqOfs + seqCount, Generous());
		for(unsigned i = 0; i < seqCount; ++i)
		{
			assignSeq(store.contigStore[seqOfs + i].seq, multiFasta[i], format);			// read Genome sequence
			assignSeqId(store.contigNameStore[seqOfs + i], multiFasta[i], format);
		}
		seqOfs += seqCount;
	}
	reserve(store.contigStore, seqOfs, Exact());
	return seqOfs > 0;
}

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences from a single file
template <typename TFSSpec, typename TFSConfig>
bool loadContigs(FragmentStore<TFSSpec, TFSConfig> &store, CharString const &fileName)
{
	StringSet<CharString> fileNames;
	appendValue(fileNames, fileName);
	return loadContigs(store, fileNames);
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
