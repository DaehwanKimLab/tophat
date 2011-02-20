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
  $Id: graph_consensus_base.h 2103 2008-05-23 07:57:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_BASE_H
#define SEQAN_HEADER_CONSENSUS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Segment Match Generation.value.Overlap_Library:
	Segment matches from overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;



//////////////////////////////////////////////////////////////////////////////
// Consensus tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Consensus Calling:
..summary:A tag that specifies how to call the consensus.
*/


/**
.Tag.Consensus Calling.value.Majority_Vote:
	A consensus based on the most common character.
*/

struct Majority_Vote_;
typedef Tag<Majority_Vote_> const Majority_Vote;

/**
.Tag.Consensus Calling.value.Bayesian:
	A consensus based on bayesian probability.
*/

struct Bayesian_;
typedef Tag<Bayesian_> const Bayesian;



//////////////////////////////////////////////////////////////////////////////
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////

struct ConsensusOptions {
public:
	// Method
	// 0: graph-based multiple sequence alignment
	// 1: realign
	int method;

	// ReAlign Method
	// 0: Needleman-Wunsch
	// 1: Gotoh
	int rmethod;

	// Bandwidth of overlap alignment
	int bandwidth;

	// Number of computed overlaps per read (at the beginning and end of a read)
	int overlaps;

	// Minimum match length of a computed overlap
	int matchlength;

	// Minimum quality (in percent identity) of a computed overlap
	int quality;

	// Window size, only relevant for insert sequencing
	// If window == 0, no insert sequencing is assumed
	int window;
	
	// Output
	// 0: seqan style
	// 1: afg output format
	int output;

	// Multi-read alignment
	bool noalign;

	// Offset all reads, so the first read starts at position 0
	bool moveToFront;

	// Include reference genome
	bool include;

	// Scoring object for overlap alignments
	Score<int> sc;

	// Various input and output files
	std::string readsfile;				// File of reads in FASTA format
	std::string afgfile;				// AMOS afg file input
	std::string outfile;				// Output file name
	
	// Initialization
	ConsensusOptions() 
	{
		sc = Score<int>(2,-6,-4,-9);
	}
};


//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TStrSpec, typename TPosPair, typename TStringSpec, typename TSpec, typename TConfig, typename TId>
inline void 
_loadContigReads(StringSet<TValue, Owner<TStrSpec> >& strSet,
				 String<TPosPair, TStringSpec>& startEndPos,
				 FragmentStore<TSpec, TConfig> const& fragStore,
				 TId const contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Sort aligned reads according to contig id
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	resize(strSet, length(fragStore.alignedReadStore));

	// Retrieve all reads, limit them to the clear range and if required reverse complement them
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TSize numRead = 0;
	TReadPos begClr = 0;
	TReadPos endClr = 0;
	TSize lenRead = 0;
	TSize offset = 0;
	for(;alignIt != alignItEnd; ++alignIt) {
		offset = _min(alignIt->beginPos, alignIt->endPos);
		getClrRange(fragStore, *alignIt, begClr, endClr);
		strSet[numRead] = infix(fragStore.readSeqStore[alignIt->readId], begClr, endClr);
		lenRead = endClr - begClr;
		if (alignIt->beginPos < alignIt->endPos) appendValue(startEndPos, TPosPair(offset, offset + lenRead), Generous());
		else {
			reverseComplementInPlace(strSet[numRead]);
			appendValue(startEndPos, TPosPair(offset + lenRead, offset), Generous());
		}
		++numRead;
	}
	resize(strSet, numRead, Exact());
}




//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize, typename TReadSlot> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigId,
				 TSize& coverage,
				 TReadSlot& slot)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TMatrix>::Type TValue;

	// Gap char
	TValue gapChar = gapValue<TValue>();

	// Sort according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	
	// Find range of the given contig
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Sort the reads according to the begin position
	sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
	TAlignIter alignItBegin = alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Get the maximum coverage and the slot for each read
	typedef String<TSize> TFirstFreePos;
	typedef typename Iterator<TFirstFreePos, Standard>::Type TPosIter;
	TFirstFreePos freePos;
	TSize pos = 0;
	TSize maxTmp = 0;
	TSize numCol = 0;
	reserve(slot, alignItEnd - alignIt);
	for(;alignIt != alignItEnd; ++alignIt) {
		TPosIter itPos = begin(freePos, Standard());
		TPosIter itPosEnd = end(freePos, Standard());
		pos = 0;
		for(;itPos != itPosEnd; ++itPos, ++pos) 
			if (*itPos < _min(alignIt->beginPos, alignIt->endPos)) break;
		if (pos + 1 > length(freePos)) resize(freePos, pos+1, Generous());
		maxTmp = _max(alignIt->beginPos, alignIt->endPos);
		freePos[pos] = maxTmp;
		if (maxTmp > numCol) numCol = maxTmp;
		appendValue(slot, pos);
	}
	coverage = length(freePos);
	clear(freePos);

	// Fill the matrix
	typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
	fill(mat, coverage * numCol, '.');
	alignIt = alignItBegin;
	TSize readPos = 0;
	TMatIter matIt = begin(mat, Standard());
	typename TFragmentStore::TReadSeq myRead;
	for(;alignIt != alignItEnd; ++alignIt, ++readPos) {
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard());
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard());
		
		// Place each read inside the matrix
		myRead = fragStore.readSeqStore[alignIt->readId];
		TSize lenRead = length(myRead);
		TSize offset = alignIt->beginPos;
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplementInPlace(myRead);
			offset = alignIt->endPos;
		}
		matIt = begin(mat, Standard());
		matIt += (slot[readPos] * numCol + offset);

		typedef typename Iterator<typename TFragmentStore::TReadSeq, Standard>::Type TReadIter;
		TReadIter seqReadIt = begin(myRead, Standard());

		// First clear range
		TSize mySeqPos = 0;
		int diff = 0;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			mySeqPos = itGaps->seqPos;
			diff = -1 * mySeqPos;
			seqReadIt += mySeqPos;
		}
		TSize clr2 = lenRead;
		TSize stop = 0;
		for(;itGaps != itGapsEnd; ++itGaps) {
			// Any clipped sequence at the end
			stop =  itGaps->seqPos;
			if (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos) > 0) 
				clr2 = stop = lenRead - (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos));
			
			for(;mySeqPos < stop; ++matIt, ++seqReadIt, ++mySeqPos) 
				*matIt = *seqReadIt;

			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i, ++matIt) 
				*matIt = gapChar;
	
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;mySeqPos < clr2; ++mySeqPos, ++seqReadIt, ++matIt) 
			*matIt = *seqReadIt;
	}
	//for(TSize row = 0; row < coverage; ++row) {
	//	for(TSize col = 0; col<numCol; ++col) {
	//		std::cout << mat[row * numCol + col];
	//	}
	//	std::cout << std::endl;
	//}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigId,
				 TSize& coverage)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	String<TSize> slot;
	return convertAlignment(fragStore, mat, contigId, coverage, slot);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	TSize coverage;
	return convertAlignment(fragStore, mat, 0, coverage);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TGappedConsensus, typename TSize> 
inline void
getGappedConsensus(FragmentStore<TSpec, TConfig>& fragStore,
				   TGappedConsensus& gappedConsensus,
				   TSize contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TGappedConsensus>::Type TValue;
	
	TValue gapChar = gapValue<TValue>();
	typedef typename Iterator<typename TFragmentStore::TContigSeq, Standard>::Type TContigIter;
	TContigIter seqContigIt = begin(fragStore.contigStore[contigId].seq, Standard());
	TContigIter seqContigItEnd = end(fragStore.contigStore[contigId].seq, Standard());
	typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor>, Standard>::Type TGapsIter;
	TGapsIter itGaps = begin(fragStore.contigStore[contigId].gaps, Standard());
	TGapsIter itGapsEnd = end(fragStore.contigStore[contigId].gaps, Standard());
	int diff = 0;
	TSize mySeqPos = 0;
	for(;itGaps != itGapsEnd; goNext(itGaps)) {
		for(;mySeqPos < itGaps->seqPos; ++seqContigIt, ++mySeqPos) 
			appendValue(gappedConsensus, *seqContigIt, Generous());
			
		for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) 
			appendValue(gappedConsensus, gapChar, Generous());
			diff = (itGaps->gapPos - itGaps->seqPos);
	}
	for(;seqContigIt != seqContigItEnd; ++seqContigIt) 
		appendValue(gappedConsensus, *seqContigIt, Generous());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TGappedConsensus, typename TSize> 
inline void
assignGappedConsensus(FragmentStore<TSpec, TConfig>& fragStore,
					  TGappedConsensus& gappedCons,
					  TSize contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TGappedConsensus>::Type TValue;
	TValue gapChar = gapValue<TValue>();

	// Update the contig
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigStoreElement;
	TContigStoreElement& contigEl = fragStore.contigStore[contigId];
	clear(contigEl.gaps);
	clear(contigEl.seq);

	// Create the sequence and the gap anchors
	typedef typename Iterator<TGappedConsensus, Standard>::Type TStringIter;
	TStringIter seqIt = begin(gappedCons, Standard());
	TStringIter seqItEnd = end(gappedCons, Standard());
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	TReadPos ungappedPos = 0;
	TReadPos gappedPos = 0;
	bool gapOpen = false;
	for(;seqIt != seqItEnd; ++seqIt, ++gappedPos) {
		if (*seqIt == gapChar) gapOpen = true;				
		else {
			if (gapOpen) {
				appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
				gapOpen = false;
			}
			Dna5Q letter = *seqIt;
			assignQualityValue(letter, 'D');
			appendValue(contigEl.seq, letter);
			++ungappedPos;
		}
	}
	if (gapOpen) 
		appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
}



//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize, typename TConfigOptions>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos,
				   TConfigOptions const& consOpt) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename Id<TOutGraph>::Type TId;

	// Initialization
	TStringSet& seqSet = stringSet(gOut);

	// Select all overlapping reads and record the diagonals of the band
	String<Pair<TId, TId> > pList;
	String<Pair<int, int> > diagList;
	if (consOpt.window == 0) selectPairs(seqSet, begEndPos, consOpt.bandwidth, pList, diagList);
	else selectPairsIndel(seqSet, begEndPos, consOpt.window, pList, diagList);

	// Set-up a sparse distance matrix
	Graph<Undirected<double> > pairGraph;
	
	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<int> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, diagList, begEndPos, consOpt.sc, consOpt.matchlength, consOpt.quality, consOpt.overlaps, matches, scores, pairGraph, Overlap_Library() );
	clear(pList);
	clear(diagList);

	// Use these segment matches for the initial alignment graph
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, consOpt.sc, ReScore() );
	clear(matches);
	clear(scores);

	// Guide Tree
	Graph<Tree<double> > guideTree;
	upgmaTree(pairGraph, guideTree);
	clear(pairGraph);

	// Triplet library extension
	graphBasedTripletLibraryExtension(g);

	// Perform a progressive alignment
	progressiveAlignment(g, guideTree, gOut);
	clear(g);
	clear(guideTree);
}

//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos) 
{
	ConsensusOptions consOpt;
	consensusAlignment(gOut, begEndPos, consOpt);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFragSpec, typename TConfig, typename TStringSet, typename TCargo, typename TSpec, typename TContigId>
inline void
updateContig(FragmentStore<TFragSpec, TConfig>& fragStore,
			 Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TContigId contigId)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<TSize, TSize> TComponentLength;
	typedef char TValue;

	// Initialization
	TStringSet& strSet = stringSet(g);
	TSize nseq = length(strSet);
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize maxCoverage = 0;
	TSize len = 0;
	String<TValue> mat;

	// Store for each read the begin position, the end position and the row in the alignment matrix
	String<TSize> readBegEndRowPos;
	resize(readBegEndRowPos, 3*nseq);

	// Strongly Connected Components, topological sort, and length of each component
	String<TSize> component;
	String<TSize> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		TSize numOfComponents = length(order);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef String<std::pair<TSize, TSize> > TComponentToRank;
		TComponentToRank compToRank;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) 
			appendValue(compToRank, std::make_pair(order[compIndex], compIndex), Generous());
		::std::sort(begin(compToRank, Standard()), end(compToRank, Standard()));

		typedef Pair<TSize, TSize> TRankPair;
		typedef String<TRankPair> TSequenceToRanks;
		TSequenceToRanks seqToRank;
		resize(seqToRank, nseq);
		typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
		TVertexIterator itVertex(g);
		for(;!atEnd(itVertex);++itVertex) {
			TVertexDescriptor vert = value(itVertex);
			TSize seq = idToPosition(strSet, sequenceId(g, vert));
			if (fragmentBegin(g, vert) == 0) 
				seqToRank[seq].i1 = ::std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), ::std::make_pair((TSize) component[vert], (TSize) 0))->second;
			if (fragmentBegin(g, vert) + fragmentLength(g, vert) == length(strSet[seq]))
				seqToRank[seq].i2 = ::std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), ::std::make_pair((TSize) component[vert], (TSize) 0))->second;
		}
		clear(compToRank);

		// Assign the sequences to rows
		String<TSize> seqToRow;
		resize(seqToRow, nseq);
		maxCoverage = 0;
		typedef String<bool> TLeftOver;
		typedef typename Iterator<TLeftOver, Standard>::Type TLeftOverIter;
		TLeftOver leftOver;
		fill(leftOver, nseq, true);
		typedef String<std::pair<TSize, TSize> > TSeqToBegin;
		typedef typename Iterator<TSeqToBegin, Standard>::Type TSeqToBeginIter;
		TSeqToBegin seqToBegin;
		TSize finishedSeq = 0;
		while(finishedSeq < nseq) {
			TLeftOverIter itL = begin(leftOver, Standard());
			TLeftOverIter itLEnd = end(leftOver, Standard());
			for(TSize pos = 0; itL != itLEnd; ++itL, ++pos) 
				if (*itL) appendValue(seqToBegin, std::make_pair((seqToRank[pos]).i1, pos), Generous());
			::std::sort(begin(seqToBegin, Standard()), end(seqToBegin, Standard()));
			
			TSize endPos = 0;
			TSeqToBeginIter itSB = begin(seqToBegin, Standard());
			TSeqToBeginIter itSBEnd = end(seqToBegin, Standard());
			for(;itSB != itSBEnd;++itSB) {
				if (endPos <= (*itSB).first) {
					TSize currentSeq = (*itSB).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver[currentSeq] = false;
					++finishedSeq;
				}	
			}
			clear(seqToBegin);
			++maxCoverage;
		}
		clear(leftOver);

		// Create the matrix
		len = 0;
		String<TSize> compOffset;
		resize(compOffset, numOfComponents);
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		fill(mat, len * maxCoverage, gapChar);

		// Fill in the segments
		typedef typename Infix<TString>::Type TInfix;
		typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
		typedef typename TGraph::TPosToVertexMap TPosToVertexMap;
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			TSize c = property(component, it->second);
			TSize row = seqToRow[idToPosition(strSet, it->first.first)];
			//if (row == 0) {
			//	std::cout << sequenceId(g, it->second) << ':' << str << ',' << strSet[sequenceId(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceId(g, it->second)]).i1 << ',' << (seqToRank[sequenceId(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str, Standard());
			TInfixIter sItEnd = end(str, Standard());
			TSize i = compOffset[c];
			for(TSize pCol = i;sIt!=sItEnd;++sIt, ++pCol, ++i) 
				mat[row * len + pCol] = *sIt;
		}
		String<bool> active;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			TSize offset = compOffset[order[compIndex]];
			TSize currentCompLength = compLength[order[compIndex]];

			clear(active);
			fill(active, maxCoverage, false);

			// Find the empty rows
			for(TSize i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) 
					active[(seqToRow[i])] = true;
			}
			
			// Substitute false gaps with special gap character
			for(TSize i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(TSize pCol = offset;pCol < offset + currentCompLength;++pCol) 
						mat[i * len + pCol] = specialGap;
				}
			}
		}

		// Get the new begin and end positions
		for(TSize i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), length(strSet[i]) - 1);
			TSize readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), 0))];
			TSize readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[3*i] = readBegin;
			readBegEndRowPos[3*i+1] = readEnd;
			readBegEndRowPos[3*i+2] = seqToRow[i];
		}

	}
	clear(component);
	clear(order);
	compLength.clear();

	
	//// Debug code
	//for(TSize row = 0; row<maxCoverage; ++row) {
	//	for(TSize col = 0; col<len; ++col) {
	//		std::cout << mat[row * len + col];			
	//	}
	//	std::cout << std::endl;
	//}

	// Create the new consensus
	typedef typename Value<TString>::Type TAlphabet;
	String<TValue> gappedCons;
	consensusCalling(mat, gappedCons, maxCoverage, TAlphabet(), Majority_Vote());

	// Assign new consensus
	assignGappedConsensus(fragStore, gappedCons, contigId);

	// Update all aligned reads
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TReadPos ungappedPos = 0;
	TReadPos gappedPos = 0;
	bool gapOpen;
	for(TSize i = 0;alignIt != alignItEnd; ++alignIt, ++i) {
		TSize lenRead = length(fragStore.readSeqStore[alignIt->readId]);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, *alignIt, begClr, endClr);
		clear(alignIt->gaps);
		ungappedPos = begClr;
		if (alignIt->beginPos > alignIt->endPos) ungappedPos = lenRead - endClr;
		if (ungappedPos != 0) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, 0));
		gappedPos = 0;
		gapOpen = false;
		for(TSize column = readBegEndRowPos[3*i]; column<readBegEndRowPos[3*i + 1]; ++column, ++gappedPos) {
			if (mat[readBegEndRowPos[3*i + 2] * len + column] == gapChar) gapOpen = true;				
			else {
				if (gapOpen) {
					appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
					gapOpen = false;
				}
				++ungappedPos;
			}
		}
		if (gapOpen) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
		if (alignIt->beginPos < alignIt->endPos) {
			if (endClr != lenRead) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - (lenRead - endClr)), Generous());
		} else {
			if (begClr != 0) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - begClr), Generous());
		}

		// Set new begin and end position
		if (alignIt->beginPos < alignIt->endPos) {
			alignIt->beginPos = readBegEndRowPos[3*i];
			alignIt->endPos = readBegEndRowPos[3*i+1];
		} else {
			alignIt->beginPos = readBegEndRowPos[3*i+1];
			alignIt->endPos = readBegEndRowPos[3*i];
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TCounters, typename TSize, typename TAlphabet>
inline void
__countLetters(String<TValue, TSpec> const& mat,
			   TCounters& counterValues,
			   TSize alignDepth,
			   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TMatrix;
	typedef typename Iterator<TMatrix, Standard>::Type TMatIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;

	// Set-up counter values
	typedef typename Value<TCounters>::Type TCounter;
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	resize(counterValues, len);
	for(TSize i=0;i<len; ++i) {
		TCounter counter;
		fill(counter, alphabetSize + 1, 0);
		counterValues[i] = counter;
	}

	// Count all 
	TMatIter matIt = begin(mat, Standard());
	TMatIter matItEnd = end(mat, Standard());
	TCounterIt countIt = begin(counterValues, Standard());
	TSize pos = 0;
	for(; matIt != matItEnd; ++matIt, ++countIt, ++pos) {
		if (pos % len == 0) countIt = begin(counterValues, Standard());
		if (*matIt != specialGap) {
			if (*matIt == gapChar) ++value(*countIt, alphabetSize);
			else ++value(*countIt, ordValue((TAlphabet) *matIt));
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TGappedCons, typename TAlignDepth, typename TAlphabet>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 TGappedCons& gappedConsensus,
				 TAlignDepth maxCoverage,
				 TAlphabet,
				 Bayesian)
{
	typedef double TProbability;
	typedef String<TProbability> TProbabilityDistribution;
	typedef String<TProbabilityDistribution> TPositionalPrDist;
	typedef typename Iterator<TPositionalPrDist, Standard>::Type TPosPrDistIter;

	typedef typename Size<String<TValue, TSpec> >::Type TSize;
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';

	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, maxCoverage, TAlphabet() );


	// Initialization
	TSize len = length(mat) / maxCoverage;
	TProbabilityDistribution backroundDist;
	fill(backroundDist, alphabetSize + 1, ((TProbability) 1 / (TProbability) (alphabetSize + 1)));
	
	// Get an initial consensus
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues, Standard());
	TCounterIt countItEnd = end(counterValues, Standard());
	TPositionalPrDist posPrDist;
	TValue c = TAlphabet();
	for(;countIt != countItEnd; ++countIt) {
		TSize max = 0;
		typedef typename Iterator<TCounter, Standard>::Type TCIt;
		TCIt cIt = begin(*countIt, Standard());
		TCIt cItEnd = end(*countIt, Standard());
		TSize pos = 0;
		for(;cIt != cItEnd; ++cIt, ++pos) {
			if (*cIt > max) {
				max = *cIt;
				c = (pos == alphabetSize) ? gapChar : (TValue) TAlphabet(pos);
			}
		}
		TProbabilityDistribution prDist;
		fill(prDist, alphabetSize + 1, 0);
		if (c == gapChar) prDist[alphabetSize] = 1;
		else prDist[ordValue((TAlphabet) c)] = 1;
		appendValue(posPrDist, prDist, Generous());
	}

	TSize run = 1;
	TProbabilityDistribution pI;
	TProbabilityDistribution pIJ;
	TProbabilityDistribution pIOld;
	TProbabilityDistribution pIJOld;
	while (run) {
		// Store the values from the last iteration
		pIOld = pI;
		pIJOld = pIJ;

		// Count all letters in the consensus
		TProbabilityDistribution nI;
		fill(nI, alphabetSize + 1, 0);
		TPosPrDistIter itPosPrDist = begin(posPrDist, Standard());
		TPosPrDistIter itPosPrDistEnd = end(posPrDist, Standard());
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) 
			for(TSize i = 0; i<(alphabetSize + 1); ++i) 
				nI[i] += (*itPosPrDist)[i];
	
		// Composition probabilities
		clear(pI);
		resize(pI, alphabetSize + 1);
		TProbability lenPosPrDist = (TProbability) length(posPrDist);
		for(TSize i = 0; i<length(pI); ++i) 
			pI[i] = nI[i] / lenPosPrDist;
		

		// Count all letters that agree / disagree with the consensus
		TProbabilityDistribution nIJ;
		fill(nIJ, (alphabetSize + 1) * (alphabetSize + 1), 0);
		typedef String<TValue, TSpec> TMatrix;
		typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
		TMatIter matIt = begin(mat, Standard());
		TMatIter matItEnd = end(mat, Standard());
		itPosPrDist = begin(posPrDist, Standard());
		TSize pos = 0;
		for(; matIt != matItEnd; ++matIt, ++itPosPrDist, ++pos) {
			if (pos % len == 0) itPosPrDist = begin(posPrDist, Standard());
			TValue c = *matIt;
			if (c != specialGap) {
				TSize fragJ = (c != gapChar) ? ordValue(TAlphabet(c)) : alphabetSize;
				for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) 
					nIJ[consI * (alphabetSize + 1) + fragJ] += (*itPosPrDist)[consI];	
			}
		}

		// Sequencing error probabilities
		clear(pIJ);
		resize(pIJ, (alphabetSize + 1) * (alphabetSize + 1));
		TProbability sumIJ = 0;
		for(TSize diag = 0; diag<(alphabetSize + 1); ++diag) sumIJ += nIJ[diag * (alphabetSize + 1) + diag];
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) 
			for(TSize fragJ = 0; fragJ<(alphabetSize + 1); ++fragJ)
				pIJ[consI * (alphabetSize + 1) + fragJ] = nIJ[consI * (alphabetSize + 1) + fragJ] / sumIJ;
	
		// Recompute positional probability distribution
		itPosPrDist = begin(posPrDist, Standard());
		TSize col = 0;
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist, ++col) {
			TProbabilityDistribution prDist;
			resize(prDist, alphabetSize + 1);
			for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
				TProbability numerator = pI[consI];
				TProbability denominator = 0;
				for(TSize allI = 0; allI<(alphabetSize + 1); ++allI) {
					TProbability denominatorSub = value(pI, allI);
					for(TSize row = 0; row < maxCoverage; ++row) {
						TValue c = mat[row * len + col];
						if (c != specialGap) {
							TSize fragJ = (c != gapChar) ? ordValue(TAlphabet(c)) : alphabetSize;
							if (allI == consI) 
								numerator *= pIJ[allI * (alphabetSize + 1) + fragJ]; 
							denominatorSub *= pIJ[allI * (alphabetSize + 1) + fragJ]; 
						}
					}
					denominator += denominatorSub;
				}
				prDist[consI] = numerator / denominator;
			}
			*itPosPrDist = prDist;
		}	

		// Check termination criterion
		TProbability eps = 0.00001;
		typedef typename Iterator<TProbabilityDistribution, Standard>::Type TProbIter;
		TProbIter pIter = begin(pIOld, Standard());
		TProbIter pIterCompare = begin(pI, Standard());
		TProbIter pIterEnd = end(pIOld, Standard());
		TSize runOld = run;
		for(;pIter != pIterEnd; ++pIter, ++pIterCompare) {
			if (*pIter > *pIterCompare) {
				if (*pIter - *pIterCompare > eps) {
					++run;
					break;
				}
			} else {
				if (*pIterCompare - *pIter > eps) {
					++run;
					break;
				}
			}
		}
		if (runOld == run) {
			pIter = begin(pIJOld, Standard());
			pIterCompare = begin(pIJ, Standard());
			pIterEnd = end(pIJOld, Standard());
			for(;pIter != pIterEnd; ++pIter, ++pIterCompare) {
				if (*pIter > *pIterCompare) {
					if (*pIter - *pIterCompare > eps) {
						++run;
						break;
					}
				} else {
					if (*pIterCompare - *pIter > eps) {
						++run;
						break;
					}
				}
			}
		}

		if (runOld == run) {
			std::cout << "Iterations: " << run << std::endl;
			run = 0;
		}
	}
	
	// Compute the most likely consensus
	TPosPrDistIter itPosPrDist = begin(posPrDist, Standard());
	TPosPrDistIter itPosPrDistEnd = end(posPrDist, Standard());
	clear(gappedConsensus);
	for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) {
		TProbability max = 0;
		TSize ind = 0;
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			if ((*itPosPrDist)[consI] > max) {
				max = (*itPosPrDist)[consI];
				ind = consI;
			}
		}
		if (ind == alphabetSize) appendValue(gappedConsensus, gapChar);
		else appendValue(gappedConsensus, TAlphabet(ind));
	}

}

//////////////////////////////////////////////////////////////////////////////


template <typename TFragSpec, typename TConfig, typename TContigId>
inline void
consensusCalling(FragmentStore<TFragSpec, TConfig>& fragStore,
				 TContigId contigId,
				 Bayesian)
{
	SEQAN_CHECKPOINT

	typedef FragmentStore<TFragSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Value<TReadSeq>::Type TAlphabet;
	typedef char TValue;

	// Convert the contig to an alignment matrix
	typedef String<TValue> TAlignMat;
	TAlignMat mat;
	TSize maxCoverage;
	convertAlignment(fragStore, mat, contigId, maxCoverage);

	// Call the consensus
	String<TValue> gappedConsensus;
	consensusCalling(mat, gappedConsensus, maxCoverage, TAlphabet(), Bayesian());

	// Assign the new consensus
	assignGappedConsensus(fragStore, gappedConsensus, contigId);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TGappedCons, typename TAlignDepth, typename TAlphabet>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 TGappedCons& gappedConsensus,
				 TAlignDepth maxCoverage,
				 TAlphabet,
				 Majority_Vote)
{
	typedef typename Size<String<TValue, TSpec> >::Type TSize;
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	
	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, maxCoverage, TAlphabet() );
	
	// Get the consensus
	typedef typename Iterator<TCounters, Standard>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues, Standard());
	TCounterIt countItEnd = end(counterValues, Standard());
	clear(gappedConsensus);
	TSize max = 0;
	TValue c = TValue();
	TSize pos = 0;
	for(;countIt != countItEnd; ++countIt) {
		max = 0;	
		typedef typename Iterator<TCounter, Standard>::Type TCIt;
		TCIt cIt = begin(*countIt, Standard());
		TCIt cItEnd = end(*countIt, Standard());
		pos = 0;
		for(;cIt != cItEnd; ++cIt, ++pos) {
			if (*cIt > max) {
				max = *cIt;
				c = (pos == alphabetSize) ? gapChar : (TValue) TAlphabet(pos);
			}
		}
		appendValue(gappedConsensus, c);
	}
}


//////////////////////////////////////////////////////////////////////////////


template <typename TFragSpec, typename TConfig, typename TContigId>
inline void
consensusCalling(FragmentStore<TFragSpec, TConfig>& fragStore,
				 TContigId contigId,
				 Majority_Vote)
{
	typedef FragmentStore<TFragSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Value<TReadSeq>::Type TAlphabet;
	typedef char TValue;

	// Convert the contig to an alignment matrix
	typedef String<TValue> TAlignMat;
	TAlignMat mat;
	TSize maxCoverage;
	convertAlignment(fragStore, mat, contigId, maxCoverage);

	// Call the consensus
	String<TValue> gappedConsensus;
	consensusCalling(mat, gappedConsensus, maxCoverage, TAlphabet(), Majority_Vote());

	// Assign the new consensus
	assignGappedConsensus(fragStore, gappedConsensus, contigId);
}


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Old proprietary FastaReadFormat
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.FastaReadFormat:
	Fasta read format to write a multi-read alignment.
*/

struct FastaReadFormat_;
typedef Tag<FastaReadFormat_> const FastaReadFormat;

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & file,
	  FragmentStore<TSpec, TConfig>& fragStore,
	  FastaReadFormat) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef char TMultiReadChar;
	TMultiReadChar gapChar = gapValue<TMultiReadChar>();

	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	for(TSize idCount = 0;contigIt != contigItEnd; ++contigIt, ++idCount) {
		// Alignment matrix
		typedef String<TMultiReadChar> TAlignMat;
		TAlignMat mat;
		TSize maxCoverage;
		String<TSize> readSlot;
		convertAlignment(fragStore, mat, idCount, maxCoverage, readSlot);
		TSize len = length(mat) / maxCoverage;
		
		// Gapped consensus sequence
		typedef String<TMultiReadChar> TGappedConsensus;
		TGappedConsensus gappedConsensus;
		getGappedConsensus(fragStore, gappedConsensus, idCount);

		// Print the alignment matrix
		String<TSize> coverage;
		fill(coverage, len, 0);
		typedef typename Iterator<TGappedConsensus, Standard>::Type TConsIter;
		TConsIter itCons = begin(gappedConsensus, Standard());
		TSize winSize = 60;
		int offset = 2;
		TSize column = 0;
		while (column<len) {
			TSize window_end = column + winSize;
			if (window_end >= len) window_end = len;
			// Position
			for(int i = 0; i<offset - 2; ++i) _streamPut(file,' ');
			_streamWrite(file,"Pos: ");
			_streamPutInt(file, column);
			_streamPut(file,'\n');
			// Ruler
			for(int i = 0; i<offset + 3; ++i) _streamPut(file,' ');
			for(TSize local_col = 1; local_col<window_end - column + 1; ++local_col) {
				if ((local_col % 10)==0) _streamPut(file, ':');
				else if ((local_col % 5)==0) _streamPut(file, '.');
				else _streamPut(file, ' ');
			}
			_streamPut(file,'\n');
			// Matrix
			for(TSize row = 0; row<maxCoverage; ++row) {
				TSize tmp = row;
				int off = 0;
				while (tmp / 10 != 0) {
					tmp /= 10;
					++off;
				}
				for(int i = 0; i<offset - off; ++i) _streamPut(file,' ');
				_streamPutInt(file, row);
				_streamPut(file,':');
				_streamPut(file,' ');
				for(TSize local_col = column; local_col<window_end; ++local_col) {
					_streamPut(file, mat[row * len + local_col]);
					if (mat[row * len + local_col] != '.') ++coverage[local_col];
				}
				_streamPut(file,'\n');
			}
			_streamPut(file,'\n');
	
			// Consensus
			for(int i = 0; i<offset; ++i) _streamPut(file,' ');
			_streamWrite(file,"C: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col, ++itCons) 
				_streamPut(file, *itCons);
			_streamPut(file,'\n');
			for(int i = 0; i<offset-1; ++i) _streamPut(file,' ');
			_streamWrite(file,">2: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col) {
				if (coverage[local_col] > 2) _streamPut(file, gappedConsensus[local_col]);
				else _streamPut(file, gapChar);
			}
			_streamPut(file,'\n');
			_streamPut(file,'\n');
			column+=winSize;
		}
		_streamPut(file,'\n');
		_streamPut(file,'\n');

		// Print all aligned reads belonging to this contig

		// Sort according to contigId
		sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	
		// Find range of the given contig
		typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
		TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());

		// Sort the reads according to the begin position
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
		TAlignIter alignItTmp = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		TAlignIter alignItTmpEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		String<std::pair<TSize, TSize> > idToPos;
		reserve(idToPos, alignItTmpEnd - alignItTmp);
		for(TSize iCount = 0; alignItTmp!=alignItTmpEnd; ++iCount, ++alignItTmp) 
			appendValue(idToPos, std::make_pair(alignItTmp->id, readSlot[iCount]));
		::std::sort(begin(idToPos, Standard()), end(idToPos, Standard()));

		// Sort the reads according to the id
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortId());
		alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());

		bool noNamesPresent = (length(fragStore.readNameStore) == 0);
		for(TSize iCount = 0;alignIt != alignItEnd; ++alignIt, ++iCount) {

			// Print all reads
			_streamWrite(file,"typ:");
			if (!noNamesPresent) {
				_streamPut(file,'R');
				_streamPutInt(file, iCount);
			} else _streamWrite(file, fragStore.readNameStore[alignIt->readId]);
			_streamPut(file,'\n');
			_streamWrite(file,"seq:");
			_streamWrite(file, fragStore.readSeqStore[alignIt->readId]);
			_streamPut(file,'\n');
			_streamWrite(file,"Pos:");
			_streamPutInt(file, alignIt->beginPos);
			_streamPut(file,',');
			_streamPutInt(file, alignIt->endPos);
			_streamPut(file,'\n');
#ifndef CELERA_OFFSET
			TSize begClr = 0;
			TSize endClr = 0;
			getClrRange(fragStore, *alignIt, begClr, endClr);
			_streamWrite(file,"clr:");
			_streamPutInt(file, begClr);
			_streamPut(file,',');
			_streamPutInt(file, endClr);
			_streamPut(file,'\n');
#endif
			std::stringstream gapCoords;
			TSize letterCount = 0;
			TSize gapCount = 0;
			for(TSize column = _min(alignIt->beginPos, alignIt->endPos); column < _max(alignIt->beginPos, alignIt->endPos); ++column) {
				if (mat[idToPos[iCount].second * len + column] == gapChar) {
					++gapCount;
					gapCoords << letterCount << ' ';
				} else ++letterCount;
			}
			_streamWrite(file,"dln:");
			_streamPutInt(file, gapCount);
			_streamPut(file,'\n');
			_streamWrite(file,"del:");
			_streamWrite(file, gapCoords.str().c_str());
			_streamPut(file,'\n');
			_streamPut(file,'\n');

		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Read simulator format: Simple fasta read file with positions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig, typename TFilePath>
inline bool 
_convertSimpleReadFile(TFile& file,
					   FragmentStore<TSpec, TConfig>& fragStore,
					   TFilePath& filePath, 
					   bool moveToFront)
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TContigPos TPos;
	typedef typename TFragmentStore::TReadSeq TReadSeq;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	

	// All maps to mirror file ids to our internal ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;


	// Parse the file and convert the internal ids
	TPos maxPos = 0;
	TPos minPos = SupremumValue<TPos>::VALUE;
	TId count = 0;
	TValue c;
	if ((!file) || (_streamEOF(file))) return false;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New read?
		if (c == '>') {
			TAlignedElement alignEl;
			TId id = count;
			TId fragId = count;
			TId repeatId = 0;
			
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);

			// Get the layout positions
			alignEl.beginPos = _parse_readNumber(file, c);
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			alignEl.endPos = _parse_readNumber(file, c);
			
			// Any attributes?
			String<char> eid;
			String<char> qlt;
			TReadSeq seq;
			if (c == '[') {
				String<char> fdIdentifier;
				while (c != ']') {
					c = _streamGet(file);
					_parse_skipWhitespace(file, c);
					clear(fdIdentifier);
					_parse_readIdentifier(file, fdIdentifier, c);
					if (fdIdentifier == "id") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
					} else if (fdIdentifier == "fragId") {
						c = _streamGet(file);
						fragId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "repeatId") {
						c = _streamGet(file);
						repeatId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(eid, c, Generous());
							c = _streamGet(file);
						}
					} else if (fdIdentifier == "qlt") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(qlt, c, Generous());
							c = _streamGet(file);
						}
					} else {
						// Jump to next attribute
						while ((c != ',') && (c != ']')) {
							c = _streamGet(file);
						}
					}
				}
			}
			_parse_skipLine(file, c);
			_parse_skipWhitespace(file, c);
			while ((!_streamEOF(file)) && (c != '>')) {
				_parse_readSequenceData(file,c, seq);
				_parse_skipWhitespace(file, c);
			}
			
			// Set quality
			typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
			typedef typename Iterator<String<char> >::Type TQualIter;
			TReadIter begIt = begin(seq, Standard() );
			TReadIter begItEnd = begin(seq, Standard() );
			if (length(qlt)) {
				TQualIter qualIt = begin(qlt);
				TQualIter qualItEnd = end(qlt);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));
			} else {
				for(;begIt != begItEnd; goNext(begIt)) assignQualityValue(value(begIt), 'D');
			}

			// Set eid if not given
			if (empty(eid)) {
				std::stringstream input;
				input << "R" << id;
				input << "-" << repeatId;
				eid = input.str().c_str();
			}

			// Insert the read
			readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
			appendRead(fragStore, seq, fragId);
			appendValue(fragStore.readNameStore, eid, Generous());

			// Insert an aligned read
			TSize readLen = length(seq);
			if (alignEl.beginPos < alignEl.endPos) {
				if (readLen != alignEl.endPos - alignEl.beginPos) {
					alignEl.endPos = alignEl.beginPos + readLen;
				}
				if (alignEl.beginPos < minPos) minPos = alignEl.beginPos;
				if (alignEl.endPos > maxPos) maxPos = alignEl.endPos;
			} else {
				if (readLen != alignEl.beginPos - alignEl.endPos) {
					alignEl.beginPos = alignEl.endPos + readLen;
				}
				if (alignEl.endPos < minPos) minPos = alignEl.endPos;
				if (alignEl.beginPos > maxPos) maxPos = alignEl.beginPos;
			}
			alignEl.readId = id;
			alignEl.pairMatchId =  fragId;
			alignEl.contigId = 0;
			alignEl.id = length(fragStore.alignedReadStore);
			appendValue(fragStore.alignedReadStore, alignEl, Generous());
			++count;
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Read contig or reference sequence
	TContigElement contigEl;
	std::string fileName = filePath + 'S';
	FILE* strmRef = fopen(fileName.c_str(), "rb");
	String<char> contigEid = "C0";
	if ((strmRef) && (!_streamEOF(strmRef))) {
		c = _streamGet(strmRef);
		while (!_streamEOF(strmRef)) {
			if (_streamEOF(strmRef)) break;
			if (c == '>') {
				clear(contigEid);
				c = _streamGet(strmRef);
				while ((c != '\r') && (c != '\n')) {
					appendValue(contigEid, c, Generous());
					c = _streamGet(strmRef);
				}
				_parse_skipLine(strmRef, c);
				_parse_skipWhitespace(strmRef, c);
				while ((!_streamEOF(strmRef)) && (c != '>')) {
					_parse_readSequenceData(strmRef,c,contigEl.seq);
					_parse_skipWhitespace(strmRef, c);
				}
			} else {
				_parse_skipLine(strmRef, c);
			}
		}
		fclose(strmRef);
	}
	if (empty(contigEl.seq)) {
		typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
		if (moveToFront) appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos - minPos), Generous());
		else appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos), Generous());
	}
	appendValue(fragStore.contigStore, contigEl, Generous());
	appendValue(fragStore.contigNameStore, contigEid, Generous());


	// Read fragments
	fileName = filePath + 'F';
	FILE* strmFrag = fopen(fileName.c_str(), "rb");
	if ((strmFrag) && (!_streamEOF(strmFrag))) {
		c = _streamGet(strmFrag);
		while (!_streamEOF(strmFrag)) {
			if (_streamEOF(strmFrag)) break;
			if (c == '>') {
				TMatePairElement matePairEl;
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmFrag, c);
			
				// Any attributes?
				std::stringstream input;
				input << "F" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmFrag);
						_parse_skipWhitespace(strmFrag, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmFrag, fdIdentifier, c);
						if (fdIdentifier == "libId") {
							c = _streamGet(strmFrag);
							matePairEl.libId = _parse_readNumber(strmFrag, c);
						} else if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmFrag);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c, Generous());
								c = _streamGet(strmFrag);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmFrag);
							}
						}
					}
				}
				_parse_skipLine(strmFrag, c);
				_parse_skipWhitespace(strmFrag, c);

				// Read the two reads belonging to this mate pair
				matePairEl.readId[0] = _parse_readNumber(strmFrag, c);
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);
				matePairEl.readId[1] = _parse_readNumber(strmFrag, c);
				_parse_skipLine(strmFrag, c);

				// Insert mate pair
				if (matePairEl.readId[0] != matePairEl.readId[1]) {
					frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl, Generous());
					appendValue(fragStore.matePairNameStore, eid, Generous());
				}
			} else {
				_parse_skipLine(strmFrag, c);
			}
		}
		fclose(strmFrag);
	}
	

	// Read libraries
	fileName = filePath + 'L';
	FILE* strmLib = fopen(fileName.c_str(), "rb");
	if ((strmLib) && (!_streamEOF(strmLib))) {
		c = _streamGet(strmLib);
		while (!_streamEOF(strmLib)) {
			if (_streamEOF(strmLib)) break;
			if (c == '>') {

				TLibraryStoreElement libEl;
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmLib, c);
			
				// Any attributes?
				std::stringstream input;
				input << "L" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmLib);
						_parse_skipWhitespace(strmLib, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmLib, fdIdentifier, c);
						if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmLib);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c, Generous());
								c = _streamGet(strmLib);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmLib);
							}
						}
					}
				}
				_parse_skipLine(strmLib, c);
				_parse_skipWhitespace(strmLib, c);

				// Read the mean and standard deviation
				libEl.mean = _parse_readNumber(strmLib, c);
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);
				libEl.std = _parse_readNumber(strmLib, c);
				_parse_skipLine(strmLib, c);

				// Insert mate pair
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl, Generous());
				appendValue(fragStore.libraryNameStore, eid, Generous());
			} else {
				_parse_skipLine(strmLib, c);
			}
		}
		fclose(strmLib);
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
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		if (moveToFront) {
			alignIt->beginPos -= minPos;
			alignIt->endPos -= minPos;
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Rudimentary write functions for CeleraFrg and Celera Cgb
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
_writeCeleraFrg(TFile& target,
				FragmentStore<TSpec, TConfig>& fragStore) 
{

	SEQAN_CHECKPOINT
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// Iterate over all aligned reads to get the clear ranges
	typedef Pair<TReadPos, TReadPos> TClrRange;
	String<TClrRange> clearStr;
	resize(clearStr, length(fragStore.readStore));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clearStr, alignIt->readId) = TClrRange(begClr, endClr);
	}

	// Write Reads
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	bool noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; goNext(readIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"act:");
		_streamPut(target, 'A');
		_streamPut(target, '\n');
		_streamWrite(target,"acc:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"src:\n");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamWrite(target, "\n.\n");
		}
		_streamWrite(target,"etm:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(value(fragStore.readSeqStore, idCount));
		TSeqIter seqItEnd = end(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 70 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 70 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		// Note: Clear range does not have to be ordered, e.g. no indication for reverse complemented reads, this is happening in cgb records
		_streamWrite(target,"clr:");
		_streamPutInt(target, (value(clearStr, idCount)).i1);
		_streamPut(target, ',');
		_streamPutInt(target, (value(clearStr, idCount)).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig>
inline void 
_writeCeleraCgb(TFile& target,
				FragmentStore<TSpec, TConfig>& fragStore) 
{
	SEQAN_CHECKPOINT
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename TFragmentStore::TReadPos TReadPos;


	// Write the first contig
	TId contigId = 0;

	// Sort the reads according to position
	sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());

	// Write Header
	_streamWrite(target,"{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\n");
	_streamWrite(target,"len:");
	_streamPutInt(target, length((value(fragStore.contigStore, contigId)).seq));
	_streamPut(target, '\n');
	_streamWrite(target,"cns:\n.\nqlt:\n.\nfor:0\n");
	_streamWrite(target,"nfr:");
	_streamPutInt(target, length(fragStore.readStore));
	_streamPut(target, '\n');

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	TSize offsetLeft = _min(alignIt->beginPos, alignIt->endPos);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (contigId != alignIt->contigId) continue;
		_streamWrite(target,"{IMP\n");
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"mid:");
		_streamPutInt(target, alignIt->readId + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"con:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"pos:");
		_streamPutInt(target, alignIt->beginPos - offsetLeft);
		_streamPut(target, ',');
		_streamPutInt(target, alignIt->endPos - offsetLeft);
		_streamPut(target, '\n');
		_streamWrite(target,"dln:0\n");
		_streamWrite(target,"del:\n");
		_streamWrite(target,"}\n");
	}
	_streamWrite(target,"}\n");
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
