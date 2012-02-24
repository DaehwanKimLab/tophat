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

#ifndef SEQAN_HEADER_STORE_ALL_H
#define SEQAN_HEADER_STORE_ALL_H

//#include <stdio.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Name Store Cache
//////////////////////////////////////////////////////////////////////////////

	template <typename TNameStore, typename TName>
	struct NameStoreLess
	{
		typedef typename Position<TNameStore>::Type TId;

		TNameStore *store;
		TName *name;
		
		NameStoreLess() {}

		NameStoreLess(TNameStore &_store, TName &_name):
			store(&_store),
			name(&_name) {}
		
		template <typename TId>
		inline bool operator() (TId a, TId b) const
		{
			if (a != maxValue(a))
			{
				if (b != maxValue(b))
					return (*store)[a] < (*store)[b];
				else
					return (*store)[a] < *name;
			} else
			{
				if (b != maxValue(b))
					return *name < (*store)[b];
				else
					return false;
			}
		}
	};

/**
.Class.NameStoreCache
..summary:Stores a mapping from names to ids.
..cat:Fragment Store
..signature:FragmentStore<>
..signature:NameStoreCache<TNameStore[, TName]>
..param.TNameStore:The name store to be cached.
...see:Class.FragmentStore
..param.TName:The name type.
...default:$Value<TNameStore>::Type$
...type:Shortcut.CharString

.Memvar.NameStoreCache#NameStoreCache
..summary:Constructor
..signature:NameStoreCache<TNameStore, TName> (nameStore)
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..class:Class.NameStoreCache
..include:seqan/store.h
*/
	
	template <typename TNameStore, typename TName = typename Value<TNameStore>::Type>
	struct NameStoreCache
	{
		typedef typename Position<TNameStore>::Type TId;
		typedef NameStoreLess<TNameStore, TName> TLess;
		typedef std::set<TId, TLess> TSet;
		
		TSet nameSet;
		TNameStore &store;
		TName name;

		NameStoreCache(TNameStore &_store):
			nameSet(TLess(_store, name)),
			store(_store)
		{
			for (unsigned i = 0; i < length(store); ++i)
				nameSet.insert(i);
		}
	};

//////////////////////////////////////////////////////////////////////////////
// refresh

/**
.Function.refresh:
..summary:Recreate a name store cache.
..cat:Fragment Store
..signature:refresh(cache)
..param.cache:A @Class.NameStoreCache@ object.
...type:Class.NameStoreCache
..see:Function.getIdByName
..include:seqan/store.h
*/
    
    template <typename TNameStore, typename TName>
    inline void
    refresh(NameStoreCache<TNameStore, TName> &cache)
    {
		cache.nameSet.clear();
		for (unsigned i = 0; i < length(cache.store); ++i)
			cache.nameSet.insert(i);
	}
		

//////////////////////////////////////////////////////////////////////////////
// getIdByName
    
/**
.Function.getIdByName:
..summary:Appends a name to a name store.
..cat:Fragment Store
..signature:getIdByName(nameStore, name, id[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be searched.
...type:Shortcut.CharString
..param.id:The resulting id.
..param.cache:A structure to efficiently retrieve the id for a given name. If ommited a brute force method is used to search.
...default:Tag.Nothing
...type:Class.NameStoreCache
..returns:$true$ if the name was found and $false$ if not.
..see:Function.getIdByName
..include:seqan/store.h
*/

    template <typename TNameStore, typename TName, typename TPos>
    inline bool 
    getIdByName(TNameStore &store, TName &name, TPos &pos)
    {
        typedef typename Iterator<TNameStore, Standard>::Type TNameStoreIter;
        
        // Iterator over read names
        for (TNameStoreIter iter = begin(store); iter != end(store); ++iter)
		{
            // if the element was found
            if (name == getValue(iter))
			{
                // set the ID
                pos = position(iter);
                // and break the loop
                return true;
            }
        }
        return false;
    }
	
    template <typename TNameStore, typename TName, typename TPos, typename TContext>
    inline bool 
    getIdByName(TNameStore &store, TName &name, TPos &pos, TContext &)
	{
		return getIdByName(store, name, pos);
	}
    
    template<typename TNameStore, typename TName, typename TPos, typename TCNameStore, typename TCName>
    inline bool 
    getIdByName(TNameStore &, TName &name, TPos &pos, NameStoreCache<TCNameStore, TCName> &context)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
		typedef typename Position<TNameStore>::Type TId;
		typedef NameStoreCache<TCNameStore, TCName> TNameStoreCache;
		typedef typename TNameStoreCache::TSet TSet;
		
		TSet const &set = context.nameSet;	
		context.name = name;
		
		typename TSet::const_iterator it = set.find(maxValue<TId>());
		if (it != set.end())
		{
			pos = *it;
			return true;
		}
		return false;
    }
	
//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendName:
..summary:Appends a name to a name store.
..cat:Fragment Store
..signature:appendName(nameStore, name[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be appended.
...type:Shortcut.CharString
..param.cache:A structure to efficiently retrieve the id for a given name. See @Function.getIdByName@.
...default:Tag.Nothing
...type:Class.NameStoreCache
..see:Function.getIdByName
..include:seqan/store.h
*/

    template <typename TNameStore, typename TName>
    inline void
    appendName(TNameStore &store, TName &name)
    {
		appendValue(store, name, Generous());
	}
	
    template <typename TNameStore, typename TName, typename TContext>
    inline void
    appendName(TNameStore &store, TName &name, TContext &)
    {
		appendName(store, name, Generous());
	}
	
    template <typename TNameStore, typename TName, typename TCNameStore, typename TCName>
    inline void
    appendName(TNameStore &store, TName &name, NameStoreCache<TCNameStore, TCName> &context)
    {
		appendValue(store, name, Generous());
		context.nameSet.insert(length(store) - 1);
    }
	

//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct FragmentStoreConfig 
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TContigFileSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TAlignedReadTagStoreSpec;
	typedef void					TAnnotationStoreElementSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Fragment Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FragmentStore
..summary:Multi-Container to store contigs, reads and a multiple alignment between them.
..cat:Fragment Store
..signature:FragmentStore<>
..signature:FragmentStore<TSpec[, TConfig]>
..param.TSpec:The specializing type.
...default:$void$
..param.TConfig:The configuration struct.
...default:$FragmentStoreConfig<TSpec>$

.Typedef.FragmentStore#TReadStore
..summary:Type of the $readStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TReadSeqStore
..summary:Type of the $readSeqStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TMatePairStore
..summary:Type of the $matePairStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TLibraryStore
..summary:Type of the $libraryStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigFileStore
..summary:Type of the $contigFileStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigStore
..summary:Type of the $contigStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignedReadStore
..summary:Type of the $alignedReadStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAnnotationStore
..summary:Type of the $annotationStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignQualityStore
..summary:Type of the $alignQualityStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAlignedReadTagStore
..summary:Type of the $alignedReadTagStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TReadNameStore
..summary:Type of the $readNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TMatePairNameStore
..summary:Type of the $matePairNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TLibraryNameStore
..summary:Type of the $libraryNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TContigNameStore
..summary:Type of the $contigNameStore$ member.
..class:Class.FragmentStore
.Typedef.FragmentStore#TAnnotationNameStore
..summary:Type of the $annotationNameStore$ member.
..class:Class.FragmentStore

.Memvar.FragmentStore#readStore
..summary:String that maps from $readId$ to $<matePairId>$.
..remarks:Value type is @Class.ReadStoreElement@.
..type:Typedef.FragmentStore#TReadStore
..class:Class.FragmentStore
.Memvar.FragmentStore#readSeqStore
..summary:String that maps from $readId$ to $readSeq$.
..type:Typedef.FragmentStore#TReadSeqStore
..class:Class.FragmentStore
.Memvar.FragmentStore#matePairStore
..summary:String that maps from $matePairId$ to $<readId[2], libId>$.
..type:Typedef.FragmentStore#TMatePairStore
..remarks:Value type is @Class.MatePairStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#libraryStore
..summary:String that maps from $libId$ to $<mean, std>$.
..type:Typedef.FragmentStore#TLibraryStore
..remarks:Value type is @Class.LibraryStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#contigFileStore
..summary:String that maps from $contigFileId$ to $<fileName, firstContigId>$.
..type:Typedef.FragmentStore#TContigFileStore
..remarks:Value type is @Class.ContigFile@.
..class:Class.FragmentStore
.Memvar.FragmentStore#contigStore
..summary:String that maps from $contigId$ to $<contigSeq, contigGaps, contigFileId>$.
..type:Typedef.FragmentStore#TContigStore
..remarks:Value type is @Class.ContigStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#alignedReadStore
..summary:String that stores $<alignId, readId, contigId, pairMatchId, beginPos, endPos, gaps>$.
..type:Typedef.FragmentStore#TAlignedReadStore
..remarks:Value type is @Class.AlignedReadStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#annotationStore
..summary:String that maps from $annoId$ to $<contigId, typeId, beginPos, endPos, parentId, lastChildId, nextSiblingId, values>$.
..type:Typedef.FragmentStore#TAnnotationStore
..class:Class.FragmentStore
.Memvar.FragmentStore#alignQualityStore
..summary:String that maps from $alignId$ to $<pairScore, score, errors>$.
..type:Typedef.FragmentStore#TAlignQualityStore
..remarks:Value type is @Class.AlignQualityStoreElement@.
..class:Class.FragmentStore
.Memvar.FragmentStore#alignedReadTagStore
..summary:String that maps from $alignId$ to $alignTag$.
..type:Typedef.FragmentStore#TAlignedReadTagStore
..class:Class.FragmentStore
.Memvar.FragmentStore#readNameStore
..summary:String that maps from $readId$ to $readName$.
..type:Typedef.FragmentStore#TReadNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#matePairNameStore
..summary:String that maps from $contigId$ to $contigName$.
..type:Typedef.FragmentStore#TMatePairNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#libraryNameStore
..summary:String that maps from $libId$ to $libName$.
..type:Typedef.FragmentStore#TLibraryNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#contigNameStore
..summary:String that maps from $contigId$ to $contigName$.
..type:Typedef.FragmentStore#TContigNameStore
..class:Class.FragmentStore
.Memvar.FragmentStore#annotationNameStore
..summary:String that maps from $annoId$ to $annoName$.
..type:Typedef.FragmentStore#TAnnotationNameStore
..class:Class.FragmentStore
..include:seqan/store.h
*/


template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
class FragmentStore
{
private:
	typedef typename TConfig::TReadStoreElementSpec			TReadStoreElementSpec;
	typedef typename TConfig::TReadSeqStoreSpec				TReadSeqStoreSpec;
	typedef typename TConfig::TMatePairStoreElementSpec		TMatePairStoreElementSpec;
	typedef typename TConfig::TLibraryStoreElementSpec		TLibraryStoreElementSpec;
	typedef typename TConfig::TContigStoreElementSpec		TContigStoreElementSpec;
	typedef typename TConfig::TContigFileSpec				TContigFileSpec;
	typedef typename TConfig::TAlignedReadStoreElementSpec	TAlignedReadStoreElementSpec;
	typedef typename TConfig::TAlignedReadTagStoreSpec		TAlignedReadTagStoreSpec;
	typedef typename TConfig::TAnnotationStoreElementSpec	TAnnotationStoreElementSpec;

public:
	typedef typename TConfig::TMean					TMean;
	typedef typename TConfig::TStd					TStd;
	typedef typename TConfig::TMappingQuality		TMappingQuality;
	
	typedef typename TConfig::TReadSeq				TReadSeq;
	typedef typename TConfig::TContigSeq			TContigSeq;

	typedef typename Position<TReadSeq>::Type		TRSeqPos_;
	typedef typename Position<TContigSeq>::Type     TCSeqPos_;
	typedef typename MakeSigned_<TRSeqPos_>::Type	TReadPos;
	typedef typename MakeSigned_<TCSeqPos_>::Type	TContigPos;
	
	typedef GapAnchor<TReadPos>						TReadGapAnchor;
	typedef GapAnchor<TContigPos>					TContigGapAnchor;
	
	typedef StringSet<CharString>					TNameStore;

	typedef AnnotationStoreElement< TContigPos, TAnnotationStoreElementSpec >	TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId								TAnnotationStoreElementId;

	typedef String< ReadStoreElement< TReadStoreElementSpec > >												TReadStore;
	typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >										TMatePairStore;
	typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >							TLibraryStore;
	typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >			TContigStore;
	typedef String< ContigFile< TContigFileSpec > >															TContigFileStore;
	typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >	TAlignedReadStore;
	typedef String< AlignQualityStoreElement< TMappingQuality >	>											TAlignQualityStore;
	typedef StringSet<CharString, TAlignedReadTagStoreSpec>													TAlignedReadTagStore;
	typedef String< TAnnotationStoreElement >																TAnnotationStore;
	typedef String< IntervalTree< TContigPos, TAnnotationStoreElementId > >									TIntervalTreeStore;
	typedef StringSet<TReadSeq, TReadSeqStoreSpec>															TReadSeqStore;
	
	typedef TNameStore																						TReadNameStore;
	typedef TNameStore																						TMatePairNameStore;
	typedef TNameStore																						TLibraryNameStore;
	typedef TNameStore																						TContigNameStore;
	typedef TNameStore																						TAnnotationNameStore;
	typedef TNameStore																						TAnnotationTypeStore;
	typedef TNameStore																						TAnnotationKeyStore;
	
	// main containers
	TReadStore			readStore;			// readId       -> matePairId
	TMatePairStore		matePairStore;		// matePairId   -> readId0, readId1, libraryId
	TLibraryStore		libraryStore;		// libraryId    -> libSizeMean, libSizeStd
	TContigStore		contigStore;		// contigId     -> contigSeq, contigGaps, contigFileId
	TContigFileStore	contigFileStore;	// contigFileId -> fileName, firstContigId
	TAlignedReadStore	alignedReadStore;	//              -> id, readId, contigId, pairMatchId (not matePairId!), beginPos, endPos, gaps
	TAnnotationStore	annotationStore;	// annoId       -> parentId, contigId, beginPos, endPos
	TIntervalTreeStore	intervalTreeStore_F;		// treeId (same as contigId)	-> intervalTree (F: forward strand)
	TIntervalTreeStore	intervalTreeStore_R;		// 						(R: reverse complement strand)

											// REMARKS: 
											// 1)
											//    beginPos <= endPos     forward strand
											//    beginPos >  endPos     backward strand (reverse complement)
											// 2) 
											//    The alignedReadStore can arbitrarily be resorted. The unique identifier id should
											//    be used to address additional information for each alignedRead in additional tables.

	// we store the read sequences in a seperate stringset to reduce the memory overhead 
	TReadSeqStore		readSeqStore;

	// extra Sam fields
	TAlignQualityStore		alignQualityStore;
	TAlignedReadTagStore	alignedReadTagStore;

	// retrieve the names of reads, mate-pairs, libraries, contigs, annotations by their ids
	TReadNameStore			readNameStore;
	TMatePairNameStore		matePairNameStore;
	TLibraryNameStore		libraryNameStore;
	TContigNameStore		contigNameStore;
	TAnnotationNameStore	annotationNameStore;
	TAnnotationTypeStore	annotationTypeStore;
	TAnnotationKeyStore		annotationKeyStore;
	
	NameStoreCache<TReadNameStore, CharString>			readNameStoreCache;
	NameStoreCache<TContigNameStore, CharString>		contigNameStoreCache;
	NameStoreCache<TAnnotationNameStore, CharString>	annotationNameStoreCache;
	NameStoreCache<TAnnotationTypeStore, CharString>	annotationTypeStoreCache;
	NameStoreCache<TAnnotationKeyStore, CharString>		annotationKeyStoreCache;

	enum {
		ANNO_ROOT,
		ANNO_DELETED,
		ANNO_GENE,
		ANNO_MRNA,
		ANNO_CDS,
		ANNO_EXON,
		ANNO_FIVE_PRIME_UTR,
		ANNO_INTRON,
		ANNO_THREE_PRIME_UTR,
		ANNO_PREDEFINED
	};
	
	FragmentStore():
		readNameStoreCache(readNameStore),
		contigNameStoreCache(contigNameStore),
		annotationNameStoreCache(annotationNameStore),
		annotationTypeStoreCache(annotationTypeStore),
		annotationKeyStoreCache(annotationKeyStore)
	{
		// ATTENTION: The order of these keywords must correspond to the order of the enums above.
		appendName(annotationTypeStore, "<root>", annotationTypeStoreCache);
		appendName(annotationTypeStore, "<deleted>", annotationTypeStoreCache);
		appendName(annotationTypeStore, "gene", annotationTypeStoreCache);
		appendName(annotationTypeStore, "mRNA", annotationTypeStoreCache);
		appendName(annotationTypeStore, "CDS", annotationTypeStoreCache);
		appendName(annotationTypeStore, "exon", annotationTypeStoreCache);
		appendName(annotationTypeStore, "five_prime_UTR", annotationTypeStoreCache);
		appendName(annotationTypeStore, "intron", annotationTypeStoreCache);
		appendName(annotationTypeStore, "three_prime_UTR", annotationTypeStoreCache);
		_storeClearAnnotations(*this);
	}
};

//////////////////////////////////////////////////////////////////////////////
    
template < typename TSpec, typename TConfig >
struct VertexDescriptor< FragmentStore<TSpec, TConfig> > 
{
	typedef FragmentStore<TSpec, TConfig>		TFragmentStore_;
	typedef typename Id<TFragmentStore_>::Type	Type;
};

//////////////////////////////////////////////////////////////////////////////

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig>, TIterSpec >::Type
begin(FragmentStore<TSpec, TConfig> &store, TIterSpec const) 
{
	return Iter<FragmentStore<TSpec, TConfig>, TIterSpec>(store);
}

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig> const, TIterSpec >::Type
begin(FragmentStore<TSpec, TConfig> const &store, TIterSpec const) 
{
	return Iter<FragmentStore<TSpec, TConfig> const, TIterSpec>(store);
}

//////////////////////////////////////////////////////////////////////////////

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig>, TIterSpec >::Type
end(FragmentStore<TSpec, TConfig> &store, TIterSpec const) 
{
	return Iter<FragmentStore<TSpec, TConfig>, TIterSpec>(store, MinimalCtor());
}

template < typename TConfig, typename TSpec, typename TIterSpec >
inline typename Iterator<FragmentStore<TSpec, TConfig> const, TIterSpec >::Type
end(FragmentStore<TSpec, TConfig> const &store, TIterSpec const) 
{
	return Iter<FragmentStore<TSpec, TConfig> const, TIterSpec>(store, MinimalCtor());
}

//////////////////////////////////////////////////////////////////////////////

template < typename TConfig, typename TSpec >
inline void
_storeClearAnnotations(FragmentStore<TSpec, TConfig> & me)
{
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	
	resize(me.annotationStore, 1);
	resize(me.annotationNameStore, 1);
	resize(me.annotationTypeStore, (unsigned)TFragmentStore::ANNO_PREDEFINED);
	
	TAnnotation root;
	root.typeId = 0;
	me.annotationStore[0] = root;
	me.annotationNameStore[0] = "<root>";
}

template < typename TConfig, typename TSpec >
inline void
_storeRemoveTempAnnoNames(FragmentStore<TSpec, TConfig> & me)
{
	typedef FragmentStore<TSpec, TConfig>						TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore			TAnnotationStore;
	typedef typename TFragmentStore::TAnnotationNameStore		TAnnotationNameStore;
	typedef typename TFragmentStore::TAnnotationTypeStore		TAnnotationTypeStore;

	typedef typename GetValue<TAnnotationNameStore>::Type		TName;
	typedef typename GetValue<TAnnotationTypeStore>::Type		TType;
	typedef typename Iterator<TAnnotationStore, Standard>::Type	TAnnoIter;
	typedef typename Position<TAnnotationStore>::Type			TPosition;
	
	TAnnoIter it = end(me.annotationStore, Standard()) - 1;
	TAnnoIter itBegin = begin(me.annotationStore, Standard());
	TPosition pos = it - itBegin;
	for (; itBegin <= it; --it, --pos)
	{
		if ((*it).typeId >= length(me.annotationTypeStore))
			continue;
		TName name = me.annotationNameStore[pos];
		TType type = me.annotationTypeStore[(*it).typeId];
		if (length(name) > length(type) + 3)
			if (prefix(name, 2) == "__" && infix(name, 2, 2 + length(type)) == type && name[2 + length(type)] == '_')
				assignValue(me.annotationNameStore, pos, "");
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TId>
inline typename GetValue<typename FragmentStore<TSpec, TConfig>::TAnnotationNameStore>::Type
getAnnoName(FragmentStore<TSpec, TConfig> & store, TId id)
{
	return store.annotationNameStore[id];
}

template <typename TSpec, typename TConfig, typename TId>
inline typename GetValue<typename FragmentStore<TSpec, TConfig>::TAnnotationTypeStore>::Type
getAnnoType(FragmentStore<TSpec, TConfig> & store, TId id)
{
	return store.annotationTypeStore[id];
}

template <typename TSpec, typename TConfig, typename TId>
inline CharString
getAnnoUniqueName(FragmentStore<TSpec, TConfig> & store, TId id)
{
	if (id < length(store.annotationNameStore))
		return getAnnoName(store, id);
	
	std::stringstream tmp;
	tmp << "__" << getAnnoType(store, store.annotationStore[id].typeId) << '_' << id;
	return tmp.str();
}


//////////////////////////////////////////////////////////////////////////////
// append functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _storeAppendRead
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
template <typename TSpec, typename TConfig, typename TId, typename TName, typename TString, typename TFlag, typename TContext>
inline void 
_storeAppendRead (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TId & readId, 
	TName & qname,
	TString & readSeq,
	TFlag & flag,
	TContext &)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

	// search for readId by name
	if (getIdByName(fragStore.readNameStore, qname, readId, fragStore.readNameStoreCache))
	{
		if ((flag & 1) == 1)
		{
			// if the read is in the store and paired
			// check the mate pair store if it is the same mate of the pair
			// assuming that only one flag 0x040 or 0x0080 is 1
			int inPair = 1 - ((flag & 0x40) >> 6);	// bit 7 is set => inPair = 0
													// else inPair = 1 (even if bits 6 and 7 are not set)
			
			TId matePairId = fragStore.readStore[readId].matePairId;
			if (matePairId != TMatePairElement::INVALID_ID)
			{
				readId = fragStore.matePairStore[matePairId].readId[inPair];
				if (readId == TMatePairElement::INVALID_ID)
				{
					// create new entry in read and read name store
					// set sequence and mate pair ID in new read store element
					readId = appendRead(fragStore, readSeq, matePairId);
					// add the identifier to the read name store
					appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
					// set the ID in the mate pair store
					fragStore.matePairStore[matePairId].readId[inPair] = readId;
					return;
				}
			}
		} else 
			return;
	}

	// if the read name is not in the store
	// create new entry in read and read name store
	readId = length(fragStore.readStore);

	// if the read is paired
	if ((flag & 1) == 1)
	{
		TMatePairElement mateElem;
		// set the first or second read ID in the mate pair element
		TId matePairId = length(fragStore.matePairStore);
		mateElem.readId[(flag & 0x80) >> 7] = readId;
		// get a new mate pair ID and add the new mate pair element
		appendValue(fragStore.matePairStore, mateElem);
		// set the new mate pair ID in the read element
		appendRead(fragStore, readSeq, matePairId);
	} 
	// if read is not paired
	else
		appendRead(fragStore, readSeq);
	
	appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
}

//////////////////////////////////////////////////////////////////////////////
// _storeAppendContig
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void 
_storeAppendContig (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TId & contigId, 
	TName & rName)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	
	if (!getIdByName(fragStore.contigNameStore, rName, contigId, fragStore.contigNameStoreCache))
	{
		// if the contig is not in the store yet
		// set the ID on the last entry after appending
		contigId = length(fragStore.contigStore);
		// append contig store
		appendName(fragStore.contigNameStore, rName, fragStore.contigNameStoreCache);
		appendValue(fragStore.contigStore, TContigElement());
//		std::cout << "added contig:" << rName << std::endl;	
	}
}

//////////////////////////////////////////////////////////////////////////////
// _annotationAppendAnnotation
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void 
_storeAppendAnnotationName (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TId & annotationId, 
	TName & annotationName)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;

	if (!empty(annotationName) && getIdByName(fragStore.annotationNameStore, annotationName, annotationId, fragStore.annotationNameStoreCache))
		return;
	// if the annotation is not in the store yet
	// set the ID on the last entry after appending
	annotationId = length(fragStore.annotationNameStore);
	// append to annotationName store
	appendName(fragStore.annotationNameStore, annotationName, fragStore.annotationNameStoreCache);
//	std::cout << "added annotation:" << annotationName << std::endl;	
}

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void 
_storeAppendType (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TId & typeId, 
	TName & annotationType)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	
	if (!getIdByName(fragStore.annotationTypeStore, annotationType, typeId, fragStore.annotationTypeStoreCache))
	{
		// if the annotation type name is not in the store yet
		// set the ID on the last entry after appending
		typeId = length(fragStore.annotationTypeStore);
		// append to annotationType store
		if (!empty(annotationType))
			appendName(fragStore.annotationTypeStore, annotationType, fragStore.annotationTypeStoreCache);
//		std::cout << "added type:" << annotationType << std::endl;	
	}
}

template <typename TSpec, typename TConfig, typename TId, typename TName>
inline void 
_storeAppendKey (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TId & keyId,
	TName & annotationKey)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	
	if (!getIdByName(fragStore.annotationKeyStore, annotationKey, keyId, fragStore.annotationKeyStoreCache))
	{
		// if the key name is not in the store yet
		// set the ID on the last entry after appending
		keyId = length(fragStore.annotationKeyStore);
		// append to annotationKey store
		if (!empty(annotationKey))
			appendName(fragStore.annotationKeyStore, annotationKey, fragStore.annotationKeyStoreCache);
//		std::cout << "added key:" << annotationKey << std::endl;	
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey, typename TValue>
inline void 
annotationAssignValueByKey (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TAnnotation & annotation,
	TKey const & key,
	TValue const & value)
{
	typedef typename TAnnotation::TValues	TValues;
	typedef typename Size<TValues>::Type	TKeyId;
	
	TKeyId keyId = 0;	
	_storeAppendKey(fragStore, keyId, key);
	if (length(annotation.values) <= keyId)
		resize(annotation.values, keyId + 1);
	assignValue(annotation.values, keyId, value);
}

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey, typename TValue>
inline bool 
annotationGetValueByKey (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TAnnotation const & annotation,
	TKey const & key,
	TValue & value)
{
	typedef typename TAnnotation::TValues	TValues;
	typedef typename Size<TValues>::Type	TKeyId;
	
	TKeyId keyId = 0;	
	if (getIdByName(fragStore.annotationKeyStore, key, keyId, fragStore.annotationKeyStoreCache))
	{
		assign(value, annotation.values[keyId]);
		return true;
	}
	return false;
}

template <typename TSpec, typename TConfig, typename TAnnotation, typename TKey>
inline CharString 
annotationGetValueByKey (
	FragmentStore<TSpec, TConfig> & fragStore, 
	TAnnotation const & annotation,
	TKey const & key)
{
	typedef typename TAnnotation::TValues	TValues;
	typedef typename Size<TValues>::Type	TKeyId;
	
	TKeyId keyId = 0;
	if (getIdByName(fragStore.annotationKeyStore, key, keyId, fragStore.annotationKeyStoreCache))
		return annotation.values[keyId];
	else
		return "";
}


//////////////////////////////////////////////////////////////////////////////
// Read Store Accessors
//////////////////////////////////////////////////////////////////////////////

/**
.Function.clearReads
..summary:Removes all reads from a fragment store.
..cat:Fragment Store
..signature:clearReads(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:This function clears the @Memvar.FragmentStore#readStore@, @Memvar.FragmentStore#readSeqStore@ and @Memvar.FragmentStore#readNameStore@.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig>
inline void
clearReads(FragmentStore<TSpec, TConfig> &me)
{
	clear(me.readStore);
	clear(me.readSeqStore);
	clear(me.readNameStore);
}

/**
.Function.appendRead:
..summary:Appends a read to a fragment store.
..cat:Fragment Store
..signature:appendRead(store, read[, matePairId])
..signature:appendRead(store, read, name[, matePairId])
..param.store:The fragment store.
...type:Class.FragmentStore
..param.read:The read sequence.
..param.name:The read name.
...type:Shortcut.CharString
..param.matePairId:Id of mate-pair this read is part of.
...default:$INVALID_ID$, which corresponds to an unmated read.
..returns:The $readId$ of the newly appended read.
..remarks:This function appends a single read to the @Memvar.FragmentStore#readStore@ and @Memvar.FragmentStore#readSeqStore@.
If name is given, it is appended to the @Memvar.FragmentStore#readNameStore@.
..see:Function.getRead
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read, 
	TId matePairId)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
	typename Value<TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readSeqStore, read, Generous());
	return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read, 
	CharString const &name,
	TId matePairId)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
	typename Value<TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readSeqStore, read, Generous());
	appendValue(me.readNameStore, name, Generous());
	return length(me.readStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read)
{
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;
    
	return appendRead(me, read, TReadStoreElement::INVALID_ID);
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read,
	CharString const &name)
{
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;
    
	return appendRead(me, read, name, TReadStoreElement::INVALID_ID);
}

/**
.Function.getRead
..summary:Returns the read with the given $readId$.
..cat:Fragment Store
..signature:getRead(store, readId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId:The read id.
..returns:The sequence of the read with id $readId$ from the @Memvar.FragmentStore#readSeqStore@.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
inline typename Value<typename FragmentStore<TSpec, TConfig>::TReadSeqStore>::Type
getRead(
	FragmentStore<TSpec, TConfig> &me, 
	TId id)
{
	return value(me.readSeqStore, id);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendAlignedRead:
..summary:Appends an aligned read entry to a fragment store.
..cat:Fragment Store
..signature:appendRead(store, readId, contigId, beginPos, endPos[, pairMatchId])
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId:The id of the read.
..param.contigId:The id of the contig.
..param.beginPos:The begin position of the alignment.
..param.endPos:The end position of the alignment.
..param.pairMatchId:Id of alignedRead pair.
...default:$INVALID_ID$, which corresponds to an unmated read.
..returns:The $alignedReadId$ of the aligned read.
..remarks:This function appends a single aligned read to the @Memvar.FragmentStore#alignedReadStore@.
Note that this really only adds a match.
To generate a global alignment out of all of these matches, use @Function.convertMatchesToGlobalAlignment@.
..see:Function.appendRead
..include:seqan/store.h
*/
template <typename TSpec, typename TConfig, typename TReadId, typename TContigId, typename TPos, typename TPairMatchId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
appendAlignedRead(
        FragmentStore<TSpec, TConfig> & store,
        TReadId const & readId,
        TContigId const & contigId,
        TPos const & beginPos,
        TPos const & endPos,
        TPairMatchId const & pairMatchId)
{
    SEQAN_CHECKPOINT;
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;

    TAlignedReadStoreElement element(length(store.alignedReadStore), readId, contigId, beginPos, endPos);
    element.pairMatchId = pairMatchId;
    appendValue(store.alignedReadStore, element);

    return back(store.alignedReadStore).id;
}

template <typename TSpec, typename TConfig, typename TReadId, typename TContigId, typename TPos>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
appendAlignedRead(
        FragmentStore<TSpec, TConfig> & store,
        TReadId const & readId,
        TContigId const & contigId,
        TPos const & beginPos,
        TPos const & endPos)
{
    SEQAN_CHECKPOINT;
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    return appendAlignedRead(store, readId, contigId, beginPos, endPos, TAlignedReadStoreElement::INVALID_ID);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.appendMatePair
..summary:Appends two paired-end reads to a fragment store.
..cat:Fragment Store
..signature:appendMatePair(store, readId1, readId2)
..signature:appendMatePair(store, readId1, readId2, name1, name2)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId1:The read sequence of the first read.
..param.readId2:The read sequence of the second read.
..param.name1:The read name of the first read.
..param.name2:The read name of the second read.
..returns:The $matePairId$ of the newly appended mate-pair.
..remarks:This function appends two reads to the @Memvar.FragmentStore#readStore@ and @Memvar.FragmentStore#readSeqStore@ 
and a mate-pair entry between both of them to the @Memvar.FragmentStore#matePairStore@.
If names are given, they are appended to the @Memvar.FragmentStore#readNameStore@.
..see:Function.appendRead
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read1, 
	TRead const &read2)
{
	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typename Value<TReadStore>::Type r;
	typename Value<TMatePairStore>::Type mp;
	r.matePairId = length(me.matePairStore);
	mp.readId[0] = length(me.readStore);
	mp.readId[1] = length(me.readStore) + 1;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readStore, r, Generous());
	appendValue(me.matePairStore, mp, Generous());
	appendValue(me.readSeqStore, read1, Generous());
	appendValue(me.readSeqStore, read2, Generous());
	return length(me.matePairStore) - 1;
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TMatePairStore>::Type
appendMatePair(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read1, 
	TRead const &read2, 
	CharString const &name1,
	CharString const &name2)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typename Value<TReadStore>::Type r;
	typename Value<TMatePairStore>::Type mp;
	r.matePairId = length(me.matePairStore);
	mp.readId[0] = length(me.readStore);
	mp.readId[1] = length(me.readStore) + 1;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readStore, r, Generous());
	appendValue(me.matePairStore, mp, Generous());
	appendValue(me.readSeqStore, read1, Generous());
	appendValue(me.readSeqStore, read2, Generous());
	appendValue(me.readNameStore, name1, Generous());
	appendValue(me.readNameStore, name2, Generous());
	return length(me.matePairStore) - 1;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.compactAlignedReads
..summary:Removes invalid aligned reads and rename $alignId$ sequentially beginning with 0.
..cat:Fragment Store
..signature:compactAlignedReads(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:The new size of the @Memvar.FragmentStore#alignedReadStore@.
..remarks:This function removes all entries from @Memvar.FragmentStore#alignedReadStore@ whose $alignId$ equals to $INVALID_ID$ as well as orphan entries
in @Memvar.FragmentStore#alignQualityStore@.
Afterwards the alignIds are renamed sequentially beginning with 0.
This function can be used to remove alignments which are selected by previously setting their id to $INVALID_ID$.
..include:seqan/store.h
*/

// 1. remove aligned reads with invalid ids
// 2. rename ids beginning with 0
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactAlignedReads(FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;

	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Size<TAlignQualityStore>::Type					TAQSize;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Iterator<TAlignQualityStore, Standard>::Type	TAlignQualityIter;
	
	sortAlignedReads(me.alignedReadStore, SortId());
	
	TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
	TAlignQualityIter itAQ = begin(me.alignQualityStore, Standard());
	TAlignQualityIter itAQbegin = itAQ;
	TAQSize aqSize = length(me.alignQualityStore);
	TId newId = 0;
	
	for (; itAR != itARend; ++itAR, ++newId)
	{
		TId id = (*itAR).id;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (id < aqSize)
		{
			*itAQ = *(itAQbegin + id);
			++itAQ;
		}
		(*itAR).id = newId;
	}
	
	resize(me.alignedReadStore, newId, Exact());
	resize(me.alignQualityStore, itAQ - itAQbegin, Exact());
	return newId;
}

/**
.Function.compactPairMatchIds
..summary:Renames $pairMatchId$ sequentially beginning with 0.
..cat:Fragment Store
..signature:compactPairMatchIds(store)
..param.store:The fragment store.
...type:Class.FragmentStore
..returns:The number of pair matches.
..remarks:This function renames the $pairMatchId$ in the @Memvar.FragmentStore#alignedReadStore@ sequentially beginning with 0.
Two read alignments can be identified to be a pair match if they have the same $pairMatchId$.
Please note that paired reads not necessarily have to mapped as a pair match, 
e.g. if they are on different contigs or have the same orientation or a wrong insert size.
..include:seqan/store.h
*/

// rename pair match ids beginning with 0, returns the number of pair matches
template <typename TSpec, typename TConfig>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
compactPairMatchIds(FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;

	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	
	sortAlignedReads(me.alignedReadStore, SortPairMatchId());
	
	TAlignedReadIter itAR = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itARend = end(me.alignedReadStore, Standard());
	if (itAR == itARend) return 0;
	
	TId lastId = (*itAR).pairMatchId;
	TId newId = 0;
	for (; itAR != itARend; ++itAR)
	{
		TId id = (*itAR).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (lastId < id)
		{
			lastId = id;
			++newId;
		}
		(*itAR).pairMatchId = newId;
	}
	return newId + 1;
}

/**
.Function.calculateInsertSizes
..summary:Calculates a string with insert sizes for each pair match.
..cat:Fragment Store
..signature:compactPairMatchIds(insertSizes, store)
..param.insertSizes:The resulting string of insert sizes.
...remarks:This string is accordingly resized and can be addressed by the $pairMatchId$.
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:This function calls @Function.compactPairMatchIds@ first and calculate the insert size for every pair match.
The insert size of a pair match is the outer distance between the two matches.
..include:seqan/store.h
*/

template <typename TLibSizeString, typename TSpec, typename TConfig>
inline void
calculateInsertSizes(TLibSizeString &insertSizes, FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename TFragmentStore::TContigPos						TGPos;

	TAlignedReadIter it = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

	resize(insertSizes, compactPairMatchIds(me), Exact());
	TId lastId = TAlignedRead::INVALID_ID;
	TGPos leftMatePos = 0;
	for (; it != itEnd; ++it)
	{
		TId id = (*it).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) break;	// we assume that invalid ids are at the end of the AlignedReadStore
		if (id != lastId)
		{
			leftMatePos = (*it).beginPos;
			lastId = id;
		} else
			insertSizes[id] = (*it).beginPos - leftMatePos;
	}
}

/**
.Function.getMateNo
..summary:Returns the mate number of read for a given $readId$.
..cat:Fragment Store
..signature:getMateNo(store, readId)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.readId:The read id.
..returns:The mate number (0..first mate, 1..second mate) of the read in its mate-pair or -1 if the read is not paired.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TId>
inline int
getMateNo(FragmentStore<TSpec, TConfig> const &me, TId readId)
{
	typedef FragmentStore<TSpec, TConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename TFragmentStore::TMatePairStore	TMatePairStore;

	typedef typename Value<TReadStore>::Type		TRead;
	typedef typename Value<TMatePairStore>::Type	TMatePair;
	
	if (readId != TRead::INVALID_ID)
	{
		TRead const &r = me.readStore[readId];
		if (r.matePairId != TRead::INVALID_ID)
		{
			TMatePair const &mp = me.matePairStore[r.matePairId];
			if (mp.readId[0] == readId) return 0;
			if (mp.readId[1] == readId) return 1;
		}
	}
	return -1;
}

/**
.Function.calculateMateIndices
..summary:Calculates a string that maps the $readId$ of a read to the $readId$ of its mate.
..cat:Fragment Store
..signature:calculateMateIndices(mateIndices, store)
..param.mateIndices:The resulting string of mate indices.
...remarks:This string is accordingly resized and can be addressed by the $readId$.
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:Entries of reads without a mate contain $INVALID_ID$.
..include:seqan/store.h
*/

// calculate index of the other mate for each pair match
template <typename TMateIndexString, typename TSpec, typename TConfig>
inline void
calculateMateIndices(TMateIndexString &mateIndices, FragmentStore<TSpec, TConfig> &me)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;

	TAlignedReadIter it = begin(me.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(me.alignedReadStore, Standard());

	for (TId idx = 0; it != itEnd; ++it, ++idx)
	{
		TId id = (*it).pairMatchId;
		if (id == TAlignedRead::INVALID_ID) continue;
		if (length(mateIndices) < 2*id + 2)
			resize(mateIndices, 2*id + 2, TAlignedRead::INVALID_ID, Generous());
		mateIndices[2*id + 1 - getMateNo(me, (*it).readId)] = idx;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Class.AlignedReadLayout
..summary:Stores a 2-dimensional visible layout of a multi-read alignment.
..cat:Fragment Store
..signature:AlignedReadLayout

.Memvar.AlignedReadLayout#contigRows
..class:Class.AlignedReadLayout
..summary:2-D multi-read layout
..remarks:Stores for a contig and a row the ids of aligned reads from left to right.
$contigRows[contigId][row]$ stores the $alignId$ of all aligned reads from left to right assigned to the same row
$row$ is the row of the alignment in the multiple sequence alignment and $contigId$ the id of the reference contig.
..include:seqan/store.h
*/
	
struct AlignedReadLayout
{
	typedef String<unsigned>	TRow;
	typedef String<TRow>		TRows;
	typedef String<TRows>		TContigRows;
	
	TContigRows contigRows;			// rows string, each row is a string of ids of alignedReads from left to right
	String<Pair<int> > mateCoords;	// coords of mate pair
};

/**
.Function.layoutAlignment
..summary:Calculates a visible layout of aligned reads.
..cat:Fragment Store
..signature:layoutAlignment(layout, store)
..param.layout:The resulting layout structure.
...type:Class.AlignedReadLayout
..param.store:The fragment store.
...type:Class.FragmentStore
..remarks:For each contig this function layouts all reads in rows from up to down reusing empty row spaces.
..see:Function.printAlignment
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig>
void layoutAlignment(AlignedReadLayout &layout, FragmentStore<TSpec, TConfig> &store)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;

	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename AlignedReadLayout::TRows						TRows;
	typedef typename AlignedReadLayout::TRow						TRow;
	typedef typename Iterator<TRows>::Type							TRowsIter;
	
	// sort matches by increasing begin positions
	sortAlignedReads(store.alignedReadStore, SortBeginPos());
	sortAlignedReads(store.alignedReadStore, SortContigId());

	clear(layout.contigRows);
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());

	for (TId id = 0; it != itEnd; ++it, ++id)
	{
		if ((*it).contigId == TAlignedRead::INVALID_ID) continue;
		if (length(layout.contigRows) <= (*it).contigId)
			resize(layout.contigRows, (*it).contigId + 1);
		
		TRowsIter lit = begin(layout.contigRows[(*it).contigId], Standard());
		TRowsIter litEnd = end(layout.contigRows[(*it).contigId], Standard());
		
		TContigPos beginPos = _min((*it).beginPos, (*it).endPos);
		
		for (; lit != litEnd; ++lit)
		{
			if (empty(*lit)) break;
			TAlignedRead &align = store.alignedReadStore[back(*lit)];
			if (_max(align.beginPos, align.endPos) < beginPos)			// maybe <= would be better
				break;													// but harder to differ between reads borders
		}
			
		if (lit == litEnd)
		{
			TRow s;
			appendValue(s, id);
			appendValue(layout.contigRows[(*it).contigId], s);
		} else
			appendValue(*lit, id);
	}
}

template <typename TStream, typename TFormatTag, typename TContigGaps, typename TReadGaps, typename TAlignedRead, typename TLine>
inline void _printRead(
	TStream &stream, 
	Tag<TFormatTag> const &format,
	AlignedReadLayout &, 
	TContigGaps &,
	TReadGaps &readGaps,
	TAlignedRead &,
	TLine)
{
	write(stream, readGaps, "", format);
}

template <typename TStream, typename TFormatTag, typename TContigGaps, typename TContigName>
inline void _printContig(
	TStream &stream,
	Tag<TFormatTag> const &format,
	AlignedReadLayout &, 
	TContigGaps &contigGaps,
	TContigName const &)
{
	write(stream, contigGaps, "", format);
}

/**
.Function.printAlignment
..summary:Prints a window of the visible layout of reads into a outstream.
..cat:Fragment Store
..signature:printAlignment(stream, format, layout, store, contigId, posBegin, posEnd, lineBegin, lineEnd)
..param.stream:A C++ outstream, e.g. std::cout.
..param.layout:A layout structure created by a previous call of @Function.layoutAlignment@.
...type:Class.AlignedReadLayout
..param.format:Output format.
...type:Tag.File Format.tag.Raw
...remarks: This tag is used for subsequent calls of @Function.write@ for contig and read gaps data structures.
..param.store:The fragment store.
...type:Class.FragmentStore
..param.contigId:The $contigId$ of the affected contig.
..param.posBegin:Window begin position in gap-space.
..param.posEnd:Window end position in gap-space.
..param.lineBegin:Begin line of the window.
..param.lineEnd:End line of the window.
..remarks:The window coordinates ($beginPos$, ...) may be chosen bigger than the layout.
The empty space is then filled with whitespaces.
..see:Function.layoutAlignment
..include:seqan/store.h
*/

template <typename TStream, typename TFormatTag, typename TSpec, typename TConfig, typename TContigId, typename TPos, typename TNum>
void printAlignment(
	TStream &stream, 
	Tag<TFormatTag> const &format,
	AlignedReadLayout &layout, 
	FragmentStore<TSpec, TConfig> &store, 
	TContigId contigId,
	TPos posBegin, TPos posEnd,
	TNum lineBegin, TNum lineEnd)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename TFragmentStore::TReadSeq						TReadSeq;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	typedef typename TContig::TContigSeq							TContigSeq;

	typedef AlignedReadLayout::TRows								TRows;
	typedef typename Value<TRows>::Type								TRow;
	typedef typename Size<TRows>::Type								TRowsSize;
	typedef typename Iterator<TRows>::Type							TRowsIter;
	typedef typename Iterator<TRow>::Type							TRowIter;

	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	typedef Gaps<CharString, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
	
	TContigGaps	contigGaps;
	if ((TId)contigId < length(store.contigStore))
	{
		set(contigGaps.data_source, store.contigStore[contigId].seq);
		set(contigGaps.data_gaps, store.contigStore[contigId].gaps);
//		TContigGaps	contigGaps(store.contigStore[contigId].seq, store.contigStore[contigId].gaps);
		setBeginPosition(contigGaps, posBegin);
		setEndPosition(contigGaps, posEnd);
		_printContig(stream, format, layout, contigGaps, store.contigNameStore[contigId]);
		stream << '\n';
	} else
		stream << '\n';
	
	if ((TId)contigId >= length(layout.contigRows))
		return;
	
	if ((TRowsSize)lineEnd > length(layout.contigRows[contigId])) lineEnd = length(layout.contigRows[contigId]);
	if ((TRowsSize)lineBegin >= (TRowsSize)lineEnd) return;

	TRowsIter lit = begin(layout.contigRows[contigId], Standard()) + lineBegin;
	TRowsIter litEnd = begin(layout.contigRows[contigId], Standard()) + lineEnd;
	TReadSeq readSeq;
	CharString readSeqString;
	setBeginPosition(contigGaps, 0);
	setEndPosition(contigGaps, _unclippedLength(contigGaps));

	for (TNum line = 1; lit < litEnd; ++lit, ++line)
	{
		TRowIter itEnd = end(*lit, Standard());
		TRowIter left = begin(*lit, Standard());
		TRowIter right = itEnd;
		TRowIter mid = itEnd;

		while (left < right)
		{
			mid = left + (right - left) / 2;
			TAlignedRead &align = store.alignedReadStore[*mid];

			if (align.contigId < (TId)contigId || (align.contigId == (TId)contigId && (TPos)_max(align.beginPos, align.endPos) <= posBegin))
				left = mid + 1;	// what we search is in the right part
			else
				right = mid;	//            ...           left part
		}
		
		TPos cursor = posBegin;
		for (; mid < itEnd; ++mid)
		{
//			if (*mid >= lastRead) continue;
			TAlignedRead &align = store.alignedReadStore[*mid];
			if (align.contigId != (TId)contigId) break;

			TReadGaps readGaps(readSeqString, align.gaps);
			TContigPos	left = align.beginPos;
			TContigPos	right = align.endPos;
			TContigPos	cBegin = _min(left, right);
			TContigPos	cEnd = _max(left, right);
			
			if ((TPos)cEnd <= posBegin) continue; // shouldn't occur
			if (posEnd <= (TPos)cBegin) break;
			
			readSeq = store.readSeqStore[align.readId];
			if (left > right)
			{
				reverseComplement(readSeq);
				readSeqString = readSeq;
				toLower(readSeqString);
			} else
				readSeqString = readSeq;
			
			if ((TPos)cBegin < posBegin)
				setBeginPosition(readGaps, posBegin - (TPos)cBegin);
			else
				for (; cursor < (TPos)cBegin; ++cursor)
					stream << ' ';
			
			if (posEnd < (TPos)cEnd)
				setEndPosition(readGaps, posEnd - (TPos)cBegin);
			
			_printRead(stream, format, layout, contigGaps, readGaps, align, line);
			cursor = cEnd;
		}
		stream << '\n';
	}
}

/**
.Function.convertMatchesToGlobalAlignment
..summary:Converts all matches to a multiple global alignment in gap-space.
..cat:Fragment Store
..signature:convertMatchesToGlobalAlignment(store, score)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.score:A score object used by @Function.globalAlignment@ in this function.
..remarks:Before calling this function all $gaps$ structures in @Memvar.FragmentStore#alignedReadStore@ and @Memvar.FragmentStore#contigStore@ must be empty, i.e. there are no gaps in the alignments.
This function iterates over entries in the @Memvar.FragmentStore#alignedReadStore@ and semi-global aligns each read to its contig segments given by begin and end position.
Gaps introduced by these pair-wise alignments are then inserted to the affected contig and reads correspondingly.
..remarks:The invariant that positions in the @Memvar.FragmentStore#alignedReadStore@ are in gap-space holds before (there were no gaps in alignments) and after calling this functions.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TScore, typename TShrinkMatches>
void convertMatchesToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TScore const & score, TShrinkMatches const &)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

	typedef typename TFragmentStore::TReadSeq						TReadSeq;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename TContig::TContigSeq							TContigSeq;
	typedef Align<TReadSeq, ArrayGaps>								TAlign;
	typedef Gaps<TReadSeq, ArrayGaps>								TGaps;

	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
	typedef typename Iterator<TContigGaps>::Type							TContigIter;
	typedef typename Iterator<TReadGaps>::Type								TReadIter;
	
	// sort matches by increasing begin positions
	sortAlignedReads(store.alignedReadStore, SortBeginPos());
	sortAlignedReads(store.alignedReadStore, SortContigId());

	TReadSeq readSeq;
	TId lastContigId = TAlignedRead::INVALID_ID;
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
	TAlignedReadIter firstOverlap = begin(store.alignedReadStore, Standard());
	for (; it != itEnd; ++it)
	{
		TContigPos	left = (*it).beginPos;
		TContigPos	right = (*it).endPos;
		TContigPos	cBegin = _min(left, right);
		TContigPos	cEnd = _max(left, right);
		TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
		TReadGaps	readGaps(readSeq, (*it).gaps);
		
		readSeq = store.readSeqStore[(*it).readId];
		if (left > right)
			reverseComplement(readSeq);
				
		// 1. Calculate pairwise alignment
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), infix(store.contigStore[(*it).contigId].seq, cBegin, cEnd));
		assignSource(row(align, 1), readSeq);
		if (IsSameType<TShrinkMatches, True>::VALUE)
		    globalAlignment(align, score, AlignConfig<true, false, false, true>(), Gotoh());
        else
		    globalAlignment(align, score);

		// 2. Skip non-overlapping matches
		cBegin = positionSeqToGap(contigGaps, cBegin);
		// Increment cBegin by the number of gaps in front of the read
		// row in the global alignment computed above.
		cBegin += beginPosition(row(align, 1));
		if (lastContigId != (*it).contigId)
		{
			firstOverlap = it;
			lastContigId = (*it).contigId;
		} else
			while (firstOverlap != it && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
				++firstOverlap;

		// 3. Iterate over alignment
		setBeginPosition(contigGaps, cBegin);
		
		TContigIter cIt = begin(contigGaps);
		TReadIter rIt = begin(readGaps);
		typename Iterator<TGaps>::Type it1 = begin(row(align, 0));
		typename Iterator<TGaps>::Type it2 = begin(row(align, 1));

        unsigned beginLocalContigGaps = beginPosition(row(align, 0));
		// Heuristic (hack) for gaps in the beginning, so the following 
		// does not happen:
		//
		// contig: XXXX------AAA
		//             CCC---AAA
		//                CCCAAA
		//
		// But instead, the C's are all aligned.
		if (beginLocalContigGaps > 0u) {
			unsigned i = beginLocalContigGaps;
			do {
				goPrevious(cIt);
				i -= 1;
				cBegin -= 1;
			} while (isGap(cIt) && i > 0u);
			if (!isGap(cIt)) {
				goNext(cIt);
				cBegin += 1;
			}
		}
        
		for (; !atEnd(cIt) && !atEnd(it1); goNext(cIt), goNext(rIt))
		{
			bool isGapContig = isGap(cIt);
			bool isGapLocalContig = (beginLocalContigGaps > 0) ? true : isGap(it1);
			if (isGapContig != isGapLocalContig)
			{
				if (isGapContig)
				{
					// *** gap in contig of the global alignment ***
					// copy exisiting contig gap
					insertGaps(rIt, 1);
					continue;
				} else
				{
					// *** gap in contig of the pairwise alignment ***
					// insert padding gaps in contig and reads
					TContigPos insPos = cIt.current.gapPos;
					insertGaps(cIt, 1);
					for (TAlignedReadIter j = firstOverlap; j != it; ++j)
					{
						TContigPos rBegin = _min((*j).beginPos, (*j).endPos);
						TContigPos rEnd = _max((*j).beginPos, (*j).endPos);
						if (rBegin < insPos && insPos < rEnd)
						{
							if (rBegin < insPos)
							{
								TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
								insertGap(gaps, insPos - rBegin);
							} else
							{
								// shift beginPos if insertion was at the front of the read
								if ((*j).beginPos < (*j).endPos)
									++(*j).beginPos;
								else
									++(*j).endPos;
							}
							// shift endPos as the alignment was elongated or shifted
							if ((*j).beginPos < (*j).endPos)
								++(*j).endPos;
							else
								++(*j).beginPos;
						} else if (insPos <= rBegin) {
                            ++(*j).endPos;
                            ++(*j).beginPos;
                        }
					}
				}
			}
			if (isGap(it2))
			{
				// *** gap in read of pairwise alignment ***
				// copy gaps from alignment
				insertGaps(rIt, 1);
			}
            if (beginLocalContigGaps == 0)
                goNext(it1);
            else
                beginLocalContigGaps -= 1;
			goNext(it2);
		}

		// store new gap-space alignment borders
		cEnd = cBegin + length(readGaps);
		if (left < right)
		{
			(*it).beginPos = cBegin;
			(*it).endPos = cEnd;
		} else
		{
			(*it).beginPos = cEnd;
			(*it).endPos = cBegin;
		}
/*		
//		if (interesting)
		{
			String<String<unsigned> > layout;
			layoutAlignment(layout, store, (*it).contigId);
			std::cout << store.readNameStore[(*it).readId] << std::endl;
			std::cout << readGaps << '\t' << cBegin << '\t' << cEnd << std::endl << std::endl;
			printAlignment(std::cout, layout, store, (*it).contigId, (int)cBegin-20, (int)cEnd+20, 0, 40, 1 + (it - begin(store.alignedReadStore, Standard())));
//			getc(stdin);
		}
*/
//		if (store.readNameStore[(*it).readId] == "read3305")
//			return;
	}
}

/**
.Function.convertPairWiseToGlobalAlignment
..summary:Converts pair-wise alignments to a multiple global alignment.
..cat:Fragment Store
..signature:convertPairWiseToGlobalAlignment(store, pairwiseContigGaps)
..param.store:The fragment store.
...type:Class.FragmentStore
..param.pairwiseContigGaps:A string of anchored contig gaps for every pairwise alignment.
..remarks:Before calling this function the $gaps$ structures in the @Memvar.FragmentStore#contigStore@ must be empty, i.e. there are no gaps in the contig.
The pairwise alignment gaps of the reads are stored in the $gaps$ structure in the @Memvar.FragmentStore#alignedReadStore@, whereas the pairwise alignment gaps of the contig are stored in the $pairwiseContigGaps$ string.
..remarks:After calling this functions all positions in the @Memvar.FragmentStore#alignedReadStore@ are in gap-space.
..include:seqan/store.h
*/

template <typename TSpec, typename TConfig, typename TContigGapsString>
void convertPairWiseToGlobalAlignment(FragmentStore<TSpec, TConfig> &store, TContigGapsString &gaps)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

	// stores
	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TContigStore					TContigStore;

	// store elements
	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TReadSeqStore>::Type						TReadSeq;
	typedef typename Value<TContigStore>::Type						TContig;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignedReadIter;
	typedef typename Id<TAlignedRead>::Type							TId;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	typedef typename TContig::TContigSeq							TContigSeq;
	typedef Align<TReadSeq, ArrayGaps>								TAlign;
	typedef Gaps<TReadSeq, ArrayGaps>								TGaps;

	// gap structures
	typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename TContig::TGapAnchors> >			TContigGapsGlobal;
	typedef Gaps</*TContigSeq*/Nothing, AnchorGaps<typename Value<TContigGapsString>::Type> >	TContigGapsPW;
	typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >			TReadGaps;
	
	// gap iterators
	typedef typename Iterator<TContigGapsGlobal>::Type								TContigGlobalIter;	
	typedef typename Iterator<TContigGapsPW>::Type									TContigPWIter;
	typedef typename Iterator<TReadGaps>::Type										TReadIter;
	
	// sort matches by increasing begin positions
	sortAlignedReads(store.alignedReadStore, SortBeginPos());
	sortAlignedReads(store.alignedReadStore, SortContigId());

	TReadSeq readSeq;
	TId lastContigId = TAlignedRead::INVALID_ID;
	TAlignedReadIter it = begin(store.alignedReadStore, Standard());
	TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
	TAlignedReadIter firstOverlap = begin(store.alignedReadStore, Standard());
	for (; it != itEnd; ++it)
	{
		TContigPos	left = (*it).beginPos;
		TContigPos	right = (*it).endPos;
		TContigPos	cBegin = _min(left, right);
		TContigPos	cEnd = _max(left, right);
		
		// 1. Initialize gap structures
		TContigGapsGlobal	contigGapsGlobal(/*store.contigStore[(*it).contigId].seq, */store.contigStore[(*it).contigId].gaps);
		TContigGapsPW		contigGapsPW(/*store.contigStore[(*it).contigId].seq, */gaps[(*it).id]);
		TReadGaps			readGaps(store.readSeqStore[(*it).readId], (*it).gaps);
		
		// 2. Skip non-overlapping matches
		cBegin = positionSeqToGap(contigGapsGlobal, cBegin);
		if (lastContigId != (*it).contigId)
		{
			firstOverlap = it;
			lastContigId = (*it).contigId;
		} else
			while (firstOverlap != it && _max((*firstOverlap).beginPos, (*firstOverlap).endPos) <= cBegin)
				++firstOverlap;

		// 3. Iterate over alignment
		setBeginPosition(contigGapsGlobal, cBegin);

		TContigGlobalIter cIt = begin(contigGapsGlobal);
		TContigPWIter pIt = begin(contigGapsPW);
		TReadIter rIt = begin(readGaps);
		
		for (; /*!atEnd(cIt) && */!atEnd(rIt); goNext(cIt), goNext(rIt))
		{
			bool isGapContig = isGap(cIt);
			bool isGapLocalContig = isGap(pIt);
			if (isGapContig != isGapLocalContig)
			{
				if (isGapContig)
				{
					// *** gap in contig of the global alignment ***
					// copy exisiting contig gap
					insertGaps(rIt, 1);
					continue;
				} else
				{
					// *** gap in contig of the pairwise alignment ***
					// insert padding gaps in contig and reads
					TContigPos insPos = cIt.current.gapPos;
					insertGaps(cIt, 1);
					for (TAlignedReadIter j = firstOverlap; j != it; ++j)
					{
                        
						TContigPos rBegin = _min((*j).beginPos, (*j).endPos);
						TContigPos rEnd = _max((*j).beginPos, (*j).endPos);
						if (rBegin < insPos && insPos < rEnd)
						{
							if (rBegin < insPos)
							{
								TReadGaps gaps(store.readSeqStore[(*j).readId], (*j).gaps);
								insertGap(gaps, insPos - rBegin);
							} else
							{
								// shift beginPos if insertion was at the front of the read
								if ((*j).beginPos < (*j).endPos)
									++(*j).beginPos;
								else
									++(*j).endPos;
							}
							// shift endPos as the alignment was elongated or shifted
							if ((*j).beginPos < (*j).endPos)
								++(*j).endPos;
							else
								++(*j).beginPos;
						}
					}
				}
			}
			goNext(pIt);
		}

		// store new gap-space alignment borders
		cEnd = cBegin + length(readGaps);
		if (left < right)
		{
			(*it).beginPos = cBegin;
			(*it).endPos = cEnd;
		} else
		{
			(*it).beginPos = cEnd;
			(*it).endPos = cBegin;
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
