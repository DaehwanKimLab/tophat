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

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	// dfs order
	struct Preorder_;
	struct Postorder_;

	template <typename TDfsOrder = Postorder_, typename THideEmptyEdges = True>
	struct VSTreeIteratorTraits {
		typedef TDfsOrder DfsOrder;
		typedef THideEmptyEdges HideEmptyEdges;
	};

/**
.Tag.Preorder:
..summary:Preorder depth-first search.
..cat:Index
..signature:Preorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
..see:Tag.Postorder
..include:seqan/index.h
*/

/**
.Tag.Postorder:
..summary:Postorder depth-first search.
..cat:Index
..signature:Postorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
..see:Tag.Preorder
..include:seqan/index.h
*/

/**
.Tag.PreorderEmptyEdges:
..summary:Preorder depth-first search in a suffix tree with leaves for every suffix.
..cat:Index
..signature:PreorderEmptyEdges
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
Empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.PostorderEmptyEdges
..see:Tag.Preorder
..include:seqan/index.h
*/

/**
.Tag.PostorderEmptyEdges:
..summary:Postorder depth-first search in a suffix tree with leaves for every suffix.
..cat:Index
..signature:PostorderEmptyEdges
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
Empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.Postorder
..include:seqan/index.h
*/

/**
.Tag.EmptyEdges:
..summary:Consider a suffix tree with leaves for every suffix.
..cat:Index
..signature:EmptyEdges
..remarks:When given as a second parameter in @Function.goDown@, empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.HideEmptyEdges
..include:seqan/index.h
*/

/**
.Tag.HideEmptyEdges:
..summary:Consider a suffix tree with no empty edges (default behaviour).
..cat:Index
..signature:HideEmptyEdges
..remarks:When given as a second parameter in @Function.goDown@, only non-empty edges are traversed.
..see:Tag.EmptyEdges
..include:seqan/index.h
*/

	// predefined iterator traits
	struct Preorder:			VSTreeIteratorTraits<Preorder_,  True> {};
	struct Postorder:			VSTreeIteratorTraits<Postorder_, True> {};
	struct PreorderEmptyEdges:	VSTreeIteratorTraits<Preorder_,  False> {};	// also iterate over
	struct PostorderEmptyEdges:	VSTreeIteratorTraits<Postorder_, False> {};	// empty edges (with $-label)

	// traits for TopDown iterators (w/o ParentLinks) for which postorder/preorder is ignored
	struct HideEmptyEdges:		VSTreeIteratorTraits<Postorder_, True> {};
	struct EmptyEdges:			VSTreeIteratorTraits<Postorder_, False> {};	// empty edges (with $-label)
	
	// MultiMems are more specialized MaxRepeats
	template <typename TSpec = void>
	struct MaxRepeats_;	// base class
	struct MultiMems_;	// subclass tag



	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;

		// top down traversal iterators
		template <typename TSpec = Preorder>
		struct TopDown {};					// starts in the suffix tree root and can go down and go right

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};			// .. can also go up

		// bottom up traversal iterators
		template <typename TSpec = Postorder>
		struct BottomUp {};					// starts in the first node of a depth-first-search and can go next

			struct	SuperMaxRepeats;					// maximal repeat and not part of a longer repeat
			struct	SuperMaxRepeatsFast;
			struct	Mums;								// Maximal Unique Match (unique in every sequence)

			typedef MaxRepeats_<void>		MaxRepeats;	// maximal repeat
			struct	MaxRepeatOccurrences;
			typedef MaxRepeats_<MultiMems_> MultiMems;	// Multiple Maximal Exact Match
			struct	MultiMemOccurences;					// i.e. maximal match over different sequences


/**
.Metafunction.GetVSTreeIteratorTraits:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:GetVSTreeIteratorTraits<TIterator>::Type
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
..include:seqan/index.h
*/

	template <typename TIterator>
	struct GetVSTreeIteratorTraits:
		DeepestSpec<TIterator> {};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSize>
	struct VertexEsa {
		Pair<TSize> range;			// current SA interval of hits (unique node identifier)
		TSize		parentRight;	// right boundary of parent node's range (allows to go right)

		VertexEsa() {}

		VertexEsa(MinimalCtor):
			range(0,0),
			parentRight(0) {}

		VertexEsa(TSize otherRangeLeft, TSize otherRangeRight, TSize otherParentRight):
			range(Pair<TSize>(otherRangeLeft, otherRangeRight)),
			parentRight(otherParentRight) {}

		VertexEsa(Pair<TSize> const &otherRange, TSize otherParentRight):
			range(otherRange),
			parentRight(otherParentRight) {}

		VertexEsa(VertexEsa const &other):
			range(other.range),
			parentRight(other.parentRight) {}
	};
	
	template <typename TSize>
	inline bool operator==(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
	{
		return a.range == b.range;
	}

	template <typename TSize>
	inline bool operator!=(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
	{
		return a.range != b.range;
	}

//////////////////////////////////////////////////////////////////////////////
///.Metafunction.VertexDescriptor.param.T.type:Spec.IndexEsa

	template < typename TText, typename TSpec >
	struct VertexDescriptor< Index<TText, IndexEsa<TSpec> > > {
		typedef typename Size< Index<TText, IndexEsa<TSpec> > >::Type TSize;
		typedef VertexEsa<TSize> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	struct ArrayGaps;

	template <typename TSource, typename TSpec>
	class Align;


//////////////////////////////////////////////////////////////////////////////
// ESA fibres

/**
.Tag.ESA Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of an @Spec.IndexEsa.ESA@ index.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of an Enhanced Suffix Array based @Spec.IndexEsa.Index@.
..cat:Index

..tag.EsaText:The original text the index should be based on.

..tag.EsaRawText:The raw text the index is really based on.
...remarks:$EsaText$ and $EsaRawText$ fibres are equal by default.
They differ if the index text is a set of strings. Then, raw text is the concatenation of all strings in this set.

..tag.EsaSA:The suffix array.
...remarks:The suffix array contains the indices of all suffices of $EsaRawText$ in lexicographical order.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of the @Metafunction.SAValue@ of $TIndex$.

..tag.EsaLcp:The lcp table.
...remarks:The lcp table contains the lcp-value of two adjacent suffices in the suffix array $EsaSA$.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.EsaChildtab:The child table.
...remarks:The child table contains structural information of the suffix tree (see Abhouelda et al.).
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.EsaBwt:The Burrows-Wheeler table.
...remarks:The Burrows-Wheeler table contains the Burrows-Wheeler transformation of $EsaRawText$.
The entries are the characters left of the corresponding suffix in the suffix array $EsaSA$.
...remarks:@Metafunction.Fibre@ returns the same type for $EsaRawText$ and for $EsaBwt$.

..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.IndexEsa
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.ESA Index Fibres

	typedef FibreText		EsaText;
	typedef FibreRawText	EsaRawText;
	typedef FibreSA		EsaSA;
	typedef FibreRawSA		EsaRawSA;
	typedef FibreSae		EsaSae;
	typedef FibreLcp		EsaLcp;
	typedef FibreLcpe		EsaLcpe;
	typedef FibreChildtab	EsaChildtab;
	typedef FibreBwt		EsaBwt;


//////////////////////////////////////////////////////////////////////////////
// ESA index

/**
.Spec.IndexEsa:
..summary:An index based on an enhanced suffix array.
..cat:Index
..general:Class.Index
..signature:Index<TText, IndexEsa<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array (see @Tag.ESA Index Fibres.EsaSA@), a lcp table (see @Tag.ESA Index Fibres.EsaLcp@), etc.
..remarks:This index can be accessed as a Suffix Tree using the @Spec.VSTree Iterator@ classes.
..include:seqan/index.h
*/

/*
	already defined in index_base.h

	template <typename TSpec = void>
	struct IndexEsa;
*/

	template < typename TText, typename TSpec >
	class Index<TText, IndexEsa<TSpec> > {
	public:
		Holder<typename Fibre<Index, EsaText>::Type>	text;
		typename Fibre<Index, EsaSA>::Type				sa;			// suffix array 
		typename Fibre<Index, EsaLcp>::Type			lcp;		// longest-common-prefix table
		typename Fibre<Index, EsaLcpe>::Type			lcpe;		// extended lcp table
		typename Fibre<Index, EsaChildtab>::Type		childtab;	// child table (tree topology)
		typename Fibre<Index, EsaBwt>::Type			bwt;		// burrows-wheeler table
		typename Cargo<Index>::Type						cargo;		// user-defined cargo

		Index() {}

		Index(Index &other):
			text(other.text),
			sa(other.sa),
			lcp(other.lcp),
			lcpe(other.lcpe),
			childtab(other.childtab),
			bwt(other.bwt),
			cargo(other.cargo) {}

		Index(Index const &other):
			text(other.text),
			sa(other.sa),
			lcp(other.lcp),
			lcpe(other.lcpe),
			childtab(other.childtab),
			bwt(other.bwt),
			cargo(other.cargo) {}

		template <typename TText_>
		Index(TText_ &_text):
			text(_text) {}

		template <typename TText_>
		Index(TText_ const &_text):
			text(_text) {}
	};

//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpec >
	void _indexRequireTopDownIteration(Index<TText, IndexEsa<TSpec> > &index) 
	{
		indexRequire(index, EsaSA());
		indexRequire(index, EsaLcp());
		indexRequire(index, EsaChildtab());
	}

	template < typename TText, typename TSpec >
	void _indexRequireBottomUpIteration(Index<TText, IndexEsa<TSpec> > &index) 
	{
		indexRequire(index, EsaSA());
		indexRequire(index, EsaLcp());
	}

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline void clear(Index<TText, IndexEsa<TSpec> > &index) {
		clear(getFibre(index, EsaSA()));
		clear(getFibre(index, EsaLcp()));
		clear(getFibre(index, EsaLcpe()));
		clear(getFibre(index, EsaChildtab()));
		clear(getFibre(index, EsaBwt()));
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	

		bool result = true;
		if ((!open(getFibre(index, EsaText()), toCString(name), openMode)) && 
			(!open(getFibre(index, EsaText()), fileName, openMode))) 
			result = false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, EsaSA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	open(getFibre(index, EsaLcp()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	open(getFibre(index, EsaChildtab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	open(getFibre(index, EsaBwt()), toCString(name), openMode);
		return result;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, EsaText()), toCString(name), openMode)) && 
			(!save(getFibre(index, EsaText()), fileName, openMode))) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, EsaSA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	save(getFibre(index, EsaLcp()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	save(getFibre(index, EsaChildtab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	save(getFibre(index, EsaBwt()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, IndexEsa<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
	}

}

#endif
