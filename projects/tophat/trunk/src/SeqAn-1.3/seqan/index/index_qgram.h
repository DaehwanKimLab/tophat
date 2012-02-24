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

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// q-gram index fibres

/**
.Tag.QGram Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.IndexQGram.q-gram Index@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a @Spec.IndexQGram.q-gram Index@.
..cat:Index

..tag.QGramText:The original text the index should be based on.

..tag.QGram_RawText:The concatenation of all text sequences.
...remarks:$QGramText$ and $QGram_RawText$ fibres are equal by default.
They differ if the index text is a set of strings. Then, raw text is the concatenation of all strings in this set.

..tag.QGramSA:The suffix array.
...remarks:Contains all occurrences of q-grams, s.t. the occurrences of a single q-gram are stored in a contiguous block (q-gram bucket).
q-grams exceeding the end of the text are ignored.
The beginning of each bucket can be determined by the q-gram directory ($QGramDir$, see below).
...remarks:It corresponds to a suffix array which is sorted by the first q-gram.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of the @Metafunction.SAValue@ of $TIndex$.

..tag.QGramDir:The directory/hash table.
...remarks:The directory contains for every possible q-gram hash value the start index of the q-gram bucket.
A q-gram bucket is a contiguous interval in the suffix array ($QGramSA$, see above).
Each suffix in this interval begins with the same q-gram.
The end index is the start index of the next bucket.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.QGramBucketMap:Maps q-gram hashes to buckets.
This fibre is used by the @Spec.OpenAddressing@ index and stores all parameters of the open addressing hash function and hash value occupancy in the QGramDir fibre.
In contrast to @Spec.OpenAddressing@, @Spec.IndexQGram@ uses a trivial 1-to-1 mapping from q-gram hash values to buckets.
For that index the fibre is of type @Tag.Nothing@.

..tag.QGramCounts:The counts array.
...remarks:Contains the numbers of occurrences per sequence of each q-gram, s.t. the numbers of the same q-gram are stored in a contiguous block (q-gram count bucket).
A bucket contains entries (seqNo,count) of sequences with at least one q-gram occurrence. q-grams exceeding the end of the text are ignored.
The beginning of each count bucket can be determined by the q-gram counts directory ($QGramCountsDir$, see below).
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of the @Metafunction.SAValue@ of $TIndex$.

..tag.QGramCountsDir:The counts directory.
...remarks:The counts directory contains for every possible q-gram hash value the start index of the q-gram count bucket.
A q-gram count bucket is a contiguous interval in the counts array ($QGramCounts$, see above).
The end index is the start index of the next bucket.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.QGramShape:The shape the index is based on.
...remarks:The q-gram index needs an underlying @Class.Shape@. This shape can be gapped or ungapped.
The number of '1's (relevant positions) in the shape determines $q$ and the size of the directory table.
...remarks:Dynamic shapes (@Spec.SimpleShape@, @Spec.GenericShape@, ...) must be initialized before the index can be used.

..tag.QGramSADir:The union of suffix array and directory.
...remarks:In most applications a q-gram index consisting of both of these table is required.
To efficiently create them at once use this tag for @Function.indexRequire@ or @Function.indexCreate@.
 
..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.IndexQGram
..include:seqan/index.h
*/

	struct FibreDir_;			// directory/hash table, contains start indices of buckets
	struct FibreSADir_;		// identifies algorithm to construct both SA and directory at once
	struct Fibre_Shape_;		// underlying shape
	struct FibreCounts_;		// counts each q-gram
	struct FibreCountsDir_;	// directory for counts buckets
	struct FibreBucketMap_;	// stores a q-gram hash value for each directory entry (-1 if empty)

	typedef Tag<FibreDir_> const		FibreDir;
	typedef Tag<FibreSADir_> const		FibreSADir;
	typedef Tag<Fibre_Shape_> const		FibreShape;
	typedef Tag<FibreCounts_> const	FibreCounts;
	typedef Tag<FibreCountsDir_> const	FibreCountsDir;
	typedef Tag<FibreBucketMap_> const	FibreBucketMap;

//////////////////////////////////////////////////////////////////////////////

	typedef FibreText		QGramText;
	typedef FibreRawText	QGram_RawText;
	typedef FibreSA		QGramSA;
	typedef FibreRawSA		QGramRawSA;
	typedef FibreDir		QGramDir;
	typedef FibreSADir		QGramSADir;
	typedef FibreShape		QGramShape;
	typedef FibreCounts	QGramCounts;
	typedef FibreCountsDir	QGramCountsDir;
	typedef FibreBucketMap	QGramBucketMap;

//////////////////////////////////////////////////////////////////////////////
// q-gram index

/**
.Spec.IndexQGram:
..summary:An index based on an array of sorted q-grams.
..cat:Index
..general:Class.Index
..signature:Index<TText, IndexQGram<TShapeSpec[, TSpec]> >
..param.TText:The text type.
...type:Class.String
..param.TShapeSpec:The @Class.Shape@ specialization type.
...note:This can be either a $TSpec$ argument (e.g. $SimpleShape$) or a complete @Class.Shape@ class (e.g. Shape<Dna, SimpleShape>).
..param.TSpec:The specializing type.
...default:Default
...type:Spec.OpenAddressing
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array sorted by the first q characters (see @Tag.QGram Index Fibres.QGramSA@) and a q-gram directory (see @Tag.QGram Index Fibres.QGramDir@).
..include:seqan/index.h
*/

	template < typename TShapeSpec, typename TSpec = Default >
	struct IndexQGram {};

	// use the index value type as shape value type
	template < typename TObject, typename TShapeSpec, typename TSpec >
	struct Fibre< Index<TObject, IndexQGram<TShapeSpec, TSpec> >, FibreShape> 
	{
		typedef Index< TObject, IndexQGram<TShapeSpec, TSpec> >	TIndex;
		typedef Shape< typename Value<TIndex>::Type, TShapeSpec >	Type;
	};

	// allow different value types for the shape
	template < typename TObject, typename TShapeValue, typename TShapeSpec, typename TSpec >
	struct Fibre< Index<TObject, IndexQGram<Shape<TShapeValue, TShapeSpec>, TSpec> >, FibreShape> 
	{
		typedef Shape<TShapeValue, TShapeSpec>	Type;
	};

	template < typename TObject, typename TShapeSpec, typename TSpec >
	struct Fibre< Index<TObject, IndexQGram<TShapeSpec, TSpec> >, FibreBucketMap>
	{
		typedef Nothing Type;
	};

#ifdef PLATFORM_WINDOWS_VS
#pragma warning( push )
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( disable: 4521 )
// Disable warning C4522 locally (multiple assignment operators).
#pragma warning( disable: 4522 )
#endif  // PLATFORM_WINDOWS_VS

	template < typename TObject, typename TShapeSpec, typename TSpec >
	class Index<TObject, IndexQGram<TShapeSpec, TSpec> > {
	public:
		typedef typename Fibre<Index, QGramText>::Type			TText;
		typedef typename Fibre<Index, QGramSA>::Type			TSA;
		typedef typename Fibre<Index, QGramDir>::Type			TDir;
		typedef typename Fibre<Index, QGramCounts>::Type		TCounts;
		typedef typename Fibre<Index, QGramCountsDir>::Type	TCountsDir;
		typedef typename Fibre<Index, QGramShape>::Type		TShape;
		typedef typename Fibre<Index, QGramBucketMap>::Type	TBucketMap;
		typedef typename Cargo<Index>::Type						TCargo;
		typedef typename Size<Index>::Type						TSize;

		Holder<TText>	text;		// underlying text
		TSA				sa;			// suffix array sorted by the first q chars
		TDir			dir;		// bucket directory
		TCounts			counts;		// counts each q-gram per sequence
		TCountsDir		countsDir;	// directory for count buckets
		TShape			shape;		// underlying shape
		TCargo			cargo;		// user-defined cargo
		TBucketMap		bucketMap;	// bucketMap table (used by open-addressing index)
		TSize			stepSize;	// store every <stepSize>'th q-gram in the index

		Index():
			stepSize(1) {}

		Index(Index &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			stepSize(1) {}

		Index(Index const &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			stepSize(1) {}

		template <typename TText_>
		Index(TText_ &_text):
			text(_text),
			stepSize(1) {}

		template <typename TText_>
		Index(TText_ const &_text):
			text(_text),
			stepSize(1) {}

		template <typename TText_, typename TShape_>
		Index(TText_ &_text, TShape_ const &_shape):
			text(_text),
			shape(_shape),
			stepSize(1) {}

		template <typename TText_, typename TShape_>
		Index(TText_ const &_text, TShape_ const &_shape):
			text(_text),
			shape(_shape),
			stepSize(1) {}
	};

#ifdef PLATFORM_WINDOWS_VS
// Reset warning state to previous values for C4521, C4522.
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS

    template < typename TText, typename TShapeSpec, typename TSpec >
    struct Value< Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, QGram_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TShapeSpec, typename TSpec >
    struct Size< Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, QGram_RawText >::Type >::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TShapeSpec, typename TSpec >
	struct DefaultIndexCreator<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA> {
        typedef Default Type;
    };

//////////////////////////////////////////////////////////////////////////////
// counts array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreCounts> {
		typedef String<
				Pair<
					typename Size< TText >::Type,
					typename Size< Index<TText, TSpec> >::Type
				>,
				typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type 
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreDir>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreDir) {
		return index.dir;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreDir>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreDir) {
		return index.dir;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreCounts>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreCounts) {
		return index.counts;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreCounts>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreCounts) {
		return index.counts;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreCountsDir>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreCountsDir) {
		return index.countsDir;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreCountsDir>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreCountsDir) {
		return index.countsDir;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBucketMap>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreBucketMap) {
		return index.bucketMap;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBucketMap>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreBucketMap) {
		return index.bucketMap;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreShape>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreShape) {
		return index.shape;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreShape>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreShape) {
		return index.shape;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexDir
..summary:Shortcut for $getFibre(.., QGramDir)$.
..cat:Index
..signature:indexDir(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:A reference to the @Tag.QGram Index Fibres.QGramDir@ fibre (q-gram directory).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreDir>::Type & 
	indexDir(Index<TText, TSpec> &index) { 
		return getFibre(index, FibreDir()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreDir>::Type & 
	indexDir(Index<TText, TSpec> const &index) { 
		return getFibre(index, FibreDir()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.dirAt
..summary:Shortcut for $value(indexDir(..), ..)$.
..cat:Index
..signature:dirAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreDir>::Type>::Type 
	dirAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreDir()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreDir>::Type>::Type 
	dirAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreDir()), i);
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexCounts
..summary:Shortcut for $getFibre(.., QGramCounts)$.
..cat:Index
..signature:indexCounts(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:A reference to the @Tag.QGram Index Fibres.QGramCounts@ fibre (counts array).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreCounts>::Type & 
	indexCounts(Index<TText, TSpec> &index) {
		return getFibre(index, FibreCounts()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreCounts>::Type & 
	indexCounts(Index<TText, TSpec> const &index) {
		return getFibre(index, FibreCounts()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexCountsDir
..summary:Shortcut for $getFibre(.., QGramCountsDir)$.
..cat:Index
..signature:indexCountsDir(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:A reference to the @Tag.QGram Index Fibres.QGramCountsDir@ fibre (counts directory).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreCountsDir>::Type & 
	indexCountsDir(Index<TText, TSpec> &index) {
		return getFibre(index, FibreCountsDir()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreCountsDir>::Type & 
	indexCountsDir(Index<TText, TSpec> const &index) {
		return getFibre(index, FibreCountsDir()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexBucketMap
..summary:Shortcut for $getFibre(.., QGramBucketMap)$.
..cat:Index
..signature:indexBucketMap(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:A reference to the @Tag.QGram Index Fibres.QGramBucketMap@ fibre (maps q-gram hashes to buckets).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBucketMap>::Type & 
	indexBucketMap(Index<TText, TSpec> &index) {
		return getFibre(index, FibreBucketMap()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBucketMap>::Type & 
	indexBucketMap(Index<TText, TSpec> const &index) {
		return getFibre(index, FibreBucketMap()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexShape
..summary:Shortcut for $getFibre(.., QGramShape)$.
..cat:Index
..signature:indexShape(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexQGram
..returns:Returns a reference to the @Class.Shape@ object of a q-gram index.
Formally, this is a reference to the @Tag.QGram Index Fibres.QGramShape@ fibre.
...type:Class.Shape
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreShape>::Type & 
	indexShape(Index<TText, TSpec> &index) { 
		return getFibre(index, FibreShape()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreShape>::Type & 
	indexShape(Index<TText, TSpec> const &index) { 
		return getFibre(index, FibreShape()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getStepSize
..summary:Return the q-gram step size used for index creation.
..cat:Index
..signature:getStepSize(index)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..returns:The step size. If $x$ is returned every $x$'th q-gram is stored in the index.
..include:seqan/index.h
*/

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline typename Size<TText>::Type 
	getStepSize(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index)
	{
		return (index.stepSize != 0)? index.stepSize: length(indexShape(index));
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.setStepSize
..summary:Change the q-gram step size used for index creation.
..cat:Index
..signature:setStepSize(index, stepSize)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.stepSize:Store every $stepSize$'th q-gram in the index.
..remarks:The default step size of a q-gram index is 1, which corresponds to all overlapping q-grams.
To take effect of changing the $stepSize$ the q-gram index should be empty or recreated. 
..remarks:A $stepSize$ of 0 corresponds to $stepSize=length(indexShape(index))$, i.e. all non-overlapping q-grams.
..see:Function.getStepSize
..include:seqan/index.h
*/

	template <typename TText, typename TShapeSpec, typename TSpec, typename TSize>
	inline void
	setStepSize(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, TSize stepSize)
	{
		index.stepSize = stepSize;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TIndex>
	inline __int64 _fullDirLength(TIndex const &index) 
	{
		typedef typename Fibre<TIndex, FibreShape>::Type	TShape;
		typedef typename Host<TShape>::Type					TTextValue;
		return _intPow((__int64)ValueSize<TTextValue>::VALUE, weight(indexShape(index))) + 1;
	}

	template <typename TIndex>
	inline __int64 _fullDir2Length(TIndex const &index) 
	{
		typedef typename Fibre<TIndex, FibreShape>::Type	TShape;
		typedef typename Host<TShape>::Type					TTextValue;
		return (_intPow(
					(__int64)ValueSize<TTextValue>::VALUE,
					weight(indexShape(index)) + 1) - 1)
				/ ((unsigned)ValueSize<TTextValue>::VALUE - 1) + 1;
	}


    //////////////////////////////////////////////////////////////////////////////
    // QGramLess
	//
	// compare two q-grams of a given text (q-grams can be smaller than q)
    template < typename TSAValue, typename TText >
	struct QGramLess_ : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin, _end;
		typename Size<TText>::Type _q;

		template <typename TSize>
		QGramLess_(TText &text, TSize q): 
			_begin(begin(text, Standard())),
			_end(end(text, Standard())),
			_q(q) {}

		// skip the first <offset> characters
		template <typename TSize1, typename TSize2>
		QGramLess_(TText &text, TSize1 q, TSize2 offset): 
			_begin(begin(text, Standard()) + offset),
			_end(end(text, Standard())),
			_q(q) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			if (a == b) return false;
			TIter itA = _begin + a;
			TIter itB = _begin + b;
			if (a <= b) {
				TIter itEnd = itB + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct QGramLess_<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q;

		template <typename TSize>
		QGramLess_(StringSet<TString, TSpec> const &text, TSize q): 
			_stringSet(text),
			_q(q) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a);
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b);
			if (suffixLength(a, _stringSet) > suffixLength(b, _stringSet)) {
				TIter _end = end(_stringSet[getValueI1(b)], Standard());
				TIter itEnd = itB + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter _end = end(_stringSet[getValueI1(a)], Standard());
				TIter itEnd = itA + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessOffset
	//
	// compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
	struct QGramLessOffset_ :
		QGramLess_<TSAValue, TText>
	{
		template <typename TSize1, typename TSize2>
		QGramLessOffset_(TText &text, TSize1 q, TSize2 offset): 
			QGramLess_<TSAValue, TText> (text, q, offset) {}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct QGramLessOffset_<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type	TIter;
		typedef typename Size<TString>::Type				TSize;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q, _offset;

		template <typename TSize1, typename TSize2>
		QGramLessOffset_(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset): 
			_stringSet(text),
			_q(q),
			_offset(offset) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TString const &sA = _stringSet[getValueI1(a)];
			TString const &sB = _stringSet[getValueI1(b)];
			TIter itA = begin(sA, Standard()) + getValueI2(a) + _offset;
			TIter itB = begin(sB, Standard()) + getValueI2(b) + _offset;
			TSize restA = length(sA) - getValueI2(a);
			TSize restB = length(sB) - getValueI2(b);
			if (restA > restB) {
				TIter itEnd;
				if (restB >= _q)
					itEnd = itB + _q;
				else
					itEnd = end(sB, Standard());
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd;
				if (restA >= _q)
					itEnd = itA + _q;
				else
					itEnd = end(sA, Standard());
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessNoCheck
	//
	// compare two q-grams of a given text (no check for q-grams smaller than q)
    template < typename TSAValue, typename TText >
	struct QGramLessNoCheck_ : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin;
		typename Size<TText>::Type _q;

		template <typename TSize>
		QGramLessNoCheck_(TText &text, TSize q): 
			_begin(begin(text, Standard())),
			_q(q) {}

		// skip the first <offset> characters
		template <typename TSize>
		QGramLessNoCheck_(TText &text, TSize q, TSize offset): 
			_begin(begin(text, Standard()) + offset),
			_q(q) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			if (a == b) return false;
			TIter itA = _begin + a;
			TIter itB = _begin + b;
			if (a <= b) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

    template < typename TSAValue, typename TString, typename TSpec >
	struct QGramLessNoCheck_<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q;

		template <typename TSize>
		QGramLessNoCheck_(StringSet<TString, TSpec> const &text, TSize q): 
			_stringSet(text),
			_q(q) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a);
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b);
			if (suffixLength(a, _stringSet) > suffixLength(b, _stringSet)) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessNoCheckOffset_
	//
	// compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
	struct QGramLessNoCheckOffset_: QGramLessNoCheck_<TSAValue, TText> 
	{
		template <typename TSize1, typename TSize2>
		QGramLessNoCheckOffset_(TText &text, TSize1 q, TSize2 offset): 
			QGramLessNoCheck_<TSAValue, TText> (text, q, offset) {}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct QGramLessNoCheckOffset_<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q, _offset;

		template <typename TSize1, typename TSize2>
		QGramLessNoCheckOffset_(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset): 
			_stringSet(text),
			_q(q),
			_offset(offset) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a) + _offset;
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b) + _offset;
			if (a <= b) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - little helpers
	//
	
	// map hashes 1:1 to directory
	template < typename THashValue >
	inline THashValue
	requestBucket(Nothing &, THashValue hash)
	{
		return hash;
	}

	template < typename THashValue >
	inline THashValue
	getBucket(Nothing const &, THashValue hash)
	{
		return hash;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 1: Clear directory
	template < typename TDir >
	inline void
	_qgramClearDir(TDir &dir, Nothing &)
	{
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);
	}

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 2: Count q-grams
	template < typename TDir, typename TBucketMap, typename TText, typename TShape, typename TStepSize >
	inline void
	_qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, TText const &text, TShape shape, TStepSize stepSize)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type						TSize;

		if (length(text) < length(shape)) return;
		TSize num_qgrams = (length(text) - length(shape)) / stepSize + 1;

		TIterator itText = begin(text, Standard());
		++dir[requestBucket(bucketMap, hash(shape, itText))];
		if (stepSize == 1)
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				++dir[requestBucket(bucketMap, hashNext(shape, itText))];
			}
		else
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				itText += stepSize;
				++dir[requestBucket(bucketMap, hash(shape, itText))];
			}
	}

	template < typename TDir, typename TBucketMap, typename TString, typename TSpec, typename TShape, typename TStepSize >
	inline void
	_qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, StringSet<TString, TSpec> const &stringSet, TShape shape, TStepSize stepSize)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type							TSize;

		if (stepSize == 1)
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = length(sequence) - length(shape) + 1;

				TIterator itText = begin(sequence, Standard());
				++dir[requestBucket(bucketMap, hash(shape, itText))];
				for(TSize i = 1; i < num_qgrams; ++i)
				{
					++itText;
					++dir[requestBucket(bucketMap, hashNext(shape, itText))];
				}
			}
		else
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

				TIterator itText = begin(sequence, Standard());
				for(TSize i = 0; i < num_qgrams; ++i)
				{
					++dir[requestBucket(bucketMap, hash(shape, itText))];
					itText += stepSize;
				}
			}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 3: Cumulative sum

	// First two entries are 0.
	// Step 4 increments the entries hash(qgram)+1 on-the-fly while filling the SA table.
	// After step 4 each entry (0..n-1) is the beginning of a qgram bucket.
	template < typename TDir, typename TWithConstraints >
	inline typename Value<TDir>::Type
	_qgramCummulativeSum(TDir &dir, TWithConstraints)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TDir>::Type						TSize;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		TSize diff = 0, diff_prev = 0, sum = 0;
		while (it != itEnd) 
		{
			if (TWithConstraints::VALUE && diff == (TSize)-1)
			{
				sum += diff_prev;
				diff_prev = 0;
				diff = *it;
				*it = (TSize)-1;								// disable bucket
			} else {
				sum += diff_prev;
				diff_prev = diff;
				diff = *it;
				*it = sum;
			}
			++it;
		}
		return sum + diff_prev;
	}

	// The first entry is 0.
	// This function is used when Step 4 is ommited.
	template < typename TDir, typename TWithConstraints >
	inline typename Value<TDir>::Type
	_qgramCummulativeSumAlt(TDir &dir, TWithConstraints const)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TDir>::Type						TSize;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		TSize diff = 0, sum = 0;
		while (it != itEnd) 
		{
			sum += diff;
			diff = *it;
			if (TWithConstraints::VALUE && diff == (TSize)-1) 
				diff = 0;										// ignore disabled buckets
			*it = sum;
			++it;
		}
		return sum + diff;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 4: Fill suffix array
	// w/o constraints
	template <
		typename TSA,
		typename TText, 
		typename TShape, 
		typename TDir, 
		typename TBucketMap, 
		typename TWithConstraints, 
		typename TStepSize >
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		TText const &text, 
		TShape shape, 
		TDir &dir, 
		TBucketMap &bucketMap, 
		TStepSize stepSize,
		TWithConstraints const)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type						TSize;

		TSize num_qgrams = length(text) - length(shape) + 1;
		TIterator itText = begin(text, Standard());

		if (TWithConstraints::VALUE) {
			TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;			// first hash
			if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = 0;						// if bucket is enabled
		} else
			sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = 0;			// first hash

		if (stepSize == 1)
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;	// next hash
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;					// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = i;	// next hash
			}
		else
			for(TSize i = 1; i < num_qgrams; i += stepSize)
			{
				itText += stepSize;
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;		// next hash (we mustn't use hashNext here)
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;					// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = i;		// next hash
			}
	}

	// multiple sequences
	template <
		typename TSA, 
		typename TString, 
		typename TSpec, 
		typename TShape, 
		typename TDir,
		typename TBucketMap,
		typename TStepSize,
		typename TWithConstraints >
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		StringSet<TString, TSpec> const &stringSet,
		TShape shape, 
		TDir &dir,
		TBucketMap &bucketMap,
		TStepSize stepSize,
		TWithConstraints const)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type							TSize;

		if (stepSize == 1)
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = length(sequence) - length(shape) + 1;

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

				TIterator itText = begin(sequence, Standard());
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;					// first hash
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;						// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;			// first hash

				for(TSize i = 1; i < num_qgrams; ++i)
				{
					++itText;
					assignValueI2(localPos, i);
					if (TWithConstraints::VALUE) {
						TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;			// next hash
						if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;					// if bucket is enabled
					} else
						sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = localPos;	// next hash
				}
			}
		else
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = length(sequence) - length(shape) + 1;

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

				TIterator itText = begin(sequence, Standard());
				for(TSize i = 1; i < num_qgrams; i += stepSize)
				{
					assignValueI2(localPos, i);
					if (TWithConstraints::VALUE) {
						TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;				// hash
						if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;					// if bucket is enabled
					} else
						sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;		// hash
					itText += stepSize;
				}
			}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Step 5: Correct disabled buckets
	template < typename TDir >
	inline void
	_qgramPostprocessBuckets(TDir &dir)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TDir>::Type						TSize;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		TSize prev = 0;
		for (; it != itEnd; ++it) 
			if (*it == (TSize)-1)
				*it = prev;
			else
				prev = *it;
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndex:
..summary:Builds a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndex(index)
..signature:createQGramIndex(sa, dir, bucketMap, text, shape, stepSize) [DEPRECATED]
..param.index:The q-gram index.
...type:Spec.IndexQGram
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..param.bucketMap:Stores the q-gram hashes for the openaddressing hash maps, see @Function.indexBucketMap@.
If bucketMap is of the type @Tag.Nothing@ the q-gram hash determines the bucket address in the index.
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.stepSize:Store every $stepSize$'th q-gram in the index.
..returns:Index contains the sorted list of qgrams. For each q-gram $dir$ contains the first position in index that corresponds to this q-gram.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The resulting tables must have appropriate size before calling this function.
..include:seqan/index.h
*/

	template < typename TIndex >
	inline bool _qgramDisableBuckets(TIndex &)
	{
		return false;	// we disable no buckets by default
	}

	template < typename TIndex >
	void createQGramIndex(TIndex &index)
	{
	SEQAN_CHECKPOINT
		typename Fibre<TIndex, QGramText>::Type const &text      = indexText(index);
		typename Fibre<TIndex, QGramSA>::Type         &sa        = indexSA(index);
		typename Fibre<TIndex, QGramDir>::Type        &dir       = indexDir(index);
		typename Fibre<TIndex, QGramShape>::Type      &shape     = indexShape(index);
		typename Fibre<TIndex, QGramBucketMap>::Type  &bucketMap = index.bucketMap;
		
		// 1. clear counters
		_qgramClearDir(dir, bucketMap);

		// 2. count q-grams
		_qgramCountQGrams(dir, bucketMap, text, shape, getStepSize(index));

		if (_qgramDisableBuckets(index))
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, True());
			
			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), True());

			// 5. correct disabled buckets
			_qgramPostprocessBuckets(dir);
		}
		else
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, False());
			
			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), False());
		} 
	}

	// DEPRECATED
	// better use createQGramIndex(index) (above)
	template <
        typename TSA,
		typename TDir,
		typename TBucketMap,
		typename TText,
		typename TShape,
		typename TStepSize >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		TBucketMap &bucketMap,
		TText const &text,
		TShape &shape,
		TStepSize stepSize)	
	{
	SEQAN_CHECKPOINT
		
		// 1. clear counters
		_qgramClearDir(dir, bucketMap);

		// 2. count q-grams
		_qgramCountQGrams(dir, bucketMap, text, shape, stepSize);

		// 3. cumulative sum
		SEQAN_DO(_qgramCummulativeSum(dir, False()) == length(sa));
		
		// 4. fill suffix array
		_qgramFillSuffixArray(sa, text, shape, dir, bucketMap, stepSize, False());
	}

	// DEPRECATED
	// better use createQGramIndex(index) (above)
	template <
        typename TSA,
		typename TDir,
		typename TBucketMap,
		typename TText,
		typename TShape,
		typename TStepSize >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		TBucketMap &bucketMap,
		TText const &text,
		TShape &shape)	
	{
		createQGramIndex(sa, dir, bucketMap, text, shape, 1);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndexSAOnly:
..summary:Builds the suffix array of a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndexSAOnly(sa, text, shape, stepSize)
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.text:The sequence.
..param.shape:The shape to be used. q is the length of this shape
...type:Class.Shape
..param.stepSize:Store every $stepSize$'th q-gram in the index.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The resulting tables must have appropriate size before calling this function.
..include:seqan/index.h
*/

	template < 
		typename TSA, 
		typename TText,
		typename TShape,
		typename TStepSize >
	void createQGramIndexSAOnly(
		TSA &sa,
		TText const &text,
		TShape &shape,
		TStepSize stepSize)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type TIter;

		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(sa, Standard());
		TIter itEnd = end(sa, Standard());
		TSize i = 0;
		for(; it != itEnd; ++it, i += stepSize)
			*it = i;

		// 2. Sort suffix array with quicksort
		TSize span = length(shape);
		if (i + span > length(text) + 1)
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				QGramLess_<typename Value<TSA>::Type, TText const>(text, span));
		else
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				QGramLessNoCheck_<typename Value<TSA>::Type, TText const>(text, span));
	}

	template < 
		typename TSA, 
		typename TString,
		typename TSpec,
		typename TShape,
		typename TStepSize >
	void createQGramIndexSAOnly(
		TSA &sa,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape,
		TStepSize stepSize)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TSA, Standard>::Type	TIter;
		typedef typename Value<TSA>::Type				TValue;
		typedef typename Size<TString>::Type			TSize;
		typedef StringSet<TString, TSpec>				TStringSet;
		
		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(sa, Standard());
		TValue pair;
		unsigned int const q1 = length(shape) - 1;
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			unsigned int const strlen = length(stringSet[seqNo]);
			if (strlen > q1) {
				assignValueI1(pair, seqNo);
				for (TSize i = 0; i < strlen - q1; ++it, i += stepSize) {
					assignValueI2(pair, i);
					*it = pair;
				}
			}
		}

		// 2. Sort suffix array with quicksort
		TSize q = length(shape);
		if (lengthSum(stringSet) == length(sa))
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				QGramLess_<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
		else
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				QGramLessNoCheck_<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
	}

	template < 
		typename TSA, 
		typename TDir,
		typename TText,
		typename TSize1,
		typename TSize2 >
	void _refineQGramIndex(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TSize1 oldQ,
		TSize2 newQ)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type		TIter;
		typedef typename Iterator<TDir, Standard>::Type		TDirIter;

		if (newQ <= (TSize2)oldQ) return;

		if (length(dir) < 2) {
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				QGramLess_<typename Value<TSA>::Type, TText const>(text, newQ));
			return;
		}

		// 1. Sort each bucket with quicksort and compare substrings s[i+oldQ..i+newQ)
		TDirIter dirIt = begin(dir, Standard());
		TDirIter dirItEnd = end(dir, Standard());
		TIter itBegin = begin(sa, Standard());
		TIter itBktBegin = itBegin + *dirIt;
		TIter itBktEnd;
		++dirIt;
		for(; dirIt != dirItEnd; ++dirIt, itBktBegin = itBktEnd) {
			itBktEnd = itBegin + *dirIt;
			if (itBktEnd - itBktBegin < 2) continue;
			::std::sort(
				itBktBegin, 
				itBktEnd, 
				QGramLessOffset_<typename Value<TSA>::Type, TText const>(text, newQ - oldQ, oldQ));
		}
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndexDirOnly:
..summary:Builds the directory of a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndexDirOnly(dir, bucketMap, text, shape, stepSize)
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..param.bucketMap:Stores the q-gram hashes for the openaddressing hash maps, see @Function.indexBucketMap@.
If bucketMap is of the type @Tag.Nothing@ the q-gram hash determines the bucket address in the index.
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.stepSize:Store every $stepSize$'th q-gram in the index.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The resulting tables must have appropriate size before calling this function.
..include:seqan/index.h
*/

	template <
		typename TDir,
		typename TBucketMap,
		typename TText,
		typename TShape,
		typename TStepSize >
	void createQGramIndexDirOnly(
		TDir &dir,
		TBucketMap &bucketMap,
		TText const &text,
		TShape &shape,
		TStepSize stepSize)
	{
	SEQAN_CHECKPOINT

		// 1. clear counters
		_qgramClearDir(dir, bucketMap);

		// 2. count q-grams
		_qgramCountQGrams(dir, bucketMap, text, shape, stepSize);

		// 3. cumulative sum (Step 4 is ommited)
		_qgramCummulativeSumAlt(dir, False());
	}

	template < 
		typename TDir,
		typename TBucketMap,
		typename TString,
		typename TSpec,
		typename TShape,
		typename TStepSize >
	void createQGramIndexDirOnly(
		TDir &dir,
		TBucketMap &bucketMap,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape,
		TStepSize stepSize)
	{
	SEQAN_CHECKPOINT

		// 1. clear counters
		_qgramClearDir(dir, bucketMap);

		// 2. count q-grams
		_qgramCountQGrams(dir, bucketMap, stringSet, shape, stepSize);

		// 3. cumulative sum (Step 4 is ommited)
		_qgramCummulativeSumAlt(dir, False());
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createCountArray:
..summary:Builds an index on a StringSet storing how often a q-gram occurs in each sequence.
..cat:Index
..signature:createCountsArray(counts, dir, bucketMap, stringSet, shape, stepSize)
..param.counts:The resulting list of pairs (seqNo,count).
..param.dir:The resulting array that indicates at which position in the count table the corresponding a certain q-gram can be found.
..param.bucketMap:Stores the q-gram hashes for the openaddressing hash maps, see @Function.indexBucketMap@.
If bucketMap is of the type @Tag.Nothing@ the q-gram hash determines the bucket address in the index.
..param.stringSet:The StringSet.
...type:Class.StringSet
..param.shape:The shape to be used.
...type:Class.Shape
..param.stepSize:Store every $stepSize$'th q-gram in the index.
..remarks:This function should not be called directly. Please use @Function.indexCreate@ or @Function.indexRequire@.
The resulting tables must have appropriate size before calling this function.
..include:seqan/index.h
*/

	template < 
		typename TCounts, 
		typename TDir,
		typename TBucketMap,
		typename TString,
		typename TSpec,
		typename TShape,
		typename TStepSize >
	void createCountsArray(
		TCounts &counts,
		TDir &dir,
		TBucketMap &bucketMap,
		StringSet<TString, TSpec> const &stringSet,
		TShape shape,
		TStepSize stepSize)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TString>::Type						TSize;

		TDir lastSeqSeen;
		resize(lastSeqSeen, length(dir));
		
		// 1. clear counters
		_qgramClearDir(dir, bucketMap);
		arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);

		// 2. count distinct sequences for each q-gram
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

			TIterator itText = begin(sequence, Standard());
			TValue bktNo = requestBucket(bucketMap, hash(shape, itText));
			lastSeqSeen[bktNo] = seqNo;
			++dir[bktNo];
			if (stepSize == 1)
				for(TSize i = 1; i < num_qgrams; ++i)
				{
					++itText;
					bktNo = requestBucket(bucketMap, hashNext(shape, itText));
					if (seqNo != lastSeqSeen[bktNo]) {
						lastSeqSeen[bktNo] = seqNo;
						++dir[bktNo];
					}
				}
			else
				for(TSize i = 1; i < num_qgrams; ++i)
				{
					itText += stepSize;
					bktNo = requestBucket(bucketMap, hash(shape, itText));
					if (seqNo != lastSeqSeen[bktNo]) {
						lastSeqSeen[bktNo] = seqNo;
						++dir[bktNo];
					}
				}
		}

		// 3. cumulative sum
		resize(counts, _qgramCummulativeSum(dir, False()));

		// 4. fill count array
		arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

			TIterator itText = begin(sequence, Standard());
			TValue bktNo = getBucket(bucketMap, hash(shape, itText));
			lastSeqSeen[bktNo] = seqNo;
			TSize pos = dir[bktNo + 1]++;
			counts[pos].i1 = seqNo;				// first hash
			counts[pos].i2 = 1;

			if (stepSize == 1)
				for(TSize i = 1; i < num_qgrams; ++i)
				{
					++itText;
					bktNo = getBucket(bucketMap, hashNext(shape, itText));
					if (seqNo == lastSeqSeen[bktNo])
						++(counts[dir[bktNo + 1] - 1].i2);
					else {
						lastSeqSeen[bktNo] = seqNo;
						pos = dir[bktNo + 1]++;
						counts[pos].i1 = seqNo;				// next hash
						counts[pos].i2 = 1;
					}
				}
			else
				for(TSize i = 1; i < num_qgrams; ++i)
				{
					itText += stepSize;
					bktNo = getBucket(bucketMap, hash(shape, itText));
					if (seqNo == lastSeqSeen[bktNo])
						++(counts[dir[bktNo + 1] - 1].i2);
					else {
						lastSeqSeen[bktNo] = seqNo;
						pos = dir[bktNo + 1]++;
						counts[pos].i1 = seqNo;				// next hash
						counts[pos].i2 = 1;
					}
				}
		}
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *one* sequence in external memory 

	// *** COMPARATORS & MAPS ***
        
    template <typename InType, typename Result = int>
    struct _qgramComp : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(InType const &a, InType const &b) const
        {
			typedef typename Value<InType, 2>::Type TQGram;
			typename Value<TQGram>::Type const *sa = a.i2.i;
            typename Value<TQGram>::Type const *sb = b.i2.i;
			typename Value<TQGram>::Type const *saEnd = sa + length(a.i2);

            for(; sa != saEnd; ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
			return posCompare(a.i1, b.i1);
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, unsigned _size, typename Result>
    struct _qgramComp< Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >, Result > :
        public ::std::binary_function<
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Result> {       
        inline Result operator()(
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &a,
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
			return posCompare(a.i1, b.i1);
        }
    };


    template <typename TValue, typename TResult = unsigned>
    struct _qgramHash : public ::std::unary_function<TValue, TResult> {
        inline TResult operator()(TValue const &a) const
        {
			typedef typename Value<TValue, 2>::Type	TQGram;
			TResult hash = 0;
			unsigned len = length(a.i2);
            for (unsigned i = 0; i < len; ++i) {
				hash *= ValueSize< typename Value<TQGram>::Type >::VALUE;
				hash += ordValue(a.i2[i]);
            }
            return hash;
        }
    };

	// TODO: replace fixed tuple size of 6 with q and add q to Shape template arguments
	template < 
		typename TSA, 
		typename TDir,
		typename TText,
		typename TValue,
		unsigned q >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		TText &text,
		Shape<TValue, UngappedShape<q> >)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef Shape<TValue, UngappedShape<q> >					TShape;
        typedef typename MakeUnsigned_<TValue>::Type				TUValue;

        // *** SPECIALIZATION ***

		typedef Pipe< TText, Source<> >								TSource;
        typedef Pipe< TSource, Caster<TUValue, CasterConvert> >		TUnsigner;
		typedef Pipe< TUnsigner, Tupler<q> >						TTupler;
						                typedef _qgramComp<TypeOf_(TTupler)> qcomp_t;
        typedef Pool< 
					TypeOf_(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, TSizeOf_(TTupler) > > 
				> TSortTuples;
										typedef _qgramHash<TypeOf_(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(text);
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead = length(sorter);
		bool first = true;

		for (leftToRead = length(sorter); leftToRead > 0; --leftToRead, ++sorter, ++itSA)
		{
			// copy occurence position
			*itSA = (*sorter).i1;
			if (first || qcomp(old_qgram, *sorter) != 0) 
			{
				old_qgram = *sorter;
				hash = qhash(old_qgram);

				SEQAN_ASSERT(old_hash < hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *multiple* sequences in external memory 

	template < 
		typename TSA, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TValue,
		unsigned q,
		typename TLimitsString >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		Shape<TValue, UngappedShape<q> >,
		TLimitsString &limits)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type	TConcat;
		typedef Shape<TValue, UngappedShape<q> >						TShape;
        typedef typename MakeUnsigned_<TValue>::Type					TUValue;
		typedef Multi<
			Tupler<q, true, Compressed>, 
			typename Value<TSA>::Type,
			typename StringSetLimits< StringSet<TString, TSpec> >::Type >		TTuplerSpec;

        // *** SPECIALIZATION ***

		typedef Pipe< TConcat, Source<> >							TSource;
        typedef Pipe< TSource, Caster<TUValue, CasterConvert> >		TUnsigner;
		typedef Pipe< TUnsigner, TTuplerSpec >						TTupler;
						                typedef _qgramComp<TypeOf_(TTupler)> qcomp_t;
        typedef Pool< 
					TypeOf_(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, TSizeOf_(TTupler) > > 
				> TSortTuples;
										typedef _qgramHash<TypeOf_(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(concat(stringSet));
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner, limits);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram = *sorter;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead;
		bool first = true;

		for (leftToRead = length(sorter); leftToRead > 0; --leftToRead, ++sorter, ++itSA)
		{
			// copy occurence position
			*itSA = (*sorter).i1;
			
			if (first || qcomp(old_qgram, *sorter) != 0) 
			{
				old_qgram = *sorter;
				hash = qhash(old_qgram);
				
				SEQAN_ASSERT_LEQ(old_hash, hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}



//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline typename Size<Index<TText, IndexQGram<TShapeSpec, TSpec> > >::Type
	_qgramQGramCount(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index)
	{
		typedef Index<TText, IndexQGram<TShapeSpec, TSpec> >	TIndex;
		typedef typename Fibre<TIndex, QGramShape>::Type		TShape;
		typedef typename Size<TIndex>::Type						TSize;

		TShape const &shape = indexShape(index);

		// count all overlapping q-grams
		TSize stepSize = getStepSize(index);
		TSize qgramCount = 0;
		for(unsigned i = 0; i < countSequences(index); ++i)
			if (sequenceLength(i, index) >= length(shape))
				qgramCount += (sequenceLength(i, index) - length(shape)) / stepSize + 1;
		return qgramCount;
	}

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline bool indexCreate(
		Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, 
		FibreSADir, 
		Default const) 
	{		
		resize(indexSA(index), _qgramQGramCount(index), Exact());
		resize(indexDir(index), _fullDirLength(index), Exact());
		createQGramIndex(index);
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexSupplied(Index<TText, TSpec> &index, FibreSADir) {
		return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreDir())));
	}

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline bool indexCreate(
		Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, 
		FibreSA, 
		Default const)
	{
		resize(indexSA(index), _qgramQGramCount(index), Exact());
		createQGramIndexSAOnly(indexSA(index), indexText(index), indexShape(index), getStepSize(index));
		return true;
	}

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline bool indexCreate(
		Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, 
		FibreCounts, 
		Default const) 
	{
		resize(indexCountsDir(index), _fullDirLength(index), Exact());
		createCountsArray(indexCounts(index), indexCountsDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));
		return true;
	}

	template <typename TText, typename TShapeSpec, typename TSpec>
	inline bool indexCreate(
		Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, 
		FibreDir, 
		Default const)
	{
		typedef Index<TText, IndexQGram<TShapeSpec, TSpec> > TIndex;

		resize(indexDir(index), _fullDirLength(index), Exact());
		createQGramIndexDirOnly(indexDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));
		return true;
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.getKmerSimilarityMatrix:
..summary:Creates a matrix storing the number of common q-grams between all pairs of sequences.
..cat:Index
..signature:getKmerSimilarityMatrix(index, distMat[, seqSet])
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.distMat:The resulting q-gram similarity matrix.
...type:Concept.Container
..param.seqSet:Contains sequence numbers if only a subset of sequences should be compared.
...type:Concept.Container
..remarks:$distMat$ will be resized to $seqCount*seqCount$, where $seqCount$ is the number of sequences in the index/in $seqSet$.
The number of common q-grams between sequence $i$ and $j$ is stored at position $i*seqCount + j$.
It sums up the minimum number of q-gram occurrences between both sequences for each q-gram.
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TDistMatrix >
	inline void getKmerSimilarityMatrix(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		TDistMatrix &distMat)
	{
		typedef Index< TObject, IndexQGram<TShapeSpec,TSpec> >	TIndex;
		typedef typename Size<TIndex>::Type						TSize;
		typedef typename Size<TDistMatrix>::Type				TSizeMat;
		typedef typename Value<TDistMatrix>::Type				TValueMat;

		typedef typename Fibre<TIndex, QGramCountsDir>::Type	TCountsDir;
		typedef typename Iterator<TCountsDir, Standard>::Type	TIterCountsDir;
		typedef typename Fibre<TIndex, QGramCounts>::Type		TCounts;
		typedef typename Iterator<TCounts, Standard>::Type		TIterCounts;

		// declare requirements
		indexRequire(index, QGramCounts());

		// initialize distance matrix
		TSizeMat seqNoLength = countSequences(index);
		clear(distMat);
		resize(distMat, seqNoLength * seqNoLength);
		arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

		TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
		TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
		TIterCounts itCountsBegin = begin(indexCounts(index), Standard());

		// for each bucket count common q-grams for each sequence pair
		TSize bucketBegin = *itCountsDir;
		for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir) 
		{
			TSize bucketEnd = *itCountsDir;

			// q-gram must occur in at least 2 different sequences
			if (bucketBegin != bucketEnd) 
			{
				TIterCounts itA = itCountsBegin + bucketBegin;
				TIterCounts itEnd = itCountsBegin + bucketEnd;
				for(; itA != itEnd; ++itA) 
				{
					TSizeMat ofs = (*itA).i1 * seqNoLength;
					TSize countA = (*itA).i2;
					TIterCounts itB = itA;

					for(; itB != itEnd; ++itB) 
					{
						TSize countB = (*itB).i2;
						if (countA < countB)
							distMat[ofs + (*itB).i1] += countA;
						else
							distMat[ofs + (*itB).i1] += countB;
					}
				}
			}
			bucketBegin = bucketEnd;
		}

		// copy upper triangle to lower triangle and scale
		for(TSizeMat row = 0; row < seqNoLength; ++row) 
		{
			TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
			for(TSizeMat col = row + 1; col < seqNoLength; ++col)
			{
				// fractional common kmer count
				TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
				TValueMat val = distMat[row * seqNoLength + col];

				// number of common q-grams / Number of possible common q-grams
				if (maxValRow < maxValCol) {
					if (maxValRow != 0)
						val /= maxValRow;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				} else {
					if (maxValCol != 0)
						val /= maxValCol;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				}
			}
		}

		// set diagonal to 1
		for(TSizeMat i = 0; i < seqNoLength; ++i)
			distMat[i * (seqNoLength + 1)] = 1;
	}


//////////////////////////////////////////////////////////////////////////////
// getKmerSimilarityMatrix 
//     for a subset of the StringSet given as a sorted string of sequence numbers

	template < 
		typename TObject, 
		typename TShapeSpec, 
		typename TSpec, 
		typename TDistMatrix, 
		typename TSeqNoString 
	>
	inline void getKmerSimilarityMatrix(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		TDistMatrix &distMat,
		TSeqNoString const &seqNo)
	{
		typedef Index< TObject, IndexQGram<TShapeSpec,TSpec> >	TIndex;
		typedef typename Size<TIndex>::Type						TSize;
		typedef typename Size<TDistMatrix>::Type				TSizeMat;
		typedef typename Value<TDistMatrix>::Type				TValueMat;

		typedef typename Fibre<TIndex, QGramCountsDir>::Type	TCountsDir;
		typedef typename Iterator<TCountsDir, Standard>::Type	TIterCountsDir;
		typedef typename Fibre<TIndex, QGramCounts>::Type		TCounts;
		typedef typename Iterator<TCounts, Standard>::Type		TIterCounts;
		typedef typename Iterator<TSeqNoString, Standard>::Type	TIterSeqNo;

		// declare requirements
		indexRequire(index, QGramCounts());

		// initialize distance matrix
		TSizeMat seqNoLength = length(seqNo);
		clear(distMat);
		resize(distMat, seqNoLength * seqNoLength);
		arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

		TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
		TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
		TIterCounts itCountsBegin = begin(indexCounts(index), Standard());
		TIterSeqNo itSetBegin = begin(seqNo, Standard());
		TIterSeqNo itSetEnd = end(seqNo, Standard());

		// for each bucket count common q-grams for each sequence pair
		TSize bucketBegin = *itCountsDir;
		for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir) 
		{
			TSize bucketEnd = *itCountsDir;

			// q-gram must occur in at least 2 different sequences
			if (bucketBegin != bucketEnd) 
			{
				TIterCounts itA = itCountsBegin + bucketBegin;
				TIterCounts itEnd = itCountsBegin + bucketEnd;
				TIterSeqNo itSetA = itSetBegin;

				while (itA != itEnd && itSetA != itSetEnd)
				{
					if ((*itA).i1 < *itSetA)
						++itA;
					else if ((*itA).i1 > *itSetA)
						++itSetA;
					else 
					{
						TSizeMat ofs = (itSetA - itSetBegin) * seqNoLength;
						TSize countA = (*itA).i2;
						TIterCounts itB = itA;
						TIterSeqNo itSetB = itSetA;

						while (itB != itEnd && itSetB != itSetEnd)
						{
							if ((*itB).i1 < *itSetB)
								++itB;
							else if ((*itB).i1 > *itSetB)
								++itSetB;
							else 
							{
								TSize countB = (*itB).i2;
								if (countA < countB)
									distMat[ofs + (itSetB - itSetBegin)] += countA;
								else
									distMat[ofs + (itSetB - itSetBegin)] += countB;
								++itB;
								++itSetB;
							}
						}
						++itA;
						++itSetA;
					}
				}
			}
			bucketBegin = bucketEnd;
		}

		// copy upper triangle to lower triangle and scale
		for(TSizeMat row = 0; row < seqNoLength; ++row) 
		{
			TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
			for(TSizeMat col = row + 1; col < seqNoLength; ++col)
			{
				// fractional common kmer count
				TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
				TValueMat val = distMat[row * seqNoLength + col];

				// number of common q-grams / Number of possible common q-grams
				if (maxValRow < maxValCol) {
					if (maxValRow != 0)
						val /= maxValRow;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				} else {
					if (maxValCol != 0)
						val /= maxValCol;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				}
			}
		}

		// set diagonal to 1
		for(TSizeMat i = 0; i < seqNoLength; ++i)
			distMat[i * (seqNoLength + 1)] = 1;
	}
	
	
//////////////////////////////////////////////////////////////////////////////
/**
.Function.range:
..signature:range(index, shape)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.shape:A shape object.
...note:The shape stores the q-gram of the last call to @Function.hash@ or @Function.hashNext@.
...type:Class.Shape
..returns:All positions where the q-gram stored in $shape$ occurs in the text (see @Tag.QGram Index Fibres.QGramText@)
are stored in a contiguous range of the suffix array.
$range$ returns begin and end position of this range.
If the type of $index$ is $TIndex$ the return type is $Pair<Size<TIndex>::Type>.
..note:The necessary index tables are built on-demand via @Function.indexRequire@ if index is not $const$.
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline Pair<typename Size< Index< TObject, IndexQGram<TShapeSpec, TSpec> > >::Type>
	range(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape)
	{
		typedef typename Size< Index< TObject, IndexQGram<TShapeSpec, TSpec> > >::Type TSize;
		typedef typename Size< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		return Pair<TSize>(indexDir(index)[bucket], indexDir(index)[bucket + 1]);
	}

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline Pair<typename Size< Index< TObject, IndexQGram<TShapeSpec, TSpec> > >::Type>
	range(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape)
	{
		indexRequire(index, QGramDir());
		return getOccurrences(const_cast<Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getOccurrence:
..signature:getOccurrence(index, shape)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.shape:A shape object.
...note:The shape stores the q-gram of the last call to @Function.hash@ or @Function.hashNext@.
...type:Class.Shape
..returns:A position where the q-gram stored in $shape$ occurs in the text (see @Tag.QGram Index Fibres.QGramText@).
If the type of $index$ is $TIndex$ the return type is $SAValue<TIndex>::Type$.
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename SAValue< Index< TObject, IndexQGram<TShapeSpec, TSpec> > >::Type
	getOccurrence(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		return saAt(indexDir(index)[getBucket(indexBucketMap(index), value(shape))], index);
	}

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename SAValue< Index< TObject, IndexQGram<TShapeSpec, TSpec> > >::Type
	getOccurrence(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramSADir());
		return getOccurrence(const_cast<Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getOccurrences:
..signature:getOccurrences(index, shape)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.shape:A shape object.
...note:The shape stores the q-gram of the last call to @Function.hash@ or @Function.hashNext@.
...type:Class.Shape
..returns:All positions where the q-gram stored in $shape$ occurs in the text (see @Tag.QGram Index Fibres.QGramText@).
If the type of $index$ is $TIndex$ the return type is $Infix<Fibre<TIndex, QGramSA>::Type const>::Type$.
..remarks:The necessary index tables are built on-demand via @Function.indexRequire@ if index is not $const$.
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	getOccurrences(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		typedef typename Size<typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		return infix(indexSA(index), indexDir(index)[bucket], indexDir(index)[bucket + 1]);
	}	

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	getOccurrences(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramSADir());
		return getOccurrences(const_cast<Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}	

//////////////////////////////////////////////////////////////////////////////
/**
.Function.countOccurrences:
..signature:countOccurrences(index, shape)
..param.index:A q-gram index.
...type:Spec.IndexQGram
..param.shape:A shape object.
...note:The shape stores the q-gram of the last call to @Function.hash@ or @Function.hashNext@.
...type:Class.Shape
..returns:The number of positions where the q-gram stored in $shape$ occurs in the text (see @Tag.QGram Index Fibres.QGramText@).
If the type of $index$ is $TIndex$ the return type is $Size<TIndex>::Type$.
..note:The necessary index tables are built on-demand via @Function.indexRequire@ if index is not $const$.
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Size< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	countOccurrences(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		typedef typename Size<typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		return indexDir(index)[bucket + 1] - indexDir(index)[bucket];
	}	

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Size< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type 
	countOccurrences(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramDir());
		return countOccurrences(const_cast<Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.countOccurrencesMultiple:
..summary:Returns the number of occurences of a q-gram for every sequence of a @Class.StringSet@ .
..signature:countOccurrencesMultiple(index, shape)
..param.index:A q-gram index of a @Class.StringSet@.
...type:Spec.IndexQGram
..param.shape:A shape object.
...note:The shape stores the q-gram of the last call to @Function.hash@ or @Function.hashNext@.
...type:Class.Shape
..returns:A sequence of @Class.Pair.pairs@ (seqNo,count), count>0.
For every @Class.StringSet@ sequence the q-gram occurs in, seqNo is the sequence number and count the number of occurrences.
If the type of $index$ is $TIndex$ the return type is $Infix<Fibre<TIndex, QGramCounts>::Type const>::Type$.
..remarks:The necessary index tables are built on-demand via @Function.indexRequire@ if index is not $const$.
..see:Function.countOccurrences
..include:seqan/index.h
*/

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreCounts>::Type const >::Type 
	countOccurrencesMultiple(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		typedef typename Size<typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreCountsDir>::Type>::Type TDirSize;
		TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
		return infix(indexCounts(index), indexCountsDir(index)[bucket], indexCountsDir(index)[bucket + 1]);
	}	

	template < typename TObject, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
	inline typename Infix< typename Fibre< Index< TObject, IndexQGram<TShapeSpec, TSpec> >, FibreCounts>::Type const >::Type 
	countOccurrencesMultiple(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index,
		Shape< TValue, TShapeSpec2 > const &shape) 
	{
		indexRequire(index, QGramCounts());
		return countOccurrencesMultiple(const_cast<Index< TObject, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
	}


//////////////////////////////////////////////////////////////////////////////
// clear

	template < typename TText, typename TShapeSpec, typename TSpec >
	inline void clear(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index)
	{
		clear(getFibre(index, QGramSA()));
		clear(getFibre(index, QGramDir()));
		clear(getFibre(index, QGramCounts()));
		clear(getFibre(index, QGramCountsDir()));
	}

//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TShapeSpec, typename TSpec >
	inline bool open(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");
		bool result = true;
		if ((!open(getFibre(index, QGramText()), toCString(name), openMode)) && 
			(!open(getFibre(index, QGramText()), fileName, openMode)))
			result = false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, QGramSA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	open(getFibre(index, QGramDir()), toCString(name), openMode);
		return result;
	}
	template < typename TObject, typename TShapeSpec, typename TSpec >
	inline bool open(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		const char *fileName)
	{
		return open(index, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TShapeSpec, typename TSpec >
	inline bool save(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");
		if ((!save(getFibre(index, QGramText()), toCString(name), openMode)) && 
			(!save(getFibre(index, QGramText()), fileName, openMode))) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, QGramSA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	save(getFibre(index, QGramDir()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TShapeSpec, typename TSpec >
	inline bool save(
		Index< TObject, IndexQGram<TShapeSpec, TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif //#ifndef SEQAN_HEADER_...
