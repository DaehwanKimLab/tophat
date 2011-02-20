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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_FORMAT_MMAP_H
#define SEQAN_HEADER_FILE_FORMAT_MMAP_H

namespace SEQAN_NAMESPACE_MAIN
{

	// define memory mapped stringset
	typedef StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >	MultiFasta;	//deprecated
	typedef StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >	MultiSeqFile;


	template <typename TValue>
	inline bool
	_isLineBreak(TValue value)
	{
		return (value == '\n' || value == '\r');
	}

	template <typename TIterator>
	inline bool
	_seekLineBreak(TIterator &it, TIterator itEnd)
	{
		while (!_isLineBreak(*it))
			if (++it == itEnd) return false;
		return true;
	}

	template <typename TIterator>
	inline bool
	_seekNonLineBreak(TIterator &it, TIterator itEnd)
	{
		if (*it == '\n')
		{
			if (++it == itEnd) return false;
			if (*it == '\r')
				if (++it == itEnd) return false;
		} else
			if (*it == '\r')
			{
				if (++it == itEnd) return false;
				if (*it == '\n')
					if (++it == itEnd) return false;
			}
		return true;
	}

	template <typename TIterator>
	inline bool
	_seekTab(TIterator& it, TIterator itEnd)
	{
		for (; it != itEnd; ++it)
			if (*it == '\t') return true;
		return false;
	}


//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta
//////////////////////////////////////////////////////////////////////////////

	// test for Fasta format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		Fasta)
	{
		return seq[0] == '>';
	}
	
	template < typename TFilename >
	inline bool
	guessFormatFromFilename(
		TFilename const & fname,
		Fasta)
	{
		typedef typename Value<TFilename>::Type									TValue;
		typedef ModifiedString<TFilename, ModView<FunctorLowcase<TValue> > >	TLowcase;
		
		typename Size<TFilename>::Type len = length(fname);
		TLowcase lowcaseFname(fname);

		if (len >= 3)
		{
			if (suffix(lowcaseFname, len - 3) == ".fa") return true;
		}
		if (len >= 4)
		{
			if (suffix(lowcaseFname, len - 4) == ".fas") return true;
			if (suffix(lowcaseFname, len - 4) == ".faa") return true;	// FASTA Amino Acid file 
			if (suffix(lowcaseFname, len - 4) == ".ffn") return true;	// FASTA nucleotide coding regions file
			if (suffix(lowcaseFname, len - 4) == ".fna") return true;	// FASTA Nucleic Acid file
			if (suffix(lowcaseFname, len - 4) == ".frn") return true;
		}
		if (len >= 6)
		{
			if (suffix(lowcaseFname, len - 6) == ".fasta") return true;
		}
		return false;
	}
	
	// split stringset into single Fasta sequences
	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		Fasta)
	{
		typedef String<TValue, MMap<TConfig> >						TString;
		typedef StringSet<TString, ConcatDirect<TDelimiter> >		TStringSet;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator itBeg = begin(me.concat, Standard());
		TIterator itEnd = end(me.concat, Standard());
		bool newLine = true;
		for (TIterator it = itBeg; it != itEnd; ++it)
		{
			TValue c = *it;
			if (newLine && c == '>')
				appendValue(me.limits, it - itBeg, Generous());
			newLine = _isLineBreak(c);
		}
		if (empty(me.limits))
			appendValue(me.limits, 0);
		appendValue(me.limits, itEnd - itBeg);
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '>')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it)
			if (!_isLineBreak(*it))
			{
				*dit = *it;
				++dit;
			}
		resize(dst, dit - begin(dst, Standard()));
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '>')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQual(
		TSeq & dst,
		TFastaSeq const &,
		Fasta)
	{
		clear(dst);
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TFastaSeq const &,
		Fasta)
	{
		clear(dst);
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fastq (Fasta extension for quality values)
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Fastq:
	FASTQ file format for sequences.
*/
struct TagFastq_;
typedef Tag<TagFastq_> const Fastq;

	// test for Fastq format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		Fastq)
	{
		return seq[0] == '@';
	}
	
	template < typename TFilename >
	inline bool
	guessFormatFromFilename(
		TFilename const & fname,
		Fastq)
	{
		typedef typename Value<TFilename>::Type									TValue;
		typedef ModifiedString<TFilename, ModView<FunctorLowcase<TValue> > >	TLowcase;
		
		typename Size<TFilename>::Type len = length(fname);
		TLowcase lowcaseFname(fname);

		if (len >= 3)
		{
			if (suffix(lowcaseFname, len - 3) == ".fq") return true;
		}
		if (len >= 6)
		{
			if (suffix(lowcaseFname, len - 6) == ".fastq") return true;
		}
		return false;
	}
	
	// split stringset into single Fasta sequences
	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		Fastq)
	{
		typedef String<TValue, MMap<TConfig> >						TString;
		typedef StringSet<TString, ConcatDirect<TDelimiter> >		TStringSet;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator itBeg = begin(me.concat, Standard());
		TIterator itEnd = end(me.concat, Standard());
		bool newLine = true;
		for (TIterator it = itBeg; it != itEnd; ++it)
		{
			if (newLine && *it == '@')
				appendValue(me.limits, it - itBeg, Generous());
			if (newLine && *it == '+')
			{
				// skip qualitity fasta id
				if (!_seekLineBreak(it, itEnd)) break;
				if (!_seekNonLineBreak(it, itEnd)) break;
				// skip qualitity values
				if (!_seekLineBreak(it, itEnd)) break;
			}
			newLine = _isLineBreak(*it);
		}
		if (empty(me.limits))
			appendValue(me.limits, 0);
		appendValue(me.limits, itEnd - itBeg);
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '@')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it) 
		{
			if (_isLineBreak(*it))
			{
				if (!_seekNonLineBreak(it, itEnd)) break;
				if (*it == '+') break;
			}
			*dit = *it;
			++dit;
		}
		resize(dst, dit - begin(dst, Standard()));
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '@')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQual(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		if (it == itEnd) return;
		if (*it == '@')
		{
			// seek quality id
			do {
				if (!_seekLineBreak(it, itEnd)) return;
				if (!_seekNonLineBreak(it, itEnd)) return;
			} while (*it != '+');

			// skip quality id
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;

			// copy sequence
			resize(dst, itEnd - it);		
			TDstIterator dit = begin(dst, Standard());
			for (; it != itEnd; ++it)
				if (!_isLineBreak(*it))
				{
					*dit = *it;
					++dit;
				}
			resize(dst, dit - begin(dst, Standard()));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it1 = itBeg;
		
		clear(dst);
		if (it1 == itEnd) return;
		if (*it1 == '@')
		{
			do {
				if (!_seekLineBreak(it1, itEnd)) return;
				if (!_seekNonLineBreak(it1, itEnd)) return;
			} while (*it1 != '+');
			TIterator it2 = it1;
			_seekLineBreak(it2, itEnd);
			assign(dst, infix(fasta, (it1 - itBeg) + 1, it2 - itBeg));
		}
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - QSeq (used by Illumina for most of their read files)
//////////////////////////////////////////////////////////////////////////////

	struct _QSeq;
	typedef Tag<_QSeq> const QSeq;

	// FIXME The following enum is more or less arbitrary since the information
	// in a QSeq file may differ depending on where they come from. Not sure if
	// this is something that needs to be fixed here, rather than in the Illu-
	// mina pipeline itself.
	//
	// Also, this enum is quite convoluted but I don't feel like spilling a lot
	// of common symbols (e.g. 'X', 'Y') into the SeqAn namespace.
	struct QSeqEntry {
		enum {
			MachineName,
			Run,
			Lane,
			Tile,
			X,
			Y,
			Index,
			Read,
			Sequence,
			Quality,
			Filter
		};
	};

	template < typename TString >
	inline bool _isQSeqFile(TString const& filename)
	{
        unsigned int const namelen = 19;
        unsigned int const pathlen = length(filename);
        if (pathlen < namelen) return false;
		::std::string str;
		assign(str, suffix(filename, pathlen - namelen));
		::std::istringstream is(str);
		unsigned int num;
		// Format (as regex): /^s_\d_\d_\d{4}_qseq.txt$/
		return is.get() == 's' && is.get() == '_' &&
			isdigit(is.get()) && is.get() == '_' &&
			isdigit(is.get()) && is.get() == '_' &&
			is >> num &&
			getline(is, str) && str == "_qseq.txt" && is.eof();
	}

//	// Needed?
//	template < typename TString >
//	inline TString _getFirstFile(
//		char const* dirname,
//		QSeq)
//	{
//		Directory dir(dirname);
//		
//		for ( ; dir; ++dir)
//			if (_isQSeqFile(*dir))
//				return *dir;
//	}
//
//	template < typename TString >
//	inline String< TString >
//	_getAllFiles(
//		char const* dirname,
//		QSeq)
//	{
//		String< TString > ret;
//		Directory dir(dirname);
//
//		for (; dir; ++dir)
//			if (_isQSeqFile(*dir))
//				appendValue(ret, *dir);
//
//		return ret;
//	}

	// test for QSeq format
	template < typename TSeq >
	inline bool
	guessFormat(
		TSeq const & seq,
		QSeq)
	{
		typedef typename Iterator<TSeq const>::Type TIter;
		TIter front = begin(seq);
		TIter const back = end(seq);
		if (!_seekTab(front, back)) return false;
		::std::string token_base;
		assign(token_base, seq);
		::std::istringstream is(token_base);
		::std::string mname;
		unsigned int numval;
		// Actual information encoded in qseq file may vary. Take a few guesses:
		//     machine name, run number, lane number, tile number
		return is >> mname >> numval >> numval >> numval;
	}
	
	template < typename TFilename >
	inline bool
	guessFormatFromFilename(
		TFilename const & fname,
		QSeq)
	{
		// QSeq files come in a variety of ways throughout the Gerald pipeline.
		// In this simplest case, "sorted.txt" is a file in a fragment genome
		// directory, each corresponding to 10MB worth of DNA.
		static CharString const standalone_name = "sorted.txt";
		typename Size<TFilename>::Type len = length(fname);
		return len >= length(standalone_name) &&
			suffix(fname, len - length(standalone_name)) == standalone_name ||
			_isQSeqFile(fname);
	}

	// split stringset into single QSeq sequences
	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		QSeq)
	{
		typedef String<TValue, MMap<TConfig> >						TString;
		typedef StringSet<TString, ConcatDirect<TDelimiter> >		TStringSet;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator const front = begin(me.concat, Standard());
		TIterator const back = end(me.concat, Standard());

		appendValue(me.limits, 0, Generous());
		for (TIterator i = front; i != back; ++i)
			if (_isLineBreak(*i))
				appendValue(me.limits, i - front, Generous());

		if (!_isLineBreak(*(back - 1))) // Ignore final line break.
			appendValue(me.limits, back - front);
	}

	template <typename TSequence, typename TSource>
	void assignQSeqEntry(
		TSequence& destination,
		TSource const& source,
		unsigned int entry
	) {
		typedef typename Iterator<TSource const>::Type TIterator;

		TIterator const front = begin(source, Standard());
		TIterator const back = end(source, Standard());

		TIterator infixStart = front;

		for (unsigned int i = QSeqEntry::MachineName; i < entry; ++i) {
			_seekTab(infixStart, back);
			++infixStart;
		}

		SEQAN_ASSERT(infixStart != back);

		TIterator infixEnd = infixStart + 1;
		_seekTab(infixEnd, back);

		assign(destination, infix(source, infixStart - front, infixEnd - front));
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
		assignQSeqEntry(dst, fasta, QSeqEntry::Sequence);
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
		// For now: just return the whole line.
		typename Position<TQSeqSeq const>::Type front = 0;
		while (_isLineBreak(fasta[front]))
			++front;
		assign(dst, infix(fasta, front, length(fasta)));
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignQual(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
		assignQSeqEntry(dst, fasta, QSeqEntry::Quality);
	}

	template <typename TSeq, typename TQSeqSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TQSeqSeq const & fasta,
		QSeq)
	{
		assignSeqId(dst, fasta, QSeq());
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - Auto-Format
//////////////////////////////////////////////////////////////////////////////

	typedef
		TagList<Fastq,
		TagList<Fasta,
		TagList<QSeq> > > 						SeqFormats;
	typedef TagSelector<SeqFormats>				AutoSeqFormat;

//____________________________________________________________________________
// guess file format

	template < typename TFileSeq >
	inline bool
	guessFormat(
		TFileSeq const &,
		TagSelector<> &format)
	{
		format.tagId = 0;
		return false;
	}
	
	template < typename TFileSeq, typename TTagList >
	inline bool
	guessFormat(
		TFileSeq const & seq,
		TagSelector<TTagList> &format)
	{
		if (guessFormat(seq, typename TTagList::Type()))
		{
			if (format.tagId == 0)
			{
				format.tagId = Length<TTagList>::VALUE;			// if tagId == 0 then store detected format
				return true;
			} else
				return format.tagId == Length<TTagList>::VALUE;	// if tagId != 0 then compare detected format with tagId
		}
		return guessFormat(seq, static_cast<typename TagSelector<TTagList>::Base &>(format));
	}
	
//____________________________________________________________________________
// guess file format from filename

	template < typename TFilename >
	inline bool
	guessFormatFromFilename(
		TFilename const &,
		TagSelector<> &format)
	{
		format.tagId = 0;
		return false;
	}
	
	template < typename TFilename, typename TTagList >
	inline bool
	guessFormatFromFilename(
		TFilename const & fname,
		TagSelector<TTagList> &format)
	{
		if (guessFormatFromFilename(fname, typename TTagList::Type()))
		{
			if (format.tagId == 0)
			{
				format.tagId = Length<TTagList>::VALUE;			// if tagId == 0 then store detected format
				return true;
			} else
				return format.tagId == Length<TTagList>::VALUE;	// if tagId != 0 then compare detected format with tagId
			return true;
		}
		return guessFormatFromFilename(fname, static_cast<typename TagSelector<TTagList>::Base &>(format));
	}

//____________________________________________________________________________
// split stringset into single sequences

	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &, 
		TagSelector<void> const &)
	{
	}
	
	template < typename TValue, typename TConfig, typename TDelimiter, typename TTagList >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		TagSelector<TTagList> const &format)
	{
		if (format.tagId == Length<TTagList>::VALUE)
			split(me, typename TTagList::Type());
		else
			split(me, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}
	
//____________________________________________________________________________
// assignSeq

	template <typename TSeq, typename TFileSeq>
	inline void
	assignSeq(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
	}

	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignSeq(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
		if (format.tagId == Length<TTagList>::VALUE)
			assignSeq(dst, seq, typename TTagList::Type());
		else
			assignSeq(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignSeqId

	template <typename TSeqId, typename TFileSeq>
	inline void
	assignSeqId(
		TSeqId &,
		TFileSeq const &,
		TagSelector<> const &)
	{
	}

	template <typename TSeqId, typename TFileSeq, typename TTagList>
	inline void
	assignSeqId(
		TSeqId & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
		if (format.tagId == Length<TTagList>::VALUE)
			assignSeqId(dst, seq, typename TTagList::Type());
		else
			assignSeqId(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignQual

	template <typename TSeq, typename TFileSeq>
	inline void
	assignQual(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
	}
	
	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignQual(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
		if (format.tagId == Length<TTagList>::VALUE)
			assignQual(dst, seq, typename TTagList::Type());
		else
			assignQual(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//____________________________________________________________________________
// assignQualId

	template <typename TSeq, typename TFileSeq>
	inline void
	assignQualId(
		TSeq &,
		TFileSeq const &,
		TagSelector<> const &)
	{
	}
	
	template <typename TSeq, typename TFileSeq, typename TTagList>
	inline void
	assignQualId(
		TSeq & dst,
		TFileSeq const & seq,
		TagSelector<TTagList> const &format)
	{
		if (format.tagId == Length<TTagList>::VALUE)
			assignQualId(dst, seq, typename TTagList::Type());
		else
			assignQualId(dst, seq, static_cast<typename TagSelector<TTagList>::Base const &>(format));
	}

//////////////////////////////////////////////////////////////////////////////
// Directory import
//////////////////////////////////////////////////////////////////////////////

	template <typename TSeqSet, typename TFilename, typename TSeqFormat>
	inline void
	appendSeqs(
		TSeqSet &seqSet,
		TFilename &dirname,
		TSeqFormat format)
	{
		typedef typename Value<TSeqSet>::Type TSeq;
		
		Directory		dir(dirname);
		TSeq			seq;
		MultiSeqFile	multiSeqFile;
		
		if (!atEnd(dir))
		{
			// dirname is path of a directory
			CharString fname = dirname;
#ifdef PLATFORM_WINDOWS
			appendValue(fname, '\\');
#else
			appendValue(fname, '/');
#endif
			size_t len = length(fname);

			for (; !atEnd(dir); goNext(dir))
			{
				assign(suffix(fname, len), value(dir));
				if (guessFormatFromFilename(fname, format))
				{
				std::cout << fname<<std::endl;
					if (!open(multiSeqFile.concat, toCString(fname), OPEN_RDONLY)) continue;

					split(multiSeqFile, format);
					unsigned seqCount = length(multiSeqFile);
					
					reserve(seqSet, length(seqSet) + seqCount, Generous());					
					for(unsigned i = 0; i < seqCount; ++i)
					{
						assignSeq(seq, multiSeqFile[i], format);
						appendValue(seqSet, seq, Generous());
					}
					close(multiSeqFile.concat);
				}
			}
		} 
		else
		{
			// dirname is path of a file
			if (guessFormatFromFilename(dirname, format))
			{
				if (!open(multiSeqFile.concat, toCString(dirname), OPEN_RDONLY)) return;

				split(multiSeqFile, format);
				unsigned seqCount = length(multiSeqFile);
				
				reserve(seqSet, length(seqSet) + seqCount, Generous());					
				for(unsigned i = 0; i < seqCount; ++i)
				{
					assignSeq(seq, multiSeqFile[i], format);
					appendValue(seqSet, seq, Generous());
				}
				close(multiSeqFile.concat);
			}
		}
	}
	
}

#endif
