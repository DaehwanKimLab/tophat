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
#include <sstream>

#ifndef SEQAN_HEADER_FILE_SAM_H
#define SEQAN_HEADER_FILE_SAM_H

namespace SEQAN_NAMESPACE_MAIN {
	
    template<typename TPropertyMap, typename TDescriptor, typename TValue>
    inline void
    assignProperty(TPropertyMap& pm,
                   TDescriptor const d,
                   TValue const val);
    
	//////////////////////////////////////////////////////////////////////////////
	// File Formats - SAM
	//////////////////////////////////////////////////////////////////////////////
	
	template<typename TSpec>
	struct TagSBAM_;
		
	/**
	 .Tag.File Format.tag.Sam:
        SAM file format for alignments.
	 ..remark:According to the http://samtools.sourceforge.net/|specification 0.1.2
	 */
	struct _SAM;
	typedef Tag<TagSBAM_<_SAM> > const SAM;
	
	/**
	 .Tag.File Format.tag.Bam:
        BAM file format for alignments.
	 ..remark:According to the http://samtools.sourceforge.net/|specification 0.1.2
	 */
	struct _BAM;
	typedef Tag<TagSBAM_<_BAM> > const BAM;
	
	/**
	 .Class.SAMFileMeta
	 ..cat:Input/Output
	 ..summary:Data structure to store information of header in the SAM of BAM file format.
	 
	 */
	class SAMFileMeta {
		
	};
	
	/**
	 .Class.SAMMeta
	 ..cat:Input/Output
	 ..summary:Data structure to store information belonging to one sequence. Needed for output in SAM and BAM file format.
	 */
	class SAMMeta {
			
	public:
		CharString qName;	
        unsigned int flag;
        CharString rName;
        // Position can be extracted from the alignment
        unsigned int mapQ;
        // CIGAR is calculated from the Alignment
        CharString MRNM;
        unsigned int mPos;
        unsigned int iSize;
        CharString qual;
        String<Triple<CharString, char, CharString> > additional;
        
	public:
		SAMMeta(){
			SEQAN_CHECKPOINT
            setData();
		}
		SAMMeta(CharString const &name){
            SEQAN_CHECKPOINT
            setData();
			qName = name;
		}
		bool operator == (SAMMeta const &other)
		{
			return qName == other.qName;
		}
	
    private:
        void setData(){
            flag = 0;
            rName = "*";
            mapQ = 0;
            MRNM = "*";
            mPos = 0;
            iSize = 0;
            qual = "*";
        }
	};
	//TODO: qName + rName + Position should be unique
    // or the adress of the align

	template<typename TFile, typename TGapsSpec, typename TIDString, typename TSAMSpec>
	inline void write(TFile & target,
					  Align<DnaString, TGapsSpec> const & source,
					  TIDString const & myId,
					  Tag<TagSBAM_<TSAMSpec> > ) {
		
		write(target, source, 0, myId, Tag<TagSBAM_<TSAMSpec> >());
	
	}
	
	//TODO: documentation
	template<typename TFile, typename TGapsSpec, typename TIDString, typename TSAMSpec>
	inline void write(TFile & target,
					  Align<DnaString, TGapsSpec> const & source,
					  const int refPosition,
					  TIDString const & myId,
					  Tag<TagSBAM_<TSAMSpec> > ) {
	
        int nrOfReads = length(rows(source));
        
        String<SAMMeta> metas;
        resize(metas, nrOfReads, Generous());
        String<unsigned int> IDs;
        resize(IDs, nrOfReads, Generous());
		
		for (int i = 0; i < refPosition; ++i){			
            assignValue(IDs, i, i);
            
			std::stringstream ss("Read_");
			ss << (i + 1);
            SAMMeta meta(ss.str());
            assignProperty(metas, i, meta);
		}
		
		for (int i = (refPosition + 1); i < nrOfReads; ++i){			
            assignValue(IDs, i, i);
            
			std::stringstream ss("Read_");
			ss << (i + 1);
            SAMMeta meta(ss.str());
            assignProperty(metas, i, meta);
		}
        write(target, source, refPosition, metas, IDs, myId, Tag<TagSBAM_<TSAMSpec> >());
         
	}
	
	//TODO: documentation
	template<typename TFile, typename TGapsSpec, typename TIDString, typename TMeta, typename TSAMSpec>
	inline void write(TFile & target,
					  Align<DnaString, TGapsSpec> const & source,
					  const int refPosition,
					  const String<TMeta> & metas,
					  const String<unsigned int> & IDs,
					  TIDString const &,
					  Tag<TagSBAM_<TSAMSpec> > ) {

		int nrOfReads = length(rows(source));

		unsigned int refID = getValue(IDs, refPosition);

		for (unsigned int i = 0; i < refPosition; ++i){
			unsigned int readID = getValue(IDs, i);
            SAMMeta meta = getProperty(metas, readID);
			_write_alignment(target, source, refPosition, i, meta, Tag<TagSBAM_<TSAMSpec> >());
		}
		
		for (unsigned int i = (refPosition + 1); i < nrOfReads; ++i){
			unsigned int readID = getValue(IDs, i);
            SAMMeta meta = getProperty(metas, readID);
			_write_alignment(target, source, refPosition, i, meta, Tag<TagSBAM_<TSAMSpec> >());
		}
	}
	
    //TODO: documentation
    // function that writes addional meta data for an alignment
	template<typename TFile>
    inline void _write_additional_meta(TFile & target,
                                       String<Triple<CharString, char, CharString> > const & meta,
                                       Tag<TagSBAM_<_SAM> > )
    {
        typedef String<Triple<CharString, char, CharString> > TList;
        typedef typename Iterator<TList >::Type TIter;
        
        TIter elem = begin(meta);
        
        while(elem != end(meta)){
            _streamPut(target, '\t');
            
            _streamWrite(target, value(elem).i1);
            _streamPut(target, ':');
            _streamPut(target, value(elem).i2);
            _streamPut(target, ':');
            _streamWrite(target, value(elem).i3);
            
            ++ elem;
        }
    }
    
    //TODO: documentation
	// function to write pair of rows is sam format alignment line
	template<typename TFile, typename TGapsSpec, typename TMeta>
	inline void _write_alignment(TFile & target,
            Align<DnaString,TGapsSpec> const & source,
            const int ref_row_pos_int,
            const int read_row_pos_int,
            const TMeta & meta,
            Tag<TagSBAM_<_SAM> > ) {
		
		typedef Align<DnaString, TGapsSpec> const TAlign;
		typedef	typename Row<TAlign>::Type TRow;
		typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
		typedef typename Position<TAlign>::Type TPosition;
		typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
		
        
		// count characters for QUAL
		int count = 0;
		
		//_streamPutInt(target, ref_row_pos_int); _streamPut(target, ','); _streamPutInt(target, read_row_pos_int);
		_streamWrite(target, meta.qName); // QNAME
		_streamPut(target, '\t');
        _streamPutInt(target, meta.flag); // FLAG
		_streamPut(target, '\t');
		_streamWrite(target, meta.rName); // RNAME
		_streamPut(target, '\t');
		_streamPutInt(target, sourceBeginPosition(row(source, ref_row_pos_int))); // POS
        _streamPut(target, '\t');
		_streamPutInt(target, meta.mapQ); // MAPQ
		_streamPut(target, '\t');
        
		// ======= Write CIGAR =======
        // get reference and read row from alignment
		TRow & ref_row = row(source, ref_row_pos_int);
		TRow & read_row = row(source, read_row_pos_int);
		
		// create iteraters in both rows at the begin position
		TIter ref_col_pos = begin(ref_row);
		TIter read_col_pos = begin(read_row);
		TIter ref_end = end(ref_row);
		
		// for leading gaps 
		int offset = beginPosition(read_row) - beginPosition(ref_row);
		// in read
		if(0 < offset){
			count += offset;
			_streamPutInt(target, offset);
			_streamPut(target, 'D');
			for(; 0 < offset; --offset){
				goNext(ref_col_pos);
			}
		} 
		// in reference
		if(0 > offset){
			_streamPutInt(target, -offset);
			_streamPut(target, 'I');
			for(; offset < 0; ++offset){
				goNext(read_col_pos);
			}
		}
		
		//TODO: anders iterieren, ohne den sonderfall der ersten schleife beachten zu muessen.
		char old_operation = 'S'; // for the start only
		char new_operation = 'M';
		int operation_count = 0;
		
		// from first common character iterate simultaneous
		for (; false and ref_col_pos != ref_end or !atEnd(read_col_pos); goNext(ref_col_pos), goNext(read_col_pos)) {
			
			// if ref is a gap
			if(isGap(ref_col_pos)) {
				// if read is a gap: Padded Gap
				if(isGap(read_col_pos)) {
					new_operation = 'P';
				}
				// if read is a Character: Insertion
				else {
					new_operation = 'I';
				}
			}
			// if ref is a character
			else {
				// if read is a gap: Deletion
				if(isGap(read_col_pos)) {
					new_operation = 'D';
				}
				// if read is a character: Match/Mismatch
				else {
					new_operation = 'M';
				}
			}
			
			// if operation changed
			if(operation_count != 0 and old_operation != new_operation) {
				// write old operation to file
				_streamPutInt(target, operation_count);
				_streamPut(target, old_operation);
				operation_count = 0;
			}
			old_operation = new_operation;
			++operation_count;
		}
		// write last operation
		_streamPutInt(target, operation_count);
		_streamPut(target, old_operation);
		
		// for ending gaps 
		offset = endPosition(read_row) - endPosition(ref_row);
		// in read
		if(0 > offset){
			count += -offset;
			_streamPutInt(target, -offset);
			_streamPut(target, 'D');
		} 
		// in reference
		if(0 < offset){
			_streamPutInt(target, offset);
			_streamPut(target, 'I');
		}
		
		//======== Write mate infos ========
        _streamPut(target, '\t');
		_streamWrite(target, meta.MRNM);	// MRNM
        _streamPut(target, '\t');
		_streamPutInt(target, meta.mPos); 		// MPOS
        _streamPut(target, '\t');
		_streamPutInt(target, meta.iSize); 		// ISIZE
        _streamPut(target, '\t');
		
		//======== Write read sequence ========
		// set iterator back on the beginning.
		read_col_pos = begin(read_row);
		// iterate over sequence
		//TODO: vielleicht besser direkt ueber den quellstring
		for (; !atEnd(read_col_pos); goNext(read_col_pos)) {
			// write only non-gap characters
			if(!isGap(read_col_pos)) {
				_streamPut(target, getValue(read_col_pos));
			}
		}
		
		//======== Write Qual ========
        _streamPut(target, '\t');
        _streamWrite(target, meta.qual);

        //======== Write addional Metadata ==
        _write_additional_meta(target, meta.additional, SAM());
        
        // end of line
		_streamPut(target, '\n');
         
	}
	
	//TODO: documentation
	template <typename TFile, typename TGapsSpec, typename TIDString, typename TSAMSpec>
	inline void _writeAlignment(
		TFile & target,
		Align<DnaString, TGapsSpec> const & source,
		int ref,
		int read,
		const CharString & refName,
		const CharString & readName,
		BAM)
	{
		//TODO: binaeres Format implementieren
	}
/*
    int _pack(CharString key){
        // String key should have the lenth 2
        // write the character values of the string in the lower 2 bytes of i
        int i = ( (value(key, 0) << 8) | (value(key, 1)) );
        return i;
    }
    
    CharString _unpack(int n){
        unsigned int last8bit = 255;
        
        CharString str;
        resize(str, 2);
        // get the last byte as char
        value(str, 1) = n & last8bit;
        // shift n by one byte so that the last byte is the next character
        n = n >> 8;
        value(str, 0) = n & last8bit;
        
        return str;
    }
*/    
    template <typename TFile, typename TGapsSpec>
    void
    read(TFile & file,
         Align<DnaString, TGapsSpec> & data,
         SAMMeta & meta,
         SAM)
    {
        while(! _streamEOF(file) && _streamGet(file) != '\n' ){
            while(_streamGet(file) != '\t'){
                
            }
            
        }
    }

	template <typename TFile>
    void
    read(TFile & file,
           String<Pair<CharString, String<Pair<CharString, CharString> > > > & data,
           SAM)
    {
        typedef typename Position<TFile>::Type TPos;
        typedef Pair<CharString, CharString> TTag;
        typedef String<TTag> TTags;
        typedef Pair<CharString, TTags> TLine;
        
        TPos posStart = _streamTellG(file);
        
        char c;
        unsigned int countLines = 0;
        
        // count line that start with '@'
        while(_streamGet(file) == '@'){
            ++countLines;
            
            while(! _streamEOF(file) ){
                c = _streamGet(file); 
                
                if (c == '\n'){
                    break;
                }
            }
        }
        
        // resize date to the number of counted lines
        resize(data, countLines);
        
        // return to first position in file
        _streamSeekG(file, posStart);
        
        // create identifier string
        CharString identifier;
        resize(identifier, 2);
        
        countLines = 0;
        while(_streamGet(file) == '@'){
            // read first two characters and safe them as identifier
            assignValue(identifier, 0, _streamGet(file));
            assignValue(identifier, 1, _streamGet(file));
            
            // skip the following tab or space
            c = _streamGet(file);
            
            // save start of tags
            TPos posTags = _streamTellG(file);
            
            // check how many tags are in this line and how long they are
            unsigned int countTagsInLine = 0;
            String<int> tagLengths;
            resize(tagLengths, 10); //TODO: correct maximal number of tags
            
            while(c != '\n' && c != '\r' && !_streamEOF(file)){
                // skip idenifier and ':'
                _streamGet(file);
                _streamGet(file);
                _streamGet(file);
                
                c = _streamGet(file);
                unsigned int countTagLength = 0;
                while(c != ' ' && c != '\t' && c != '\n' && c != '\r' && !_streamEOF(file)){
                    ++countTagLength;
                    c = _streamGet(file);
                }
                assignValue(tagLengths, countTagsInLine, countTagLength);
                
                ++countTagsInLine;
            }
            
            _streamSeekG(file, posTags);
            
            TTags tags;
            resize(tags, countTagsInLine);
            
            for(int i = 0; i < countTagsInLine; ++i){
                CharString tagIdentifier;
                resize(tagIdentifier, 2);
                
                // read first two characters and safe them as tag identifier
                assignValue(tagIdentifier, 0, _streamGet(file));
                assignValue(tagIdentifier, 1, _streamGet(file));
                
                // skip ':'
                _streamGet(file);
                
                CharString tagValue;
                resize(tagValue, value(tagLengths, i));
                
                for(int j = 0; j < value(tagLengths, i); ++j){
                    assignValue(tagValue, j, _streamGet(file));
                }
                
                TTag tag(tagIdentifier, tagValue);
                assignValue(tags, i, tag);
                
                // skip space or tab
                c = _streamGet(file);
            }
            
            assignValue(data, countLines, TLine(identifier, tags));
            
            ++countLines;
            
        }
    }
    
    template<typename TFile>
    inline void write_header_line(TFile & target,
                                  const CharString & cat, 
                                  const String<Pair<CharString, CharString> > & meta,
                                  Tag<TagSBAM_<_SAM> > )
    {
        typedef String<Pair<CharString, CharString> > TList;
        typedef typename Iterator<TList >::Type TIter;
        
        TIter elem = begin(meta);
        TIter lastElem = end(meta); --lastElem;
        while(elem != end(meta)){
            
            _streamWrite(target, value(elem).i1);
            _streamPut(target, ':');
            _streamWrite(target, value(elem).i2);
            
            if(elem != lastElem){
                _streamPut(target, ' ');
            }
            
            ++elem;
        }
    }
    
    template<typename TFile>
    inline void write(TFile & target,
                                  const String<Pair<CharString, String<Pair<CharString, CharString> > > > & meta,
                                  const Tag<TagSBAM_<_SAM> > & sam)
    {
        typedef String<Pair<CharString, String<Pair<CharString, CharString> > > > TList;
        typedef typename Iterator<TList>::Type TIter;
        
        TIter line = begin(meta);
        while(line != end(meta)){
            _streamPut(target, '@');
            _streamWrite(target, value(line).i1);
            _streamPut(target, '\t');
            
            write_header_line(target, value(line).i1, value(line).i2, sam);
            
            _streamPut(target, '\n');
            
            ++line;
        }
        
    }
    
    
//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
