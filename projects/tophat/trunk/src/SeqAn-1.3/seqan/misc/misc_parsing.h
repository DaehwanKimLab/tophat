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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_MISC_PARSING_H
#define SEQAN_HEADER_MISC_PARSING_H

#include <cmath>


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General parsing funtions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline void 
_parseSkipLine(TFile& file, TChar& c)
{
	if (c == '\n') {
		c = _streamGet(file);
		return;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n') break;
	}
	c = _streamGet(file);
}
//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TChar>
inline void 
_parseSkipWhitespace(TFile& file, TChar& c)
{
	if ((unsigned) c > 32) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((unsigned) c > 32) break;
	}
}

template<typename TFile, typename TChar>
inline void 
_parseSkipSpace(TFile& file, TChar& c)
{
	if (c != '\t' && c != ' ') return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c != '\t' && c != ' ') break;
	}
}


/**
.Internal._parseSkipUntilChar:
..summary:Skip to the next ocurrence of x in file.
..cat:Miscenalleous
..signature:_parseSkipUntilChar(file, x, c)
..param.file:The file to read from.
..param.x:The character to skip to.
..param.c:Parser state character.
 */
template<typename TFile, typename TChar>
inline void 
_parseSkipUntilChar(TFile& file, const TChar &x, TChar& c)
{
	if (c == x) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == x) break;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parseIsDigit(TChar const c)
{
	return (((unsigned) c >  47) && ((unsigned) c <  58));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parseIsLetter(TChar const c)
{
	return ( (((unsigned) c > 64) && ((unsigned) c < 91)) || (((unsigned) c > 96) && ((unsigned) c < 123)) );
}

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): The name of this function is WRONG.
template<typename TChar>
inline bool
_parseIsAlphanumericChar(TChar const c)
{
	return ((_parseIsDigit(c)) || (_parseIsLetter(c)) || (c == '_') || (c == '.') || (c == '-') || (c == '|') || (c == '/') || (c == ':'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline int
_parseReadNumber(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsDigit(c)) break;
		append(str, c);
	}
 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline double
_parseReadDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parseReadIdentifier(TFile & file, TChar& c)
{
	// Read identifier
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsAlphanumericChar(c)) break;
		append(str, c);
	}
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline char
_parseReadChar(TFile & file, TChar& c)
{
    char result = c;
    if (!_streamEOF(file))
        c = _streamGet(file);
    return result;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TChar>
inline void
_parseReadIdentifier(TFile & file, TString& str, TChar& c)
{
	// Read identifier
	append(str, c, Generous());
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsAlphanumericChar(c)) break;
		append(str, c, Generous());
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parseReadWord(TFile & file, TChar& c)
{
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsLetter(c)) break;
		append(str, c);
	}
	return str;
}


// parse word up to a maximum length
template<typename TFile, typename TChar, typename TSize>
inline String<char>
_parseReadWord(TFile & file, TChar& c, TSize max_len)
{
	// Read word
	String<char> str(c);
	--max_len;
	TSize i = 0;
	while (!_streamEOF(file) ) {
		c = _streamGet(file);
		if (!_parseIsLetter(c) || i >= max_len) break;
		append(str, c);
		++i;
	}
	return str;
}




//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar>
inline String<char>
_parseReadFilepath(TFile& file, TChar& c)
{
	String<char> str(c);
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return str;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
		append(str, c);
	}
	typename Iterator<String<char>,Rooted >::Type str_it = end(str);	
	while(str_it != begin(str)) {
		--str_it;
		if(*str_it != ' ' && *str_it != '\t'){
		++str_it;
		break;
		}
	}
	resize(str,position(str_it));
	return str;
}


//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar>
inline String<char>
_parseReadWordUntilWhitespace(TFile& file, TChar& c)
{
	String<char> str(c);
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return str;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
		append(str, c);
	}
	return str;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TString>
inline void
_parseReadSequenceData(TFile & file,
						TChar & c,
						TString& str)
{
	SEQAN_CHECKPOINT

	append(str, c);

	// Read sequence
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parseIsLetter(c)) break;
		else append(str, c);
	}
}



template<typename TFile, typename TChar>
inline void 
_parseSkipBlanks(TFile& file, TChar& c)
{
	if ((c != ' ') && (c != '\t')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c != ' ') && (c != '\t')) break;
	}
}

template<typename TFile, typename TChar>
inline void 
_parseSkipLine2(TFile& file, TChar& c)
{
	if (c != '\n' && c != '\r')
		while (!_streamEOF(file)) {
			c = _streamGet(file);
			if (c == '\n' || c == '\r') break;
		}
	if (!_streamEOF(file))
		c = _streamGet(file);
}



//////////////////////////////////////////////////////////////////////////////
template<typename TFile, typename TChar>
inline double
_parseReadEValue(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT

	// Read number
	String<char> str(c);
	bool e = false;
	double val1 = 0;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if(!e && c == 'e'){
			e = true;
			val1 = atof(toCString(str));
			c = _streamGet(file);
			resize(str,0);
		}
		if (!_parseIsDigit(c) && c != '.' && c != '-' && c != '+') break;
		append(str, c);
	}
	if(e)
	{
		return val1 * pow((double)10.0,(double)atof(toCString(str)));
	}	
 	else 
		return (double)atof(toCString(str));
}





/////////////////////////////////////////////////////////////////////////////////
// read floating point value
template<typename TFile, typename TChar>
inline float
_parseReadFloat(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c != '.' && c != ',' && !_parseIsDigit(c)) break;
		append(str, c);
	}
 	return atof(toCString(str));
}




/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with character x (skip whitespaces)
// zeigt am ende darauf!!!
template<typename TFile, typename TChar>
inline bool
_parseUntilBeginLine(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	_parseSkipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		_parseSkipLine(file, c);
		_parseSkipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with word
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	_parseSkipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parseReadWord(file,c,len))
				break;
		_parseSkipLine(file, c);
		_parseSkipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with word (parse no more than num_lines lines)
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len, TSize num_lines)
{
SEQAN_CHECKPOINT
	_parseSkipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	TSize i = 0;
	bool found = false;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parseReadWord(file,c,len))
			{
				found = true;
				break;
			}
		if(i >= num_lines)
			break;
		++i;
		_parseSkipLine(file, c);
		_parseSkipWhitespace(file,c);
	}
	if(!_streamEOF(file) && found) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with one of the characters in string x (skip whitespaces)
//zeigt am ende darauf!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntilBeginLineOneOf(TFile & file, TChar& c, String<TChar> & x, TSize len)
{
SEQAN_CHECKPOINT
	_parseSkipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	bool found = false;
	while (!_streamEOF(file)){
		for(int i = 0; i < len; ++i)
			if(c == x[i]) 
			{
				found = true;
				break;
			}
		if(found) break;
		_parseSkipLine(file, c);
		_parseSkipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until c == x
//zeigt am ende darauf!
template<typename TFile, typename TChar>
inline bool
_parseUntil(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}



/////////////////////////////////////////////////////////////////////////////////
//parse until word
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntil(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parseReadWord(file,c,len))
				break;
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until c == x or new line
//zeigt am ende darauf!
template<typename TFile, typename TChar>
inline bool
_parseLineUntil(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		if (c == '\n' || c == '\r')
		{
			_streamSeekG(file,pos);
			c = c_before;
			return false;
		}
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse this line until word
//zeigt am ende hinter wort if true, oder auf ende der zeile
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseLineUntil(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
		{	if(word == _parseReadWord(file,c,len))
				break;
		}
		else if (c == '\n' || c == '\r')
			{
				_streamSeekG(file,pos);
				c = c_before;
				return false;
			}
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}






}

#endif

