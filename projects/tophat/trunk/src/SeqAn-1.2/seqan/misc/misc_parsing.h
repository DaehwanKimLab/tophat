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

#ifndef SEQAN_HEADER_MISC_PARSING_H
#define SEQAN_HEADER_MISC_PARSING_H



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
_parse_skipLine(TFile& file, TChar& c)
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
_parse_skipWhitespace(TFile& file, TChar& c)
{
	if ((unsigned) c > 32) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((unsigned) c > 32) break;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isDigit(TChar const c)
{
	return (((unsigned) c >  47) && ((unsigned) c <  58));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isLetter(TChar const c)
{
	return ( (((unsigned) c > 64) && ((unsigned) c < 91)) || (((unsigned) c > 96) && ((unsigned) c < 123)) );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isAlphanumericChar(TChar const c)
{
	return ((_parse_isDigit(c)) || (_parse_isLetter(c)) || (c == '_') || (c == '.') || (c == '-') || (c == '|'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline int
_parse_readNumber(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline double
_parse_readDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parse_readIdentifier(TFile & file, TChar& c)
{
	// Read identifier
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c);
	}
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TChar>
inline void
_parse_readIdentifier(TFile & file, TString& str, TChar& c)
{
	// Read identifier
	append(str, c, Generous());
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c, Generous());
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parse_readWord(TFile & file, TChar& c)
{
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
	return str;
}


// parse word up to a maximum length
template<typename TFile, typename TChar, typename TSize>
inline String<char>
_parse_readWord(TFile & file, TChar& c, TSize max_len)
{
	// Read word
	String<char> str(c);
	--max_len;
	TSize i = 0;
	while (!_streamEOF(file) ) {
		c = _streamGet(file);
		if (!_parse_isLetter(c) || i >= max_len) break;
		append(str, c);
		++i;
	}
	return str;
}




//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar>
inline String<char>
_parse_readFilepath(TFile& file, TChar& c)
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
_parse_readWordUntilWhitespace(TFile& file, TChar& c)
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
_parse_readSequenceData(TFile & file,
						TChar & c,
						TString& str)
{
	SEQAN_CHECKPOINT

	append(str, c);

	// Read sequence
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		else append(str, c);
	}
}



template<typename TFile, typename TChar>
inline void 
_parse_skipBlanks(TFile& file, TChar& c)
{
	if ((c != ' ') && (c != '\t')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c != ' ') && (c != '\t')) break;
	}
}

template<typename TFile, typename TChar>
inline void 
_parse_skipLine2(TFile& file, TChar& c)
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
_parse_readEValue(TFile & file, TChar& c)
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
		if (!_parse_isDigit(c) && c != '.' && c != '-' && c != '+') break;
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
_parse_readFloat(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c != '.' && c != ',' && !_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atof(toCString(str));
}




/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with character x (skip whitespaces)
// zeigt am ende darauf!!!
template<typename TFile, typename TChar>
inline bool
_parse_untilBeginLine(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
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
_parse_untilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
				break;
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
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
_parse_untilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len, TSize num_lines)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	TSize i = 0;
	bool found = false;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
			{
				found = true;
				break;
			}
		if(i >= num_lines)
			break;
		++i;
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
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
_parse_untilBeginLineOneOf(TFile & file, TChar& c, String<TChar> & x, TSize len)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
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
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
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
_parse_until(TFile & file, TChar& c, TChar x)
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
_parse_until(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
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
_parse_lineUntil(TFile & file, TChar& c, TChar x)
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
_parse_lineUntil(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
		{	if(word == _parse_readWord(file,c,len))
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

