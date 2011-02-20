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
 $Id: blast_parsing.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_PARSING_H
#define SEQAN_HEADER_BLAST_PARSING_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Internal Blast Parsing Functions
// remark: uses a lot of functions from graph_utility_parsing.h
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
// read alignment string (letters and gaps)
// steht am ende dahinter 
template<typename TFile, typename TChar>
inline String<char>
_parse_readAlignmentString(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c) && !(c=='-')) break;
		append(str, c);
	}
	return str;
}





////////////////////////////////////////////////////////////////////////////
// parses query and database name (file should be pointing to the beginning of the file)
template<typename TFile, typename TChar>
typename Position<TFile>::Type
_parse_readQueryAndDBName(TFile & file,
						  TChar & c,
						  String<char> & query_name,
						  String<char> & db_name)
{
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	TChar c_before = c;
	TPosition act_pos = _streamTellG(file);
	TPosition query_pos,db_pos;

	//String<char> delim = "DQ";
	//_parse_untilBeginLineOneOf(file,c,delim,2);


	String<char> query = "Query";
	if(_parse_untilBeginLine(file,c,query,6))
		query_pos = _streamTellG(file);
	else
		return act_pos;
	_streamSeekG(file,act_pos);
	c = c_before;
	String<char> database = "Database";
	if(_parse_untilBeginLine(file,c,database,8))
		db_pos = _streamTellG(file);
	else
		return act_pos;
	
	
	//getQueryName
	_streamSeekG(file,query_pos);
	_parse_skipWhitespace(file,c);
	c = _streamGet(file);
	_parse_skipWhitespace(file,c);
	query_name = _parse_readWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		query_name += _parse_readWord(file, c);
	
	//getDBName
	_streamSeekG(file,db_pos);
	c = _streamGet(file);
	_parse_skipWhitespace(file,c);
	db_name = _parse_readWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		db_name += _parse_readWord(file, c);
	_parse_skipWhitespace(file,c);
		
	c = _streamGet(file);

	return _streamTellG(file); 
	
}








}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
