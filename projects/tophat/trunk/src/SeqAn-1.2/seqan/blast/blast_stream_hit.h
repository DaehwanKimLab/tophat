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
 $Id: blast_stream_hit.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_STREAM_HIT_H
#define SEQAN_HEADER_BLAST_STREAM_HIT_H


namespace SEQAN_NAMESPACE_MAIN
{


////////////////////////////////////////////////////////////////////////////////////////////
//  Blast Hit storing only one hsp at a time
////////////////////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
class BlastHit<TBlastHsp, StreamReport<TFile> > 
{
	public:
		typedef typename Position<TFile>::Type TPosition;

		String<char> name;
		unsigned int length; //length of whole sequence  
		TBlastHsp act_hsp;
		TPosition begin_pos, first_hsp_pos;
		
		BlastReport<TBlastHsp,StreamReport<TFile> >* data_host;

	
		BlastHit()
		{
		}

		~BlastHit()
		{
		}

};




//parse BlastHit
template<typename TFile, typename TChar, typename TBlastSpec>
inline typename Position<TFile>::Type
_parseBlastHit(TFile & file,
			TChar & c, 
			BlastHit<TBlastSpec,StreamReport<TFile> > & hit)
{
	typedef typename Position<TFile>::Type TPosition;
	typedef BlastHit<TBlastSpec,StreamReport<TFile> > TBlastHit;
	typedef typename Hsp<TBlastHit>::Type TBlastHsp;

	String<char> pword;
	int pint;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	if(_parse_untilBeginLine(file,c,'>'))
	{
		hit.begin_pos = _streamTellG(file);
		c = _streamGet(file);
		pword = _parse_readWord(file, c);
		while (!_streamEOF(file) && c != '\n' && c != '\r')
			pword += _parse_readWord(file, c);
		if(pword[length(pword)-1] == ' ')
			resize(pword,length(pword)-1);
		hit.name = pword;
		_parse_skipWhitespace(file,c);
		String<char> search = "Length";
		if(_parse_untilBeginLine(file,c,search,6))
		{
			_parse_skipWhitespace(file,c);
			if(c == '=')
				c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			pint = _parse_readNumber(file, c);
			hit.length = pint;
		}
//		TPosition temp = _streamTellG(file);
		//foreach Hsp
		//if(_parse_untilBeginLine(file,c,'S') && _parse_readWord(file,c)=="Score")
		search = "Score";
		if(_parse_untilBeginLine(file,c,search,5))
		{
			hit.first_hsp_pos = _streamTellG(file);
			//if(_parse_untilBeginLine(file,c,'>'))
			//	return _streamTellG(file);
		}

		
	}//end hit
	_streamSeekG(file,act_pos);
	c = '>';
	return act_pos;
}







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
