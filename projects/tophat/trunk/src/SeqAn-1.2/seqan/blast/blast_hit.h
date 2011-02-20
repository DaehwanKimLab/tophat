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
 $Id: blast_hit.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_HIT_H
#define SEQAN_HEADER_BLAST_HIT_H


namespace SEQAN_NAMESPACE_MAIN
{

template<typename TBlastHsp, typename TStoreSpec>
class BlastHit;


/**
.Class.BlastHit:
..cat:Blast
..summary:Object for storing Blast hits. 
..signature:BlastHit<TBlastHsp, TSpec>  
..param.TBlastHsp:The type of HSPs that are stored.
..param.TSpec:The specializing type.
...type:Spec.StreamReport
...type:Spec.StoreReport
..remarks:Use Metafunction.Hit to get the BlastHit type used in a BlastReport object.
*/
template<typename TBlastHsp, typename TSpec>
class BlastHit<TBlastHsp, StoreReport<TSpec> > 
{
	public:
		String<char> name;
		unsigned int length; //length of whole sequence 

		String<TBlastHsp> hsps;
		
		BlastHit()
		{
		SEQAN_CHECKPOINT
		}

		//BlastHit(String<char> name, unsigned int len)
		//{
		//	name = name;
		//	length = length;
		//	clear(hsps);
		//}

		BlastHit(BlastHit const& other)
		{
		SEQAN_CHECKPOINT
			assign(hsps,other.hsps);
			name = other.name;
			length = other.length;
		}

		BlastHit & operator = (BlastHit const & other)
		{
		SEQAN_CHECKPOINT
			assign(hsps,other.hsps);
			name = other.name;
			length = other.length;
			return *this;
		}

		~BlastHit()
		{
		}

};



template<typename TBlastHsp, typename TSpec>
inline void
clear(BlastHit<TBlastHsp, StoreReport<TSpec> >& blastHit)
{
SEQAN_CHECKPOINT
	
	for(unsigned int i = 0; i < length(blastHit.hsps); ++i)
		clear(blastHit.hsps[i]);
	resize(blastHit.hsps,0);
	resize(blastHit.name,0);
}


template<typename TBlastHsp, typename TStoreSpec>
inline String<char> &
name(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.name;
}

template<typename TBlastHsp, typename TStoreSpec>
inline String<char> 
getName(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.name;
}

template<typename TBlastHsp, typename TStoreSpec>
inline unsigned int &
length(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.length;
}

template<typename TBlastHsp, typename TStoreSpec>
inline unsigned int 
getLength(BlastHit<TBlastHsp, TStoreSpec>& blastHit)
{
SEQAN_CHECKPOINT
	return blastHit.length;
}

// for StoreReport only
template<typename TBlastHsp, typename TSpec>
inline unsigned int 
numHsps(BlastHit<TBlastHsp, StoreReport<TSpec> >& blastHit)
{
SEQAN_CHECKPOINT
	return length(blastHit.hsps);
}		


/////////////////////////////////////////////////////////////////////

//parse BlastHit
template<typename TFile, typename TChar, typename TBlastHit>
inline typename Position<TFile>::Type
_parseBlastHit(TFile & file,
			TChar & c, 
			TBlastHit & hit)
{
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Hsp<TBlastHit>::Type TBlastHsp;

	String<char> pword;
	int pint;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	if(_parse_untilBeginLine(file,c,'>'))
	{
		start_pos = _streamTellG(file);
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
			//c = _streamGet(file);
			bool in_hit = true;
			TPosition act_hsp_pos,next_hsp_pos;
			while(in_hit){
				act_hsp_pos = _streamTellG(file);
				TBlastHsp hsp;
				next_hsp_pos = _parseBlastHsp(file,c,hsp);
				//resize(hit.hsps, length(hit.hsps)+1);
				//hit.hsps[length(hit.hsps)-1] = hsp;
				appendValue(hit.hsps,hsp);
				//append(hit.hsps,hsp);
				if(next_hsp_pos == act_hsp_pos)
					in_hit = false;
				if(next_hsp_pos == (TPosition)0)
					return (TPosition) 0;
			}
			_streamSeekG(file,next_hsp_pos);
			c = _streamGet(file);
			if(_parse_untilBeginLine(file,c,'>'))
				return _streamTellG(file);
		}

		
	}//end hit
	_streamSeekG(file,act_pos);
	return act_pos;
}






////////////////////////// MetaFunctions ////////////////////////////////

template<typename TBlastHsp, typename TStoreSpec>
struct Value<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Value<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hit<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hit<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef BlastHit<TBlastHsp,TStoreSpec> Type;
};



template<typename TBlastHsp, typename TStoreSpec>
struct Hsp<BlastHit<TBlastHsp, TStoreSpec> > 
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TStoreSpec>
struct Hsp<BlastHit<TBlastHsp, TStoreSpec> const> 
{
	typedef TBlastHsp Type;
};


/////////////komisch///////////

	template<typename TBlastHsp, typename TStoreSpec>
	struct Hsp<String<BlastHit<TBlastHsp, TStoreSpec> > > 
	{
		typedef TBlastHsp Type;
	};

	template<typename TBlastHsp, typename TStoreSpec>
	struct Hsp<String<BlastHit<TBlastHsp, TStoreSpec> const> > 
	{
		typedef TBlastHsp Type;
	};

	template<typename TBlastSpec, typename TStoreSpec>
	struct Hsp<String<BlastHsp<TBlastSpec, TStoreSpec> > > 
	{
		typedef BlastHsp<TBlastSpec, TStoreSpec> Type;
	};

	template<typename TBlastSpec, typename TStoreSpec>
	struct Hsp<String<BlastHsp<TBlastSpec, TStoreSpec> const> > 
	{
		typedef BlastHsp<TBlastSpec, TStoreSpec> Type;
	};




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
