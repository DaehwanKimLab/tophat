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
 $Id: blast_stream_hsp_iterator.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_STREAM_HSP_ITERATOR_H
#define SEQAN_HEADER_BLAST_STREAM_HSP_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Hsp Iterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



////**
//.Spec.HspIterator:
//..cat:Blast
//..summary:Hsp iterator for @Class.BlastHit@.
//..signature:Iterator<TBlastHit, HspIterator>
//..param.TBlastHit:A Blast hit.
//...type:Class.BlastHit
//..general:Class.Iter
//*/
template<typename TBlastHsp, typename TFile>
class Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > 
{
public:

	typedef BlastHit<TBlastHsp,StreamReport<TFile> > TBlastHit;
	typedef typename Position<TFile>::Type TPosition;

	TBlastHsp data_hsp;
	TBlastHit* data_host;
	TPosition data_pos, data_next_pos, data_hsp_begin_pos;
	bool data_at_end;


	Iter()	
	{
	data_at_end = false;
	}
	
	Iter(TBlastHit & blast) 
	{
	SEQAN_CHECKPOINT
		data_host = &blast; 
		data_pos = blast.first_hsp_pos;
		data_next_pos = data_pos;
		data_hsp_begin_pos = (TPosition) 0;
		data_at_end = false;
	}

	Iter(Iter const& other): 
		data_host(other.data_host), 
		data_pos(other.data_pos), 
		data_next_pos(other.data_next_pos), 
		data_at_end(other.data_at_end),
		data_hsp(other.data_hsp),
		data_hsp_begin_pos(other.data_hsp_begin_pos) 
	{
	SEQAN_CHECKPOINT
	}

	~Iter() 
	{
	SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & other) 
	{
		SEQAN_CHECKPOINT
		if (this == &other) return *this;
		data_host = other.data_host;
		data_pos = other.data_pos;
		data_at_end = other.data_at_end;
		data_hsp = other.data_hsp;
		data_next_pos = other.data_next_pos;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Blast StreamHspIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHsp, typename TFile>
struct Iterator<BlastHit<TBlastHsp,StreamReport<TFile> >, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastHit<TBlastHsp,StreamReport<TFile> > const, HspIterator>
{	
	typedef Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> >, HspIterator>
{	
	typedef Iter<typename Hit<BlastReport<TBlastHsp,StreamReport<TFile> > >::Type, StreamBlastIterator<HspIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> > const, HspIterator>
{	
	typedef Iter<typename Hit<BlastReport<TBlastHsp,StreamReport<TFile> > >::Type, StreamBlastIterator<HspIterator> > Type;
};


//////////////////////////////////////////////////////////////////////////////
template<typename TBlastHsp, typename TFile>
struct Value<Iter<BlastHit<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};

template<typename TBlastHsp, typename TFile>
struct Value<Iter<BlastHit<TBlastHsp,StreamReport<TFile> > const, StreamBlastIterator<HspIterator> > >
{
	typedef TBlastHsp Type;
};









///.Metafunction.Host.param.T.type:Class.BlastHit
template<typename TBlastHit>
struct Host<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{	
	typedef TBlastHit Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
struct Reference<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type& Type;
};

template<typename TBlastHit>
struct Reference<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
struct GetValue<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type Type;
};

template<typename TBlastHit>
struct GetValue<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >
{
	typedef typename Value<Iter<TBlastHit const, StreamBlastIterator<HspIterator> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast StreamBlastIterator<HspIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit, typename TFile>
inline bool
atBegin(TFile &,
		Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == it.data_host->first_hsp_pos);	
}

//////////////////////////////////////////////////////////////////////////////



template<typename TBlastHit, typename TFile>
inline void
goBegin(TFile &,
		Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_host->first_hsp_pos;
	it.data_next_pos = it.data_pos;
	it.data_at_end = false;
}



//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit, typename TFile>
inline void
goNext(TFile & file,
	   Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(file,it)) 
	{
		if(it.data_pos == it.data_next_pos)
			getNextHspFilePos(file,it);
		if(it.data_pos == it.data_next_pos)
			it.data_at_end = true;
		else
            it.data_pos = it.data_next_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////


//template<typename TBlastHit>
//inline void
//goPrevious(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	if (!atBegin(it)) --it.data_pos;
//}
//

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
inline bool
operator ==(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it1,
			Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos && it1.data_host==it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHit>
inline bool
operator !=(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it1,
			Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos || it1.data_host!=it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline typename GetValue<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type
getValue(TFile & file,
		 Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != it.data_hsp_begin_pos)
	{
		_streamSeekG(file,it.data_pos);
		(it.data_host->data_host)->act_c = ' ';
		it.data_hsp_begin_pos = it.data_pos;
		typename Position<TFile>::Type pot_next_pos = _parseBlastHsp(file,(it.data_host->data_host)->act_c,it.data_hsp);
		if(pot_next_pos > it.data_pos)
			it.data_next_pos = pot_next_pos;
		//if(_parseBlastHit(it.data_host->strm,it.data_host->act_c,it.data_hsp) == it.data_pos)
		//	it.data_at_end = true;
	}
	return it.data_hsp;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline typename Reference<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type
value(TFile & file,
	  Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != it.data_hsp_begin_pos)
		it.data_hsp = getValue(file,it);
	return it.data_hsp;
}


//////////////////////////////////////////////////////////////////////////////
//
//
//template<typename TBlastHit>
//inline typename Host< typename Iterator< typename Host< Iter<TBlastHit, StreamBlastIterator<HspIterator> >::Type>::Type, StreamBlastIterator<HitIterator> >::Type >::Type const&
//hostReport(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	return *(it.data_host->data_host);
//} 

//////////////////////////////////////////////////////////////////////////////



/**
.Function.hostHit:
..cat:Blast
..summary:The BlastHit this iterator is working on.
..signature:hostHit(it)
..param.it:An iterator.
...type:Spec.HspIterator
..returns:A pointer to the host BlastHit.
*/
template<typename TBlastHit>
inline typename Host<Iter<TBlastHit, StreamBlastIterator<HspIterator> > >::Type const&
hostHit(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastHit, typename TFile>
inline bool
atEnd(TFile &,
	  Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
{
SEQAN_CHECKPOINT
//	return (it.data_last_pos != it.data_pos);	
	return it.data_at_end;	
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastHit>
//inline void
//goEnd(Iter<TBlastHit, StreamBlastIterator<HspIterator> >& it)
//{
//	it.data_pos = doof;
//}


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
inline void
getNextHspFilePos(TFile & file,
				  Iter<BlastHit<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HspIterator> >& it)
{
	typedef typename Position<TFile>::Type TPosition;

	_streamSeekG(file,it.data_pos);
	char c = 'e';
	(it.data_host->data_host)->act_c = c;

	_parse_skipWhitespace(file,c);
	_parse_skipLine(file,c);

	TPosition next_event_pos;
	bool last_hit = true;
	String<char> delim = ">";
	if(_parse_untilBeginLine(file,c,'>'))
	{
		last_hit = false;
		next_event_pos = _streamTellG(file);
		if((it.data_host->data_host)->next_report && (next_event_pos > (it.data_host->data_host)->next_report_pos))
			next_event_pos = (it.data_host->data_host)->next_report_pos;
	}
	_streamSeekG(file,it.data_pos);
	c = 'e';

	String<char> search = "Score";
	if(_parse_untilBeginLine(file,c,search,5))
	{
		if(!last_hit && (_streamTellG(file) > next_event_pos))
	        _streamSeekG(file,it.data_pos);
		else
            it.data_next_pos = _streamTellG(file);
	}//end hsp
	else
        _streamSeekG(file,it.data_pos);

}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
