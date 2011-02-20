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
 $Id: blast_stream_report.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_STREAM_REPORT_H
#define SEQAN_HEADER_BLAST_STREAM_REPORT_H


namespace SEQAN_NAMESPACE_MAIN
{


//TODO macht noch nicht so richtig sinn mit dem stringSet 
//+vielleicht w�r ne map von fasta id -> stringset id gut
template<typename TBlastHsp, typename TFile>
class BlastReport<TBlastHsp, StreamReport<TFile> > 
{
	public:
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > TBlastHit;
		typedef typename Position<TFile>::Type TPosition;
	
		String<char> query_name;
		String<char> db_name;
	
	
		TPosition first_hit_pos;
		char act_c; 
		bool hits_found;
		bool next_report;
		TPosition next_report_pos;


		BlastReport()
		{
		SEQAN_CHECKPOINT
			next_report = true;
			next_report_pos = 0;
		}

		BlastReport(BlastReport const& other)
		{
		SEQAN_CHECKPOINT

			query_name = other.query_name;
			db_name = other.db_name;
			act_c = other.act_c;
			hits_found = other.hits_found;
			first_hit_pos = other.first_hit_pos;
			next_report = other.next_report;
			next_report_pos = other.next_report_pos;
		}

		//BlastReport(TFile file)
		//{
		//	read(file,*this, Blast());
		//}

		~BlastReport()
		{
		SEQAN_CHECKPOINT
		}


};





//read a  blast report and set the filestream to the first hit position (first '>')
template<typename TBlastHsp, typename TFile>
void 
read(TFile & file,
	 BlastReport<TBlastHsp, StreamReport<TFile> >& blastObj,	 
	 Tag<TagBlast_>) 
{
SEQAN_CHECKPOINT
 


	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	TValue c = blastObj.act_c;
	if(blastObj.next_report)
		_streamSeekG(file,blastObj.next_report_pos);

	blastObj.next_report = false;

	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	String<char> query_name, db_name;
	
	//get query and database names
	_parse_readQueryAndDBName(file,c,query_name,db_name);
	blastObj.query_name = query_name;
	blastObj.db_name = db_name;

	TPosition after_dbquery_pos = _streamTellG(file);
	TValue c_before = c;

	blastObj.hits_found = false;
	TPosition next_event_pos = after_dbquery_pos;

	String<char> delim = "Reference";
	if(_parse_untilBeginLine(file,c,delim,9))
	{
		blastObj.next_report_pos = _streamTellG(file);
		blastObj.next_report = true;
	}

	_streamSeekG(file,after_dbquery_pos);
	if(_parse_untilBeginLine(file,c,'>'))
	{
		next_event_pos = _streamTellG(file);
		if(!blastObj.next_report || next_event_pos < blastObj.next_report_pos)
		{
			blastObj.hits_found = true;
			blastObj.first_hit_pos = next_event_pos;
		}
		else
			next_event_pos = after_dbquery_pos;
	}
	//get some more values aber erst sp�ter
	//_readParameters(file,c,blastObj) ;

	if(blastObj.hits_found)
	{
		_streamSeekG(file,next_event_pos);
		c = '>';
	}
	else
	{
		if(blastObj.next_report)
		{
			c = ':';
			_streamSeekG(file,blastObj.next_report_pos);
		}
		else
		{
			_streamSeekG(file,next_event_pos);
			c = c_before;
		}
	}

	blastObj.act_c = c;

}



//////////////////// Metafunctions /////////////////////////////

template<typename TBlastHsp, typename TFile>
struct Value<BlastReport<TBlastHsp, StreamReport<TFile> > > 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Value<BlastReport<TBlastHsp, StreamReport<TFile> > const> 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};


template<typename TBlastHsp, typename TFile>
struct Hit<BlastReport<TBlastHsp, StreamReport<TFile> > > 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Hit<BlastReport<TBlastHsp, StreamReport<TFile> > const> 
{
	typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
};


/////////////////////////////////////////////
//todo

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastReport<TBlastHsp, StreamReport<TFile> > > > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastReport<TBlastHsp, StreamReport<TFile> > const> > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastHit<TBlastHsp, StreamReport<TFile> > > > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

	template<typename TBlastHsp, typename TFile>
	struct Hit<String<BlastHit<TBlastHsp, StreamReport<TFile> > const> > 
	{
		typedef BlastHit<TBlastHsp,StreamReport<TFile> > Type;
	};

///////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
