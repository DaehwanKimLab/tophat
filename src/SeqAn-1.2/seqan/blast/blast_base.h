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
 $Id: blast_base.h 4714 2009-08-24 13:38:02Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/
#ifndef SEQAN_HEADER_BLAST_BASE_H
#define SEQAN_HEADER_BLAST_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Blast Report types


struct FullInfo;
//...remarks:BasicInfo stores begin and end positions on query and database sequence, as well as the alignment. FullInfo stores additional information such as score, e-value...

struct BasicInfo;



/**
.Spec.StoreReport:
..cat:Blast
..general:Class.BlastReport
..summary:BlastReport specialization. Parses a Blast report and stores all hits and HSPs.
..signature:BlastReport<TBlastHsp,StoreReport<TSpec> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TSpec:The specializing type.
...default:BasicInfo
...type:Spec.BasicInfo
...type:Spec.FullInfo 
..include:blast.h
*/
//
//...remarks:BasicInfo only stores query name, database name and a String of all hits found. FullInfo also stores the following 
//parameters: lambda, k, h, gapped_lambda, gapped_k, gapped_h, gap_open, gap_extension; String<char> matrix; double min_expect;
//
template<typename TInfoSpec = BasicInfo>
struct StoreReport;		//stores the whole report


/**
.Spec.StreamReport:
..cat:Blast
..general:Class.BlastReport
..summary:BlastReport specialization that works on a file stream (parses hits/HSPs when iterating over them).
..signature:BlastReport<TBlastHsp,StreamReport<TFile> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TFile:The type of the stream.
...default:std::fstream
..include:blast.h
*/
template<typename TFile = std::fstream>    //works on a stream
struct StreamReport;




//////////////////////////////////////////////////////////////////////////////
//Blast Meta functions


/**
.Metafunction.Hit:
..summary:Blast Hit type of a Blast object.
..signature:Hsp<T>::Type
..param.T:A Blast report object.
...type:Class.BlastReport
..returns.param.Type:BlastHit type.
*/
template<typename T>
struct Hit;

/**
.Metafunction.Hsp:
..summary:Blast HSP type of a Blast object.
..signature:Hsp<T>::Type
..param.T:A Blast object.
...type:Class.BlastReport
...type:Class.BlastHit
..returns.param.Type:BlastHsp type.
*/
template<typename T>
struct Hsp;


//////////////////////////////////////////////////////////////////////////////
// Blast Tag

struct TagBlast_;
typedef Tag<TagBlast_> const Blast;


//////////////////////////////////////////////////////////////////////////////
// Blat Tag 
//already defined in seeds/seedHandlingTags.h (49)
//moved to basic_tag
/*
struct TagBlat_;
typedef Tag<TagBlat_> const Blat;
*/
//////////////////////////////////////////////////////////////////////////////



/**
.Spec.BlastN:
..cat:Blast
..general:Class.BlastHsp
..summary:For BlastN Blast reports.
..signature:BlastN
..include:blast.h
*/


/**
.Spec.BlastP:
..cat:Blast
..general:Class.BlastHsp
..summary:For BlastP Blast reports.
..signature:BlastP
..include:blast.h
*/



struct TagBlastN_;
struct TagMegaBlast_;
struct TagBlastP_;
struct TagBlastX_;
struct TagTBlastN_;
struct TagTBlastX_;

template<typename TSpec = TagBlastN_>
class NucleotideBlast{
public:
	NucleotideBlast(){}
	~NucleotideBlast(){}
};

template<typename TSpec = TagBlastP_>
class ProteinBlast{
public:
	ProteinBlast(){}
	~ProteinBlast(){}
};

typedef NucleotideBlast<TagBlastN_> BlastN;
typedef NucleotideBlast<TagMegaBlast_> MegaBlast;

typedef ProteinBlast<TagBlastP_> BlastP;
typedef ProteinBlast<TagBlastX_> BlastX;
typedef ProteinBlast<TagTBlastN_> TBlastN;
typedef ProteinBlast<TagTBlastX_> TBlastX;



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
