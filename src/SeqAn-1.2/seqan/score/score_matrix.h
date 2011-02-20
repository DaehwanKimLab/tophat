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
  $Id: score_matrix.h 1351 2007-11-30 15:00:38Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SCORE_MATRIX_H
#define SEQAN_HEADER_SCORE_MATRIX_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSequenceValue, typename TSpec>
struct _ScoringMatrixData;

//////////////////////////////////////////////////////////////////////////////

template <typename TSequenceValue = AminoAcid, typename TSpec = Default>
struct ScoreMatrix;


/**
.Tag.File Format.tag.Fasta:
	FASTA file format for sequences.
*/
struct TagScoreMatrixFile_;
typedef Tag<TagScoreMatrixFile_> const ScoreMatrixFile;

//////////////////////////////////////////////////////////////////////////////


/**
.Spec.Score Matrix:
..cat:Scoring
..summary:A general scoring matrix.
..general:Class.Score
..signature:Score<TValue, ScoreMatrix<TSequenceValue, TSpec> >
..param.TValue:Type of the score values.
...default:$int$
..param.TSequenceValue:Type of alphabet underlying the matrix.
...default:$AminoAcid$
*/

template <typename TValue, typename TSequenceValue, typename TSpec>
class Score<TValue, ScoreMatrix<TSequenceValue, TSpec> >
{
public:
	enum
	{
		VALUE_SIZE = ValueSize<TSequenceValue>::VALUE,
		TAB_SIZE = VALUE_SIZE * VALUE_SIZE
	};

	TValue data_tab[TAB_SIZE];
	TValue data_gap_extend;
	TValue data_gap_open;

public:
	Score(TValue _gap_extend = -1):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_extend)
	{
		setDefaultScoreMatrix(*this, TSpec());
	}
	Score(TValue _gap_extend, TValue _gap_open):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_open)
	{
		setDefaultScoreMatrix(*this, TSpec());
	}
	template <typename TString>
	Score(TString const & filename, TValue _gap_extend = -1):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_extend)
	{
		loadScoreMatrix(*this, filename);
	}
	template <typename TString>
	Score(TString const & filename, TValue _gap_extend, TValue _gap_open):
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_open)
	{
		loadScoreMatrix(*this, filename);
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2>
inline TValue
score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc,
	  TVal1 val1,
	  TVal2 val2)
{
	typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
	unsigned int i = (TSequenceValue) val1; // conversion TVal1 => TSequenceValue => integral
	unsigned int j = (TSequenceValue) val2; // conversion TVal2 => TSequenceValue => integral
	return sc.data_tab[i * TScore::VALUE_SIZE + j];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2, typename T>
inline void
setScore(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
	  TVal1 val1,
	  TVal2 val2,
	  T score)
{
	typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
	unsigned int i = (TSequenceValue) val1; // conversion TVal1 => TSequenceValue => integral
	unsigned int j = (TSequenceValue) val2; // conversion TVal2 => TSequenceValue => integral
	sc.data_tab[i * TScore::VALUE_SIZE + j] = score;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSequenceValue, typename TSpec, typename TTag>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
					  TTag)
{
	typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
	TValue const * tab = _ScoringMatrixData<TValue, TSequenceValue, TTag>::getData();
	arrayCopy(tab, tab + TScore::TAB_SIZE, sc.data_tab);
}
template <typename TValue, typename TSequenceValue, typename TSpec>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
					  Default)
{
	typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
	arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());
}

//////////////////////////////////////////////////////////////////////////////

//helper: wrapper for sscanf
inline void
_sscanfValue(char * buf, unsigned int & val)
{
	std::sscanf(buf, "%u", & val);
}
inline void
_sscanfValue(char * buf, int & val)
{
	std::sscanf(buf, "%i", & val);
}
inline void
_sscanfValue(char * buf, float & val)
{
	std::sscanf(buf, "%f", & val);
}
inline void
_sscanfValue(char * buf, double & val)
{
	float f;
	std::sscanf(buf, "%f", & f);
	val = f;
}

//____________________________________________________________________________

template <typename TFile, typename TMeta>
void
readMeta(TFile & fl,
		 TMeta & meta,
		 ScoreMatrixFile)
{
	clear(meta);
	if (_streamEOF(fl)) return;

	typedef typename Value<TMeta>::Type TValue;
	TValue c = _streamGet(fl);

	while (!_streamEOF(fl) && (c == '#'))
	{
		c = _streamGet(fl);
		_stream_appendLine(fl, meta, c);
		appendValue(meta, '\n');
	}
	_streamUnget(fl);
}

//____________________________________________________________________________

template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
void
read(TFile & fl,
	 Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
	 ScoreMatrixFile)
{
	typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
	typedef typename Value<TFile>::Type TFileValue;

	//clears the matrix
	arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());


	//start reading
	if (_streamEOF(fl)) return;

	TFileValue c = _streamGet(fl);
	String<TFileValue> s;


	//search for alphabet line 
	//this is the first line that does not start with '#'
	//e.g.: "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *"
	do
	{
		clear(s);
		_stream_appendLine(fl, s, c);
	} while (!_streamEOF(fl) && (empty(s) || (s[0] == '#')));

	if (_streamEOF(fl)) return;


	//build table to map the alphabet line to the TSequenceValue values
	typedef typename Iterator<String<TFileValue>, Standard>::Type TIterator;
	TIterator it = begin(s);
	TIterator it_end = end(s);
	String<unsigned int> mapping_;
	String<unsigned int> column_;
	for (unsigned int i = 0; it != it_end; ++i)
	{
		if ((*it) != ' ')
		{
			unsigned int pos = (TSequenceValue) *it;		//conversion TFileValue => TSequenceValue => integral
			appendValue(mapping_, pos);
			appendValue(column_, i);						//this marks the end of the column
		}
		++it;
	}

	//read the matrix itself
	while (!_streamEOF(fl))
	{
		clear(s);
		_stream_appendLine(fl, s, c);

		if (empty(s) || (s[0] == '#')) continue; //skip empty lines and comments

		//read first character = alphabet column
		unsigned int row = (TSequenceValue) s[0]; //conversion TFileValue => TSequenceValue => integral
		unsigned int offset = row * TScore::VALUE_SIZE;

		//read rest of line
		unsigned int right;
		unsigned int left = 0;

		TFileValue buf[100]; //100 is enough, believe me!
		buf[99] = 0;

		for (unsigned int i = 0; i < length(column_); ++i)
		{//read column i

			//scan cell into buffer
			right = column_[i];
			TFileValue * it;
			for (it = buf + 99; it >= buf; --right)
			{
				if (right <= left) break;
				if (s[right] == ' ') break;
				--it;
				*it = s[right];
			}

			//parse buffer
			TValue val;
			_sscanfValue(it, val);

			sc.data_tab[offset + mapping_[i]] = val;

			left = column_[i];
		}
	}
}

template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
inline void
read(TFile & fl,
	 Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc)
{
	read(fl, sc, ScoreMatrixFile());
}

//____________________________________________________________________________

template <typename TValue, typename TSequenceValue, typename TSpec, typename TString>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
				TString & filename)
{
		FILE * fl;
		_streamOpen(fl, filename);
		read(fl, sc);
		_streamClose(fl);
}
template <typename TValue, typename TSequenceValue, typename TSpec, typename TString, typename TMeta>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
				TString & filename,
				TMeta & meta)
{
		FILE * fl;
		_streamOpen(fl, filename);
		readMeta(fl, meta, ScoreMatrixFile());
		read(fl, sc, ScoreMatrixFile());
		_streamClose(fl);
}

//____________________________________________________________________________

//wrapper for sprintf
inline void
_sprintfValue(char * buf, unsigned int val)
{
	std::sprintf(buf, "%u", val);
}
inline void
_sprintfValue(char * buf, int val)
{
	std::sprintf(buf, "%d", val);
}
inline void
_sprintfValue(char * buf, float val)
{
	double d = val;
	std::sprintf(buf, "%G", d);
}
inline void
_sprintfValue(char * buf, double val)
{
	std::sprintf(buf, "%G", val);
}

//____________________________________________________________________________

template <typename TSequenceValue, typename TFile, typename TValue, typename TMeta>
void
_writeScoringMatrix(TFile & fl,
					TValue * tab,
					TMeta & meta)
{
	typedef typename Value<TFile>::Type TFileValue;

	enum
	{
		VALUE_SIZE = ValueSize<TSequenceValue>::VALUE,
		TAB_SIZE = VALUE_SIZE * VALUE_SIZE
	};

	typedef typename Value<TFile>::Type TFileValue;

	//write meta
	if (!empty(meta))
	{
		bool line_begin = true;
		for (unsigned int i = 0; i < length(meta); ++i)
		{
			if (line_begin)
			{//escape each line with a starting #
				_streamPut(fl, '#');
				line_begin = false;
			}
			if (meta[i] == '\r') continue;
			if (meta[i] == '\n') line_begin = true;
			_streamPut(fl, meta[i]);
		}
		if (!line_begin) _streamPut(fl, '\n');
	}

	//determine column width
	unsigned int col_width = 1;
	char buf[100]; //100 is enough, believe me!
	for (unsigned int i = 0; i < TAB_SIZE; ++i)
	{
		_sprintfValue(buf, tab[i]);
		unsigned int cell_width = std::strlen(buf);
		if (cell_width > col_width) col_width = cell_width; //compute maximum
	}
	
	++col_width; //increase col_width for an additional blank 

	//write alphabet line
	_streamPut(fl, ' '); // a blank for alphabet column
	for (unsigned int j = 0; j < VALUE_SIZE; ++j)
	{
		TFileValue val = (TSequenceValue) j; //conversion integral => TSequenceValue => TFileValue
		//leading blanks for column j
		for (unsigned int k = 1; k < col_width; ++k) _streamPut(fl, ' ');
		_streamPut(fl, val);
	}
	_streamPut(fl, '\n');

	//write rest of matrix
	for (unsigned int i = 0; i < VALUE_SIZE; ++i)
	{
		//write alphabet column cell
		TFileValue val = (TSequenceValue) i; //conversion integral => TSequenceValue => TFileValue
		_streamPut(fl, val);

		//write rest of line i
		unsigned int offset = i * VALUE_SIZE;
		for (unsigned int j = 0; j < VALUE_SIZE; ++j)
		{
			_sprintfValue(buf, tab[offset + j]);
			unsigned int len = strlen(buf);

			//leading blanks
			for (unsigned int k = 0; k < col_width - len; ++k) _streamPut(fl, ' ');

			//write cell
			for (unsigned int k = 0; k < len; ++k) _streamPut(fl, buf[k]);
//			_streamPut(fl, ',');
		}
		_streamPut(fl, '\n');
	}
}

//____________________________________________________________________________

template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec, typename TMeta>
inline void
write(TFile & fl,
	  Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc,
	  TMeta & meta)
{
	_writeScoringMatrix<TSequenceValue>(fl, sc.data_tab, meta);
}

template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
inline void
write(TFile & fl,
	  Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc)
{
	write(fl, sc, "");
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
