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
  $Id: hirschberg_set.h 3743 2009-03-23 15:50:31Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_HIRSCHBERG_SET_H
#define SEQAN_HEADER_HIRSCHBERG_SET_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

class _HirschbergSet
{

public:
	int x1,x2,y1,y2;
	int score;

public:
	_HirschbergSet()
		: x1(0),x2(0),y1(0),y2(0)
	{
SEQAN_CHECKPOINT
	}

	_HirschbergSet(int a1,int a2,int b1,int b2,int sc)
		: x1(a1),x2(a2),y1(b1),y2(b2),score(sc)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT(a1 <= b1);
		SEQAN_ASSERT(a2 <= b2);
	}

	_HirschbergSet & 
	operator = (_HirschbergSet const & other_)
	{
SEQAN_CHECKPOINT
		x1 = other_.x1;
		x2 = other_.x2;
		y1 = other_.y1;
		y2 = other_.y2;
		score = other_.score;
		return *this;
	}

};

////////////////////////////////////////////////////////////////////////////////////////////////
// Accessor Methods
////////////////////////////////////////////////////////////////////////////////////////////////

	// Sequence 1


inline 
int& 
_begin1(_HirschbergSet & me) {
	SEQAN_CHECKPOINT
	return me.x1;
}


inline 
int const& 
_begin1(_HirschbergSet const & me) {
	SEQAN_CHECKPOINT
	return me.x1;
}

inline 
void
_setBegin1(_HirschbergSet & me, int const & new_begin) {
	SEQAN_CHECKPOINT
	me.x1 = new_begin;
}

inline 
int&
_end1(_HirschbergSet & me) {
	SEQAN_CHECKPOINT
	return me.x2;
}

inline 
int const& 
_end1(_HirschbergSet const & me) {
	SEQAN_CHECKPOINT
	return me.x2;
}

inline 
void
_setEnd1(_HirschbergSet & me, int const & new_end) {
	SEQAN_CHECKPOINT
	me.x2 = new_end;
}

// Sequence 2

inline int&
_begin2(_HirschbergSet & me) {
	SEQAN_CHECKPOINT
	return me.y1;
}

inline 
int const&
_begin2(_HirschbergSet const & me) {
	SEQAN_CHECKPOINT
	return me.y1;
}

inline 
void
_setBegin2(_HirschbergSet & me, int const & new_begin) {
	SEQAN_CHECKPOINT
	me.y1 = new_begin;
}

inline 
int&
_end2(_HirschbergSet & me) {
	SEQAN_CHECKPOINT
	return me.y2;
}

inline 
int const&
_end2(_HirschbergSet const & me) {
	SEQAN_CHECKPOINT
	return me.y2;
}

inline 
void
_setEnd2(_HirschbergSet & me, int const & new_end) {
	SEQAN_CHECKPOINT
	me.y2 = new_end;
}

// //////////////////////////////////////////////////////////////////////////////////////////////
// Score Methods
// //////////////////////////////////////////////////////////////////////////////////////////////

inline 
int&
_score(_HirschbergSet & me)	{
	SEQAN_CHECKPOINT
	return me.score;
}

inline 
int const&
_score(_HirschbergSet const & me) {
	SEQAN_CHECKPOINT
	return me.score;
}

inline 
void
_setScore(_HirschbergSet & me,int new_score) {
	SEQAN_CHECKPOINT
	me.score = new_score;
}

// //////////////////////////////////////////////////////////////////////////////////////////////
//  Debug Methods
//		functions are only used for debugging or verbose output, therefore they
//      are only active in SEQAN_DEBUG
// //////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SEQAN_DEBUG
	
inline 
void
print(_HirschbergSet const & me) {
	std::cout << me.x1 << " " << me.x2 << "\t" << me.y1 << " " << me.y2 << std::endl;
}
#endif



inline bool
operator == (_HirschbergSet const & lhs, 
		_HirschbergSet const & rhs)
{
	return ((_begin1(lhs) == _begin1(rhs)) && (_end1(lhs) == _end1(rhs)) && (_begin2(lhs) == _begin2(rhs)) && (_end2(lhs) == _end2(rhs)));
}

}

#endif // #ifndef SEQAN_HEADER_HIRSCHBERG_SET_H
