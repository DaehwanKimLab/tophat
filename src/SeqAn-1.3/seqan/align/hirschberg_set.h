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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

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
