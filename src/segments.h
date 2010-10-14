/*
 *  segments.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 2/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <vector>
#include <map>
#include "bwt_map.h"

enum eREAD
  {
    READ_DONTCARE = 0,
    READ_LEFT,
    READ_RIGHT
  };

struct RefSeg
{
RefSeg() : ref_id(0), points_left(false), antisense(false), read(READ_DONTCARE), left(0), right(0) {}
RefSeg(uint32_t i, bool p, bool antisense, eREAD read, int l, int r) :
  ref_id(i), points_left(p), antisense(antisense), read(read), left(l), right(r) {}
  
  bool operator<(const RefSeg& rhs) const
  {
    if (ref_id != rhs.ref_id)
      return ref_id < rhs.ref_id;
    if (left != rhs.left)
      return left < rhs.left;
    if (right != rhs.right)
      return right < rhs.right;
    return false;
  }
  
  bool operator==(const RefSeg& rhs) const
  {
    return (ref_id == rhs.ref_id &&
	    left == rhs.left &&
	    right == rhs.right && 
	    points_left == rhs.points_left &&
	    antisense == rhs.antisense &&
	    read == rhs.read);
  }
  
  uint32_t ref_id;
  bool points_left;
  bool antisense;
  eREAD read;
  int left;
  int right;
};	
