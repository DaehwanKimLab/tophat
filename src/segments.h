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

enum ePOINT_DIR
  {
    POINT_DIR_DONTCARE = 0,
    POINT_DIR_LEFT,
    POINT_DIR_RIGHT,
    POINT_DIR_BOTH
  };

struct RefSeg
{
RefSeg() :
  ref_id(0),
    points_where(POINT_DIR_DONTCARE),
    antisense(false),
    read(READ_DONTCARE),
    left(0),
    right(0),
    support_read("")
  {}
  
RefSeg(uint32_t i, ePOINT_DIR p, bool antisense, eREAD read, int l, int r, const string& support_read = "") :
  ref_id(i),
    points_where(p),
    antisense(antisense),
    read(read),
    left(l),
    right(r),
    support_read(support_read)
  {}
  
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
	    points_where == rhs.points_where &&
	    antisense == rhs.antisense &&
	    read == rhs.read &&
	    support_read == rhs.support_read);
  }
  
  uint32_t ref_id;
  ePOINT_DIR points_where;
  bool antisense;
  eREAD read;
  int left;
  int right;
  
  string support_read;
};	
