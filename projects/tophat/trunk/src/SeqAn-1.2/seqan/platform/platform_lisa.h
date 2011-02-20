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
  $Id$
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef PLATFORM_GCC
  #define PLATFORM_GCC
#endif

// should be set before including anything
#define _FILE_OFFSET_BITS 64
#include <unistd.h>

#define finline __inline__

// default 64bit type
typedef __int64_t __int64;
typedef __uint64_t __uint64;

//define SEQAN_SWITCH_USE_FORWARDS t generated forwards 
#define SEQAN_SWITCH_USE_FORWARDS

#include "platform_generated_forwards.h"

#ifndef SEQAN_HEADER_PLATFORM_GENERATED_FORWARDS_H
#error To use the SeqAn library you first have to execute 'make forwards' in the root directory
#endif
