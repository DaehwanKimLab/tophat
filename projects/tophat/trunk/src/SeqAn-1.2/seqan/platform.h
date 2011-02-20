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

#ifndef SEQAN_PLATFORM_H
#define SEQAN_PLATFORM_H

#ifdef __MINGW32__
#error MINGW
	#include "platform/platform_mingw.h"
#elif _MSC_VER
	#include "platform/platform_windows.h"
#elif __SUNPRO_C
	#include "platform/platform_solaris.h"
#else
	#include "platform/platform_gcc.h"
#endif

#endif
