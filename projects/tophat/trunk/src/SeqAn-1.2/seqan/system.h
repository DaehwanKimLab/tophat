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
  $Id: system.h 3160 2008-12-14 22:48:33Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SYSTEM_H
#define SEQAN_HEADER_SYSTEM_H

//____________________________________________________________________________
// prerequisites

#include <seqan/file.h>

#include <cstdio>
#include <ctime>
#include <string>
#include <iostream>

#ifdef PLATFORM_WINDOWS
# include <windows.h>
#else //#ifdef PLATFORM_WINDOWS
# include <cstdlib>
# include <climits>
# include <pthread.h>
# include <errno.h>
# include <semaphore.h>
# include <aio.h>
# include <sys/mman.h>

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

#ifndef O_DIRECT
#define O_DIRECT 0
#endif

#endif //#ifdef PLATFORM_WINDOWS

#ifdef SEQAN_SWITCH_USE_FORWARDS
# include <seqan/system/system_manual_forwards.h>
# ifndef PLATFORM_WINDOWS
#  include <seqan/system/file_manual_forwards.h>
# endif
#endif

//____________________________________________________________________________
// multi-threading

#include <seqan/system/system_base.h>
#include <seqan/system/system_mutex.h>
#include <seqan/system/system_sema.h>
#include <seqan/system/system_event.h>
#include <seqan/system/system_thread.h>

//____________________________________________________________________________
// synchronous and asynchronous files

#include <seqan/system/file_sync.h>
#include <seqan/system/file_async.h>
#include <seqan/system/file_directory.h>

#endif //#ifndef SEQAN_HEADER_...
