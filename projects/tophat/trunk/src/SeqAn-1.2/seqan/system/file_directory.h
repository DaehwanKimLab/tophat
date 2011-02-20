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

#ifndef SEQAN_HEADER_FILE_DIRECTORY_H
#define SEQAN_HEADER_FILE_DIRECTORY_H

#ifdef PLATFORM_WINDOWS
# include <io.h>
#else
# include <dirent.h>
#endif


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    class Directory
    {
	protected:
	
		intptr_t			handle;
		struct _finddata_t	entry;
		bool				_atEnd;

		friend inline bool open(Directory &dir, char const *dirName);
		friend inline bool close(Directory &dir);
		friend inline char const * value(Directory &dir);
		friend inline char const * value(Directory const &dir);
		friend inline Directory & goNext(Directory &dir);
		friend inline bool atEnd(Directory &dir);
		friend inline bool atEnd(Directory const &dir);

	public:

		Directory()
		{
			handle = 0;
			_atEnd = true;
		}
		
		Directory(char const *dirName)
		{
			open(*this, dirName);
		}
		
		~Directory()
		{
			close(*this);
		}
		
		inline char const * operator* () const
		{
			return value(*this);
		}

		inline Directory & operator++ ()
		{
			return goNext(*this);
		}
		
		inline operator bool () const
		{
			return !_atEnd;
		}
	};

//////////////////////////////////////////////////////////////////////////////	

	inline bool
	open(Directory &dir, char const *dirName)
	{
		CharString selection = dirName;
		append(selection, "\\*");
		dir._atEnd = ((dir.handle = _findfirst(toCString(selection), &dir.entry)) == -1L);
		return !dir._atEnd;
	}

	inline bool
	close(Directory &dir)
	{
		int result = 0;
		if (dir.handle)
			result = _findclose(dir.handle);

		dir._atEnd = true;
		dir.handle = 0;
		return result == 0;
	}

	inline char const *
	value(Directory &dir)
	{
		return dir.entry.name;
	}

	inline char const *
	value(Directory const &dir)
	{
		return dir.entry.name;
	}

	inline Directory &
	goNext(Directory &dir)
	{
		dir._atEnd = (_findnext(dir.handle, &dir.entry) != 0);
		return dir;
	}

	inline bool
	atEnd(Directory &dir)
	{
		return dir._atEnd;
	}
	
	inline bool
	atEnd(Directory const &dir)
	{
		return dir._atEnd;
	}
	
//////////////////////////////////////////////////////////////////////////////	
	
	
#else


//////////////////////////////////////////////////////////////////////////////	

    class Directory
    {
	protected:
	
		DIR		*handle;
		dirent	*it;
		
		friend inline bool open(Directory &dir, char const *dirName);
		friend inline bool close(Directory &dir);
		friend inline char const * value(Directory &dir);
		friend inline char const * value(Directory const &dir);
		friend inline Directory & goBegin(Directory &dir);
		friend inline Directory & goNext(Directory &dir);
		friend inline bool atEnd(Directory &dir);
		friend inline bool atEnd(Directory const &dir);

	public:

		Directory()
		{
			handle = NULL;
			it = NULL;
		}
		
		Directory(char const *dirName)
		{
			open(*this, dirName);
		}
		
		~Directory()
		{
			close(*this);
		}
		
		inline char const * operator* () const
		{
			return value(*this);
		}

		inline Directory & operator++ ()
		{
			return goNext(*this);
		}
		
		inline operator bool () const
		{
			return !atEnd(*this);
		}
	};
	
//////////////////////////////////////////////////////////////////////////////	
	
	inline bool
	open(Directory &dir, char const *dirName)
	{
		if ((dir.handle = opendir(dirName)) != NULL)
		{
			goNext(dir);
			return true;
		}
		dir.it = NULL;
		return false;
	}

	inline bool
	close(Directory &dir)
	{
		int result = 0;
		if (dir.handle != NULL)
			result = closedir(dir.handle);

		dir.handle = NULL;
		dir.it = NULL;
		return result == 0;
	}

	inline char const *
	value(Directory &dir)
	{
		return dir.it->d_name;
	}

	inline char const *
	value(Directory const &dir)
	{
		return dir.it->d_name;
	}

	inline Directory &
	goBegin(Directory &dir)
	{
		rewinddir(dir.handle);
		return goNext(dir);
	}

	inline Directory &
	goNext(Directory &dir)
	{
		dir.it = readdir(dir.handle);
		return dir;
	}

	inline bool
	atEnd(Directory &dir)
	{
		return dir.it == NULL;
	}
	
	inline bool
	atEnd(Directory const &dir)
	{
		return dir.it == NULL;
	}
	
#endif

}

#endif
