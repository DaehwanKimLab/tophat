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
  $Id: string_mmap.h 4782 2009-09-09 22:15:32Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STRING_MMAP_H
#define SEQAN_HEADER_STRING_MMAP_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    template < typename TConfig = ExternalConfig<> >
    struct MMap;
	
	
	//////////////////////////////////////////////////////////////////////////////
    // Memory Mapped String
    //////////////////////////////////////////////////////////////////////////////

    template < typename TValue,
               typename TConfig >
	class String<TValue, MMap<TConfig> >
	{
	public:

        typedef typename TConfig::TFile		TFile;
        typedef typename TConfig::TSize		TSize;

		TValue				*data_begin;
		TValue				*data_end;
		TSize				data_capacity;

		TFile				file;
		int					_openMode;
        bool                _temporary, _ownFile;

#ifdef PLATFORM_WINDOWS
        HANDLE				handle;
#endif

		String(TSize size = 0):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
            _temporary = true;
            _ownFile = false;

			resize(*this, size);
        }

		String(TFile &_file):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
			open(*this, _file);
        }

		String(const char *fileName, int openMode = DefaultOpenMode<TFile>::VALUE):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
			open(*this, fileName, openMode);
        }

		template <typename TSource>
		String & operator =(TSource const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}
		String & operator =(String const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}

		~String() 
		{
			close(*this);
		}

//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<String>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<String const>::Type 
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

//____________________________________________________________________________

        inline operator bool() 
   	{
            return file;
        }

//____________________________________________________________________________

};

   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
	begin(String<TValue, MMap<TConfig> > & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}
   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
	begin(String<TValue, MMap<TConfig> > const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_begin;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
   inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
	end(String<TValue, MMap<TConfig> > & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}
   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
	end(String<TValue, MMap<TConfig> > const & me,
		Standard)
	{
SEQAN_CHECKPOINT
		return me.data_end;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
	inline typename Size<String<TValue, MMap<TConfig> > >::Type
	capacity(String<TValue, MMap<TConfig> > & me) 
	{
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

   template < typename TValue, typename TConfig >
	inline typename Size<String<TValue, MMap<TConfig> > >::Type
	capacity(String<TValue, MMap<TConfig> > const & me) 
	{
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
	inline void 
	_setLength(
		String<TValue, MMap<TConfig> > & me, 
		size_t new_length)
	{
SEQAN_CHECKPOINT
		me.data_end = me.data_begin + new_length;
	}

//____________________________________________________________________________

	template < typename TValue, typename TConfig >
   inline void 
	_setCapacity(
		String<TValue, MMap<TConfig> > & me, 
		size_t new_capacity)
	{
SEQAN_CHECKPOINT
		me.data_capacity = new_capacity;
	}


    //////////////////////////////////////////////////////////////////////////////
    // meta-function interface

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, MMap<TConfig> > >
    {
        typedef size_t Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, MMap<TConfig> > >
    {
		typedef typename _MakeSigned<size_t>::Type Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct DefaultOverflowExplicit<String<TValue, MMap<TConfig> > >
	{
		typedef Generous Type;
	};

    template < typename TValue, typename TConfig >
	struct DefaultOverflowImplicit<String<TValue, MMap<TConfig> > >
	{
		typedef Generous Type;
	};
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct IsContiguous< String<TValue, MMap<TConfig> > >
	{
		typedef True Type;
		enum { VALUE = true };
	};

    template < typename TValue, typename TConfig >
	struct AllowsFastRandomAccess< String<TValue, MMap<TConfig> > >
	{
		typedef False Type;
		enum { VALUE = false };
	};


	//////////////////////////////////////////////////////////////////////////////
    // global interface

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	inline void 
	waitForAll(String<TValue, MMap<TConfig> > &)
	{
	}
	
#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES MMapStringDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, MMap<TConfig> > &) 
	{
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, MMap<TConfig> > &)
	{
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline void 
	flushAndFree(String<TValue, MMap<TConfig> > &)
	{
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, size_t new_capacity) 
	{
		if (new_capacity > 0) 
		{
			resize(me.file, new_capacity * sizeof(TValue));
			DWORD prot = 0;
			DWORD access = 0;
			if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
			{
				prot = PAGE_READONLY;
				access = FILE_MAP_READ;
			} else {
				prot = PAGE_READWRITE;
				access = FILE_MAP_ALL_ACCESS;
			}
            LARGE_INTEGER largeSize;
			largeSize.QuadPart = new_capacity;
			largeSize.QuadPart *= sizeof(TValue);

			me.handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, largeSize.HighPart, largeSize.LowPart, NULL);
			if (me.handle == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}

			void *addr = MapViewOfFile(me.handle, access, 0, 0, 0);	
			if (addr == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
		bool result = true;
		if (me.data_begin) 
		{
			if (UnmapViewOfFile(me.data_begin) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "UnmapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			if (CloseHandle(me.handle) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CloseHandle failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return result;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin) 
		{
			if (new_capacity > 0) 
			{
				TSize seq_length = length(me);
				
				DWORD prot = 0;
				DWORD access = 0;
				if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
				{
					prot = PAGE_READONLY;
					access = FILE_MAP_READ;
				} else {
					prot = PAGE_READWRITE;
					access = FILE_MAP_ALL_ACCESS;
				}
				LARGE_INTEGER largeSize;
				largeSize.QuadPart = new_capacity;
				largeSize.QuadPart *= sizeof(TValue);

				bool result = true;
				result &= (UnmapViewOfFile(me.data_begin) != 0);
				result &= (CloseHandle(me.handle) != 0);

				HANDLE handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, largeSize.HighPart, largeSize.LowPart, NULL);
				if (handle == NULL)
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
				#endif
					return false;
				}

				void *addr = MapViewOfFile(handle, access, 0, 0, 0);	
				if (addr == NULL)
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
				#endif
					return false;
				}

				if (capacity(me) > new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

				me.handle = handle;
				me.data_begin = (TValue*) addr;
				_setLength(me, seq_length);
				_setCapacity(me, new_capacity);
				return true;
			} else
				return _unmap(me);
		} else
			return _allocateStorage(me, new_capacity) != NULL;
	}

#else

    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, MMap<TConfig> > &me) 
	{
		msync(me.data_begin, length(me) * sizeof(TValue), MS_SYNC);
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, MMap<TConfig> > &me)
	{
		msync(me.data_begin, capacity(me) * sizeof(TValue), MS_INVALIDATE);
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline void 
	flushAndFree(String<TValue, MMap<TConfig> > &me)
	{
		madvise(me.data_begin, capacity(me) * sizeof(TValue), MADV_DONTNEED);
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, size_t new_capacity) 
	{
		if (new_capacity > 0) 
		{
			resize(me.file, new_capacity * sizeof(TValue));
			int prot = 0;
			if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
			if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
			void *addr = mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
			
			if (addr == MAP_FAILED)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "mmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
		if (me.data_begin) 
		{
			int error = munmap(me.data_begin, capacity(me) * sizeof(TValue));
			if (error != 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "munmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return true;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin)
		{
			if (new_capacity > 0) 
			{
				TSize seq_length = length(me);
				
				if (capacity(me) < new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

#ifdef MREMAP_MAYMOVE
				void *addr = mremap(me.data_begin, capacity(me) * sizeof(TValue), new_capacity * sizeof(TValue), MREMAP_MAYMOVE);
#else
				// for BSD systems without mremap(..) like Mac OS X ...
				int prot = 0;
				if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
				if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
	//			void *addr = mmap(me.data_begin, new_capacity * sizeof(TValue), prot, MAP_SHARED | MAP_FIXED, me.file.handle, 0);
				munmap(me.data_begin, capacity(me) * sizeof(TValue));
				void *addr = mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
#endif

				if (addr == MAP_FAILED) 
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "mremap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
				#endif
					return false;
				}

				if (capacity(me) > new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

				me.data_begin = (TValue*) addr;
				_setLength(me, seq_length);
				_setCapacity(me, new_capacity);
				return true;
			} else
				return _unmap(me);
		} else
			return _map(me, new_capacity);
	}

#endif

	template < typename TValue, typename TConfig >
    inline void 
    clear(String<TValue, MMap<TConfig> > &me) 
	{
		cancel(me);
		_unmap(me);
		resize(me.file, 0);
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig, typename TSize >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type * 
    _allocateStorage(String<TValue, MMap<TConfig> > &me, TSize new_capacity) 
	{
		TSize size = _computeSize4Capacity(me, new_capacity);
		_map(me, size);
		return NULL;
	}

	template < typename TValue, typename TConfig >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type * 
    _reallocateStorage(
		String<TValue, MMap<TConfig> > &me,
		typename Size< String<TValue, MMap<TConfig> > >::Type new_capacity) 
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;
		TSize size = _computeSize4Capacity(me, new_capacity);
		_remap(me, size);
		return NULL;
	}

	template < typename TValue, typename TConfig, typename TSize >
    inline void
    _deallocateStorage(String<TValue, MMap<TConfig> > &/*me*/, TValue * /*ptr*/, TSize /*capacity*/)
	{
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName, int openMode) 
	{
		close(me);
		me._temporary = false;
				
		if ((me._ownFile = open(me.file, fileName, openMode))) 
		{
			me._openMode = openMode;
			return _map(me, (size_t)(size(me.file) / sizeof(TValue)));
		}

		return false;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName) 
	{
		typedef typename String<TValue, MMap<TConfig> >::TFile	TFile;
		return open(me, fileName, DefaultOpenMode<TFile>::VALUE);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, typename TConfig::TFile file) 
	{
		close(me);
		me.file = file;
        me._temporary = false;
        me._ownFile = false;

		if (me.file) 
		{
			me._openMode = OPEN_RDWR;
			return _map(me, size(me.file) / sizeof(TValue));
		}

		return false;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    openTemp(String<TValue, MMap<TConfig> > &me) 
	{
		close(me);
        me._temporary = true;
		me._openMode = OPEN_RDWR;
		
		return me._ownFile = openTemp(me.file);
    }
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
	inline void _ensureFileIsOpen(String<TValue, MMap<TConfig> > &me) 
	{
		if (!me.file)
		{
			me._temporary = true;
			if (!(me._ownFile = openTemp(me.file)))
				::std::cerr << "Memory Mapped String couldn't open temporary file" << ::std::endl;
		}
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, const char *fileName, int openMode) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, const char *fileName) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, typename TConfig::TFile file) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    close(String<TValue, MMap<TConfig> > &me) 
	{
		if (me.file) 
		{
			// close associated file
			if (me._temporary) 
			{
				me._temporary = false;
				cancel(me);
			}

			_unmap(me);

			if (me._ownFile) 
			{
				me._ownFile = false;
				return close(me.file);
			}
		}
		return true;
    }

//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
