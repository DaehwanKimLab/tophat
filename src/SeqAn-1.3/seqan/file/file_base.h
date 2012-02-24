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

#ifndef SEQAN_HEADER_FILE_BASE_H
#define SEQAN_HEADER_FILE_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	// To override the system's default temporary directory use the following:
	//#define SEQAN_DEFAULT_TMPDIR "/var/tmp"

	// To use direct I/O access define SEQAN_DIRECTIO (not completely tested yet)
	//#define SEQAN_DIRECTIO


/**
.Spec.Sync:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous input/output access.
..signature:File<Sync<> >
..remarks:This class suports pseudo-asynchronous access methods, i.e. the methods to initiate a I/O request return after request completion.
..include:seqan/file.h
*/

	template <typename TSpec = void>
    struct Sync;

/**
.Spec.Async:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous and asynchronous input/output access.
..signature:File<Async<> >
..include:seqan/file.h
*/

	template <typename TSpec = void>
    struct Async;


/**
.Class.File:
..cat:Input/Output
..summary:Represents a file.
..signature:File<TSpec>
..param.TSpec:The specializing type.
...default:$Async<>$, see @Spec.Async@.
..include:seqan/file.h
*/

	template <typename TSpec = Async<> >
    class File;

/**
.Spec.Chained:
..cat:Files
..general:Class.File
..summary:Splits a large file into a chain of smaller files.
..signature:File<Chained<FileSize, TFile> >
..param.FileSize:The maximal split file size in byte.
...default:2^31-1 (~2GB)
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a chain of $TFile$ files, whose file sizes are at most $FileSize$ bytes.
Chained Files should be used for file systems or $TFile$ types that don't support large files (e.g. FAT32, C-style FILE*).
..remarks:The chain can be used as if it were one contiguous file.
..include:seqan/file.h
*/

	// chained file's default filesize is 2gb-1byte (fat16 filesize limitation)
	template < __int64 FileSize_ = ~(((__int64)1) << 63), typename TFile = File<> >
	struct Chained;

/**
.Spec.Striped:
..cat:Files
..general:Class.File
..summary:Stripes a file across multiple files.
..signature:File<Chained<FileCount, TFile> >
..param.FileCount:The number of files used for striping.
...default:2
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a software striping without redundance (see RAID0) to accelerate I/O access when using more than one disks.
..remarks:Striped files should only be used in @Class.Pool@s or external Strings as they only support block operations and no random accesses.
..include:seqan/file.h
*/

	template < unsigned FileCount_ = 2, typename TFile = File<> >
	struct Striped;

    enum FileOpenMode {
        OPEN_RDONLY     = 1,
        OPEN_WRONLY     = 2,
        OPEN_RDWR       = 3,
        OPEN_MASK       = 3,
        OPEN_CREATE     = 4,
        OPEN_APPEND     = 8,
        OPEN_ASYNC      = 16,
		OPEN_TEMPORARY	= 32,
		OPEN_QUIET		= 128
    };

	template <typename T>
	struct DefaultOpenMode {
		enum { VALUE = OPEN_RDWR | OPEN_CREATE | OPEN_APPEND };
	};

	template <typename T>
	struct DefaultOpenTempMode {
		enum { VALUE = OPEN_RDWR | OPEN_CREATE };
	};

    enum FileSeekMode {
        SEEK_BEGIN   = 0,
        SEEK_CURRENT = 1
#ifndef SEEK_END
      , SEEK_END     = 2
#endif
    };


    //////////////////////////////////////////////////////////////////////////////
    // result type of asynch. functions
    // you have to call release(AsyncRequest<T>) after a finished *event based* transfer
	struct AsyncDummyRequest {};

/**
.Class.AsyncRequest:
..cat:Input/Output
..summary:Associated with an asynchronous I/O request.
..signature:AsyncRequest<TFile>
..param.TFile:A File type.
..remarks:This structure is used to identify asynchronous requests after their initiation.
..include:seqan/file.h
*/

    template < typename T >
    struct AsyncRequest
    {
        typedef AsyncDummyRequest Type;
    };
/*
    //////////////////////////////////////////////////////////////////////////////
    // event to represent asynchronous transfers
    // you can wait for it or test it
    template < typename T >
    struct aEvent
    {
        typedef DummyEvent Type;
    };

    ////////////////////////////////////////////////////////////////////////////////
    // callback hint parameter type
    // hint lets you recognize the finished asynch. transfer in your own callback routine
    template < typename T >
    struct aHint
    {
        typedef void Type;
    };

    //////////////////////////////////////////////////////////////////////////////
    // callback function interface
    template < typename T >
    struct aCallback
    {
        typedef void Type(aHint<T> *);
    };

    //////////////////////////////////////////////////////////////////////////////
    // file queue interface
    template < typename T >
    struct aQueue
    {
        typedef Nothing Type;
    };
*/

    //////////////////////////////////////////////////////////////////////////////
    // generic open/close interface

/**
.Function.open:
..summary:Opens a file.
..cat:Input/Output
..signature:open(file, fileName[, openMode])
..param.file:A File object.
...type:Class.File
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName, int openMode) 
	{
        return me.open(fileName, openMode);
    }

    template < typename TSpec >
    inline bool open(File<TSpec> &me, const char *fileName) 
	{
		return open(me, fileName, DefaultOpenMode<File<TSpec> >::VALUE);
    }

/**
.Function.openTemp:
..summary:Opens a temporary file.
..cat:Input/Output
..signature:openTemp(file)
..param.file:A File object.
...type:Class.File
..remarks:After closing this file will automatically be deleted.
..remarks:The openmode (see @Function.open@) is $OPEN_RDWR | OPEN_CREATE$.
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me) 
	{
        return me.openTemp();
    }

    template < typename TSpec >
    inline bool openTemp(File<TSpec> &me, int openMode) 
	{
        return me.openTemp(openMode);
    }

    template < typename File >
    inline void reopen(File &, int) 
	{
	}
    
/**
.Function.close:
..cat:Input/Output
..summary:Closes a file.
..signature:close(file)
..param.file:A File object.
...type:Class.File
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline bool close(File<TSpec> & me) 
	{
        return me.close();
    }

    template < typename TSpec >
    inline unsigned sectorSize(File<TSpec> const & /*me*/) 
	{
        return 4096;
    }


    //////////////////////////////////////////////////////////////////////////////
    // generic read(At)/write(At) interface

/**
.Function.read:
..cat:Input/Output
..summary:Loads records from a file.
..signature:read(file, memPtr, count)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..returns:A $bool$ which is $true$ on success.
..remarks:The records are read from the position pointed by the current file pointer (see @Function.seek@).
..include:seqan/file.h
*/

	template < typename TSpec, typename TValue, typename TSize >
    inline bool read(File<TSpec> & me, TValue *memPtr, TSize const count) 
	{
		return me.read(memPtr, count * sizeof(TValue));
    }
    
/**
.Function.write:
..cat:Input/Output
..summary:Saves records to a file.
..signature:write(file, memPtr, count)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..returns:A $bool$ which is $true$ on success.
..remarks:The records are written at the position pointed by the current file pointer (see @Function.seek@).
..include:seqan/file.h
*/

	template < typename TSpec, typename TValue, typename TSize >
    inline bool write(File<TSpec> & me, TValue const *memPtr, TSize const count) 
	{
		return me.write(memPtr, count * sizeof(TValue));
    }

/**
.Function.readAt:
..summary:Loads records from a specific position in a file.
..cat:Input/Output
..signature:readAt(file, memPtr, count, fileOfs)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool readAt(TFile & me, TValue *memPtr, TSize const count, TPos const fileOfs) 
	{
		typedef typename Position<TFile>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
		return read(me, memPtr, count);
    }
    
/**
.Function.writeAt:
..summary:Saves records to a specific position in a file.
..cat:Input/Output
..signature:writeAt(file, memPtr, count, fileOfs)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TFile, typename TValue, typename TSize, typename TPos >
    inline bool writeAt(TFile & me, TValue const *memPtr, TSize const count, TPos const fileOfs) 
	{
		typedef typename Position<TFile>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
		return write(me, memPtr, count);
    }



    //////////////////////////////////////////////////////////////////////////////
    // generic seek/tell/size/resize interface

/**
.Function.seek:
..summary:Changes the current file pointer.
..cat:Input/Output
..signature:seek(file, fileOfs[, origin])
..param.file:A File object.
...type:Class.File
..param.fileOfs:A file offset measured in bytes relative to $origin$.
..param.origin:Selects the origin from where to calculate the new position.
...default:$SEEK_BEGIN$
...remarks:For $SEEK_BEGIN$, $SEEK_CURRENT$, or $SEEK_END$ the origin is the beginning, the current pointer, or the end of the file.
..returns:The new file position measured in bytes from the beginning.
..include:seqan/file.h
*/

	template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs, int origin) 
	{
		typedef typename Position< File<TSpec> >::Type TFilePos;
		TFilePos newOfs = me.seek(fileOfs, origin);
        #ifdef SEQAN_DEBUG_OR_TEST_
			if (origin == SEEK_BEGIN && newOfs != (TFilePos)fileOfs) {
				::std::cerr << "seek returned " << ::std::hex << newOfs << " instead of " << fileOfs << ::std::dec << ::std::endl;
			}
        #endif
        return newOfs;
    }
    
	template < typename TSpec, typename TPos >
    inline typename Position< File<TSpec> >::Type seek(File<TSpec> &me, TPos const fileOfs) 
	{
		return seek(me, fileOfs, SEEK_BEGIN);
	}
/**
.Function.tell:
..summary:Gets the current file pointer.
..cat:Input/Output
..signature:tell(file)
..param.file:A File object.
...type:Class.File
..returns:The current file position measured in bytes from the beginning.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline typename Position< File<TSpec> >::Type tell(File<TSpec> &me) 
	{
        return me.tell();
    }

/**
.Function.rewind:
..summary:Sets the current file pointer to the beginning.
..cat:Input/Output
..signature:rewind(file)
..param.file:A File object.
...type:Class.File
..remarks:Calls @Function.seek@$(file, 0)$ by default.
..include:seqan/file.h
*/

    template < typename File >
    inline void rewind(File &me) 
	{
		seek(me, 0);
    }
    
/**
.Function.size:
..summary:Gets the file size.
..cat:Input/Output
..signature:size(file)
..param.file:A File object.
...type:Class.File
..returns:The file size measured in bytes.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline typename Size<File<TSpec> >::Type size(File<TSpec> &me) 
	{
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        typename Size<File<TSpec> >::Type result = seek(me, 0, SEEK_END);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

/**
.Function.resize:
..cat:Input/Output
..signature:resize(file, new_length)
..param.file:A File object.
...type:Class.File
..param.new_length:The new file size measured in bytes.
..include:seqan/file.h
*/

    template < typename TSpec, typename TSize >
    inline void resize(File<TSpec> &me, TSize new_length) 
	{
        typename Size<File<TSpec> >::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        setEof(me);
        seek(me, old_pos, SEEK_BEGIN);
    }

/**
.Function.setEof:
..summary:Sets the file end to the current pointer.
..cat:Input/Output
..signature:setEof(file)
..param.file:A File object.
...type:Class.File
..include:seqan/file.h
*/

    template < typename TSpec >
    inline bool setEof(File<TSpec> &/*me*/) 
	{ 
		return true; 
	}


    //////////////////////////////////////////////////////////////////////
    // Pseudo asynchronous Methods
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // callback based read/write
/*
    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncRead(File & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        result = read(me, memPtr, count);
        cb(hint);
        return NULL;
    }
    
    template < typename File, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncWrite(File & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        write(me, memPtr, count);
        cb(hint);
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        readAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }
    
    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename AsyncRequest<File>::Type
    asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        result = writeAt(me, memPtr, count, fileOfs);
        cb(hint);
        return NULL;
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncRead(File & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        read(me, memPtr, count);
        event.signal();
        return NULL;
    }
    
    template < typename File, typename TValue, typename TSize,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncWrite(File & me, TValue const *memPtr, TSize const count,
        aEvent &event)
    {
        write(me, memPtr, count);
        event.signal();
        return NULL;
    }

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        readAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }
    
    template < typename File, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename AsyncRequest<File>::Type
    asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        writeAt(me, memPtr, count, fileOfs);
        event.signal();
        return NULL;
    }
*/

    //////////////////////////////////////////////////////////////////////
    // queue-less request based pseudo asychronous read/write

/**
.Function.asyncReadAt:
..summary:Asynchronously loads records from a specific position in a file.
..cat:Input/Output
..signature:asyncReadAt(file, memPtr, count, fileOfs, request)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first destination record in memory.
..param.count:The amount of records to be read.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..param.request:Reference to a structure that will be associated with this asynchronous request.
...type:Class.AsyncRequest
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename AsyncRequest >
    inline bool 
	asyncReadAt(File & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        AsyncRequest &)
    {
        return readAt(me, memPtr, count, fileOfs);
    }
    
/**
.Function.asyncWriteAt:
..summary:Asynchronously saves records to a specific position in a file.
..cat:Input/Output
..signature:asyncWriteAt(file, memPtr, count, fileOfs, request)
..param.file:A File object.
...type:Class.File
..param.memPtr:A pointer to the first source record in memory.
..param.count:The amount of records to be written.
..param.fileOfs:The absolute file position in bytes measured from the beginning.
..param.request:Reference to a structure that will be associated with this asynchronous request.
...type:Class.AsyncRequest
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename File, typename TValue, typename TSize, typename TPos,
               typename AsyncRequest >
    inline bool
	asyncWriteAt(File & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        AsyncRequest &)
    {
        return writeAt(me, memPtr, count, fileOfs);
    }

	
	//////////////////////////////////////////////////////////////////////
    // pseudo queue specific functions

/**
.Function.flush:
..summary:Waits for all open requests to complete.
..cat:Input/Output
..signature:flush(file)
..param.file:A File object.
...type:Class.File
..remarks:$flush$ returns after all pending requests are completed.
..include:seqan/file.h
*/

    template < typename TSpec >
    inline void flush(File<TSpec> &) 
	{
	}

/**
.Function.waitFor:
..summary:Waits for an asynchronous request to complete.
..cat:Input/Output
..signature:waitFor(request[, timeout_millis])
..param.request:Reference to an AsyncRequest object.
...type:Class.AsyncRequest
..param.timeout_millis:Timout value in milliseconds.
...remarks:A value of 0 can be used to test for completion without waiting.
...default:Infinity.
..returns:A $bool$ which is $true$ on completion and $false$ on timeout.
..remarks:$waitFor$ suspends the calling process until $request$ is completed or after $timeout_millis$ milliseconds.
..include:seqan/file.h
*/

    inline bool waitFor(AsyncDummyRequest &) 
	{ 
		return true; 
	}

	template < typename TTime >
    inline bool waitFor(AsyncDummyRequest &, TTime) 
	{ 
		return true; 
	}

	// deprecated
	template < typename TSpec, typename AsyncRequest >
    inline void release(File<TSpec> &, AsyncRequest &) 
	{
	}

/**
.Function.cancel:
..summary:Cancels an asynchronous request.
..cat:Input/Output
..signature:cancel(file, request)
..param.file:A File object.
...type:Class.File
..param.request:Reference to an AsyncRequest object.
...type:Class.AsyncRequest
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

    template < typename TSpec, typename AsyncRequest >
    inline bool cancel(File<TSpec> &, AsyncRequest &) 
	{
		return true; 
	}


	// little helpers

	template <typename T1, typename T2> inline
	T1 enclosingBlocks(T1 _size, T2 _blockSize) 
	{
		return (_size + _blockSize - 1) / _blockSize;
	}

	template <typename T1, typename T2> inline
	T1 alignSize(T1 _size, T2 _aligning) 
	{
        if (_size < _aligning)
            return _aligning;
        else
		    return (_size / _aligning) * (T1)_aligning;
	}

}

#endif
