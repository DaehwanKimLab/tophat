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

#ifndef SEQAN_HEADER_FILE_CSTYLE_H
#define SEQAN_HEADER_FILE_CSTYLE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/*    template <>
    struct Value< FILE* >
    {
	    typedef unsigned char Type;
    };
*/
/*  already defined in "seqan/cstream.h" ...

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
/*
    template <>
    struct Size< FILE* >
    {
	    typedef size_t Type;
    };
*/
/*

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
    template <>
    struct Difference< FILE* >
    {
	    typedef long Type;
    };


	inline const char * 
	_getCStyleOpenMode(int openMode) 
	{
		switch (openMode & OPEN_MASK) {
            case OPEN_WRONLY:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w";
                    else
                        return "r+";
                else
                    return "a";
            case OPEN_RDWR:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w+";
                    else
                        return "r+";
                else
                    return "a+";
			default:
		        return "r";
		}
    }

    inline bool 
	open(FILE* &me, const char *fileName, int openMode) 
	{
		SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
        return (me = fopen(fileName, _getCStyleOpenMode(openMode))) != NULL;
    }

    inline bool 
	open(FILE* &me, const char *fileName) 
	{
		return open(me, fileName, DefaultOpenMode<FILE*>::VALUE);
	}

    inline bool 
	openTemp(FILE* &me) 
	{
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return (me = tmpfile()) != NULL;
    }

    inline bool 
	close(FILE* me) 
	{
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return fclose(me) == 0;
    }

    inline unsigned 
	sectorSize(FILE* const &) 
	{
        return 4096;
    }

    template < typename TPos >
    inline Size<FILE*>::Type 
	seek(FILE* me, TPos const fileOfs, int origin) 
	{
        fseek(me, fileOfs, origin);
		return ftell(me);
    }
    template < typename TPos >
    inline Size<FILE*>::Type 
	seek(FILE* me, TPos const fileOfs) 
	{
		return seek(me, fileOfs, SEEK_BEGIN);
    }

    inline Size<FILE*>::Type 
	tell(FILE* me) 
	{
		return ftell(me);
    }

    template < typename TValue, typename TSize >
    inline bool 
	read(FILE* me, TValue *memPtr, TSize const count) 
	{
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize >
    inline bool 
	write(FILE* me, TValue const *memPtr, TSize const count) 
	{
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize, typename TPos >
    inline bool 
	readAt(FILE* me, TValue *memPtr, TSize const count, TPos const fileOfs) 
	{
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }
    
    template < typename TValue, typename TSize, typename TPos >
    inline bool 
	writeAt(FILE* me, TValue const *memPtr, TSize const count, TPos const fileOfs) 
	{
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == (size_t)count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    inline Size<FILE*>::Type 
	size(FILE* me) 
	{
        Size<FILE*>::Type old_pos = tell(me);
        Size<FILE*>::Type result = 0;
        if (seek(me, 0, SEEK_END) == 0)
            result = tell(me);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

    template < typename TSize >
    inline void 
	resize(FILE* me, TSize new_length) 
	{
        Size<FILE*>::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        seek(me, old_pos, SEEK_BEGIN);
    }

	inline bool 
	flush(FILE*) 
	{
		return true; 
	}

    template < typename AsyncRequest >
	inline void 
	release(FILE*, AsyncRequest &) 
	{
	}

    template < typename AsyncRequest >
    inline bool 
	cancel(FILE*, AsyncRequest &) 
	{
		return true; 
	}

}

#endif
