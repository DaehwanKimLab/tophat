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
  $Id: file_page.h 3080 2008-11-19 09:33:16Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_PAGE_H
#define SEQAN_HEADER_FILE_PAGE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


	//////////////////////////////////////////////////////////////////////////////
	// base class for memory buffers

	template < typename TValue >
	struct SimpleBuffer 
	{
		typedef	typename Size<SimpleBuffer>::Type		TSize;
		typedef	typename Iterator<SimpleBuffer>::Type	TIterator;

		TIterator	begin;      // the beginning of the buffer
        TIterator	end;        // end of valid data
        TSize		pageSize;   // size of allocated memory

        SimpleBuffer():
            begin(NULL),
            end(NULL) {}

        SimpleBuffer(TIterator _begin, TIterator _end):
            begin(_begin),
            end(_end) {}

        SimpleBuffer(TIterator _begin, TSize _size):
            begin(_begin),
            end(_begin + _size) {}

        SimpleBuffer(TSize _pageSize):
            begin(NULL),
            end(NULL),
            pageSize(_pageSize) {}

        inline TValue       & operator[](TSize i)       { return begin[i]; }
        inline TValue const & operator[](TSize i) const { return begin[i]; }
	};


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template < typename TValue >
    struct Value< SimpleBuffer<TValue> >	{ typedef TValue Type; };

    template < typename TValue >
	struct Size< SimpleBuffer<TValue> >		{ typedef size_t Type; };

    template < typename TValue >
    struct Iterator< SimpleBuffer<TValue>, Standard >		{ typedef TValue *Type; };

    template < typename TValue >
    struct Iterator< SimpleBuffer<TValue> const, Standard >	{ typedef TValue const *Type; };


	//////////////////////////////////////////////////////////////////////////////
	// global interface

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    pageSize(SimpleBuffer<TValue> &me) {
        return me.pageSize;
    }

    template < typename TValue, typename TSize >
    inline void setPageSize(SimpleBuffer<TValue> &me, TSize size) {
        me.pageSize = size;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    size(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    length(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue, typename TSize >
    inline void resize(SimpleBuffer<TValue> &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TSize, typename T >
	inline void allocPage(SimpleBuffer<TValue> &pf, TSize size, T const & me) {
        setPageSize(pf, size);
        allocate(me, pf.begin, pageSize(pf));
        resize(pf, size);
	}

	template < typename TValue, typename T > inline
	void freePage(SimpleBuffer<TValue> &pf, T const & me) {
		deallocate(me, pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
        setPageSize(pf, 0);
	}

	template < typename TValue >
	inline TValue* begin(SimpleBuffer<TValue> &pf, Standard) {
		return pf.begin;
	}

	template < typename TValue >
	inline TValue const * begin(SimpleBuffer<TValue> const &pf, Standard) {
		return pf.begin;
	}

	template < typename TValue >
	inline TValue * end(SimpleBuffer<TValue> &pf, Standard) {
		return pf.end;
	}

	template < typename TValue >
	inline TValue const * end(SimpleBuffer<TValue> const &pf, Standard) {
		return pf.end;
	}



    //////////////////////////////////////////////////////////////////////////////
	// a bucket is a structure to represent a small window of a page
    // used by algorithms which need a global view of all pages (merge sort, mapper)

	template < typename TValue >
    struct PageBucket
	{
		typedef	typename Iterator<PageBucket>::Type TIterator;

        unsigned    pageOfs;                // begin of bucket window with relation to page begin
        TIterator	begin, cur, end;        // begin/end of buckets memory buffer and a pointer
    };

    template < typename TValue >
    struct PageBucketExtended : public PageBucket< TValue > {
		int     	pageNo;		            // related page (needed by merger sort)
    };

	template < typename TValue >
    ::std::ostream& operator<<(::std::ostream &out, const PageBucketExtended<TValue> &pb) {
        for(TValue *cur = pb.begin; cur != pb.end; cur++)
            out << *cur << " ";
        return out;
    }


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template < typename TValue >
    struct Value< PageBucket<TValue> >		{ typedef TValue Type; };

    template < typename TValue >
    struct Size< PageBucket<TValue> >		{ typedef size_t Type; };

    template < typename TValue >
    struct Iterator< PageBucket<TValue>, Standard >			{ typedef TValue *Type; };

    template < typename TValue >
    struct Iterator< PageBucket<TValue> const, Standard >	{ typedef TValue const *Type; };


	//////////////////////////////////////////////////////////////////////////////

	template < typename TValue, typename TFile, typename TSpec >
    struct PageFrame;    

	//////////////////////////////////////////////////////////////////////////////

	template < typename TValue, typename TFile, typename TSpec >
    struct Value< PageFrame< TValue, TFile, TSpec > > {
        typedef TValue Type;
    };

	template < typename TValue, typename TFile, typename TSpec >
    struct Size< PageFrame< TValue, TFile, TSpec > > {
        typedef size_t Type;
    };

	template < typename TValue, typename TFile, typename TSpec >
    struct Iterator< PageFrame< TValue, TFile, TSpec >, Standard > {
        typedef TValue *Type;
    };

	template < typename TValue, typename TFile, typename TSpec >
    struct Iterator< PageFrame< TValue, TFile, TSpec > const, Standard> {
        typedef TValue const *Type;
    };


	//////////////////////////////////////////////////////////////////////////////
	// page frame of dynamic size

    template < typename TSpec = void >
    struct Dynamic;

    template < typename TValue, typename TFile >
	struct PageFrame< TValue, TFile, Dynamic<> >: public SimpleBuffer< TValue >
	{
		typedef TFile							File;
		typedef SimpleBuffer<TValue>	        TBase;
        typedef typename aRequest<TFile>::Type  aRequest;

        enum Status		{ READY, READING, WRITING };

		bool			dirty;		// data needs to be written to disk before freeing
		unsigned   		pageNo;		// maps frames to pages (reverse vector mapper)
        aRequest        request;    // request structure of the async io process
		Status			status;
        PageFrame       *next;      // next buffer in a chained list

        PageFrame():
			TBase(),
            dirty(false),
            pageNo(-1),
			status(READY),
			next(NULL) {}
    };


	//////////////////////////////////////////////////////////////////////////////
	// page frame of static size

    template < unsigned _PageSize >
    struct Fixed;

    typedef ::std::list<int>		PageLRUList;    // least recently usage list
	typedef PageLRUList::iterator	PageLRUEntry;

    template < typename TValue,
               typename TFile,
               unsigned _PageSize >
	struct PageFrame<TValue, TFile, Fixed<_PageSize> >
	{
		typedef TFile								File;
        typedef typename aRequest<TFile>::Type		aRequest;
		typedef	typename Size<PageFrame>::Type		TSize;
		typedef	typename Iterator<PageFrame>::Type	TIterator;

		enum            { PageSize = _PageSize };
		enum Status		{ READY, READING, WRITING };
		enum DataStatus	{ ON_DISK = -1, UNINITIALIZED = -2 };
		enum Priority	{ NORMAL_LEVEL = 0, PREFETCH_LEVEL = 1, ITERATOR_LEVEL = 2, PERMANENT_LEVEL = 3 };

		bool			dirty;		// data needs to be written to disk before freeing
		int     		pageNo;		// maps frames to pages (reverse vector mapper)
		TIterator		begin;	    // start address of page memory
        aRequest        request;    // request structure of the async io process
		Status			status;
		DataStatus		dataStatus;
		PageLRUEntry	lruEntry;   // priority based lru
        Priority        priority;

		PageFrame():
            dirty(false),
            pageNo(-1),
			begin(NULL),
			status(READY),
			dataStatus(UNINITIALIZED),
            priority(NORMAL_LEVEL) {}

		inline TValue       & operator[](TSize i)       { return begin[i]; }
        inline TValue const & operator[](TSize i) const { return begin[i]; }
	};


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template < typename TValue, typename TFile, unsigned _PageSize >
    struct Iterator< PageFrame< TValue, TFile, Fixed<_PageSize> >, Standard >
    {
        typedef VolatilePtr<TValue> Type;
    };

	template < typename TValue, typename TFile, unsigned _PageSize >
    struct Iterator< PageFrame< TValue, TFile, Fixed<_PageSize> > const, Standard>
    {
        typedef VolatilePtr<TValue> const Type;
    };


	//////////////////////////////////////////////////////////////////////////////
	// global interface

	template < typename TValue, typename TFile, typename TSpec, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Dynamic<TSpec> > &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    size(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    length(PageFrame<TValue, TFile, Fixed<_PageSize> > const &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    pageSize(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/, TSize /*size*/) {}


	//////////////////////////////////////////////////////////////////////////////
	// various page frame methods

	template < typename TValue, typename TFile, typename TSpec >
    ::std::ostream& operator<<(::std::ostream &out, const PageFrame<TValue, TFile, TSpec > &pf) 
	{
        out << "PageFrame @ " << pf.pageNo;
        if (pf.dirty)
            out << " DIRTY";
        else
            out << " CLEAN";

        switch (pf.status) {
			case PageFrame<TValue, TFile, TSpec >::READY:
                out << " READY";
                break;
			case PageFrame<TValue, TFile, TSpec >::READING:
                out << " READING";
                break;
			case PageFrame<TValue, TFile, TSpec >::WRITING:
                out << " WRITING";
        }

        if (pf.dataStatus == pf.ON_DISK)
            out << " ON_DISK";
        else
            out << " UNITIALIZED";

        out << " Prio:" << pf.priority;
        out << " Buffer:" << (TValue*)pf.begin;

        return out;
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void allocPage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) 
	{
		TValue* tmp = NULL;
		allocate(me, tmp, pageSize(pf));
		pf.begin = tmp;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "allocPage: " << ::std::hex << tmp << ::std::dec << ::std::endl;
		#endif
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void freePage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) 
	{
		#ifdef SEQAN_VVERBOSE
			if ((TValue*)pf.begin)
				::std::cerr << "freePage:  " << ::std::hex << (TValue*)pf.begin << ::std::dec << ::std::endl;
		#endif
        nukeCopies(pf.begin);
		deallocate(me, (TValue*)pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
//        resize(pf, pageSize(pf));
		return areadAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.status = pf.WRITING;
//        resize(pf, pageSize(pf));
		return awriteAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize> inline
    bool readLastPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file, TSize size) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return readAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize > inline
	bool writeLastPage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file, TSize size) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf) 
	{
		if ((pf.status != pf.READY) && waitFor(pf.request)) 
		{
			pf.status = pf.READY;
			pf.dirty = false;
            return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec, typename TTime > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf, TTime timeOut) 
	{
		if ((pf.status != pf.READY) && waitFor(pf.request, timeOut)) {
			pf.status = pf.READY;
			pf.dirty = false;
			return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool cancel(PageFrame<TValue, TFile, TSpec> &pf, TFile &file) 
	{
        waitFor(pf, 0);
		if (pf.status != pf.READY) 
		{
            if (!cancel(file, pf.request)) return false;
            pf.status = pf.READY;
        }
        return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page based read/write methods used by Pool classes

    template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) 
	{
        if (size(pf) == pageSize(pf))
            return readPage(pf.pageNo, pf, file);
        else
            return readLastPage(pf.pageNo, pf, file, size(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) 
	{
        if (size(pf) == pageSize(pf))
            return writePage(pf, pf.pageNo, file);
        else
            return writeLastPage(pf, pf.pageNo, file, size(pf));
	}

	template < typename TValue, typename TFile > inline
	unsigned readBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, unsigned dataSize, TFile &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = _min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readBucket:  " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)pageNo * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << readSize << ::std::endl;
		#endif
        if (readSize && readAt(file, b.begin, readSize, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, typename TFile > inline
	bool writeBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, TFile &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)pageNo * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << b.cur - b.begin << ::std::endl;
		#endif
        if ((b.cur == b.begin) || writeAt(file, b.begin, b.cur - b.begin, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writeBucket(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, unsigned &pageOfs, TFile &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << pf.begin;
			::std::cerr << " from page " << ::std::dec << pf.pageNo << " at " << (pos_t)pf.pageNo * (pos_t)pageSize(pf) + pageOfs;
			::std::cerr << " size " << size(pf) << ::std::endl;
		#endif
        if (pf.end == pf.begin) return true;
        if (awriteAt(file, pf.begin, size(pf), (pos_t)pf.pageNo * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}


	//////////////////////////////////////////////////////////////////////////////

    template < typename TPageFrame >
    struct PageChain 
	{
        TPageFrame          *first, *last;
        unsigned            frames, maxFrames;
        
        PageChain(unsigned _maxFrames = UINT_MAX):
            first(NULL),
            last(NULL),
            frames(0),
            maxFrames(_maxFrames)
        {
            for(unsigned i = 0; i < _maxFrames; ++i)
                pushBack();
        }
        
        ~PageChain()
        {
            while (first)
                popFront();
        }
        
        inline TPageFrame & operator[](int k) 
		{
            TPageFrame *p = first;
            for (; k > 0; --k)
                p = p->next;
            return *p;
        }
        
        inline TPageFrame const & operator[](int k) const 
		{
            TPageFrame *p = first;
            for (; k > 0; --k)
                p = p->next;
            return *p;
        }

        inline TPageFrame * getReadyPage() 
		{
            if (!first || (frames < maxFrames && waitFor(*first, 0)))
                return pushBack();
            else {
                waitFor(*first);
                return firstToEnd();
            }
        }

        template < typename TFile >
        inline void cancelAll(TFile &file) 
		{
            TPageFrame *p = first;
            for (; p != NULL; p = p->next)
                cancel(*p, file);
        }

        inline void waitForAll() 
		{
            TPageFrame *p = first;
            for (; p != NULL; p = p->next)
                waitFor(*p);
        }

    private:

        inline TPageFrame * firstToEnd() 
		{
            last->next = first;
            last = first;
            first = first->next;
            last->next = NULL;
            return last;
        }

        inline TPageFrame * pushBack() 
		{
            TPageFrame * p = new TPageFrame();
            if (p) {
                if (last)
                    last->next = p;
                else
                    first = p;
                last = p;
                ++frames;
            }
            return p;
        }

        inline TPageFrame * popFront() 
		{
            TPageFrame *p = first;
            if (p) {
                first = first->next;
                if (!first) last = NULL;
                --frames;
                delete p;
            }
            return p;
        }
    };

	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

    template < typename TPageFrame >    
	struct Value< PageChain<TPageFrame> > 
	{
		typedef TPageFrame Type;
	};

    template < typename TPageFrame >    
	struct Size< PageChain<TPageFrame> > 
	{
		typedef unsigned Type;
	};

    template < typename TPageFrame >    
	struct Iterator< PageChain<TPageFrame> > 
	{
		typedef TPageFrame *Type;
	};

    template < typename TPageFrame >    
	struct Iterator< PageChain<TPageFrame> const > 
	{
		typedef TPageFrame const *Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// page container with lru mechanism
	// the lru strategy uses different priorities
	// the page with the least priority is used first
	// 0..random access pages
	// 1..forward iterator pages
	// 2..quasi permanent pages

    template < typename TPageFrame,
               unsigned FRAMES,
               unsigned PRIORITY_LEVELS = TPageFrame::PERMANENT_LEVEL + 1 >
	struct PageContainer
	{
		typedef String<TPageFrame>					TPages;
		typedef typename Position<TPages>::Type		TPos;

		enum { PriorityLevels = PRIORITY_LEVELS };

		TPages			pages;
        PageLRUList		*lruList;

		PageContainer()
		{
			lruList = new PageLRUList[PRIORITY_LEVELS];
			resize(pages, FRAMES, Exact());
            for(TPos i = 0; i < FRAMES; ++i)
			    pages[i].lruEntry = lruList[0].insert(lruList[0].end(), i);
        }

		~PageContainer()
		{
			delete[] lruList;
		}

        inline TPageFrame       & operator[](TPos i)       { return pages[i]; }
        inline TPageFrame const & operator[](TPos i) const { return pages[i]; }

		inline void push_back() 
		{
			TPos last = endPosition(pages);
			resize(pages, last + 1);
			pages[last].lruEntry = lruList[0].insert(lruList[0].end(), last);
		}

		inline void erase(int frameNo) 
		{
			lruList[pages[frameNo].priority].erase(pages[frameNo].lruEntry);
            erase(pages, frameNo);
		}

        inline void rename(int frameNo) 
		{
            *(pages[frameNo].lruEntry) = frameNo;
        }

		inline void pop_back() 
		{
			lruList[back(pages).priority].erase(back(pages).lruEntry);
            erase(pages, endPosition(pages) - 1);
		}


		//////////////////////////////////////////////////////////////////////////////
		// lru strategy interface

		inline void upgrade(const TPageFrame &pf) 
		{
			lruList[pf.priority].splice(lruList[pf.priority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].begin();
		}

		inline void downgrade(const TPageFrame &pf) 
		{
			lruList[pf.priority].splice(lruList[pf.priority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].end();
			--pf.lruEntry;
		}

		inline void upgrade(TPageFrame &pf, int newPriority) 
		{
			lruList[newPriority].splice(lruList[newPriority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].begin();
			pf.priority = static_cast<typename TPageFrame::Priority> (newPriority);
		}

		inline void downgrade(TPageFrame &pf, int newPriority) 
		{
			lruList[newPriority].splice(lruList[newPriority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].end();
			--pf.lruEntry;
			pf.priority = static_cast<typename TPageFrame::Priority> (newPriority);
		}

		inline void _dump() 
		{
			for(unsigned i = 0; i < PRIORITY_LEVELS; ++i) {
                ::std::cerr << "|";
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
                    TPageFrame &pf = pages[*I];
                    ::std::cerr << pf.pageNo;
                    if (pf.dirty) ::std::cerr << "*";
                    else          ::std::cerr << " ";
                    if (pf.status == TPageFrame::READY) ::std::cerr << "  ";
                    else                                ::std::cerr << ". ";
				};
            }
            ::std::cerr << ::std::endl;
		}

        // Function is a functor which is called with a PageFrame object,
        // that is dirty or not READY (in an IO transfer)
		template <class Function>
		inline int mru(Function _Func, unsigned maxLevel = PRIORITY_LEVELS - 1) 
		{
			for(unsigned i = 0; i <= maxLevel; ++i) {
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
					TPageFrame& pf = pages[*I];
					if (pf.status == TPageFrame::READY && !pf.dirty)
						return *I;
					else
						if (_Func(pf)) return *I;
				};
            }
			#ifdef SEQAN_VVERBOSE
				::std::cerr << "ALL PAGES DIRTY OR IN USE (try to use const iterators) :-(" << ::std::endl;
			#endif
			return -1;
		}

		inline int mruDirty() 
		{
            for(unsigned i = 0; i < PRIORITY_LEVELS; ++i)
                if (!lruList[i].empty())
                    return lruList[i].back();
			return -1;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// meta-function interface

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> > 
	{
		typedef String<TPageFrame> Type;
	};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const > 
	{
		typedef String<TPageFrame> const Type;
	};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Value< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> > 
	{
		typedef TPageFrame Type;
	};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Size< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Size< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Position< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Position< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >:
		public Iterator< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type> {};

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	struct Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const >:
		public Iterator< typename Host< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const>::Type> {};


	//////////////////////////////////////////////////////////////////////////////
	// global interface

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS, typename TExpand >
	inline void reserve(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		unsigned _Count,
		Tag<TExpand> const expand)
	{
		reserve(pageCont.pages, _Count, expand);
	}

    template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline void resize(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont, 
		unsigned _Count) 
	{
		unsigned _Size = length(pageCont.pages);
		if (_Size < _Count) {
			reserve(pageCont.pages, _Count);
            for(unsigned i = _Size; i < _Count; ++i)
                pageCont.push_back();
		} else 
			if (_Size > _Count)
				for(unsigned i = _Count; i < _Size; ++i)
					pageCont.pop_back();
	}

    template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline typename Size< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> >::Type
	length(PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont) 
	{
		return length(pageCont.pages);
	}

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS>, Standard >::Type
	begin(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		Standard const) 
	{
		return begin(pageCont.pages, Standard());
	}

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const, Standard >::Type
	begin(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont,
		Standard const) 
	{
		return begin(pageCont.pages, Standard());
	}

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS>, Standard >::Type
	end(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> &pageCont,
		Standard const) 
	{
		return end(pageCont.pages, Standard());
	}

	template < typename TPageFrame, unsigned FRAMES, unsigned PRIORITY_LEVELS >
	inline typename Iterator< PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const, Standard >::Type
	end(
		PageContainer<TPageFrame, FRAMES, PRIORITY_LEVELS> const &pageCont,
		Standard const) 
	{
		return end(pageCont.pages, Standard());
	}




	//////////////////////////////////////////////////////////////////////////////

    template < typename TValue, typename TSize, typename T, class Function >
    inline bool equiDistantDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
			::std::cerr << "equiDistantDistribution: _pages is null!" << ::std::endl;
            return false;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        allocPage(_clusterBuffer, _bufferSize, me);
        PageBucketExtended<TValue> pb;
        pb.begin = _clusterBuffer.begin;

        unsigned clusterSize = _bufferSize / pages;
        if (lastPageSize > 0 && clusterSize >= lastPageSize) {
            // last page bucket would get more memory than page would need
            // --> exclude from equi-size distribution
            if (--pages) {
                _bufferSize -= lastPageSize;
                clusterSize = _bufferSize / pages;
            }
        }

        if (pages) {
            unsigned remainder = _bufferSize % pages;
            for(unsigned i = 0, numerator = 0; i < pages; ++i) {
                pb.end = pb.begin + clusterSize;
                if ((numerator += remainder) >= pages) {    // simple bresenham for distribution
                    numerator -= pages;
                    ++pb.end;
                }
                pb.cur = pb.begin;
                pb.pageOfs = 0;
			    _Func(pb);
                pb.begin = pb.end;
            }
        }

        if (pages < _pages) {
            pb.end = pb.begin + lastPageSize;
            pb.cur = pb.begin;
            pb.pageOfs = 0;
			_Func(pb);
        }

        return true;
    }

    template < typename TValue, typename TSize, typename T, class Function >
    inline unsigned equiDistantAlignedDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned aligning, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
			::std::cerr << "equiDistantDistribution: _pages is null!" << ::std::endl;
            return 0;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantAlignedDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        unsigned clusterSize = _bufferSize / pages;
        unsigned aclusterSize = (clusterSize / aligning) * aligning;
        if (clusterSize - aclusterSize > aligning / 2)
            aclusterSize += aligning;

		if (aclusterSize != 0) {

			if (lastPageSize > 0 && aclusterSize > lastPageSize) {
				// last page bucket would get more memory than page would need
				// --> exclude from equi-size distribution
				--pages;
				allocPage(_clusterBuffer, aclusterSize * pages + lastPageSize, me);
			} else
				allocPage(_clusterBuffer, aclusterSize * pages, me);

			PageBucketExtended<TValue> pb;
			pb.begin = _clusterBuffer.begin;

			if (pages) {
				for(unsigned i = 0; i < pages; ++i) {
					pb.end = pb.begin + aclusterSize;
					pb.cur = pb.begin;
					pb.pageOfs = 0;
					_Func(pb);
					pb.begin = pb.end;
				}
			}

			if (pages < _pages) {
				pb.end = pb.begin + lastPageSize;
				pb.cur = pb.begin;
				pb.pageOfs = 0;
				_Func(pb);
			}
		}

        return aclusterSize;
    }

}

#endif
