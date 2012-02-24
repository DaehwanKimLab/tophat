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

// THIS FILE IS CURRENTLY NOT USED IN SEQAN

#ifndef SEQAN_HEADER_BASIC_COUNTED_PTR_H
#define SEQAN_HEADER_BASIC_COUNTED_PTR_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// counted pointer

	template < typename Type >
	struct CountedPtr
	{
		typedef CountedPtr		Self_;
		typedef CountedPtr*	    SelfPtr_;
		typedef CountedPtr&	    SelfRef_;

		typedef Type&			reference;
		typedef const Type&		const_reference;
		typedef Type*			pointer;

        explicit CountedPtr(pointer p = 0): // allocate a new counter
            itsCounter(0)
        {
            if (p) itsCounter = new counter(p);
        }

        CountedPtr(const Self_& r) throw() {
            acquire(r.itsCounter);
        }

        ~CountedPtr() {
            release();
        }

        CountedPtr& operator=(const Self_& r)
        {
            if (this != &r) {
                release();
                acquire(r.itsCounter);
            }
            return *this;
        }

        reference operator*() const throw() {
            return *itsCounter->ptr;
        }

        pointer operator->() const throw() {
            return itsCounter->ptr;
        }

        pointer get() const throw() {
            return itsCounter ? itsCounter->ptr : 0;
        }

        bool unique() const throw() {
            return (itsCounter ? itsCounter->count == 1 : true);
        }

		inline operator pointer () const {
            return get();
		}

    private:

        struct counter {
            pointer     ptr;
            unsigned    count;
            counter(pointer p = 0, unsigned c = 1):
                ptr(p),
                count(c) { }
        }* itsCounter;

        void acquire(counter* c) throw()
        { // increment the count
            itsCounter = c;
            if (c) ++c->count;
        }

        void release()
        { // decrement the count, delete if it is 0
            if (itsCounter) {
                if (--itsCounter->count == 0) {
                    delete itsCounter->ptr;
                    delete itsCounter;
                }
                itsCounter = 0;
            }
        }
    };

}

#endif
