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

#ifndef SEQAN_HEADER_BASIC_AGGREGATES_H
#define SEQAN_HEADER_BASIC_AGGREGATES_H

namespace SEQAN_NAMESPACE_MAIN
{

//____________________________________________________________________________

    struct Compressed_;
	typedef Tag<Compressed_> Compressed;

	template <unsigned BITSIZE1 = 16, unsigned BITSIZE2 = 16>
	struct BitCompressed;

/**
.Class.Pair:
..cat:Aggregates
..summary:Stores two arbitrary objects.
..signature:Pair<T1[, T2[, TSpec]]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
...default:$T1$
..param.TSpec:The specializing type.
...default:$void$, no compression (faster access).
.Memfunc.Pair#Pair:
..class:Class.Pair
..summary:Constructor
..signature:Pair<T1, T2[, TSpec]> ()	
..signature:Pair<T1, T2[, TSpec]> (pair)
..signature:Pair<T1, T2[, TSpec]> (i1, i2)
..param.pair:Other Pair object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
.Memvar.Pair#i1:
..class:Class.Pair
..summary:T1 object
.Memvar.Pair#i2:
..class:Class.Pair
..summary:T2 object
..include:seqan/basic.h
*/

	// standard storage 
	template <typename T1_, typename T2_ = T1_, typename TSpec = void>
    struct Pair {
        typedef T1_ T1;
        typedef T2_ T2;
	    T1_ i1;
	    T2_ i2;
		inline Pair() {}
		inline Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		inline Pair(T1_ const &_i1, T2_ const &_i2): i1(_i1), i2(_i2) {}

		template <typename T1__, typename T2__, typename TSpec__>
		inline Pair(Pair<T1__, T2__, TSpec__> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)) {}
    };

/**
.Spec.Packed Pair:
..cat:Aggregates
..general:Class.Pair
..summary:Stores two arbitrary objects. Saves memory by disabling memory alignment.
..signature:Pair<T1, T2, Compressed>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..notes:Useful for external storage.
..remarks:Memory access could be slower. Direct access to members by pointers is not allowed on all platforms.
..include:seqan/basic.h
.Memfunc.Pair#Pair.class:Spec.Packed Pair
.Memvar.Pair#i1.class:Spec.Packed Pair
.Memvar.Pair#i2.class:Spec.Packed Pair
*/

	// unaligned and unpadded storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename T1_, typename T2_>
    struct Pair<T1_, T2_, Compressed> {
        typedef T1_ T1;
        typedef T2_ T2;
        T1_ i1;
        T2_ i2;
		inline Pair() {}
		inline Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		inline Pair(T1_ const &_i1, T2_ const &_i2): i1(_i1), i2(_i2) {}

		template <typename T1__, typename T2__, typename TSpec__>
		inline Pair(Pair<T1__, T2__, TSpec__> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)) {}
	}
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif


/**
.Spec.Bit Compressed Pair:
..cat:Aggregates
..general:Class.Pair
..summary:Stores two arbitrary objects. Saves memory by packing bits with bit fields.
..signature:Pair<T1, T2, BitCompressed<BITSIZE1, BITSIZE2> >
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.BITSIZE1:Number of bits to store $T1$.
..param.BITSIZE2:Number of bits to store $T2$.
..notes:Useful for external storage.
..remarks:Memory access could be slower. Direct access to members by pointers is not allowed.
..include:seqan/basic.h
.Memfunc.Pair#Pair.class:Spec.Bit Compressed Pair
.Memvar.Pair#i1.class:Spec.Bit Compressed Pair
.Memvar.Pair#i2.class:Spec.Bit Compressed Pair
*/

    template <typename T1_, typename T2_, unsigned BITSIZE1, unsigned BITSIZE2>
    struct Pair<T1_, T2_, BitCompressed<BITSIZE1, BITSIZE2> > {
        typedef T1_ T1;
        typedef T2_ T2;
	    T1_ i1:BITSIZE1;
	    T2_ i2:BITSIZE2;
		inline Pair() {}
		inline Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		inline Pair(T1_ const &_i1, T2_ const &_i2): i1(_i1), i2(_i2) {}

		template <typename T1__, typename T2__, typename TSpec__>
		inline Pair(Pair<T1__, T2__, TSpec__> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)) {}
	};



    template <typename T1_, typename T2_, typename TSpec>
	std::ostream& operator<<(std::ostream &out, Pair<T1_,T2_,TSpec> const &p) {
		out << "< " << getValueI1(p) << " , " << getValueI2(p) << " >";
		return out;
	}

	template <typename T1, typename T2, typename TSpec>
	struct Value< Pair<T1, T2, TSpec>, 1 > {
		typedef T1 Type;
	};

	template <typename T1, typename T2, typename TSpec>
	struct Value< Pair<T1, T2, TSpec>, 2 > {
		typedef T2 Type;
	};

	template <typename T1, typename T2, typename TSpec>
	struct Spec< Pair<T1, T2, TSpec> > {
		typedef TSpec Type;
	};


//____________________________________________________________________________

	template <typename TKey, typename TObject, typename TSpec>
	struct Key< Pair<TKey, TObject, TSpec> > 
	{
		typedef TKey Type;
	};

	template <typename TKey, typename TCargo, typename TSpec>
	struct Cargo< Pair<TKey, TCargo, TSpec> > 
	{
		typedef TCargo Type;
	};
//____________________________________________________________________________

/**
.Class.Triple:
..cat:Aggregates
..summary:Stores three arbitrary objects.
..signature:Triple<T1[, T2[, T3[, TSpec]]]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
...default:$T1$
..param.T3:The type of the third object.
...default:$T2$
..param.TSpec:The specializing type.
...default:$void$, no compression (faster access).
.Memfunc.Triple#Triple:
..class:Class.Triple
..summary:Constructor
..signature:Triple<T1, T2, T3[, TSpec]> ()
..signature:Triple<T1, T2, T3[, TSpec]> (triple)
..signature:Triple<T1, T2, T3[, TSpec]> (i1, i2, i3)
..param.triple:Other Triple object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
..param.i3:T3 object.
.Memvar.Triple#i1:
..class:Class.Triple
..summary:T1 object
.Memvar.Triple#i2:
..class:Class.Triple
..summary:T2 object
.Memvar.Triple#i3:
..class:Class.Triple
..summary:T3 object
..include:seqan/basic.h
*/

	// standard storage 
	template <typename T1_, typename T2_ = T1_, typename T3_ = T1_, typename TSpec = void>
    struct Triple {
        typedef T1_ T1;
        typedef T2_ T2;
        typedef T3_ T3;
        T1_ i1;
        T2_ i2;
        T3_ i3;
		inline Triple() {}
		inline Triple(Triple const &_p):
			i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		inline Triple(T1_ const &_i1, T2_ const &_i2, T3_ const &_i3):
			i1(_i1), i2(_i2), i3(_i3) {}

		template <typename T1__, typename T2__, typename T3__, typename TSpec__>
		inline Triple(Triple<T1__, T2__, T3__, TSpec__> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}

        inline bool
        operator==(Triple const & other) const
        {
            return i1 == other.i1 && i2 == other.i2 && i3 == other.i3;
        }
        
        inline bool
        operator<(Triple const & other) const
        {
            if (i1 < other.i1)
                return true;
            if (i1 == other.i1 && i2 < other.i2)
                return true;
            if (i1 == other.i1 && i2 == other.i2 && i3 < other.i3)
                return true;
            return false;
        }
	};
	
/**
.Spec.Packed Triple:
..cat:Aggregates
..general:Class.Triple
..summary:Stores three arbitrary objects. Saves memory by disabling memory alignment.
..signature:Triple<T1, T2, T3, Compressed>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.T3:The type of the third object.
..notes:Useful for external storage.
..remarks:Memory access could be slower. Direct access to members by pointers is not allowed on all platforms.
..include:seqan/basic.h
.Memfunc.Triple#Triple.class:Spec.Packed Triple
.Memvar.Triple#i1.class:Spec.Packed Triple
.Memvar.Triple#i2.class:Spec.Packed Triple
.Memvar.Triple#i3.class:Spec.Packed Triple
*/

	// unaligned and unpadded storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename T1_, typename T2_, typename T3_>
    struct Triple<T1_, T2_, T3_, Compressed> {
        typedef T1_ T1;
        typedef T2_ T2;
        typedef T3_ T3;
        T1_ i1;
        T2_ i2;
        T3_ i3;
		inline Triple() {}
		inline Triple(Triple const &_p):
			i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		inline Triple(T1_ const &_i1, T2_ const &_i2, T3_ const &_i3):
			i1(_i1), i2(_i2), i3(_i3) {}

		template <typename T1__, typename T2__, typename T3__, typename TSpec__>
		inline Triple(Triple<T1__, T2__, T3__, TSpec__> const &_p):
			i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}
	}
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif

	template <typename T1_, typename T2_, typename T3_, typename TSpec>
	std::ostream& operator<<(std::ostream &out, Triple<T1_,T2_,T3_,TSpec> const &t) {
		out << "< " << getValueI1(t) << " , " << getValueI2(t) << " , " << getValueI3(t) << " >";
		return out;
	}

	template <typename T1, typename T2, typename T3, typename TSpec>
	struct Value< Triple<T1, T2, T3, TSpec>, 1 > {
		typedef T1 Type;
	};

	template <typename T1, typename T2, typename T3, typename TSpec>
	struct Value< Triple<T1, T2, T3, TSpec>, 2 > {
		typedef T2 Type;
	};

	template <typename T1, typename T2, typename T3, typename TSpec>
	struct Value< Triple<T1, T2, T3, TSpec>, 3 > {
		typedef T3 Type;
	};

	template <typename T1, typename T2, typename T3, typename TSpec>
	struct Spec< Triple<T1, T2, T3, TSpec> > {
		typedef TSpec Type;
	};


//____________________________________________________________________________

/**
.Class.Tuple:
..cat:Aggregates
..summary:A plain fixed-length string.
..signature:Tuple<T, SIZE[, TSpec]>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.SIZE:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..param.TSpec:The specializing type.
...default:$void$, no compression (faster access).
..include:seqan/basic.h
*/

	// standard storage 
	template <typename T_, unsigned _size, typename TSpec = void>
    struct Tuple {
        typedef T_ T;
        enum { size = _size };
        T_ i[_size];

		template <typename TPos>
        inline T_& operator[](TPos k) {
            SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
			SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(size));
            return i[k];
        }
		template <typename TPos>
        inline const T_& operator[](TPos k) const {
            SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
			SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(size));
            return i[k];
        }
		inline T_* operator&() { return i; }
		inline const T_* operator&() const { return i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename TPos, typename tmpS>
		inline tmpS const assignValueAt(TPos k, tmpS const source) {
			return i[k] = source;
		}
    };


    template < unsigned char _size >
	struct BitVector_ {
        typedef typename BitVector_<_size + 1>::Type Type;
    };

    template <> struct BitVector_<8> { typedef unsigned char Type; };
    template <> struct BitVector_<16> { typedef unsigned short Type; };
    template <> struct BitVector_<32> { typedef unsigned long Type; };
    template <> struct BitVector_<64> { typedef __uint64 Type; };
    template <> struct BitVector_<255> { typedef __uint64 Type; };

/**
.Spec.Bit Packed Tuple:
..cat:Aggregates
..general:Class.Tuple
..summary:A plain fixed-length string. Saves memory by packing bits.
..signature:Tuple<T, SIZE, Compressed>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.SIZE:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..notes:The characters are stored as a bit sequence in an ordinal type (char, ..., __int64).
..remarks:Only useful for small alphabets and small tuple sizes (|Sigma|^size <= 2^64) as for @Spec.Dna@ or @Spec.AminoAcid@ m-grams)
..see:Spec.Sampler
..include:seqan/basic.h
*/

	// bit-compressed storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
    template <typename T_, unsigned _size>
    struct Tuple<T_, _size, Compressed> {
        typedef T_ T;
        enum { size = _size };
        enum { bitSize = BitsPerValue<T_>::VALUE };
        enum { bitMask = (1 << bitSize) - 1 };
        enum { mask = (1 << (size * bitSize)) - 1 };
        typedef typename BitVector_< bitSize * size >::Type CT;
        
        CT i;
/*
		inline Tuple() {
			SEQAN_ASSERT(bitSize * size <= sizeof(CT) * 8);
		}
*/
		template <typename TPos>
        inline const T_ operator[](TPos k) const {
            SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
			SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(size));
            return (i >> (size - 1 - k) * bitSize) & bitMask;
        }
		template <unsigned size__>
		inline Tuple operator=(Tuple<T_, size__, Compressed> const &_right) {
			i = _right.i;
			return *this;
		}
		template <typename TShiftSize>
        inline CT operator<<=(TShiftSize shift) {
            return i = (i << (shift * bitSize)) & mask;
        }
		template <typename TShiftSize>
        inline CT operator<<(TShiftSize shift) const {
            return (i << (shift * bitSize)) & mask;
        }
		template <typename TShiftSize>
        inline CT operator>>=(TShiftSize shift) {
            return i = (i >> (shift * bitSize));
        }
		template <typename TShiftSize>
        inline CT operator>>(TShiftSize shift) const {
            return i >> (shift * bitSize);
        }
        template <typename T>
        inline void operator|=(T const &t) {
            i |= t;
        }
        template <typename T, typename TSpec>
        inline void operator|=(SimpleType<T, TSpec> const &t) {
            i |= t.value;
        }
		inline CT* operator&() { return &i; }
		inline const CT* operator&() const { return &i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename TPos, typename tmpS>
		inline tmpS const assignValueAt(TPos k, tmpS const source) {
			typedef Tuple<T_, _size, Compressed> Tup;
			typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * bitSize);
			i = (i & ~mask) | ((CT)ordValue(source) << ((_size - 1 - k) * bitSize));
			return source;
		}
    }
#ifndef PLATFORM_WINDOWS
	__attribute__((packed))
#endif
	;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif


//////////////////////////////////////////////////////////////////////////////
// length

    template <typename T_, unsigned _size, typename TSpec>
	inline unsigned length(Tuple<T_, _size, TSpec> const &) { return _size; }

	///.Metafunction.LENGTH.param.T.type:Class.Tuple
    template <typename T_, unsigned _size, typename TSpec>
	struct LENGTH< Tuple<T_, _size, TSpec> >
	{
		enum { VALUE = _size };
	};

//////////////////////////////////////////////////////////////////////////////
// assignValueAt

    template <typename TObject, typename TPos, typename TSource>
    inline TSource & 
	assignValueAt(TObject &me, TPos k, TSource &source) {
        assign(value(me, k), source);
		return source;
    }

    template <typename TObject, typename TPos, typename TSource>
    inline TSource const & 
	assignValueAt(TObject &me, TPos k, TSource const &source) {
        assign(value(me, k), source);
		return source;
    }

    template <typename T_, unsigned _size, typename tmpS, typename TPos>
    inline tmpS const assignValueAt(Tuple<T_, _size, void> &me, TPos k, tmpS const source) {
        return me.i[k] = source;
    }

    template <typename T_, unsigned _size, typename tmpS, typename TPos>
    inline tmpS const assignValueAt(Tuple<T_, _size, Compressed> &me, TPos k, tmpS const source) {
        typedef Tuple<T_, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source << ((_size - 1 - k) * me.bitSize);
        return source;
    }

    template <typename T_, typename tmpS, typename Spec_, unsigned _size, typename TPos>
    inline SimpleType<tmpS, Spec_> const & assignValueAt(Tuple<T_, _size, Compressed> &me, TPos k, SimpleType<tmpS, Spec_> const &source) {
        typedef Tuple<T_, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source.value << ((_size - 1 - k) * me.bitSize);
        return source;
    }

//////////////////////////////////////////////////////////////////////////////
// clear

	template <typename T_, unsigned _size, typename TSpec>
	inline void clear(Tuple<T_, _size, TSpec> &me) {
        memset<sizeof(me.i), 0>(&(me.i));
	}
    template <typename T_, unsigned _size>
	inline void clear(Tuple<T_, _size, Compressed> &me) {
		me.i = 0; 
	}

//////////////////////////////////////////////////////////////////////////////
// optimized compares

	template <typename T_, unsigned _sizeL, unsigned _sizeR>
	inline bool operator<(Tuple<T_, _sizeL, Compressed> const &_left, Tuple<T_, _sizeR, Compressed> const &_right) {
		return _left.i < _right.i;
	}
	template <typename T_, unsigned _sizeL, unsigned _sizeR>
	inline bool operator>(Tuple<T_, _sizeL, Compressed> const &_left, Tuple<T_, _sizeR, Compressed> const &_right) {
		return _left.i > _right.i;
	}
	template <typename T_, unsigned _sizeL, unsigned _sizeR>
	inline bool operator==(Tuple<T_, _sizeL, Compressed> const &_left, Tuple<T_, _sizeR, Compressed> const &_right) {
		return _left.i == _right.i;
	}
	template <typename T_, unsigned _sizeL, unsigned _sizeR>
	inline bool operator!=(Tuple<T_, _sizeL, Compressed> const &_left, Tuple<T_, _sizeR, Compressed> const &_right) {
		return _left.i != _right.i;
	}

//////////////////////////////////////////////////////////////////////////////
// optimized shifts

    struct TupleShiftLeftWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg[I-1] = arg[I];
        }
    };

    struct TupleShiftRightWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg[I] = arg[I-1];
        }
    };

	template <typename T_, unsigned _size, typename TSpec>
	inline void shiftLeft(Tuple<T_, _size, TSpec> &me) {
		Loop<TupleShiftLeftWorker_, _size - 1>::run(me);
	}

	template <typename T_, unsigned _size, typename TSpec>
	inline void shiftRight(Tuple<T_, _size, TSpec> &me) {
		LoopReverse<TupleShiftRightWorker_, _size - 1>::run(me);
	}

	template <typename T_, unsigned _size>
	inline void shiftLeft(Tuple<T_, _size, Compressed> &me) {
		me<<=1;
	}

	template <typename T_, unsigned _size>
	inline void shiftRight(Tuple<T_, _size, Compressed> &me) {
		me>>=1;
	}

//////////////////////////////////////////////////////////////////////////////
// standard output

	template <typename T_, unsigned _size, typename TSpec>
	std::ostream& operator<<(std::ostream& out, Tuple<T_,_size,TSpec> const &a) {
		out << "[";
		if (a.size > 0)
			out << a[0];
		for(unsigned j = 1; j < a.size; ++j)
			out << " " << a[j];
		out << "]";
		return out;
	}

	template <typename T_, unsigned _size, typename TSpec>
	struct Value< Tuple<T_, _size, TSpec> > {
		typedef T_ Type;
	};

	template <typename T_, unsigned _size, typename TSpec>
	struct Spec< Tuple<T_, _size, TSpec> > {
		typedef TSpec Type;
	};

//////////////////////////////////////////////////////////////////////////////
// getValueIx

	template <typename T1, typename T2, typename TSpec>
	inline T1 getValueI1(Pair<T1, T2, TSpec> const &pair) {
		return pair.i1;
	}

	template <typename T1, typename T2, typename TSpec>
	inline T2 getValueI2(Pair<T1, T2, TSpec> const &pair) {
		return pair.i2;
	}

//____________________________________________________________________________

	template <typename T1, typename T2, typename T3, typename TSpec>
	inline T1 getValueI1(Triple<T1, T2, T3, TSpec> const &triple) {
		return triple.i1;
	}

	template <typename T1, typename T2, typename T3, typename TSpec>
	inline T2 getValueI2(Triple<T1, T2, T3, TSpec> const &triple) {
		return triple.i2;
	}

	template <typename T1, typename T2, typename T3, typename TSpec>
	inline T3 getValueI3(Triple<T1, T2, T3, TSpec> const &triple) {
		return triple.i3;
	}

//////////////////////////////////////////////////////////////////////////////
// assignValueIx

	template <typename T1, typename T2, typename TSpec, typename T>
	inline void assignValueI1(Pair<T1, T2, TSpec> &pair, T const &_i) {
		pair.i1 = _i;
	}

	template <typename T1, typename T2, typename TSpec, typename T>
	inline void assignValueI2(Pair<T1, T2, TSpec> &pair, T const &_i) {
		pair.i2 = _i;
	}

//____________________________________________________________________________

	template <typename T1, typename T2, typename T3, typename TSpec, typename T>
	inline T const assignValueI1(Triple<T1, T2, T3, TSpec> &triple, T const &_i) {
		return triple.i1 = _i;
	}

	template <typename T1, typename T2, typename T3, typename TSpec, typename T>
	inline T const assignValueI2(Triple<T1, T2, T3, TSpec> &triple, T const &_i) {
		return triple.i2 = _i;
	}

	template <typename T1, typename T2, typename T3, typename TSpec, typename T>
	inline T const assignValueI3(Triple<T1, T2, T3, TSpec> &triple, T const &_i) {
		return triple.i3 = _i;
	}

//////////////////////////////////////////////////////////////////////////////
// operator ==/!= for pairs and triples

	template <typename L1, typename L2, typename LCompression, typename R1, typename R2, typename RCompression>
	inline bool operator==(Pair<L1, L2, LCompression> const &_left, Pair<R1, R2, RCompression> const &_right) {
		return _left.i1 == _right.i1 && _left.i2 == _right.i2;
	}
	template <typename L1, typename L2, typename LCompression, typename R1, typename R2, typename RCompression>
	inline bool operator!=(Pair<L1, L2, LCompression> const &_left, Pair<R1, R2, RCompression> const &_right) {
		return _left.i1 != _right.i1 || _left.i2 != _right.i2;
	}

//____________________________________________________________________________

	template <
		typename L1, typename L2, typename L3, typename LCompression, 
		typename R1, typename R2, typename R3, typename RCompression>
	inline bool operator==(Triple<L1, L2, L3, LCompression> const &_left, Triple<R1, R2, R3, RCompression> const &_right) {
		return _left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 == _right.i3;
	}
	template <
		typename L1, typename L2, typename L3, typename LCompression, 
		typename R1, typename R2, typename R3, typename RCompression>
	inline bool operator!=(Triple<L1, L2, L3, LCompression> const &_left, Triple<R1, R2, R3, RCompression> const &_right) {
		return _left.i1 != _right.i1 || _left.i2 != _right.i2 || _left.i3 != _right.i3;
	}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
