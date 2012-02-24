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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Tags, declarations and generic code for the String class and its
// specializations.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_ARRAY_BASE_H
#define SEQAN_HEADER_SEQUENCE_ARRAY_BASE_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////


template <typename TSpec = void>
struct Alloc {};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/**
.Class.String:
..cat:Sequences
..summary:General purpose container for sequences.
..signature:String<TValue, TSpec>
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Alloc<>$, see @Spec.Alloc String@.
..implements:Concept.Container
..include:seqan/sequence.h
*/

template <typename TValue, typename TSpec = Alloc<> >
class String;




//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.String

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value<String<TValue, TSpec> >
{
    typedef TValue Type;
};
template <typename TValue, typename TSpec>
struct Value<String<TValue, TSpec> const >:
    public Value<String<TValue, TSpec> >
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.String

template <typename TValue, typename TSpec>
struct Spec<String<TValue, TSpec> >
{
    typedef TSpec Type;
};
template <typename TValue, typename TSpec>
struct Spec<String<TValue, TSpec> const>:
    public Spec<String<TValue, TSpec> >
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IsSequence.param.T.type:Class.String

template <typename TValue, typename TSpec>
struct IsSequence<String<TValue, TSpec> > {
    typedef True Type;
    enum { VALUE = true };
};

//////////////////////////////////////////////////////////////////////////////
// Returns a Class that can be used to store a temporary copy of a String

template <typename T>
struct TempCopy_
{
    typedef typename Value<T>::Type TValue;
    typedef typename RemoveConst_<TValue>::Type TValue_NotConst;
    typedef String<TValue_NotConst, Alloc<> > Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

///.Function.id.param.object.type:Class.String
///.Function.empty.param.object.type:Class.String
///.Function.capacity.param.object.type:Class.String


//////////////////////////////////////////////////////////////////////////////
// swap
//////////////////////////////////////////////////////////////////////////////

// forwards for swap

template<typename TAlphabet, typename TSpec> class String;
template<typename TString, typename TSpec> class StringSet;

/**
.Function.swap:
..summary:Swaps the contents of two values.
..cat:Content Manipulation
..signature:swap(left, right)
..param.left:The first value.
...type:Class.String
...type:Class.StringSet
..param.right:The second value.
...type:Class.String
...type:Class.StringSet
..remarks:The function swaps the values of variables left and right.
This is equivalent to using move three times with a temporary variable.

Note that this function has the same name as the STL function $std::swap$ but is in a different namespace.
Argument Dependent Lookup (ADL, aka Koenig lookup) will take care that the right $swap$ function is called from STL $sort$, for example.
We only specialize it for Class.String and Class.StringSet.
..see:Function.move
..include:seqan/sequence.h
*/

template <typename TAlphabet, typename TSpec>
inline void
swap(String<TAlphabet, TSpec> & left,
     String<TAlphabet, TSpec> & right)
{
    SEQAN_CHECKPOINT;

    typedef String<TAlphabet, TSpec> TString;

    TString tmp(left, Move());
    move(left, right);
    move(right, tmp);
}

// TODO(holtgrew): This actually belongs into sequence_multiple.h.
template <typename TString, typename TSpec>
inline void
swap(StringSet<TString, TSpec> & left,
     StringSet<TString, TSpec> & right)
{
    SEQAN_CHECKPOINT;
    typedef StringSet<TString, TSpec> TStringSet;

    TStringSet tmp(left, Move());
    move(left, right);
    move(right, tmp);
}


//////////////////////////////////////////////////////////////////////////////
//shareResources

///.Function.shareResources.param.sequence1, sequence2.type:Class.String

template <typename TValue, typename TSpec>
inline bool
shareResources(String<TValue, TSpec> const & obj1,
                    TValue const & obj2)
{
SEQAN_CHECKPOINT
    return (begin(obj1) >= &obj2) && (end(obj1) <= &obj2);
}

template <typename TValue, typename TSpec>
inline bool
shareResources(TValue const & obj1,
                    String<TValue, TSpec> const & obj2)
{
SEQAN_CHECKPOINT
    return (begin(obj2) >= &obj1) && (end(obj2) <= &obj1);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.begin.param.object.type:Class.String
///.Function.end.param.object.type:Class.String
//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.container.type:Class.String

template <typename TValue, typename TSpec, typename TPos>
inline typename Reference< String<TValue, TSpec> >::Type
value(String<TValue, TSpec> & me,
      TPos pos)
{
    typedef typename Position< String<TValue, TSpec> >::Type TStringPos;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LT_MSG(static_cast<TStringPos>(pos), static_cast<TStringPos>(length(me)), "Trying to access an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

template <typename TValue, typename TSpec, typename TPos>
inline typename Reference< String<TValue, TSpec> const >::Type
value(String<TValue, TSpec> const & me,
      TPos pos)
{
    typedef typename Position< String<TValue, TSpec> const >::Type TStringPos;
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LT_MSG(static_cast<TStringPos>(pos), static_cast<TStringPos>(length(me)), "Trying to access an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.String

template <typename TValue, typename TSpec>
inline typename Size< String<TValue, TSpec> const>::Type
length(String<TValue, TSpec> const & me)
{
SEQAN_CHECKPOINT
    return end(me, Standard()) - begin(me, Standard());
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.String

template <typename TValue, typename TSpec>
inline bool
empty(String<TValue, TSpec> const & me)
{
    SEQAN_CHECKPOINT
    return end(me, Standard()) == begin(me, Standard());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.clear:
..cat:Containers
..summary:Resets an object.
..signature:clear(object)
..param.object:The object that will be resetted.
...type:Class.String
..remarks:$object$ is set to a state that is equivalent to a default constructed object of the same type.
..remarks:If $object$ is a container, then all elements are removed from this container.
The length is set to 0.
The capacity can be changed, depending on the implementation.
..see:Function.resize
..see:Function.length
..include:seqan/sequence.h
*/

template <typename TValue, typename TSpec>
inline void
clear(String<TValue, TSpec> & me)
{
SEQAN_CHECKPOINT
    arrayDestruct(begin(me, Standard()), end(me, Standard()));
    _setLength(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct ClearSpaceStringBase_
{
};

//____________________________________________________________________________

template <>
struct ClearSpaceStringBase_<Insist>
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
SEQAN_CHECKPOINT
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        if (limit < size)
        {
SEQAN_CHECKPOINT
            size = limit;
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
SEQAN_CHECKPOINT
        typename Size<T>::Type new_length = length(seq) + size - (end - start);
        arrayClearSpace(begin(seq, Standard()) + start, length(seq) - start, end - start, size);
        _setLength(seq, new_length);
        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Value<T>::Type * seq_buffer = begin(seq);
        typename Size<T>::Type seq_length = length(seq);

        if (limit > start + size)
        {
SEQAN_CHECKPOINT
            typename Size<T>::Type removed_size = end - start;
            typename Size<T>::Type new_length = seq_length - removed_size + size;
            if (limit < new_length)
            {
SEQAN_CHECKPOINT
                arrayDestruct(seq_buffer + limit, seq_buffer + new_length);
                seq_length -= new_length - limit;
            }
            arrayClearSpace(seq_buffer + start, seq_length - start, end - start, size);
            _setLength(seq, new_length);
            return size;
        }
        else
        {
SEQAN_CHECKPOINT
            arrayDestruct(seq_buffer + start, seq_buffer + seq_length);
            _setLength(seq, limit);
            if (limit > start) return limit - start;
            else return 0;
        }
    }
/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};


//____________________________________________________________________________


template <>
struct ClearSpaceStringBase_<Limit>
{

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
SEQAN_CHECKPOINT
        return _clearSpace(seq, size, capacity(seq), Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type seq_capacity = capacity(seq);
        if (limit > seq_capacity)
        {
SEQAN_CHECKPOINT
            limit = seq_capacity;
        }
        return _clearSpace(seq, size, limit, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
SEQAN_CHECKPOINT
        return _clearSpace(seq, size, start, end, capacity(seq), Insist());
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type seq_capacity = capacity(seq);
        if (limit > seq_capacity)
        {
SEQAN_CHECKPOINT
            limit = seq_capacity;
        }
        return _clearSpace(seq, size, start, end, limit, Insist());
    }

/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};

//____________________________________________________________________________

template <typename TExpand>
struct ClearSpaceExpandStringBase_
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, size, TExpand());
        if (old_array)
        {
SEQAN_CHECKPOINT
            _deallocateStorage(seq, old_array, old_capacity);
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        if (limit < size)
        {
SEQAN_CHECKPOINT
            size = limit;
        }
        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, size, limit, TExpand());
        if (old_array)
        {
            _deallocateStorage(seq, old_array, old_capacity);
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
        typename Size<T>::Type old_length = length(seq);
        typename Size<T>::Type removed_size = end - start;
        typename Size<T>::Type new_length = old_length - removed_size + size;

        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, new_length, TExpand());
        typename Value<T>::Type * seq_array = begin(seq);

        if (old_array)
        {
SEQAN_CHECKPOINT
            arrayConstructMove(old_array, old_array + start, seq_array);
            arrayConstructMove(old_array + end, old_array + old_length, seq_array + start + size);
            _deallocateStorage(seq, old_array, old_capacity);
        }
        else
        {
SEQAN_CHECKPOINT
            arrayClearSpace(seq_array + start, old_length - start, removed_size, size);
        }

        _setLength(seq, new_length);

        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type old_length = length(seq);
        typename Size<T>::Type removed_size = end - start;
        typename Size<T>::Type need_length = old_length - removed_size + size;

        typename Size<T>::Type new_length = need_length;
        typename Size<T>::Type length_to_copy = old_length;
        if (limit < need_length)
        {
SEQAN_CHECKPOINT
            new_length = limit;
            length_to_copy = new_length - size + removed_size;
        }

        bool keep_second_part = (new_length > start + size);

        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, new_length, limit, TExpand());
        typename Value<T>::Type * seq_array = begin(seq);

        if (old_array)
        {//new buffer allocated
            typename Size<T>::Type keep_start_length = (start > new_length) ? new_length : start;
            arrayConstructMove(old_array, old_array + keep_start_length, seq_array);
            if (keep_second_part)
            {
                arrayConstructMove(old_array + end, old_array + length_to_copy, seq_array + start + size);
            }
            _deallocateStorage(seq, old_array, old_capacity);
        }
        else
        {
            if (keep_second_part)
            {
                arrayClearSpace(seq_array + start, length_to_copy - start, end - start, size);
                if (length_to_copy < old_length)
                {
                    arrayDestruct(seq_array + length_to_copy, seq_array + old_length);
                }
            }
            else
            {
                arrayDestruct(seq_array + start, seq_array + old_length);
            }
        }

        _setLength(seq, new_length);

        if (keep_second_part) return size;
        else if (new_length > start) return new_length - start;
        else return 0;
    }

/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, TExpand());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, TExpand());
    }
*/
};

//____________________________________________________________________________

template <>
struct ClearSpaceStringBase_<Exact>:
    ClearSpaceExpandStringBase_<Exact>
{
};

//____________________________________________________________________________

template <>
struct ClearSpaceStringBase_<Generous>:
    ClearSpaceExpandStringBase_<Generous>
{
};

//____________________________________________________________________________
/**
.Internal._clearSpace:
..cat:Functions
..summary:Makes space in container
..signature:_clearSpace(object, size [, pos_begin, pos_end] [, limit], resize_tag)
..param.object:The container.
..param.size:Length of the freed space.
..param.pos_begin:Position of the first item in $object$ that is to be destroyed. (optional)
..param.pos_end:Position behind the last item in $object$ that is to be destroyed. (optional)
...remarks:If $pos_end == pos_begin$, no item in $object$ will be destroyed.
..param.limit:Maximal length $object$ can get after this operation. (optional)
..param.resize_tag:Strategy that is applied if $object$ has not enough capacity to store the complete content.
..returns:The number of free characters.
...remarks:Depeding on the @Tag.Overflow Strategy.overflow strategy@ specified by $resize_tag$,
this could be $size$ or less than $size$ if $object$ has not enough @Function.capacity@.
..remarks:This function is similar to @Function.resizeSpace@ and @Function.fillSpace@.
The main difference is that $_clearSpace$ does not construct objects in the new created space.
*/
template<typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
        TSize size,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> const>::_clearSpace_(me, size);
}

template<typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
        TSize size,
        TSize limit,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> const>::_clearSpace_(me, size, limit);
}

template<typename TValue, typename TSpec, typename TSize, typename TPosition, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename TSpec, typename TSize, typename TPosition, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            TSize limit,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> const>::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.resizeSpace.param.object.type:Class.String

template<typename TValue, typename TSpec, typename TSize, typename TBeginPosition, typename TEndPosition, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resizeSpace(String<TValue, TSpec> & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
    typedef typename Size< String<TValue, TSpec> >::Type TSize_;
    typedef typename Position<String<TValue, TSpec> >::Type TPos_;
    TSize_ ret_ = _clearSpace(
        me,
        static_cast<TSize_>(size),
        static_cast<TPos_>(pos_begin),
        static_cast<TPos_>(pos_end),
        tag);
    arrayConstruct(iter(me, pos_begin), iter(me, pos_begin) + ret_);
    return ret_;
}

template<typename TValue, typename TSpec, typename TSize, typename TBeginPosition, typename TEndPosition, typename TLimit, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resizeSpace(String<TValue, TSpec> & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            TLimit limit,
            Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
    typedef typename Size< String<TValue, TSpec> >::Type TSize_;
    typedef typename Position<String<TValue, TSpec> >::Type TPos_;
    TSize_ ret_ = _clearSpace(
        me,
        static_cast<TSize_>(size),
        static_cast<TPos_>(pos_begin),
        static_cast<TPos_>(pos_end),
        static_cast<TSize_>(limit),
        tag);
    arrayConstruct(iter(me, pos_begin), iter(me, pos_begin) + ret_);
    return ret_;
}

//////////////////////////////////////////////////////////////////////////////
// assignValue (2)
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TValue2>
inline void
assignValue(String<TValue, TSpec> & /*me*/,
            TValue2 const & /*_value*/)
{
//  assign(me, toString(_value)); ???TODO
}

//???TODO: moveValue (2)

//////////////////////////////////////////////////////////////////////////////
// assign
//////////////////////////////////////////////////////////////////////////////

///.Function.assign.param.target.type:Class.String
///.Function.assign.param.source.type:Class.String

//overload of binary version for strings:

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
      TSource & source)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
      TSource const & source)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AssignString_
{
    template <typename TTarget, typename TSource>
    static inline void
    assign_(
        TTarget & target,
        TSource & source)
    {
        if (empty(source) && empty(target))
            return;  // Do nothing if both source and target are empty.
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()));
        }
        else
        {
SEQAN_CHECKPOINT
            if ((void *) &target == (void *) &source) return;

            typename TempCopy_<TSource>::Type temp(source, length(source));
            assign(target, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    assign_(
        TTarget & target,
        TSource & source,
        typename Size<TTarget>::Type limit)
    {
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type part_length = _clearSpace(target, typename Size<TTarget>::Type(length(source)), limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()));
        }
        else
        {
SEQAN_CHECKPOINT
            if ((void *) &target == (void *) &source) return;

            typename Size<TTarget>::Type source_length = length(source);
            if (source_length > limit) source_length = limit;

            typename TempCopy_<TSource>::Type temp(source, source_length);
            assign(target, temp, TExpand());
        }
    }
};

//////////////////////////////////////////////////////////////////////////////

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       Tag<TExpand> const)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AssignString_<Tag<TExpand> const>::assign_(target, source);
}
template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       TSize limit,
       Tag<TExpand> const)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AssignString_<Tag<TExpand> const>::assign_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSourceValue const * source,
       Tag<TExpand> const)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AssignString_<Tag<TExpand> const>::assign_(target, source);
}
template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSourceValue const * source,
       TSize limit,
       Tag<TExpand> const)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AssignString_<Tag<TExpand> const>::assign_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// move
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
//implementation of move for contiguous sequences
//note: there is a problem, if sizeof(TSourceValue) and sizeof(TTargetValue) are not a multiple
//  of each other, since in this case the correct size cannot be determined afterwards
//  when calling the deallocate function.
//  ???TODO
template <typename TTarget, typename TSource>
void
_moveContiguous(TTarget & target,
                TSource & source)
{
    typedef typename Value<TSource>::Type TSourceValue;
    typedef typename Value<TTarget>::Type TTargetValue;

    clear(target);

    typename Iterator<TSource, Standard>::Type source_begin = begin(source, Standard());
    typename Iterator<TTarget, Standard>::Type target_begin = (typename Iterator<TTarget, Standard>::Type) begin(source, Standard());

    typename Size<TTarget>::Type size = sizeof(TSourceValue) * capacity(source);
    if (size >=  sizeof(TTargetValue))
    {
SEQAN_CHECKPOINT
        if (sizeof(TSourceValue) <= 2) ++size; //regard the "end of string termination" case
        typename Size<TTarget>::Type target_capacity = size / sizeof(TTargetValue);
        if (sizeof(TTargetValue) <= 2) --target_capacity; //regard the "end of string termination" case

        typename Size<TTarget>::Type target_length = length(source);
        if (target_length > target_capacity)
        {
            target_length = target_capacity;
        }

        if (sizeof(TSourceValue) >= sizeof(TTargetValue))
        {
            arrayMoveForward(source_begin, source_begin + target_length, target_begin);
        }
        else
        {
            arrayMoveBackward(source_begin, source_begin + target_length, target_begin);
        }

        _setBegin(target, target_begin);
        _setLength(target, target_length);
        _setCapacity(target, target_capacity);

        _setBegin(source, 0);
        _setLength(source, 0);
        _setCapacity(source, 0);
    }
    else
    {
        clear(source);
    }
}
//____________________________________________________________________________

//overload of binary version for strings:

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource & source)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource const & source)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

//____________________________________________________________________________

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TTag>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource & source,
     Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    assign(target, source, tag);
}

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TTag>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource const & source,
     Tag<TTag> const tag)
{
SEQAN_CHECKPOINT
    assign(target, source, tag);
}

//////////////////////////////////////////////////////////////////////////////
// valueConstructMove:
// it is usually better for strings to default construct and move instead of
// copy construct strings

/*
template <typename TIterator, typename TValue, typename TSpec>
inline void
valueConstructMove(TIterator it,
                   String<TValue, TSpec> const & value)
{
    valueConstruct(it);
    move(*it, value);
}
*/

//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

///.Function.append.param.target.type:Class.String
///.Function.append.param.source.type:Class.String

template <typename TExpand>
struct AppendString_
{
    template <typename TTarget, typename TSource>
    static inline void
    append_(TTarget & target,
            TSource & source)
    {
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type target_length = length(target);
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), target_length, target_length, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + target_length);
        }
        else
        {
SEQAN_CHECKPOINT
            typename TempCopy_<TSource>::Type temp(source, length(source));
            append(target, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    append_(TTarget & target,
            TSource & source,
            typename Size<TTarget>::Type limit)
    {
        typename Iterator<TTarget, Standard>::Type target_begin = begin(target, Standard());
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type target_length = length(target);
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), target_length, target_length, limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + target_length);
        }
        else
        {
            typename Size<TTarget>::Type target_length = length(target);
            if (target_length >= limit)
            {
SEQAN_CHECKPOINT
                arrayDestruct(target_begin + limit, target_begin + target_length);
                _setLength(target, limit);
            }
            else
            {
SEQAN_CHECKPOINT
                limit -= target_length;
                typename Size<TTarget>::Type source_length = length(source) ;
                if (source_length > limit) source_length = limit;

                typename TempCopy_<TSource>::Type temp(source, source_length);
                append(target, temp, TExpand());
            }
        }
    }
};

//////////////////////////////////////////////////////////////////////////////

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AppendString_<Tag<TExpand> const>::append_(target, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AppendString_<Tag<TExpand> const>::append_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSourceValue * source,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AppendString_<Tag<TExpand> const>::append_(target, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSourceValue * source,
       typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
       Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    AppendString_<Tag<TExpand> const>::append_(target, source, limit);
}



//////////////////////////////////////////////////////////////////////////////
// appendValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct AppendValueToString_
{
    template <typename T, typename TValue>
    static inline void
    appendValue_(T & me,
                TValue & _value)
    {
SEQAN_CHECKPOINT
        typename Position<T>::Type me_length = length(me);
        if (capacity(me) <= me_length)
        {
            typename Value<T>::Type temp_copy(_value); //temp copy because resize could invalidate _value
            typename Size<T>::Type new_length = resize(me, me_length + 1, TExpand());
            if (me_length < new_length)
            {
                valueConstruct(begin(me) + me_length, temp_copy); //??? this should be valueMoveConstruct
            }
        }
        else
        {
            valueConstruct(begin(me, Standard()) + me_length, _value);
            _setLength(me, me_length + 1);
        }
    }
};

//____________________________________________________________________________

template <typename TTargetValue, typename TTargetSpec, typename TValue, typename TExpand>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    AppendValueToString_<Tag<TExpand> const>::appendValue_(me, _value);
}

//////////////////////////////////////////////////////////////////////////////
// insertValue
//////////////////////////////////////////////////////////////////////////////
/**
.Function.insertValue:
..include:seqan/sequence.h
*/

template <typename TExpand>
struct InsertValueToString_
{
    template <typename T, typename TPosition, typename TValue>
    static inline void
    insertValue_(T & me,
                TPosition pos,
                TValue & _value)
    {
SEQAN_CHECKPOINT
        typename Value<T>::Type temp_copy = _value; //temp copy because resizeSpace could invalidate _value
        resizeSpace(me, 1, pos, pos, TExpand());
        if ((typename Size<T>::Type) pos < length(me))
        {
            moveValue(me, pos, temp_copy);
        }
    }
};

//____________________________________________________________________________

template <typename TTargetValue, typename TTargetSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(String<TTargetValue, TTargetSpec> & me,
            TPosition pos,
            TValue const & _value,
            Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    InsertValueToString_<Tag<TExpand> const>::insertValue_(me, pos, _value);
}

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

///.Function.replace.param.target.type:Class.String
///.Function.replace.param.source.type:Class.String

template <typename TExpand>
struct ReplaceString_
{
    template <typename TTarget, typename TSource>
    static inline void
    replace_(TTarget & target,
             typename Size<TTarget>::Type pos_begin,
             typename Size<TTarget>::Type pos_end,
             TSource & source)
    {
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + pos_begin);
        }
        else
        {
SEQAN_CHECKPOINT
            typename TempCopy_<TSource>::Type temp(source, length(source));
            replace(target, pos_begin, pos_end, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    replace_(TTarget & target,
             typename Size<TTarget>::Type pos_begin,
             typename Size<TTarget>::Type pos_end,
             TSource & source,
             typename Size<TTarget>::Type limit)
    {
        if (!id(source) || !shareResources(target, source))
        {
SEQAN_CHECKPOINT
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + pos_begin);
        }
        else
        {
            if (pos_begin >= limit)
            {
SEQAN_CHECKPOINT
                arrayDestruct(begin(target) + limit, end(target));
                _setLength(target, limit);
            }
            else
            {
SEQAN_CHECKPOINT
                limit -= pos_begin;
                typename Size<TTarget>::Type source_length = length(source) ;
                if (source_length > limit) source_length = limit;

                typename TempCopy_<TSource>::Type temp(source, source_length);
                replace(target, pos_begin, pos_end, temp, limit, TExpand());
            }
        }
    }

};

//////////////////////////////////////////////////////////////////////////////

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    ReplaceString_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    ReplaceString_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSourceValue const * source,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    ReplaceString_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSourceValue const * source,
        typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
        Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    typedef String<TTargetValue, TTargetSpec> TTarget;
    ReplaceString_<Tag<TExpand> const>::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Internal._reallocateStorage:
..cat:Functions
..summary:Allocates a new buffer if needed.
..signature:_reallocateStorage(object, new_capacity, resize_tag)
..param.object:A container for which the buffer is reallocated.
...type:Class.String
..param.new_capacity:The capacity $object$ will get after reallocating the buffer.
..param.resize_tag:Strategy that is used for changing the capacity.
..returns:Returns the old buffer, if a new buffer has been allocated, $0$ otherwise.
..remarks:This function only allocates a new buffer if the current capacity is less then $new_capacity$.
A new buffer is not filled with any content, all copy operations must be done by the caller.
..remarks:If $object$ never had a buffer, or the buffer is not changed by the function,
the returned pointer is 0.
*/

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity)
{
SEQAN_CHECKPOINT
    return _allocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    Exact)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me)) return 0;
    else
    {
SEQAN_CHECKPOINT
        return _reallocateStorage(me, new_capacity);
    }
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    TSize2 limit,
    Exact)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me)) return 0;
    else
    {
SEQAN_CHECKPOINT
        if (new_capacity > limit) new_capacity = limit;
        return _reallocateStorage(me, new_capacity);
    }
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    Generous)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me))
    {
        return 0;
    }
    else
    {
SEQAN_CHECKPOINT
        new_capacity = computeGenerousCapacity(me, new_capacity);
        return _reallocateStorage(me, new_capacity);
    }
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    TSize2 limit,
    Generous)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me)) return 0;
    else
    {
SEQAN_CHECKPOINT
        new_capacity = computeGenerousCapacity(me, new_capacity);
        if (new_capacity > limit) new_capacity = limit;
        return _reallocateStorage(me, new_capacity);
    }
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    Insist)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    TSize2,
    Insist)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    Limit)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    TSize2,
    Limit)
{
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.reserve.param.object.type:Class.String

template <typename TValue, typename TSpec, typename TSize_>
inline void
_reserveStorage(
    String<TValue, TSpec> & /*seq*/,
    TSize_ /*new_capacity*/,
    Insist)
{
    // do nothing
}

template <typename TValue, typename TSpec, typename TSize_>
inline void
_reserveStorage(
    String<TValue, TSpec> & /*seq*/,
    TSize_ /*new_capacity*/,
    Limit)
{
    // do nothing
}

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline void
_reserveStorage(
    String<TValue, TSpec> & seq,
    TSize_ new_capacity,
    Tag<TExpand> const tag)
{
    typedef typename Size< String<TValue, TSpec> >::Type TSize;

    TSize old_capacity = capacity(seq);

//  if (old_capacity == (TSize)new_capacity) return;
//  if (!IsSameType<TExpand,TagExact_>::VALUE && old_capacity > (TSize)new_capacity) return;
    if (old_capacity >= (TSize)new_capacity) return;

    TSize seq_length = length(seq);
    typename Value< String<TValue, TSpec> >::Type * old_array = _reallocateStorage(seq, new_capacity, tag);
    if (old_array)
    {//buffer was replaced, destruct old buffer
//      arrayConstruct(begin(seq, Standard()), begin(seq, Standard()) + seq_length);
//      arrayMoveForward(old_array, old_array + seq_length, begin(seq, Standard()));
        arrayConstructCopy(old_array, old_array + seq_length, begin(seq, Standard()));
        arrayDestruct(old_array, old_array + seq_length);
        _deallocateStorage(seq, old_array, old_capacity);
        _setLength(seq, seq_length);
    }
    else if (!old_capacity)
    {//new buffer created and the string had no buffer yet
        _setLength(seq, seq_length);
    }
}

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
reserve(
    String<TValue, TSpec> & seq,
    TSize_ new_capacity,
    Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
    _reserveStorage(seq, new_capacity, tag);
    return _capacityReturned(seq, new_capacity, tag);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.resize.param.object.type:Class.String

template <typename TExpand>
struct _Resize_String
{
    template <typename T>
    static inline typename Size<T>::Type
    resize_(
        T & me,
        typename Size<T>::Type new_length)
    {
        typedef typename Size<T>::Type TSize;
        TSize me_length = length(me);
        if (new_length < me_length)
        {
SEQAN_CHECKPOINT
            arrayDestruct(begin(me, Standard()) + new_length, begin(me, Standard()) + me_length);
        }
        else
        {
            typename Size<T>::Type me_capacity = capacity(me);
            if (new_length > me_capacity)
            {
SEQAN_CHECKPOINT
                TSize new_capacity = reserve(me, new_length, TExpand());
                if (new_capacity < new_length)
                {
                    new_length = new_capacity;
                }
            }
            if (new_length > me_length)
            {
SEQAN_CHECKPOINT
                arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length);
            }
        }

        _setLength(me, new_length);
        return new_length;
    }

    template <typename T, typename TValue>
    static inline typename Size<T>::Type
    resize_(
        T & me,
        typename Size<T>::Type new_length,
        TValue const & val)
    {
        typedef typename Size<T>::Type TSize;
        TSize me_length = length(me);
        if (new_length < me_length)
        {
SEQAN_CHECKPOINT
            arrayDestruct(begin(me, Standard()) + new_length, begin(me, Standard()) + me_length);
        }
        else
        {
            TSize me_capacity = capacity(me);
            if (new_length > me_capacity)
            {
SEQAN_CHECKPOINT
                TValue tempCopy = val;  // reserve could invalidate val
                TSize new_capacity = reserve(me, new_length, TExpand());
                if (new_capacity < new_length)
                {
                    new_length = new_capacity;
                }
                arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length, tempCopy);
            } else
                if (new_length > me_length)
                {
SEQAN_CHECKPOINT
                    arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length, val);
                }
        }

        _setLength(me, new_length);
        return new_length;
    }
};

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resize(
    String<TValue, TSpec> & me,
    TSize new_length,
    Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return _Resize_String<Tag<TExpand> const>::resize_(me, new_length);
}

template <typename TValue, typename TSpec, typename TSize, typename TValue2, typename TExpand>
inline TSize
resize(String<TValue, TSpec> & me,
     TSize new_length,
     TValue2 const & val,
     Tag<TExpand> const)
{
SEQAN_CHECKPOINT
    return _Resize_String<Tag<TExpand> const>::resize_(me, new_length, val);
}

//////////////////////////////////////////////////////////////////////////////

/*
template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
String<TLeftValue, TLeftSpec> const &
operator += (String<TLeftValue, TLeftSpec> & left,
             String<TRightValue, TRightSpec> const & right)
{
    append(left, right);
    return left;
}

template <typename TLeftValue, typename TLeftSpec, typename TRight>
String<TLeftValue, TLeftSpec> const &
operator += (String<TLeftValue, TLeftSpec> & left,
             TRight const & right)
{
    append(left, right);
    return left;
}

template <typename TLeft, typename TRightValue, typename TRightSpec>
TLeft const &
operator += (TLeft & left,
             String<TRightValue, TRightSpec> const & right)
{
    append(left, right);
    return left;
}
*/

template <typename TLeftValue, typename TLeftSpec, typename TRight >
String<TLeftValue, TLeftSpec> const &
operator += (String<TLeftValue, TLeftSpec> & left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    append(left, right);
    return left;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight >
inline bool
operator == (String<TLeftValue, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<String<TLeftValue, TLeftSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator == (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<String<TRightValue, TRightSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight >
inline bool
operator !=(String<TLeftValue, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<String<TLeftValue, TLeftSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator != (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<String<TRightValue, TRightSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator < (String<TLeftValue, TLeftSpec> const & left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator < (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator <= (String<TLeftValue, TLeftSpec> const & left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator <= (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator > (String<TLeftValue, TLeftSpec> const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator > (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator >= (String<TLeftValue, TLeftSpec> const & left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator >= (TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator << (TStream & target,
             String<TValue, TSpec> const & source)
{
SEQAN_CHECKPOINT
    write(target, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator >> (TStream & source,
             String<TValue, TSpec> & target)
{
SEQAN_CHECKPOINT
    read(source, target);
    return source;
}


//////////////////////////////////////////////////////////////////////////////

/*
//////////////////////////////////////////////////////////////////////////////
// Converter
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource, typename TSpec>
struct Value<String<Convert<TTarget, TSource>, TSpec> >:
    Value< String<TSource, TSpec> >
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TAccessor, typename TConverter>
struct AccessorConverter_
{
    typedef T Type;
};
template <typename T, typename TAccessor, typename TConverter>
struct AccessorConverter_<T, TAccessor &, typename TConverter &>
{
    typedef T & Type;
};

template <typename TTarget, typename TSource, typename TSpec>
struct GetValue<String<Converter<TTarget, TSource>, TSpec> >:
    AccessorConverter_<
        TTarget,
        typename GetValue<String<TSource, TSpec> >::Type,
        typename Converter<TTarget, TSource>::Type
    >
{
};
*/

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
