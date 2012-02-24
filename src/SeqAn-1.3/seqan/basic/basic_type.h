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

#ifndef SEQAN_HEADER_BASIC_TYPE_H
#define SEQAN_HEADER_BASIC_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Value:
..summary:Type of the items in the container. 
..signature:Value<T>::Type
..param.T:Type for which the value type is determined.
..returns.param.Type:Value type of $T$.
..remarks.text:The value type of a container $T$ is the type of the elements in $T$.
    For example, the value type of a sequence of $int$ is $int$.
..example.code:Value<String<char> >::Type c; //c has type char
..include:seqan/basic.h
*/

template <typename T, const int i = 0>
struct Value
{
	typedef T Type;
};
/*
template <typename T>
struct Value<T const>
{
	typedef T Type;
};
*/
//____________________________________________________________________________

/**
.Metafunction.GetValue:
..summary:Type for reading values. 
..signature:GetValue<T>::Type
..param.T:Type of container that holds a value.
..returns.param.Type:GetValue type of $T$.
..remarks.text:Depending on $T$, the $GetValue$-type can either be $Value<T>::Type &$ or $Value<T>::Type$.
..text:$GetValue$ is the return type of @Function.getValue@ that allows a (read-only) access to objects.
Do not confuse it with @Function.value@ that returns a @Metafunction.Reference.reference@ to the value.
..see:Metafunction.Value
..see:Function.getValue
..include:seqan/basic.h
*/
template <typename T>
struct GetValue
{
	typedef typename Value<T>::Type const & Type;
};
template <typename T>
struct GetValue<T const>:
	public GetValue<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Reference:
..summary:Reference type. 
..signature:Reference<T>::Type
..param.T:A Type.
..returns.param.Type:Either $T &$ or a proxy object @Class.Proxy@ for $T$.
..see:Metafunction.Value
..see:Metafunction.GetValue
..include:seqan/basic.h
*/
template <typename T>
struct Reference
{
	typedef typename Value<T>::Type & Type;
};
template <typename T>
struct Reference<T const>
{
	typedef typename Value<T>::Type const & Type;
};

//____________________________________________________________________________


/**
.Metafunction.Size:
..summary:Type of an object that is suitable to hold size information.
..signature:Size<T>::Type
..param.T:Type for which the size type is determined.
..returns.param.Type:Size type of $T$.
..remarks.text:In most cases this type is $size_t$.
..include:seqan/basic.h
*/
template <typename T>
struct Size
{
	typedef size_t Type;
};
template <typename T>
struct Size<T const>:
	Size<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Difference:
..summary:Type of an object that stores the difference between two iterators.
..signature:Difference<T>::Type
..param.T:Type for which the difference type is determined.
...type:Class.Iter
..returns.param.Type:Difference type of $T$.
..remarks.text:In most cases this type is $ptrdiff_t$.
..see:Metafunction.Size
..include:seqan/basic.h
*/
template <typename T>
struct Difference
{
	typedef ptrdiff_t Type;
};
template <typename T>
struct Difference<T const>:
	Difference<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Position:
..summary:Type of an object that represents a position in a container.
..signature:Position<T>::Type
..param.T:Type for which the position type is determined.
...type:Class.Iter
...type:Class.String
..returns.param.Type:Position type of $T$.
..see:Metafunction.Iterator
..include:seqan/basic.h
*/
template <typename T>
struct Position
{
	typedef typename Size<T>::Type Type;
};
template <typename T>
struct Position<T const>:
	Position<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Host:
..summary:Type of the object a given object depends on.
..signature:Host<T>::Type
..param.T:Type for which the host type is determined.
..returns.param.Type:Host type of $T$.
..include:seqan/basic.h
*/
template <typename T>
struct Host
{
	typedef T Type;
};

//____________________________________________________________________________

/**
.Metafunction.Spec:
..summary:The spec of a class. 
..signature:Spec<T>::Type
..param.T:Type for which the spec is determined.
..returns.param.Type:Spec of $T$.
..remarks:The spec of a SeqAn type is the class that is used in template subclassing 
 to specify the specialization. 
 For example, the spec of $String<char, Alloc<> >$ is $Alloc<>$.
..include:seqan/basic.h
*/

// default case
template <typename T>
struct Spec {
	typedef void Type;
};


// one argument case
template <template <typename> class T, typename TSpec>
struct Spec< T<TSpec> > {
	typedef TSpec Type;
};

template <typename T>
struct Spec<T const>:
	public Spec<T> {};

//____________________________________________________________________________

/**
.Metafunction.DeepestSpec:
..summary:The deepest spec of a class with nested template arguments.
..signature:DeepestSpec<T>::Type
..param.T:Type for which the deepest spec is determined.
..returns.param.Type:Deepest spec of $T$.
..remarks:The spec of a SeqAn type is the innermost class that is used in nested subclassing.
 For example, the deepest spec of $Iter<..., VSTree<BottomUp<Mums> > >$ is $Mums$.
..include:seqan/basic.h
*/

// default case
template <typename T>
struct DeepestSpec {
	typedef T Type;
};

// recursion for 1 argument
template <
	template <typename> class T, 
	typename T1 >
struct DeepestSpec< T<T1> > {
	typedef typename 
		If<
			IsSameType<T1, void>::VALUE,										// is T1 void?
			T<T1>,															// yes, end of recursion
			typename DeepestSpec< typename Spec< T<T1> >::Type >::Type		// no,  recurse
		>::Type Type;
};

// recursion for 2 arguments
template <
	template <typename, typename> class T, 
	typename T1, typename T2 >
struct DeepestSpec< T<T1,T2> >:
	DeepestSpec< typename Spec< T<T1,T2> >::Type > {};

// recursion for 3 arguments
template <
	template <typename, typename, typename> class T, 
	typename T1, typename T2, typename T3 >
struct DeepestSpec< T<T1,T2,T3> >:
	DeepestSpec< typename Spec< T<T1,T2,T3> >::Type > {};

// recursion for 4 arguments
template <
	template <typename, typename, typename, typename> class T, 
	typename T1, typename T2, typename T3, typename T4 >
struct DeepestSpec< T<T1,T2,T3,T4> >:
	DeepestSpec< typename Spec< T<T1,T2,T3,T4> >::Type > {};

// recursion for 5 arguments
template <
	template <typename, typename, typename, typename, typename> class T, 
	typename T1, typename T2, typename T3, typename T4, typename T5 >
struct DeepestSpec< T<T1,T2,T3,T4,T5> >:
	DeepestSpec< typename Spec< T<T1,T2,T3,T4,T5> >::Type > {};

template <typename T>
struct DeepestSpec<T const>:
	public DeepestSpec<T> {};

//____________________________________________________________________________

/**
.Metafunction.Cargo:
..summary:Type of additional data stored in an object. 
..signature:Cargo<T>::Type
..param.T:Type for which the cargo tyoe is determined.
..returns.param.Type:Cargo of $T$.
..remarks:The definition of Cargo allows the addition of user specific data to existing data structures.
..include:seqan/basic.h
*/

template <typename T>
struct Cargo {
	typedef Nothing Type;
};
template <typename T>
struct Cargo<T const> {
	typedef typename Cargo<T>::Type const Type;
};

//____________________________________________________________________________

/**
.Metafunction.VertexDescriptor:
..summary:Type of an object that represents a vertex descriptor.
..signature:VertexDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use ids as vertex descriptors.
..returns.param.Type:VertexDescriptor type.
..remarks.text:The vertex descriptor is a unique handle to a vertex in a graph.
It is used in various graph functions, e.g., to add edges, to create OutEdge Iterators or to remove a vertex.
It is also used to attach properties to vertices.
..example.code:VertexDescriptor<Graph<> >::Type vD; //vD is a vertex descriptor
..include:seqan/basic.h
*/

template <typename T>
struct VertexDescriptor {
	typedef void* Type;
};
template <typename T>
struct VertexDescriptor<T const>:
	public VertexDescriptor<T> {};


//____________________________________________________________________________

	
/**
.Metafunction.Id:
..summary:Type of an object that represents an id.
..signature:Id<T>::Type
..param.T:Type for which a suitable id type is determined.
..returns.param.Type:Id type.
..remarks.text:The id type of a container is the type that is used to uniquely identify its elements.
In most cases this type is unsigned int.
..example.code:Id<Graph<> >::Type id; //id has type unsigned int
..include:seqan/basic.h
*/
template<typename T>
struct Id {
	typedef unsigned int Type;
};

//____________________________________________________________________________

template<typename T>
struct Id<T const> {
	typedef unsigned int Type;
};

//____________________________________________________________________________

/**
.Metafunction.Key:
..summary:Key type of a key to cargo mapping.
..signature:Key<T>::Type
..param.T:Type for which a key type is determined.
..returns.param.Type:Key type.
...default:The type $T$ itself.
..include:seqan/basic.h
*/
template< typename T >
struct Key
{
	typedef T Type;
};

template <typename T>
struct Key<T const>:
	Key<T> {};

//____________________________________________________________________________

/*VERALTET
.Metafunction.Object:
..summary:Object type of a key to object mapping.
..signature:Object<T>::Type
..param.T:Type for which a object type is determined.
..returns.param.Type:Object type.
..include:seqan/basic.h
*/

template<typename T>
struct Object; 

template <typename T>
struct Object<T const>:
	Object<T> {};


//____________________________________________________________________________

/**
.Metafunction.Source
..include:seqan/basic.h
*/

template < typename TSpec = void >
struct Source
{
	typedef TSpec Type;
};

template <typename T>
struct Source<T const>:
	Source<T>
{
};

//____________________________________________________________________________

/**
.Internal.Parameter_:
..cat:Metafunctions
..summary:Type for function parameters and return values.
..signature:Parameter_<T>::Type
..param.T:A type.
..returns.param.Type:The parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $Parameter_<T>::Type$ is $T$, 
otherwise $Parameter_<T>::Type$ is $T &$.
*/
template <typename T>
struct Parameter_
{
	typedef T & Type;
};

template <typename T>
struct Parameter_<T *>
{
	typedef T * Type;
};
template <typename T, size_t I>
struct Parameter_<T [I]>
{
	typedef T * Type;
};


/**
.Internal._toParameter:
..cat:Functions
..summary:Transforms pointers to parameter types.
..signature:_toParameter<T>(pointer)
..param.pointer:A pointer.
..param.T:A Type.
...text:$object$ is transformed into the parameter type of $T$ that is given by @Internal.Parameter_@.
...note:This type must be explicitely specified.
..returns:To $TParameter$ transformed $object$.
..see:Internal.Parameter_
*/
template <typename T>
typename Parameter_<T>::Type
_toParameter(T * _object)
{
SEQAN_CHECKPOINT
	return * _object;
}
template <typename T>
typename Parameter_<T>::Type
_toParameter(T _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//____________________________________________________________________________

/**
.Internal.ConstParameter_:
..cat:Metafunctions
..summary:Type for constant function parameters and return values.
..signature:ConstParameter_<T>::Type
..param.T:A type.
..returns.param.Type:The const parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $Parameter_<T>::Type$ is a pointer to a const array, 
otherwise $Parameter_<T>::Type$ is $T const &$.
..see:Internal.Parameter_
*/
template <typename T>
struct ConstParameter_
{
	typedef T const & Type;
};
template <typename T>
struct ConstParameter_<T const>:
	public ConstParameter_<T> {};

template <typename T>
struct ConstParameter_<T *>
{
	typedef T const * Type;
};
template <typename T>
struct ConstParameter_<T const *>
{
	typedef T const * Type;
};

template <typename T, size_t I>
struct ConstParameter_<T [I]>
{
	typedef T const * Type;
};
template <typename T, size_t I>
struct ConstParameter_<T const [I]>
{
	typedef T const * Type;
};

//____________________________________________________________________________

/**
.Internal.Pointer_:
..cat:Metafunctions
..summary:The associated pointer type.
..signature:Pointer_<T>::Type
..param.T:A type.
..returns.param.Type:A pointer type for $T$.
...text:if $T$ is already a pointer type, then $Pointer_<T>::Type$ is $T$,
otherwise $Pointer_<T>::Type$ is $T *$.
..see:Internal.Parameter_
..see:Internal._toParameter
*/
template <typename T>
struct Pointer_
{
	typedef T * Type;
};

template <typename T>
struct Pointer_<T *>
{
	typedef T * Type;
};
template <typename T>
struct Pointer_<T * const>
{
	typedef T * Type;
};

template <typename T, size_t I>
struct Pointer_<T [I]>
{
	typedef T * Type;
};

//non const version of Pointer_ for return values

template <typename T>
struct NonConstPointer_:
	Pointer_<T>
{
};
template <typename T>
struct NonConstPointer_<T * const>
{
	typedef T * Type;
};

/**
.Internal._toPointer:
..cat:Functions
..summary:Transforms types into pointers.
..signature:_toPointer(object)
..param.object:An object.
..returns:$object$, transformed to a pointer. 
...text:The type of the returned pointer is given by @Internal.Pointer_@.
..see:Internal.Pointer_
*/
template <typename T>
typename NonConstPointer_<T>::Type
_toPointer(T & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}
template <typename T>
typename NonConstPointer_<T const>::Type
_toPointer(T const & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}

template <typename T>
typename NonConstPointer_<T *>::Type
_toPointer(T * _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//____________________________________________________________________________


/**
.Metafunction.LENGTH:
..summary:Number of elements in a fixed-size container.
..signature:LENGTH<T>::VALUE
..param.T:Type for which the number of elements is determined.
..returns.param.VALUE:Number of elements.
..remarks.text:The default return value is 1 for dynamic-size containers.
..include:seqan/basic.h
*/
template <typename T>
struct LENGTH
{
	enum { VALUE = 1 };
};
template <typename T>
struct LENGTH<T const>:
	LENGTH<T>
{
};

/**
.Metafunction.WEIGHT:
..summary:Number of relevant positions in a shape.
..signature:WEIGHT<T>::Type
..param.T:Shape type for which the number of relevant positions is determined.
...type:Class.Shape
..returns.param.VALUE:Number of relevant positions.
..remarks.text:The default return value is the result of the @Metafunction.LENGTH@ function.
For gapped shapes this is the number of '1's.
..include:seqan/basic.h
*/
template <typename T>
struct WEIGHT:
	LENGTH<T>
{
};
template <typename T>
struct WEIGHT<T const>:
	WEIGHT<T>
{
};

//////////////////////////////////////////////////////////////////////////////

//Iterator: see basic_iterator.h

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.IsIntegral:
..summary:Tests for a type to be of integral value.
..signature:IsIntegral<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is a simple type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..include:seqan/basic.h
 */
template <typename T>
struct IsIntegral
{
	typedef
		// Implicitely signed.
		typename If< IsSameType<T, char>::VALUE,  True,
		typename If< IsSameType<T, char>::VALUE,  True,
		typename If< IsSameType<T, short>::VALUE, True,
		typename If< IsSameType<T, int>::VALUE,   True,
		typename If< IsSameType<T, long>::VALUE,  True,
		typename If< IsSameType<T, __int64>::VALUE,      True,
		// Explicitely signed.
		typename If< IsSameType<T, signed char>::VALUE,    True,
		typename If< IsSameType<T, signed char>::VALUE,    True,
		typename If< IsSameType<T, signed short>::VALUE,   True,
		typename If< IsSameType<T, signed int>::VALUE,     True,
		typename If< IsSameType<T, signed long>::VALUE,    True,
		// Explicitely unsigned.
		typename If< IsSameType<T, unsigned char>::VALUE,  True,
		typename If< IsSameType<T, unsigned char>::VALUE,  True,
		typename If< IsSameType<T, unsigned short>::VALUE, True,
		typename If< IsSameType<T, unsigned int>::VALUE,   True,
		typename If< IsSameType<T, unsigned long>::VALUE,  True,
		typename If< IsSameType<T, __uint64>::VALUE,       True,
		False
		>::Type>::Type>::Type>::Type>::Type>::Type
		>::Type>::Type>::Type>::Type>::Type
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct IsIntegral<T const> {
	typedef typename IsIntegral<T>::Type Type;
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
