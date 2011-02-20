#ifndef SEQAN_STRING_JOURNAL_FORWARDS_H
#define SEQAN_STRING_JOURNAL_FORWARDS_H

struct Node;
struct tree_visitor;

namespace seqan{
   
   struct Sloppy;
   struct Strict;
   
   //TODO: Move SloppySpec into different file?
   template < typename TValue, typename TSloppySpec >
   struct SloppyValue;

   template < typename TValue >
   struct SloppyValue< TValue, Sloppy > {
      typedef TValue & Type;
   };

   template < typename TValue >
   struct SloppyValue< TValue, Strict > {
      typedef TValue const & Type;
   };

   template< typename TValue, typename TSloppySpec >
   struct journal_iterator_proxy{      
      typedef typename SloppyValue< TValue, TSloppySpec >::Type Type;
      journal_iterator_proxy( TValue & value ) : m_value( value ){};
      Type m_value;
   };
   
   template< typename TValue, typename TSloppySpec >
   typename SloppyValue< TValue, TSloppySpec >::Type operator*( journal_iterator_proxy< TValue, TSloppySpec > jtp ){
      return jtp.m_value;
   }

   template< typename TValue, typename TSpec = Alloc<>, typename TStringSpec = Alloc<>, typename TSloppySpec = Strict >
   class Journal;
   
   template < typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   class String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >;

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   class jiter;
}

#endif // ndef(SEQAN_STRING_JOURNAL_FORWARDS_H)
