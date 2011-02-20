#ifndef SEQAN_STRING_JOURNAL_INTERFACE_H
#define SEQAN_STRING_JOURNAL_INTERFACE_H

namespace seqan{

   /*
   template< typename TValue, typename TSpec, typename TStringSpec >
   TValue const & String< TValue, Journal< TValue, TSpec, TStringSpec, Strict > >::operator[]( size_t position ){
      return this->getjournal().get( position );
   }

   template< typename TValue, typename TSpec, typename TStringSpec >
   TValue & String< TValue, Journal< TValue, TSpec, TStringSpec, Sloppy > >::operator[]( size_t position ){
      return this->getjournal().get_r( position );
   }
   */
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TPosition >
   TValue const & value( String< TValue, Journal< TValue, TSpec, TStringSpec, Strict > > const & string, TPosition position ) {
      return string.getjournal().get( position );
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TPosition >
   TValue & value( String< TValue, Journal< TValue, TSpec, TStringSpec, Sloppy > > & string, TPosition position ) {
      return string.getjournal().get( position );
   }
   /*
   template< typename TValue, typename TSpec, typename TStringSpec >
   void assignValue( String< Journal< TValue, TSpec, TStringSpec, Sloppy > > & string, size_t position, TValue const & value ){
      string.getjournal().set( position, value );
   }*/
   

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard > {
      typedef jiter< TValue, TSpec, TStringSpec, TSloppySpec > Type;
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const, Standard > {
      typedef jiter< TValue, TSpec, TStringSpec, TSloppySpec > Type;
   };
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Value< jiter< TValue, TSpec, TStringSpec, TSloppySpec > > {
      typedef TValue Type;
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Value< jiter< TValue, TSpec, TStringSpec, TSloppySpec > const > {
      typedef TValue Type;
   };
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Reference< String<TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > > {
      typedef TValue const & Type;
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct Reference< String<TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const > {
      typedef TValue const & Type;
   };

   //////////////////////////////////////////////////////////////////////////////

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct IsContiguous<String<TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > > {
      enum { VALUE = false };
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct IsContiguous<String<TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const> {
      enum { VALUE = false };
   };
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct DefaultIteratorSpec< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > > {
      typedef Standard Type;
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct DefaultIteratorSpec< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const > {
      typedef Standard Type;
   };
   
   //////////////////////////////////////////////////////////////////////////////
   // suffix array type

   template < typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   struct SAValue< Index< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, TSpec> >{
      typedef typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type Type;
   };

	template < typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
	struct Fibre< Index< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >, Fibre_SA> {
		typedef String<
			typename SAValue< Index< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, TSpec> >::Type,
			Journal< typename SAValue< Index< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, TSpec> >::Type, TSpec, TStringSpec, Sloppy > 
		> Type;		
	};
	
	template < typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
	struct Fibre< Index< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >, Fibre_LCP> {
		typedef String<
			typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type,
			Journal< typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type, TSpec, TStringSpec, Sloppy > 
		> Type;		
	};
	
   //////////////////////////////////////////////////////////////////////////////

   template< typename TValue, typename TSpec, typename TStringSpec, typename TString, typename TSloppySpec >
   void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TString &insert_string ){
      journal_string.insert( position, insert_string );
   }
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TValue value ){
      String<TValue> tmpstr;              // create temporary string
      tmpstr += value;       // save in temporary string
      journal_string.insert( position, begin( tmpstr ), 1 ); // insert one TValue from temporary string :)
   }


   template< typename TValue, typename TSpec, typename TStringSpec, typename TIterator, typename TSloppySpec >
   void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TIterator it, size_t number ){
      journal_string.insert( position, it, number );
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSource >
   void append( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, TSource const& source ){
      target.append( begin( source ), length( source ) );
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSource >
   void append( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, TSource & source ){
      target.append( begin( source ), length( source ) );
   }
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TString, typename TSloppySpec >
   void replace( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TString &insert_string ){
      journal_string.getjournal().replace( position, begin( insert_string ), length( insert_string ) );
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TIterator, typename TSloppySpec >
   void replace( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TIterator it, size_t number ){
      journal_string.getjournal().replace( position, it, number );
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition >
   void erase( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TPosition position ){
      journal_string.del( position, 1 );
   }
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition >
   void erase( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TPosition position, TPosition position_end ){
      journal_string.del( position, position_end - position );
   }
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void flatten( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string ){
      journal_string.flatten();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec>
   inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type resizeSpace( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, 
	      		size_t size, 
			      size_t pos_begin, 
			      size_t pos_end,
               size_t limit = 0 )
   {
      SEQAN_CHECKPOINT
	   return me.resizeSpace( size, pos_begin, pos_end, limit );
   }
   
   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition, typename TExpand>
   inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type 
   resizeSpace( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, 
			   typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type size, 
			   TPosition pos_begin, 
			   TPosition pos_end, 
			   Tag<TExpand> const)
	{
		SEQAN_CHECKPOINT
	   return me.resizeSpace( size, pos_begin, pos_end );
	}
   
   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TLength, typename TExpand >
   inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type resize( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, TLength new_length, Tag<TExpand> const ){
      std::cout << "Resizing JournalString: " << &me << " to length:" << new_length << std::endl;
      if( length( me ) > (unsigned int)new_length ){
         erase( me, (unsigned int)new_length, length( me ) - (unsigned int)new_length );
      }else{
         return me.getjournal().resize_me( new_length );
      }
      return new_length;
   }
   
   template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPos >
   inline typename IndexOperatorValue< TValue, TSloppySpec >::Type getValue( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, TPos pos ){
      return me[pos];
   }
   
   template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPos >
   inline TValue const & getValue( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const & me, TPos pos ){
      return me[pos];
   }
   
   template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSize, typename TExpand> 
   inline typename Size< String<TValue, Journal < TSpec, TStringSpec, TSloppySpec > > >::Type reserve( String<TValue, Journal < TSpec, TStringSpec, TSloppySpec > > & me, TSize new_capacity, Tag<TExpand> const tag){
      return length( me );
   }
   
   template<typename TValue, typename TSpec, typename TStringSpec, typename TSourceSpec, typename TSloppySpec >
   void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, TSourceSpec > &source ){
      target.assign_string( source );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &source ){
      target.assign_journal( source );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TTargetSpec, typename TSloppySpec >
   void assign( String< TValue, TTargetSpec > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &source ){
      assign( target, source.get_coherent_string() );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSourceSpec, typename TSloppySpec >
   void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, TSourceSpec > const &source ){
      target.assign_string( source );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const & source ){
      target.assign_journal( source );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TTargetSpec, typename TSloppySpec >
   void assign( String< TValue, TTargetSpec > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &source ){
      assign( target, source.get_coherent_string() );
   }

   /*template<typename TValue, typename TSpec, typename TStringSpec, typename TTagSpec>
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec > > >::Type begin( String< TValue, Journal< TValue, TSpec, TStringSpec > > const &me, Tag<TTagSpec> const ){
      return me.it_begin();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TTagSpec>
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec > > >::Type end( String< TValue, Journal< TValue, TSpec, TStringSpec > > const &me, Tag<TTagSpec> const ){
      return me.it_end();
   }*/


   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type begin( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, Standard ){
      return me.it_begin();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type end( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, Standard ){
      return me.it_end();
   }
   
   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type begin( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me, Standard ){
      return me.it_begin();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type end( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me, Standard ){
      return me.it_end();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   size_t length( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me ){
      return me.length();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void const * id( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me ){
      return me.getjournal().get_id();
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void _setLength( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > * me, size_t new_length )
   {
   SEQAN_CHECKPOINT
   	me.del( new_length, length( me ) - new_length );
   }

   template<typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   void _setLength( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, size_t new_length )
   {
   SEQAN_CHECKPOINT
   	me.del( new_length, length( me ) - new_length );
   }
   
   /*template<typename TValue, typename TSpec, typename TStringSpec >
   inline void _allocateStorage(seqan::String<TValue, seqan::Journal<TValue, TSpec, TStringSpec, Sloppy> > & me, size_t & size){
      resize( me, size, Exact() );
   }*/
}

#endif // ndef(SEQAN_STRING_JOURNAL_INTERFACE_H)
