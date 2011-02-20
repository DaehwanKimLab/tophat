#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

#ifndef SEQAN_STRING_JOURNAL_H
#define SEQAN_STRING_JOURNAL_H

namespace seqan{

   template <typename TValue, typename TSloppySpec>
   struct IndexOperatorValue {
      typedef TValue & Type;
   };
   
   template <typename TValue>
   struct IndexOperatorValue<TValue, Strict> {
      typedef TValue const & Type;
   };

   template < typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   class String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > {
   public:
      //typedef TValue* TIter;

      String(){ };
      
      template< typename _TText >
      inline String( _TText &_text ) : m_journal( _text )
      {
         SEQAN_CHECKPOINT //TODO: check why this doesn't work with fibres?
      }

      /*String( String< TValue, TSpec > &other ) : m_journal( other )
      {
         SEQAN_CHECKPOINT
      }*/
      
      template< typename TSize >
      inline typename IndexOperatorValue< TValue, TSloppySpec >::Type operator[]( TSize pos ){
         return m_journal.get( pos );      
      }
      
      template< typename TSize >
      inline TValue const & operator[]( TSize pos ) const{
         return m_journal.get( pos );
      }

      template< typename TString >
      inline void insert( size_t position, TString &insert_string ){
         m_journal.insert( position, begin( insert_string ), seqan::length( insert_string ) );
      }
      
      template< typename TIterator >
      inline void insert( size_t position, TIterator array_start, size_t number ){
         m_journal.insert( position, array_start, number );
      }
      
      inline void inorder( tree_visitor &v, size_t node_pos ){
         m_journal.inorder( v, node_pos );
      }
#ifndef NDEBUG
      inline void inorder_dbg( tree_visitor &v, size_t node_pos ){
         m_journal.inorder_dbg( v, node_pos );
      }
#endif
      inline void print_sequence(){
         m_journal.print_sequence();
      }

      inline void print_nodes_inorder(){
         m_journal.print_nodes_inorder();
      }

      inline void del( size_t position, size_t number ){
         m_journal.del( position, number );
      }

      inline void insert( size_t position, typename Iterator< String<TValue> >::Type it, size_t number ){
         m_journal.insert( position, it, number );
      }

      template< typename TIterator >
      inline void append( TIterator array_start, size_t number ){
         m_journal.append( array_start, number );
      }

      inline void nodes_inorder( std::vector<Node> &the_vector, size_t node_pos = 0 ){
         m_journal.nodes_inorder( the_vector, node_pos );
      }

      inline void print_nodes( size_t node_pos = 0 ){
         m_journal.print_nodes( node_pos );
      }

      inline size_t resizeSpace( size_t size, size_t pos_begin, size_t pos_end, size_t limit = 0 ){
         assert( pos_begin <= pos_end );
         if( limit != 0 && m_journal.length() + size > limit ){
            size = limit - m_journal.length() ;
         }
         
         if( pos_end - pos_begin < size ){
            String<TValue> tmpstr;
            for( size_t i = 0; i < size - pos_end + pos_begin; ++i ){
               tmpstr += (TValue)0;
            }
            m_journal.insert( pos_begin, begin( tmpstr ), size - pos_end + pos_begin );
         }else{
            m_journal.del( pos_begin + size, pos_end - size );
         }
         return size;
      }

      inline void assign_string( String< TValue, TSpec > const &source ){
         m_journal.assign_string( source );
      }

      inline void assign_journal( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const & source ){
         m_journal = source.getjournal();
      };

      inline void assign_string( String< TValue, TSpec > &source ){
         m_journal.assign_string( source );
      }

      inline void assign_journal( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & source ){
         m_journal = source.getjournal();
      };

      inline Journal< TValue, TSpec, TStringSpec, TSloppySpec> & getjournal(){
         return m_journal;
      };

      inline Journal< TValue, TSpec, TStringSpec, TSloppySpec > getjournal() const{
         return m_journal;
      };
      
      inline void flatten(){
         m_journal.flatten();
      }

      inline String< TValue > get_coherent_string(){
         m_journal.flatten();
         return value( m_journal.get_holder() );
      }

      inline typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it_begin( ) const{
         return m_journal.it_begin();
      }

      inline typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it_end( ) const{
         return m_journal.it_end();
      }

      inline size_t length() const{
         return m_journal.length();
      }

   private:
      Journal< TValue, TSpec, TStringSpec, TSloppySpec > m_journal;
   };

}

#endif // ndef(SEQAN_STRING_JOURNAL_H)
