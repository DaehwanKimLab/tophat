#ifndef SEQAN_ITERATOR_JOURNAL_H
#define SEQAN_ITERATOR_JOURNAL_H

namespace seqan{

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec = Strict >
   class jiter{

   public:
      inline jiter() :
         m_journal( 0 ),
         m_recalc( 0 ),
         m_it_tree( 0 ),
         m_it_outer( 0 ),
         m_it_inner( 0 )
      {
      }

      inline jiter( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const & other ) :
         m_journal( other.journal() ),
         m_recalc( other.recalc() ),
         m_it_tree( other.it_tree() ),
         m_it_outer( other.it_outer() ),
         m_it_inner( other.it_inner() )
      {
      }

      inline jiter( Journal< TValue, TSpec, TStringSpec, TSloppySpec > const *journal, typename Iterator< String< Node, TStringSpec > >::Type it )
       : m_journal( journal )
       , m_it_tree( it )
       {
         m_it_outer = m_journal->get_outer_begin() + ( *m_it_tree ).index;
         m_it_inner = m_journal->get_inner_begin() + ( *m_it_tree ).index;
         m_recalc = ( *m_it_tree ).length;
      }

      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > operator=( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const &other ){
         m_journal = other.journal();
         m_recalc = other.recalc();
         m_it_tree = other.it_tree();
         m_it_outer = other.it_outer();
         m_it_inner = other.it_inner();
         return *this;
      }

#ifndef NDEBUG
      inline void debug(){
         (*m_it_tree).print_debug();
         std::cout << std::endl;
      }
#endif
      inline TValue & operator*() const{
         if( ( *m_it_tree ).is_internal ){
            journal_iterator_proxy< TValue, TSloppySpec > jip( *m_it_inner );
            return *jip;
         }else{
            journal_iterator_proxy< TValue, TSloppySpec > jip( *m_it_outer );
            return *jip;
         }
      }

      inline operator void * () {
        	return m_it_inner;
      }

      inline operator TValue * () {
         return &( **this );
      }

      inline void operator+=( size_t by ){
         size_t new_pos = abs_pos() + by;
         /*std::cout << "JIter+=(" << by << "):\n" << "> current position: " << abs_pos() << "\n> desired position: " << new_pos << std::endl;
         std::cout << "> current interval:" << std::endl;*/
         while( new_pos >= (*m_it_tree).position + (*m_it_tree).length && j_goNext( m_it_tree ) ){}
         m_recalc = (*m_it_tree).length - ( new_pos - (*m_it_tree).position );
         if( (*m_it_tree).is_internal ){
            m_it_inner = m_journal->get_inner_begin() + ( ( *m_it_tree ).index + ( new_pos - (*m_it_tree).position ));
         }else{
            m_it_outer = m_journal->get_outer_begin() + ( ( *m_it_tree ).index + ( new_pos - (*m_it_tree).position ));
         }
      }

      inline void operator-=( size_t by ){
         while( by != 0 ){
            operator--();
            --by;
         }
      }

      inline bool operator!=( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const &other ) const{
         return !operator==( other );
      }

      inline bool operator==( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const &other ) const{
         return ( other.journal() == m_journal ) && ( other.abs_pos() == this->abs_pos() );
      }

//      template< typename TIntegral >
      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > operator-( int by ) const{
         jiter new_iter = *this;
         new_iter -= by;
         return new_iter;
      }

//      template< typename TIntegral >
      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > operator+( int by ) const{
         jiter new_iter = *this;
         new_iter += by;
         return new_iter;
      }

      inline size_t operator-( jiter<TValue, TSpec, TStringSpec, TSloppySpec> & other ) const{
         return this->abs_pos() - other.abs_pos();
      }

      inline size_t operator-( jiter<TValue, TSpec, TStringSpec, TSloppySpec> const & other ) const{
         return this->abs_pos() - other.abs_pos();
      }

      inline bool operator>(jiter<TValue, TSpec, TStringSpec, TSloppySpec> const & other ) const{
         return this->abs_pos() > other.abs_pos();
      }

      inline bool operator<(jiter<TValue, TSpec, TStringSpec, TSloppySpec> const & other ) const{
         return this->abs_pos() < other.abs_pos();
      }

      inline bool operator>=(jiter<TValue, TSpec, TStringSpec, TSloppySpec> const & other ) const{
         return this->abs_pos() >= other.abs_pos();
      }

      inline bool operator<=(jiter<TValue, TSpec, TStringSpec, TSloppySpec> const & other ) const{
         return this->abs_pos() <= other.abs_pos();
      }

      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > & operator++(){
         if( --m_recalc <= 0 ){
            next_node_forward();
         }else{
            add_one();
         }
         return *this;
      }

      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > operator++(int){
         jiter< TValue, TSpec, TStringSpec, TSloppySpec > tmp = *this;
         if( --m_recalc <= 0 ){
            next_node_forward();
         }else{
            add_one();
         }
         return tmp;
      }

      inline void next_node_forward(){
         if( j_goNext( m_it_tree ) ){
            if( (*m_it_tree).length == 0 ){
               next_node_forward();
            }else{
               //if( (*m_it_tree).is_internal ){
                  m_it_inner = m_journal->get_inner_begin() + ( *m_it_tree ).index;
               //}else{
                  m_it_outer = m_journal->get_outer_begin() + ( *m_it_tree ).index;
               //}
               m_recalc = ( *m_it_tree ).length;
            }
         }else{
            while( j_goDownRight( m_it_tree ) ) { }; //TODO: this is prolly unneccessary
         }
      }

      inline jiter< TValue, TSpec, TStringSpec, TSloppySpec > & operator--(){
         if( ++m_recalc > ( *m_it_tree ).length ){
            next_node_backward();
         }else{
            minus_one();
         }
         return *this;
      }

      inline void next_node_backward(){
         if( j_goPrev( m_it_tree ) ){
            if( (*m_it_tree).length == 0 ){
               next_node_backward();
            }else{
               //if( (*m_it_tree).is_internal ){
                  m_it_inner = m_journal->get_inner_begin() + ( *m_it_tree ).index + ( *m_it_tree ).length - 1;
               //}else{
                  m_it_outer = m_journal->get_outer_begin() + ( *m_it_tree ).index + ( *m_it_tree ).length - 1;
               //}
               m_recalc = 1;
            }
         }else{
            while( j_goDownLeft( m_it_tree ) ){ };
         }
      }

      inline void add_one( ){
//         if( (*m_it_tree).is_internal ){
            ++m_it_inner;
 //        }else{
            ++m_it_outer;
   //      }
      }

      inline void minus_one( ){
     //    if( (*m_it_tree).is_internal ){
            --m_it_inner;
       //  }else{
            --m_it_outer;
         //}
      }

      inline bool is_journal( Journal< TValue > const &jrn ) const{
         return jrn == *m_journal;
      }

      inline bool at_node( typename Iterator< String< Node, TStringSpec > >::Type node_iterator ) const{
         return m_it_tree == node_iterator;
      }

      inline bool is_outer_iterator( typename Iterator< String< TValue, TSpec > >::Type const it ) const{
         return m_it_outer == it;
      }

      // BEGIN Getters

      inline Journal< TValue, TSpec, TStringSpec, TSloppySpec > const *journal() const{
         return m_journal;
      }

      inline size_t recalc() const{
         return m_recalc;
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type it_tree() const{
         return m_it_tree;
      }

      inline typename Iterator< String< TValue, TSpec > >::Type it_outer() const{
         return m_it_outer;
      }

      inline typename Iterator< String< TValue, TStringSpec > >::Type it_inner() const{
         return m_it_inner;
      }

      inline int abs_pos() const{
         return (*m_it_tree).position + (*m_it_tree).length - m_recalc;
      }

      // END Getters

   private:
      Journal<TValue, TSpec, TStringSpec, TSloppySpec> const *m_journal;
      int m_recalc;
      typename Iterator< String< Node, TStringSpec > >::Type m_it_tree;
      typename Iterator< String< TValue, TSpec > >::Type m_it_outer;
      typename Iterator< String< TValue, TStringSpec > >::Type m_it_inner;
   };

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   inline TValue const & value( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const & me ){
      return *me;
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   inline TValue const & value( jiter< TValue, TSpec, TStringSpec, TSloppySpec > & me ){
      return *me;
   }

   template< typename TIter >
   inline bool j_goNext( TIter &it ){
      DEBUG_OUT("Processing Node-Position: " << (*it).position);
      if ( j_goDownRight( it ) ){
         while( j_goDownLeft( it ) ) { };
      }else{
         size_t i;
         do{
            i = j_goUp( it );
#ifndef NDEBUG
               std::cout << "      i -> " << i << std::endl;
#endif
            if( i == 0 ) return false;
         }while ( i != 1 );
      }
      return true;
   }

   template< typename TIter >
   inline bool j_goDownRight( TIter &it ){
      if( (*it).rightChild != 0){
#ifndef NDEBUG
            std::cout << "   -> moving to right Child!";
#endif
            it += ( (*it).rightChild - (*it).tree_index );
#ifndef NDEBUG
            std::cout << " Now at position: " << (*it).position << std::endl;
#endif
            return true;
      }
      return false;
   }

   template< typename TIter >
   inline bool j_goDownLeft( TIter &it ){
      if( (*it).leftChild != 0){
#ifndef NDEBUG
            std::cout << "   -> moving to left Child!";
#endif
            it += ( (*it).leftChild - (*it).tree_index );
#ifndef NDEBUG
            std::cout << " Now at position: " << (*it).position << std::endl;
#endif
            return true;
      }
      return false;
   }

   template< typename TIter >
   inline size_t j_goUp( TIter &it ){
      if( (*it).parent != (*it).tree_index ){
#ifndef NDEBUG
            std::cout << "   -> moving up!";
#endif
         size_t retval = (*it).tree_index;
         it += ( (*it).parent - (*it).tree_index );
         if( (*it).rightChild == retval ){
            retval = 2;
         }else{
            retval = 1;
         }
#ifndef NDEBUG
            std::cout << " Now at position: " << (*it).position << std::endl;
#endif
         return retval;
      }else{
         return 0;
      }
   }

   template< typename TIter >
   inline bool j_goPrev( TIter &it ){
#ifndef NDEBUG
         std::cout << "Processing Node-Position: " << (*it).position << std::endl;
#endif
      if ( j_goDownLeft( it ) ){
         while( j_goDownRight( it ) ) {
         };
      }else{
         size_t i;
         do{
            i = j_goUp( it );
#ifndef NDEBUG
               std::cout << "      i -> " << i << std::endl;
#endif
            if( i == 0 ) return false;
         }while ( i != 2 );
      }
      return true;
   }
}

#endif // ndef(SEQAN_ITERATOR_JOURNAL_H)
