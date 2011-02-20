#ifndef SEQAN_JOURNAL_H
#define SEQAN_JOURNAL_H

#include <iostream>
#include <vector>
#include <cassert>

#include <seqan/basic.h>
#include <seqan/sequence.h>

namespace seqan {

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec >
   class Journal{
   public:

      inline Journal(): m_holder( "" ){ //TODO: check for behavior
         assert( !empty(m_holder) );
         appendValue( m_tree, Node( false, 0, 0, 0, 0, 0, 0, 0 ) );
         m_length = 0;
      }

      template< typename TText >
      inline Journal( TText & underlying_string ):
         m_holder( underlying_string )
      {
         assert( !empty(m_holder) );
         appendValue( m_tree, Node( false, 0, 0, 0, 0, 0, 1, seqan::length( underlying_string ), new Operation( seqan::length( underlying_string ) ) ) );
         appendValue( m_tree, Node( true, seqan::length( underlying_string ), seqan::length( m_insertion_string ), 1, 0, 0, 0, 1, new Operation( 1 ) ) ); //Dummy trail-node to catch iterator overflow
         m_length = seqan::length( value( m_holder ) );

      }

      ~Journal(){
      }

      inline size_t resize_me( size_t new_length ){
         flatten();
         resize( value( m_holder ), new_length );
         (*begin( m_tree )).length = new_length;
         m_length = new_length;
         return m_length;
      }

      inline void assign_string( String< TValue, TSpec > &underlying_string ){
         Holder<String< TValue, TSpec > > holder;
         assign( m_holder, holder );
         setValue( m_holder, underlying_string );
         assert( !empty(m_holder) );
         assert( dependent(m_holder) );
         Node root( false, 0, 0, 0, 0, 0, 0, seqan::length( underlying_string ), new Operation( seqan::length( underlying_string ) ) );
         String<Node, TStringSpec> tree;
         tree += root;
         assign( m_tree, tree );
         m_length = seqan::length( value( m_holder ) );
      }

      inline TValue const & operator[]( size_t position ) const{
         return get( position );
      }

      inline TValue const & get( size_t position ) const{
         typename Iterator< String< Node, TStringSpec > >::Type it = find_node( position );
         if( !(*it).is_internal ){
            return value( m_holder )[ (*it).index + (*it).offset( position ) ];
         }else{
            return m_insertion_string[ (*it).index + (*it).offset( position ) ];
         }
      }

      inline TValue & get( size_t position ){
         typename Iterator< String< Node, TStringSpec > >::Type it = find_node( position );
         if( !(*it).is_internal ){
            return value( m_holder )[ (*it).index + (*it).offset( position ) ];
         }else{
            return value( m_insertion_string, (*it).index + (*it).offset( position ));
         }
      }

      inline void set( size_t position, TValue const & value ){
         typename Iterator< String< Node, TStringSpec > >::Type it = find_node( position );
         if( !(*it).is_internal ){
            value( m_holder )[ (*it).index + (*it).offset( position ) ] = value;
         }else{
            m_insertion_string[ (*it).index + (*it).offset( position ) ] = value;
         }
      }

      inline void print_sequence(){
         std::cout << std::endl << "Content:" << std::endl;
         for( size_t i = 0; i < m_length; ++i ){
            std::cout << get( i ) << " ";
         }
         std::cout << std::endl;
      }

      inline void print_nodes( size_t node_pos = 0 ){

         if( m_tree[ node_pos ].leftChild != 0 )
            print_nodes( m_tree[ node_pos ].leftChild );

         std::cout << get_by_node( node_pos );

         if( m_tree[ node_pos ].rightChild != 0 )
            print_nodes( m_tree[ node_pos ].rightChild );
      }

      inline String< TValue, TSpec> get_by_node( size_t node_pos ){
         Node node = m_tree[node_pos];
         String< TValue, TSpec> output = "";
         if( node.op()->deletion() ){
            return output;
         }
         resize( output, node.length, Exact() );
         if( node.is_internal ){
            typename Iterator< String< TValue, TStringSpec > >::Type it = seqan::begin( m_insertion_string, Standard() ) + node.index;
            for( size_t i = 0; i < node.length; ++i ){
               output[i] = *it;
               ++it;
            }
         }else{
            typename Iterator< String< TValue, TSpec > >::Type it = seqan::begin( value( m_holder ), Standard() ) + node.index;
            for( size_t i = 0; i < node.length; ++i ){
               output[i] = *it;
               ++it;
            }
         }
         return output;
      }

      inline Node find_node( size_t position, Node parent ) const{
         if( parent.in_interval( position) ){
            return parent;
         }else if( parent.left_of( position ) ){
            return find_node( position, m_tree[parent.leftChild] );
         }else{
            return find_node( position, m_tree[parent.rightChild] );
         }
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type find_node( size_t position ) const{
         typename Iterator< String< Node, TStringSpec > >::Type it_tree = begin( m_tree, Standard() );
         if( (*it_tree).left_of( position ) ){
            while( !(*it_tree).in_interval( position ) && j_goPrev( it_tree ) ){}
         }else{
            while( !(*it_tree).in_interval( position ) && j_goNext( it_tree ) ){}
         }
         return it_tree;
      }


      inline Node get_node( size_t index ){
         return m_tree[index];
      }

      inline Operation* find_operation( size_t pos ){
         if( pos >= m_length ){
            return new Operation();
         }
         return find_operation( pos, begin( m_tree, Standard() ) );
      }

      inline Operation* find_operation( size_t pos, typename Iterator< String< Node, TStringSpec > >::Type it_tree ){
         if( (*it_tree).left_of( pos ) ){
            while( !(*it_tree).influenced( pos ) && j_goPrev( it_tree ) ){};
         }else{
            while( !(*it_tree).influenced( pos ) && j_goNext( it_tree ) ){};
         }
         return (*it_tree).operation;
      }

      inline size_t find_node_index( size_t position ) const{
         return (*find_node( position )).tree_index;
      }

      inline void inorder( tree_visitor &v, size_t node_pos ){
         v.pre( m_tree[node_pos] );
         if( m_tree[node_pos].leftChild != 0 )
            inorder( v, m_tree[node_pos].leftChild );
         v.acc( m_tree[node_pos] );
         if( m_tree[node_pos].rightChild != 0 )
            inorder( v, m_tree[node_pos].rightChild );
         v.post( m_tree[node_pos] );
      }

#ifndef NDEBUG
      inline void inorder_dbg( tree_visitor &v, size_t node_pos ){
         std::vector<Node> vec;
         for( size_t i = 0; i < seqan::length( m_tree ); ++i ){
            vec.push_back( m_tree[i] );
         }
         v.pre( m_tree[node_pos] );
         if( m_tree[node_pos].leftChild != 0 )
            inorder_dbg( v, m_tree[node_pos].leftChild );
         v.acc( m_tree[node_pos] );
         std::cout << get_by_node( node_pos );
         std::cout << " | " << m_tree[node_pos].position << " | " << m_tree[node_pos].length << " >";
         if( m_tree[node_pos].rightChild != 0 )
            inorder_dbg( v, m_tree[node_pos].rightChild );
         v.post( m_tree[node_pos] );
      }
#endif
      inline void up_and_right( tree_visitor &v, size_t node_pos ){
         v.pre( m_tree[node_pos] );

         if( m_tree[node_pos].parent != 0 )
            up_and_right( v, m_tree[node_pos].parent );

         if( m_tree[m_tree[node_pos].parent].leftChild == node_pos ){
            v.acc( m_tree[m_tree[node_pos].parent] );

         if( m_tree[m_tree[node_pos].parent].rightChild != 0 )
            inorder( v, m_tree[m_tree[node_pos].parent].rightChild );
         }

         v.post( m_tree[node_pos] );
      }

      template< typename TIterator >
      inline void append( TIterator array_start, size_t number ){
         typename Iterator< String<Node, TStringSpec> >::Type it_tree = begin( m_tree, Standard() );

         while( j_goUp( it_tree ) ){}; // go to Root Node
         while( j_goDownRight( it_tree ) ){}; // go to last Node (dummy block)
         Node & dummy = *it_tree;
         j_goUp( it_tree ); //last data-containing node
         Node & parent = *it_tree;

         Node rightChild( true, m_length, seqan::length( m_insertion_string ), seqan::length( m_tree ) , parent.tree_index, 0, dummy.tree_index , number, new Insertion( number ) );

         dummy.parent = seqan::length( m_tree ); //dummy block has appended node as parent
         dummy.index += number; //have dumy point to new end of m_insertion_string

         parent.rightChild = seqan::length( m_tree );
         appendValue( m_tree, rightChild, Generous() );
         for( size_t i = 0; i < number; ++i, ++array_start ){
            appendValue( m_insertion_string, *array_start, Generous() );
         }
         m_length += number;
      }

      template< typename TIterator >
      inline void insert( size_t position, TIterator array_start, size_t number ){
         size_t parent_index = find_node_index( position );
         Node parent = m_tree[ parent_index ];
         size_t left_index = seqan::length( m_tree );
         size_t right_index = seqan::length( m_tree ) + 1;
         Node leftChild( parent.is_internal, parent.position, parent.index, left_index , parent_index, parent.leftChild, 0, parent.offset( position ), parent.operation->copy() );
         Node rightChild( parent.is_internal, position, parent.index + parent.offset( position ), right_index , parent_index, 0, parent.rightChild, parent.length - parent.offset( position ), parent.operation->copy() );
         leftChild.synchronize_operation();
         rightChild.synchronize_operation();
         appendValue( m_tree, leftChild, Generous() );
         parent.is_internal = true;
         if( parent.leftChild != 0 )
            m_tree[parent.leftChild].parent = left_index;
         parent.leftChild = seqan::length( m_tree ) - 1;
         appendValue( m_tree, rightChild );
         if( parent.rightChild != 0 )
            m_tree[parent.rightChild].parent = right_index;
         parent.rightChild = seqan::length( m_tree ) - 1;
         parent.position = position;
         parent.length = number;
         parent.index = seqan::length( m_insertion_string );
         delete parent.operation; //TODO: is this necessary, possible memory leak?
         parent.operation = new Insertion( number );
         for( size_t i = 0; i < number; ++i, ++array_start ){
            appendValue( m_insertion_string, *array_start, Generous() );
         }
         m_tree[ parent_index ] = parent;
         inorder_offset io(number);
         inorder( io, right_index );
         up_and_right( io, parent_index );
         m_length += number;

         typename Iterator< String<Node, TStringSpec> >::Type it_tree = begin( m_tree, Standard() );

         while( j_goUp( it_tree ) ){}; // go to Root Node
         while( j_goDownRight( it_tree ) ){}; // go to last Node (dummy block)
         (*it_tree).index += number; //have dumy point to new end of m_insertion_string

      }

      inline void nodes_inorder( std::vector<Node> &the_vector, size_t node_pos = 0 ){
         Node const& node = m_tree[node_pos];
         if( node.leftChild != 0 ){
            nodes_inorder( the_vector, node.leftChild );
         }
         the_vector.push_back( node );
         if( node.rightChild != 0 ){
            nodes_inorder( the_vector, node.rightChild );
         }
      }

      inline void print_nodes_inorder(){
         typename Iterator< String< Node, TStringSpec > >::Type it = get_first_node();
         DEBUG_OUT(">>");
#ifdef NDEBUG
            (*it).print_info();
#else
            std::cout << ">> " << (*it).position;
#endif
         while( j_goNext( it ) ){
#ifdef NDEBUG
            (*it).print_info();
#else
            std::cout << " ," << (*it).position ;
#endif
         }
         std::cout << " !" << std::endl;
      }

      inline void flatten(){
         String< TValue, TSpec > & the_string = value( m_holder );
         resize( the_string, seqan::length(the_string) + seqan::length(m_insertion_string) );
         std::vector< Node > vec;
         nodes_inorder( vec );
         std::vector< Node >::reverse_iterator vit;
         typename Iterator< String< TValue, TSpec > >::Type sit;
         sit = begin( the_string, Standard() ) + seqan::length( the_string );
         for( vit = vec.rbegin() - 1; vit < vec.rend(); ++vit ){
            Node node = *vit;
            if( node.op()->deletion() ){
               continue;
            }
            sit -= node.length;
            if( node.is_internal ){
               for( size_t i = 1; i <= node.length; ++i ){
                  value(sit + node.length - i) = m_insertion_string[ node.index + node.length - i ];
               }
            }else{
               for( size_t i = 1; i <= node.length; ++i ){
                  value(sit + node.length - i) = the_string[ node.index + node.length - i ];
               }
            }
         }
         clear( m_tree );
         seqan::erase( the_string, (size_t)0, seqan::length(the_string) - m_length );
         Node root( false, 0, 0, 0, 0, 0, 0, seqan::length( the_string ), new Operation( m_length ) );
         m_tree += root;
         appendValue( m_tree, Node( true, seqan::length( value( m_holder ) ), seqan::length( m_insertion_string ), 1, 0, 0, 0, 1, new Operation( 1 ) ) ); //Dummy trail-node to catch iterator overflow
         vec.clear();
         clear( m_insertion_string );
      }

      inline void del( size_t position, int number ){
#ifndef NDEBUG
         std::cout << "del< " << position << ", " << number << " >\tJournal" << std::endl;
#endif
         assert( position + number <= m_length );
         size_t parent_index = find_node_index( position );
         Node parent = m_tree[ parent_index ];
#ifndef NDEBUG
         std::cout << "deletion-parent: " << parent.operation->info() << std::endl;
#endif
         int remaining_length = number - ( parent.position + parent.length - position );
         int overlap_correction = 0;

         if( remaining_length > 0 ){
            number -= remaining_length;
         }
         size_t left_index = seqan::length( m_tree );
         size_t right_index = seqan::length( m_tree ) + 1;

         Node leftChild( parent.is_internal, parent.position, parent.index, left_index, parent_index, parent.leftChild, 0, parent.offset( position ), parent.operation->copy() );
         Node rightChild( parent.is_internal, position + number, parent.index + parent.offset( position + number ), right_index, parent_index, 0, parent.rightChild, parent.length - parent.offset( position + number ), parent.operation->copy() );
         leftChild.synchronize_operation();
         rightChild.synchronize_operation();
         if( parent.op()->insertion() ){
            overlap_correction = (int)parent.length - (int)leftChild.length - (int)rightChild.length;
         }

         appendValue( m_tree, leftChild );
         if( parent.leftChild != 0 )
            m_tree[parent.leftChild].parent = left_index;
         parent.is_internal = true;
         parent.leftChild = left_index;
         appendValue( m_tree, rightChild );
         if( parent.rightChild != 0 )
            m_tree[parent.rightChild].parent = right_index;
         parent.rightChild = right_index;
         parent.position = position;
         parent.length = 0; //TODO: implement if neccessary: -number implies node.length to be int;
         parent.index = 0;
         inorder_offset io( -number );
         inorder( io, right_index );
         up_and_right( io, parent_index );
         m_length -= number;
         delete parent.operation; //TODO: is this necessary, possible memory leak?
         parent.operation = new Deletion( overlap_correction - number );
         m_tree[ parent_index ] = parent;
#ifndef NDEBUG
         std::cout << "new-deletion-parent: " << parent.operation->info() << std::endl;
#endif
         if( remaining_length > 0 ){
            del( position, remaining_length );
         }
      }

      template< typename TIterator >
      inline void replace( size_t position, TIterator array_start, int number ){
         del( position, number );
         insert( position, array_start, number );
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type get_first_node() const{
         return begin( m_tree, Standard() ) + find_node_index( 0 );
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type get_last_node() const{
         return begin( m_tree, Standard() ) + find_node_index( m_length - 1 );
      }

      inline typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it_begin() const{
         typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it( this, get_first_node() );
         return it;
      }

      inline typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it_end() const{
         typename Iterator<String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type it( this, get_last_node() );
         it += (*get_last_node()).length;
         it.add_one();
         return it;
      }

      inline typename Iterator< String< TValue, TSpec > >::Type get_outer_begin() const{
         return begin( value( m_holder ), Standard() );
      }

      inline typename Iterator< String< TValue, TStringSpec > >::Type get_inner_begin() const{
         return begin( m_insertion_string, Standard() );
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type get_tree_begin() const{
         return begin( m_tree, Standard() );
      }

      inline typename Iterator< String< Node, TStringSpec > >::Type get_zero_node_iterator() const{
         return begin( m_tree, Standard() ) + find_node_index( 0 );
      }

      inline String< TValue, TSpec > get_outer() const{
         return value( m_holder );
      }

      inline String< TValue, TStringSpec > get_inner() const{
         return m_insertion_string;
      }

      inline String< Node, TStringSpec > get_tree() const{
         return m_tree;
      }

      inline Holder< String< TValue, TSpec > > get_holder() const{
         return m_holder;
      }

      inline size_t length() const{
         return m_length;
      }

      inline bool operator!=( Journal< TValue, TSpec, TStringSpec, TSloppySpec > const &other ) const{
         return !operator==( other );
      }

      inline bool operator==( Journal< TValue, TSpec, TStringSpec, TSloppySpec > const &other ) const{
         return get_tree() == other.get_tree() && get_inner() == other.get_inner() && get_outer() == other.get_outer();
      }

      inline void const * get_id(){
         return id( value( m_holder ) );
      }

   private:
      Holder< String<TValue, TSpec> > m_holder;
      String<TValue, TStringSpec> m_insertion_string;
      String<Node, TStringSpec> m_tree;
      size_t m_length;
   };
} // namespace seqan
#endif // ndef(SEQAN_JOURNAL_H)
