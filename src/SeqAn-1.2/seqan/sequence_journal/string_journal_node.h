#ifndef SEQAN_STRUCT_NODE_H
#define SEQAN_STRUCT_NODE_H

#include <cstddef>
#include <iostream>

struct Node{

   inline Node() : is_internal( true ), position( 0 ), index( 0 ), tree_index( 0 ), parent( 0 ), leftChild( 0 ), rightChild( 0 ), length( 0 ), operation( new Operation){};

   inline Node( bool is_internal, std::size_t position, std::size_t index, std::size_t tree_index, std::size_t parent, std::size_t leftChild, std::size_t rightChild, std::size_t length, Operation* op = new Operation )
   :  is_internal( is_internal ),
      position( position ),
		index( index ),
		tree_index( tree_index ),
      parent( parent ),
      leftChild( leftChild ),
      rightChild( rightChild ),
      length( length ),
      operation( op ){};

   inline Node( Node const & other ) :
      is_internal( other.is_internal ),
      position( other.position ),
      index( other.index ),
	   tree_index( other.tree_index ),
      parent( other.parent ),
      leftChild( other.leftChild ),
      rightChild( other.rightChild ),
      length( other.length ),
      operation( other.operation->copy() ) { }

   inline ~Node(){
      delete operation;
   }

   inline Node& operator = ( Node const& other ) {
      this->position = other.position;
      this->is_internal = other.is_internal;
      this->length = other.length;
      this->index = other.index;
      this->parent = other.parent;
      this->leftChild = other.leftChild;
      this->rightChild = other.rightChild;
      this->operation = other.operation->copy();
      return *this;
   }

   inline bool in_interval( std::size_t pos ) const{
      return ( pos >= position && pos < position + length );
   }

   inline bool influenced( std::size_t pos ) const{
      return ( pos >= position && pos < position + abs( operation->by() ) );
   }

   inline bool left_of( std::size_t pos ) const{
      return ( pos < position );
   }

   inline bool operator==( Node &other ){
      return (
               ( this->position == other.position )       &&
               ( this->is_internal == other.is_internal ) &&
               ( this->length == other.length )           &&
               ( this->index == other.index )             &&
               ( this->parent == other.parent )           &&
               ( this->leftChild == other.leftChild )     &&
               ( this->rightChild == other.rightChild )
              );
   }

   inline bool operator<( Node &other ){
      return this->position < other.position;
   }

   inline bool operator>( Node &other ){
      return this->position > other.position;
   }

   inline void synchronize_operation(){
      operation->set_by( length );
   }

   inline std::size_t offset( std::size_t pos ){
      return pos - position;
   }

   inline Operation* op(){
      return operation;
   }

   inline void print_debug(){
	   std::cout << ( is_internal ? "Internal " : "External " ) << "Node( #" << (int)this << " ): Position: " << position << " Length: " << length << " Indices: " << parent << "|" << tree_index << "|" << leftChild << "|" << rightChild << " !";
   }

    inline void print_info(){
        std::cout << ( is_internal ? " <Int " : " <Ext " ) << position << " | " << length << " - " << operation->info() << " " << operation->by() << " >" << std::endl;
    }

   bool is_internal;
   std::size_t position;
   std::size_t index;
   std::size_t tree_index;
   std::size_t parent;
   std::size_t leftChild;
   std::size_t rightChild;
   std::size_t length;
   Operation* operation;
};

namespace seqan {

   struct tree_visitor{

      virtual void pre(Node const& ) { }
      virtual void acc(Node & ) { }
      virtual void post(Node const& ) { }

       virtual ~tree_visitor() = 0;
   };

   tree_visitor::~tree_visitor() { }

   struct inorder_print : tree_visitor{
      inline void acc(Node & n) { std::cout << "< Node( " << n.length <<  " ) " << n.position << "| "; }
      inline void post(Node const& n) { std::cout << " | " << n.position + n.length << " >"; }
   };

   struct inorder_dbg_print : tree_visitor{
      //void acc(Node & n) { std::cout << "< " << (n.is_internal ? "i" : "e") << "|" << n.parent << "|" << n.tree_index << "|" << n.leftChild << "|" << n.rightChild << "| "; }
      inline void acc(Node & n) { std::cout << "< " << n.operation->info() << " | " << n.operation->by() << " | "; }
   };

   struct inorder_offset : tree_visitor{
      inorder_offset( size_t by ) : m_offset( by ){}

      inline void acc(Node & n) {
         n.position += m_offset;
      }

   private:
      size_t m_offset;
   };
}

#endif // ndef(SEQAN_STRUCT_NODE_H)
