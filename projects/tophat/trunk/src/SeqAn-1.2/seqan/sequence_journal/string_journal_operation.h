#include <string>
#include <sstream>

struct Operation{
   Operation( ) : m_by( 0 ) {};

   Operation( int by ) : m_by( by ) {};

   virtual std::string info() const{
      return "Operation";
   }
   
   virtual bool insertion() const {
      return false;
   }
   
   virtual bool deletion() const {
      return false;
   }
   
   virtual ~Operation(){}
   
   virtual Operation* copy() const {
      return new Operation( *this );
   };
   
   int by() const {
      return m_by;
   }
   
   void set_by( int to ){
      m_by = to;
   }
   
private:
   int m_by;
   
   
   Operation& operator =( Operation const& other );
   
protected:
   Operation( Operation const & other ) : m_by(other.m_by) {}
   
};

struct Insertion : Operation{
public:
   Insertion( int by ) : Operation( by ) {};
   
   ~Insertion(){};
   
   std::string info() const{
      return "Insertion";
   }
   
   bool insertion() const {
      return true;
   }
   
   Operation* copy() const {
      return new Insertion( *this );
   };
};

struct Deletion : Operation{
public:
   Deletion( int by ) : Operation( by ) {};
   
   ~Deletion(){};
   
   std::string info() const {
      return "Deletion";
   }
   
   bool deletion() const {
      return true;
   }
   
   Operation* copy() const {
      return new Deletion( *this );
   };
};

struct Entry{
public:

   Entry() : m_operation(), m_position(0) {};
   Entry( Entry const & other ) : m_operation( other.m_operation->copy() ), m_position( other.m_position ) {};
   Entry( Operation* operation, size_t position ) : m_operation( operation ), m_position( position ) {};
   
   ~Entry() { delete m_operation; }
   
   Entry& operator = ( Entry const& other ) {
      this->m_operation = other.m_operation->copy();
      this->m_position = other.m_position;
      return *this;
   }

   size_t position() const {
      return m_position;
   }

   int count() const {
      return m_operation->by();
   }

   void print() const {
      std::cout << this->to_s() << std::endl;
   }
   
   bool deletion() const {
      return m_operation->deletion();
   }
   
   bool insertion() const {
      return m_operation->insertion();
   }
   
   std::string to_s() const{
      std::stringstream r;
      r << "JournalEntry\t-\t" << m_operation->info() << "\tat: " << m_position << "\tcount: " << count();
      return r.str();
   }


private:
   Operation* m_operation;
   size_t m_position;
};

