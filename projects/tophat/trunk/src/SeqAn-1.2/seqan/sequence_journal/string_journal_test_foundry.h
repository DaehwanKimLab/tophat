#include <cstdlib>

   unsigned int randomNumber( int upper_bound )
   {
      int value = ( (float)rand() / RAND_MAX )*( upper_bound + 1 );
      if( value > upper_bound )
      {
         value = 0;
      }
      return (unsigned int)floor(value);
   }

   template< typename TSpec >
   void generate_random_string( unsigned int length, seqan::String< char, TSpec > & the_string, seqan::String< char > const & alphabet_string ){
      seqan::resize( the_string, length );
      typename seqan::Iterator< seqan::String< char, TSpec > >::Type it = seqan::begin( the_string );
      typename seqan::Iterator< seqan::String< char, TSpec > >::Type the_end = seqan::end( the_string );
      while( it != the_end ){
         *it = alphabet_string[ randomNumber( seqan::length( alphabet_string ) - 1 ) ];
         std::cout << *it << " ";
         ++it;
      }
      std::cout << std::endl;
      return; //TODO: best practice, is that useful / necessary
   }

template< typename TString >
class JournalTest{

public:
   JournalTest( TString& string ): the_holder( string ) {
      srand( time(NULL) );
   };
   
   void random_insertions( unsigned int number, unsigned int length, seqan::String< char > const & alphabet_string ){
      if( number > 0 ){
         seqan::String< char > insertion_string;
         seqan::resize( insertion_string, length );
         for( unsigned int i = 0; i < length; ++i ){
            insertion_string[i] = alphabet_string[ randomNumber( seqan::length( alphabet_string ) - 1 ) ];
         }
         seqan::insert( randomNumber( seqan::length( seqan::value( the_holder ) ) - 1 ), seqan::value( the_holder ), insertion_string ); //insert teh shit into teh string
         random_insertions( --number, length, alphabet_string );
      }
      return; //TODO: best practice, is that useful / necessary
   }

   void random_deletions( unsigned int number, int length = 10 ){
      if( number > 0 ){
         int position = randomNumber( seqan::length( seqan::value( the_holder ) ) - length );
         seqan::erase( seqan::value( the_holder ), position, position + length ); //erase some stuff from teh string
         random_deletions( --number, length );
      }
      return; //TODO: best practice, is that useful / necessary
   }

private:
   seqan::Holder< TString > the_holder;
};

