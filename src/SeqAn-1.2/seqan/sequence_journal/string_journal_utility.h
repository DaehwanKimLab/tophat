namespace seqan{

   template< typename TString > //TODO: rewrite to encompass general Container functionality
   struct suffix_compare_functor{

      suffix_compare_functor( TString & string ) : m_string( string ){}

      template< typename TPos1, typename TPos2 >
      bool operator()( TPos1 & first, TPos2 & second ){
         return suffix( m_string, first ) < suffix( m_string, second );
      }

   private:
      TString & m_string;
   };

    template< typename TIteratorA, typename TIteratorB >
    inline bool _suffix_bigger( TIteratorA it_a, TIteratorB it_b ){
        //std::cout << "Bigger?" << std::endl;
        while( *it_a == *it_b ){
            //std::cout << *it_a << "==" << *it_b << std::endl;
            ++it_a;
            ++it_b;
        }
        //std::cout << *it_a << "> " << *it_b << " ?" << std::endl;
        return *it_a > *it_b;
    }

    template< typename TIteratorA, typename TIteratorB >
    inline size_t _lcp_length( TIteratorA it_a, TIteratorB it_b ){
        size_t lcplength = 0;
        //std::cout << "LCP ";
        while( *it_a == *it_b ){
            ++it_a;
            ++it_b;
            ++lcplength;
            //std::cout << '.';
        }
        //std::cout << lcplength << std::endl;
        return lcplength;
    }

   inline int get_shift( String< Pair< int, int > > const & limits, int position ){   //TODO: implement binary search variant
      int shift = 0;
      for( unsigned int i = 0; i < length( limits ); ++i ){
         if( limits[i].i1 <= position ){
            shift = limits[i].i2;
         }else{
            break;
         }
      }
      return shift;
   }

   template< typename TString, typename TPos >
   inline int get_shift_b( TString const & limits, TPos position ){
      typename Iterator< TString >::Type it_upper = begin( limits );
      typename Iterator< TString >::Type it_pivot = begin( limits ) + (unsigned int)floor( length( limits ) / 2.0 ); //TODO: check integer division behavior
      typename Iterator< TString >::Type it_lower = end( limits );
      while( true ){
         if( it_pivot->i1 <= position ){
            it_upper = it_pivot;
         }else{
            it_lower = it_pivot;
         }
         if( it_upper == it_pivot || it_lower == it_pivot ){
            return it_pivot->i2;
         }
      }
   }

   bool deleted( String< Pair< int, int > > const & deletions, int position ){   //TODO: implement binary search variant
      for( unsigned int i = 0; i < length( deletions ); ++i ){
         if( position >= deletions[i].i1 && position < deletions[i].i1 + deletions[i].i2 ){
            return true;
         }
      }
      return false;
   }

   template< typename TString, typename TPos > //TString needs to be String< Pair< int, int > > for this to work
   inline bool deleted_b( TString const & deletions, TPos position ){
      if( length( deletions ) == 0 ){  //No deletions in this String
         return false;
      }

      typename Iterator< TString >::Type it_upper = begin( deletions );
      typename Iterator< TString >::Type it_pivot = begin( deletions ) + (unsigned int)floor( length( deletions ) / 2.0 ); //TODO: check integer division behavior
      typename Iterator< TString >::Type it_lower = end( deletions );

      bool found = false;

      while( !found ){
         if( position >= it_pivot->i1 ){
            if( position < it_pivot->i1 + it_pivot->i2 ){
               found = true;
               return found;
            }else{
               it_upper = it_pivot;
               if( it_lower == it_pivot || it_upper == it_pivot ){
                  return found;
               }
            }
         }else{
            it_lower = it_pivot; //TODO: check why this is so nonsensical
            if( it_lower == it_pivot || it_upper == it_pivot ){
               return found;
            }
         }
         it_pivot = it_lower + (int)floor( ( it_lower - it_upper ) / 2.0 ); //TODO: check integer division behavior
      }
      return false;
   }

   template< typename TValue1, typename TValue2, typename TTag >
   void print_pairs( String< Pair< TValue1, TValue2, TTag > > const & string ){
      for( unsigned int i = 0; i < length( string ); ++i ){
         std::cout << "< " << string[i].i1 << ", " << string[i].i2 << " > ";
      }
      std::cout << std::endl;
   }

   template< typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TIndex, typename TPos >  //TString needs to be String< TValue, Journal< TValue, ... > >
   inline void generate_shifts_and_deletions( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & string, String< TIndex > & indices, String< Pair< TPos, TPos > > & shifts, String< Pair< TPos, TPos > > & deletions ){

      shifts += Pair< TPos, TPos >( 0, 0 );

      Iterator< String< Node, Alloc<> > >::Type it_tree = string.getjournal().get_tree_begin();

      while( j_goDownLeft( it_tree ) ) {}; //go to first node; TODO:check why making string a const-Ref causes segfault and messes up get_tree_begin's returned Node

      do{
         if( (*it_tree).operation->by() < 0 && (*it_tree).operation->deletion() ){

            Pair< TPos, TPos > p( (*it_tree).position - back( shifts ).i2 + abs( (*it_tree).operation->by() ), back( shifts ).i2 + (*it_tree).operation->by() );
            if( p.i1 > back( shifts ).i1 ){
               shifts += p;
            }else if( p.i1 == back( shifts ).i1 ){
               back( shifts ).i2 = p.i2;
            }
            deletions += Pair<TPos, TPos>( p.i1 + (*it_tree).operation->by(), abs( (*it_tree).operation->by() ) );

         }else if( (*it_tree).operation->by() > 0 && (*it_tree).operation->insertion() ){

            Pair< TPos, TPos > p( (*it_tree).position - back( shifts ).i2, back( shifts ).i2 + (*it_tree).operation->by() );
            if( p.i1 > back( shifts ).i1 ){
               shifts += p;
            }else if( p.i1 == back( shifts ).i1 ){
               back( shifts ).i2 = p.i2;
            }
            for( int i = 0; i < (*it_tree).operation->by(); ++i ){
               indices += (*it_tree).position + i; //TODO: adapt this to (global)pairs as position
            }
         }
      }while( j_goNext( it_tree ) ); //iterate over the nodes inorder
   }

    template< typename TIndex, typename TString, typename TSpec >
    inline void synchronize_index( TIndex & index, StringSet< TString, TSpec > & stringset ){

        typedef typename Position< StringSet< TString, TSpec > >::Type TPos;
        typedef String< typename SAValue< TIndex >::Type, Journal< typename SAValue< TIndex >::Type, Alloc<>, Alloc<>, Sloppy > > TSA;
        typedef String< size_t, Journal< size_t, Alloc<>, Alloc<>, Sloppy > > TLCP;

        assert( indexSupplied( index, ESA_SA() ) && "No SuffixArray supplied for this index. No synchronization necessary!");
        assert( indexSupplied( index, ESA_LCP() ) && "No LCP-Table supplied for this index. No synchronization possible!\nUse 'indexCreate( index, ESA_SA() )' to recreate the SuffixArray!" );

        _refreshStringSetLimits( stringset ); //refreshing limits

        indexRequire( index, ESA_LCP() );

        StringSet< String< Pair< TPos, TPos > > > shiftset;
        StringSet< String< Pair< TPos, TPos > > > deletionset;
        //String< Journal< char, Alloc<>, Alloc<>, Sloppy > > journalset;
        String< Pair< unsigned int, unsigned int > > indices;

        resize( shiftset, length( stringset ) );
        resize( deletionset, length( stringset ) );
        //resize( journalset, length( stringset ) );


        for( unsigned int i = 0; i < length( stringset ); ++i ){
            String< unsigned int > tmpindex;
            generate_shifts_and_deletions( stringset[i], tmpindex, shiftset[i], deletionset[i] );
            for( unsigned int j = 0; j < length( tmpindex ); ++j ){
                indices += Pair< unsigned int, unsigned int >( i, tmpindex[j] );
            }
            //journalset[i] = stringset[i].getjournal();
            /*print_pairs( shiftset[i] );
            print_pairs( deletionset[i] );*/
        }

        String< size_t > dels; //Stores the positions of indices that need to be removed
        reserve( dels, 200 ); //TODO: Guesswork for avg number of recomputes dependent on Tree Size needed!!!

        TSA fibre_sa = indexSA( index );

        TLCP fibre_lcp = indexLCP( index );

        assert( length( fibre_sa ) == length( fibre_lcp ) && "Length mismatch in LCP-Table / SuffixArray, successfull synchronization not possible!\nUse 'indexCreate( index, ESA_SA() )' to recreate the SuffixArray!" );

        // BEGIN: Debug Output
		/*std::cout << "Suffix array state initial:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << "\t>" << i << "\t: " << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }*/
		// END: Debug Output+

		/*std::cout << "\nIndices:" << std::endl;
        for( size_t i = 0; i < length( indices ); ++i ){
         std::cout << indices[i].i1 << ", " << indices[i].i2 << "\t: " << suffix( stringset, indices[i] ) << "..." << std::endl;
        }*/

        unsigned int last_lcp = 0;
        size_t len = length( fibre_sa ); //minimize computational effort of loop
        Operation* op;                   //dito

        if( deleted_b( deletionset[ fibre_sa[0].i1 ], fibre_sa[0].i2 ) ){ //TODO: use std::binary_search instead
            appendValue( dels, 0, Generous() ); //Is fibre_sa[i] a Pair?
        }else{
            fibre_sa[0].i2 += get_shift( shiftset[ fibre_sa[0].i1 ], fibre_sa[0].i2 ); //TODO: use std::binary_search instead
            for( unsigned int j = 0; j <= _max( value( fibre_lcp, 0 ), last_lcp ); ++j ){
                op = stringset[ fibre_sa[0].i1 ].getjournal().find_operation( fibre_sa[0].i2 + j );
                if( op->insertion() || op->deletion() && op->by() != 0 ){
                    appendValue( indices, fibre_sa[0], Generous() ); // i-th Suffix is influenced by an operation and needs recalculation
                    appendValue( dels, 0, Generous() );
                }
            }
        }

        last_lcp = fibre_lcp[0];

        for( size_t i = 1; i < len; ++i ){
          if( deleted_b( deletionset[ fibre_sa[i].i1 ], fibre_sa[i].i2 ) ){ //TODO: use std::binary_search instead
             appendValue( dels, i, Generous() ); //Is fibre_sa[i] a Pair?
             last_lcp = fibre_lcp[i];
             /*fibre_lcp[i - 1] = _min( fibre_lcp[i - 1], last_lcp ); //TODO: localize fibre_lcp[i-1] for performance*/
             continue;
          }

          //std::cout << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t-- " << get_shift( shiftset[ fibre_sa[i].i1 ], fibre_sa[i].i2 ) << " ->\t";
          fibre_sa[i].i2 += get_shift( shiftset[ fibre_sa[i].i1 ], fibre_sa[i].i2 ); //TODO: use std::binary_search instead
          //std::cout << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << std::endl;

          for( unsigned int j = 0; j <= _max( fibre_lcp[i], last_lcp ); ++j ){
             op = stringset[ fibre_sa[i].i1 ].getjournal().find_operation( fibre_sa[i].i2 + j );
             if( op->insertion() || op->deletion() && op->by() != 0 ){
                appendValue( indices, fibre_sa[i], Generous() ); // i-th Suffix is influenced by an operation and needs recalculation
                last_lcp = fibre_lcp[i];
                fibre_lcp[i - 1] = _min( fibre_lcp[i - 1], last_lcp );
                appendValue( dels, i, Generous() );
                break;
             }
          }
          last_lcp = fibre_lcp[i];
        }

        /*std::cout << "\nSuffix array after shifts:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }

        std::cout << "\nDels:" << std::endl;
        for( size_t i = 0; i < length( dels ); ++i ){
         std::cout << dels[i] << "\t| "<< fibre_sa[dels[i]].i1 << ", " << fibre_sa[dels[i]].i2 << "\t: " << suffix( stringset, fibre_sa[dels[i]] ) << "...\t" << fibre_lcp[dels[i]] << std::endl;
        }*/

        suffix_compare_functor< StringSet< TString, TSpec > > cmp( stringset );

        std::sort( begin( indices ), end( indices ), cmp );

        //indices = prefix( indices, std::unique( begin( indices ), end( indices ) ) );


        /*std::cout << "\nIndices:" << std::endl;
        for( size_t i = 0; i < length( indices ); ++i ){
         std::cout << indices[i].i1 << ", " << indices[i].i2 << "\t: " << suffix( stringset, indices[i] ) << "..." << std::endl;
        }*/

        // BEGIN: Debug Output
		/*std::cout << "Suffix array state before:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << "\t>" << i << "\t: " << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }*/
		// END: Debug Output

		typename Iterator< String< size_t > >::Type it_dels = begin( dels );
        typename Iterator< String< size_t > >::Type it_dels_end = end( dels );

        size_t pos_begin = *it_dels;
        size_t pos_end = *it_dels;
		size_t correctional_factor = 0;
		size_t local_min;

		// BEGIN: Deleting modified and removed suffixes block-wise
        while( it_dels != it_dels_end ){
            if( *(++it_dels) != ++pos_end ){
				//std::cout << "Deleting " << pos_begin << " to " << pos_end << " with cf " << correctional_factor << std::endl;
                erase( fibre_sa, pos_begin - correctional_factor, pos_end - correctional_factor );
				// TODO: std::min for intervals?
				local_min = fibre_lcp[pos_begin - correctional_factor - 1];
				for( size_t i = pos_begin - correctional_factor; i < pos_end - correctional_factor; ++i ){
					if( fibre_lcp[i] < local_min ){
						local_min = fibre_lcp[i];
					}
				}
				fibre_lcp[pos_begin - correctional_factor - 1] = local_min;
                erase( fibre_lcp, pos_begin - correctional_factor, pos_end - correctional_factor );
                correctional_factor += pos_end - pos_begin;
                pos_begin = *it_dels;
                pos_end = *it_dels;
            }
        }
		// END: Deleting modified and removed suffixes block-wise

		/*std::cout << "\nSuffix array after deletions:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }

        std::cout << "\nIndices:" << std::endl;
        for( size_t i = 0; i < length( indices ); ++i ){
         std::cout << indices[i].i1 << ", " << indices[i].i2 << std::endl; //"\t: " << infix( stringset, indices[i] ) << "..." << std::endl;
        }*/

        flatten( fibre_sa );
        flatten( fibre_lcp );

		typename Iterator< String< Pair< unsigned int, unsigned int > > >::Type it_index = begin( indices );
		typename Iterator< String< Pair< unsigned int, unsigned int > > >::Type it_index_tmp = begin( indices );
		typename Iterator< String< Pair< unsigned int, unsigned int > > >::Type it_index_end = end( indices );
        String< size_t > templcp;
        resize( templcp, length( indices ) );
        len = length( fibre_sa );
        typename Iterator< TString >::Type it_text = begin( stringset[(*it_index).i1] ) + (*it_index).i2;

		// BEGIN: Inserting modified indices block-wise
		for( size_t i = 0; i < len; ++i ){
            if( _suffix_bigger( begin( stringset[ fibre_sa[i].i1 ] ) + fibre_sa[i].i2, it_text ) ){
			    it_index_tmp = it_index;

			    while( _suffix_bigger( begin( stringset[ fibre_sa[i].i1 ] ) + fibre_sa[i].i2, begin( stringset[(*it_index).i1] ) + (*it_index).i2 ) && ++it_index < it_index_end ){}

                //std::cout << "Inserting into SA of length: " << length(fibre_sa) << " at i = " << i << " length: " << it_index - it_index_tmp << std::endl;

			    insert( i, fibre_sa, it_index_tmp, it_index - it_index_tmp );
			    insert( i, fibre_lcp, begin(templcp), it_index - it_index_tmp );

			    for( size_t j = i - 1; j < i + it_index - it_index_tmp; ++j ){
			        //std::cout << "\t< " << fibre_sa[j].i1 << ", " << fibre_sa[j].i2 << " >\t< " << fibre_sa[j+1].i1 << ", " << fibre_sa[j+1].i2 << " >" << std::endl;
                    fibre_lcp[j] = _lcp_length( begin( stringset[fibre_sa[j].i1] ) + fibre_sa[j].i2, begin( stringset[fibre_sa[j + 1].i1] ) + fibre_sa[j + 1].i2 );
                }

                if( it_index >= it_index_end ){
                    break;
                }

                /*std::cout << "\nSuffix array intermediate:" << std::endl;
                    for( size_t i = 0; i < length( fibre_sa ); ++i ){
                    std::cout << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
                }*/

                it_text = begin( stringset[(*it_index).i1] ) + (*it_index).i2;
                len = length( fibre_sa );
            }
        }

        if( it_index < it_index_end ){
            size_t i = length( fibre_sa ) - 1;
            append( fibre_sa, suffix( indices, it_index - begin( indices ) ) );
            append( fibre_lcp, suffix( templcp, it_index - begin( indices ) ) );
            len = length( fibre_sa ) - 1;
            while( i < len ){
                fibre_lcp[i] = lcpLength( suffix( stringset, fibre_sa[i] ), suffix( stringset, fibre_sa[i + 1] ) );
                ++i;
            }
            fibre_lcp[++len] = 0;
        }


		// END: Inserting modified indices block-wise

		// BEGIN: Debug Output
		/*std::cout << "Suffix array state after:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << "\t>" << i << "\t: " << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }*/
		// END: Debug Output

        /*
        resize( fibre_sa, len + remaining_indices );
        resize( fibre_lcp, len + remaining_indices );

        typename Iterator< TSA >::Type it_sa = end( fibre_sa ) - 1;
        typename Iterator< TLCP >::Type it_lcp = end( fibre_lcp ) - 1;
        *it_lcp = 0; //Last LCP Value is always Zero
        --len; //position of last data containing field

        //BEGIN: Loop Unrolling for last position
        while( std::binary_search( begin( dels ), end( dels ), len ) ){
            --len;
        }
        if( suffix( stringset, fibre_sa[ len ] ) <= suffix( stringset, *it_index ) ){
            *it_sa = *it_index;
            --it_index;
            --remaining_indices;
        }else{
            *it_sa = fibre_sa[len];
            *it_lcp = fibre_lcp[len];
        }
        --it_sa;
        --it_lcp;
        --len;
        //END: Loop Unrolling for last position

        bool recompute = false;

        for( int i = len; i >= 0; --i ){
          if( !std::binary_search( begin( dels ), end( dels ), i ) ){
             while( remaining_indices > 0 && suffix( stringset, fibre_sa[ i ] ) <= suffix( stringset, *it_index ) ){
                *it_sa = *it_index;
                *it_lcp = lcpLength( suffix( stringset, *it_sa ), suffix( stringset, *(it_sa + 1) ) );
                --it_index;
                --it_sa;
                --it_lcp;
                --remaining_indices;
                recompute = true;
             }
             *it_sa = fibre_sa[i];
             if( recompute ){
               *it_lcp = lcpLength( suffix( stringset, *it_sa ), suffix( stringset, *(it_sa + 1) ) );
               recompute = false;
             }else{
               *it_lcp = fibre_lcp[i];
             }
             --it_sa;
             --it_lcp;
          }else{
            *it_lcp = fibre_lcp[i+1];
          }
        }

        erase( fibre_sa, (size_t)0, length( dels ) );
        erase( fibre_lcp, (size_t)0, length( dels ) );

        */

        // Cleanup: flatten the Journals to write back to fibres
        /*
        for( unsigned int i = 0; i < length( stringset ); ++i ){
          flatten( stringset[i] );
        }
        */

        /*std::cout << "Suffix array:" << std::endl;
        for( size_t i = 0; i < length( fibre_sa ); ++i ){
         std::cout << "\t> " << i << "\t: " << fibre_sa[i].i1 << ", " << fibre_sa[i].i2 << "\t: " << suffix( stringset, fibre_sa[i] ) << "...\t" << fibre_lcp[i] << std::endl;
        }

        std::cout << stringset[0] << "(" << length( stringset[0] ) << ")" << std::endl;//<< stringset[1] << "(" << length( stringset[1] ) << ")"  << " " << stringset[2] << "(" << length( stringset[2] ) << ")"  << std::endl;

        std::cout << "The End!" << std::endl;*/
    }

   template< typename TIndex, typename TString >
   inline void synchronize_index( TIndex & index, TString & string ){
      typedef typename Position< TString >::Type TPos;
      String< Pair< TPos, TPos > > shifts;
      String< Pair< TPos, TPos > > deletions;
      String< TPos > indices;

      generate_shifts_and_deletions( string, indices, shifts, deletions );

      String< typename SAValue< TIndex >::Type, Journal< typename SAValue< TIndex >::Type, Alloc<>, Alloc<>, Sloppy > > fibre_sa( indexSA( index ) );

      indexRequire( index, ESA_LCP() );

      String< size_t, Journal< size_t, Alloc<>, Alloc<>, Sloppy > > fibre_lcp( indexLCP( index ) );

      unsigned int last_lcp = 0;
#ifndef NDEBUG
      std::cout << "Shifts: ";
      print_pairs( shifts );
      std::cout << "Deletions: ";
      print_pairs( deletions );
#endif
      if( deleted_b( deletions, fibre_sa[0] ) ){
         erase( fibre_sa, 0 );
         erase( fibre_lcp, 0 );
      }else{
         value( fibre_sa, 0 ) += get_shift_b( shifts, fibre_sa[0] );
         Operation* op;
         unsigned int num = _max( fibre_sa[0], last_lcp );
         for( unsigned int j = 0; j <= num; ++j ){
            op = string.getjournal().find_operation( fibre_sa[0] + j );
            if( op->insertion() || op->deletion() && op->by() != 0 ){
               indices += value( fibre_sa, 0 ); // i-th Suffix is influenced by an operation and needs recalculation
               erase( fibre_sa, 0 );
               erase( fibre_lcp, 0 );
            }
         }
      }

      for( unsigned int i = 1; i < length( fibre_sa ); ++i ){
#ifndef NDEBUG
         std::cout << "SA[" << i << "] == " << fibre_sa[i] << "\t";
#endif
         if( deleted_b( deletions, fibre_sa[i] ) ){
            erase( fibre_sa, i );
            erase( fibre_lcp, i );
            if( i > 0 ){
               value( fibre_lcp, i - 1 ) = _min( value( fibre_lcp, i ), last_lcp );
            }
            --i; //Now position i is former position i+1 so we need to step back 1 in the loop
#ifndef NDEBUG
            std::cout << "deleted!" << std::endl;
#endif
            continue;
         }
#ifndef NDEBUG
         std::cout << "shifted by: " << get_shift( shifts, fibre_sa[i] ) << "\t";
#endif
         value( fibre_sa, i ) += get_shift_b( shifts, fibre_sa[i] );
         //assert( get_shift( shifts, fibre_sa[i] ) == get_shift_b( shifts, fibre_sa[i] ) && "Shifts Binary not werkin!" );

         Operation* op;// = string.getjournal().find_operation( fibre_sa[i] );

         for( unsigned int j = 0; j <= _max( value( fibre_lcp, i ), last_lcp ); ++j ){
            op = string.getjournal().find_operation( value( fibre_sa, i ) + j );
            if( op->insertion() || op->deletion() && op->by() != 0 ){
               indices += value( fibre_sa, i ); // i-th Suffix is influenced by an operation and needs recalculation
               if( i > 0 ){
                  value( fibre_lcp, i - 1 ) = _min( value( fibre_lcp, i ), last_lcp );
               }
               erase( fibre_sa, i );
               erase( fibre_lcp, i );
               if( i > 0 ){
                  last_lcp = fibre_lcp[i];
               }
               --i; //Now position i is former position i+1 so we need to step back 1 in the loop
#ifndef NDEBUG
               std::cout << "influenced! Break!" << std::endl;
#endif
               break;
            }
         }
         last_lcp = value( fibre_lcp, i );
#ifndef NDEBUG
         std::cout << std::endl;
#endif
      }

      // Now sorting the new indices lexicographically

      suffix_compare_functor< TString > cmp( string );

      std::sort( begin( indices ), end( indices ), cmp );

      typename Iterator< String< size_t > >::Type it_indices_end = std::unique( begin( indices ), end( indices ) ); //remove duplicates; TODO: check where they come from and if not inserting by checking is more feasible

      indices = prefix( indices, it_indices_end );
#ifndef NDEBUG
      std::cout << "Sorted indices( SA and LCP):" << std::endl;
      for( unsigned int i = 0; i < length( indices ); ++i ){
         std::cout << "> " << indices[i] << " " << infix( string, indices[i], 12 ) << "\t" << ( ( i == length( indices ) - 1) ? 0 : lcpLength( suffix( string, indices[i] ), suffix( string, indices[i + 1] ) ) ) << std::endl;
      }
      fibre_sa.getjournal().print_sequence();
      std::cout << "##################################" << std::endl;
      fibre_lcp.getjournal().print_sequence();
      std::cout << "##################################" << std::endl;
#endif
      /****************************************************************************************************************

       State:
         - fibre_sa and fibre_lcp contain only suffix-indices that need no recomputation, or adding
         - deletions have been applied and the corresponding suffixes are removed
         - both fibres are still in sync with each other and in lexicographical order
         - indices contains the indices of suffixes that need adding or recomputation and is ordered lexicographically

       ****************************************************************************************************************/

      // Now merging both fibres and recomputing lcp-values

      Iterator< String< size_t> >::Type it_index = begin( indices );

      for( unsigned int j = 0; j < length( fibre_sa ); ++j ){
         if( suffix( string, fibre_sa[ j ] ) > suffix( string, *it_index ) ){

            insert( j, fibre_sa, it_index, 1 ); // insert this index at the corresponding position
            //recompute lcp value vor predecessor and successor
            insert( j, fibre_lcp, lcpLength( suffix( string, fibre_sa[ j ] ), suffix( string, fibre_sa[ j + 1 ] ) ) );
            value( fibre_lcp, j-1 ) = lcpLength( suffix( string, fibre_sa[ j - 1 ] ), suffix( string, fibre_sa[ j ] ) );

            if( ++it_index == end( indices ) ){
               break;      /* all indices were merged and no further processing is required
                              fibre_sa is in order and insync with fibre_lcp */
            }
         }
      }

      if( it_index != end( indices ) ){ //some indices left to be appended to fibre_sa
         value(fibre_lcp, length( fibre_sa ) - 1) = lcpLength( suffix( string, fibre_sa[ length( fibre_sa ) - 1 ] ), suffix( string, *it_index ) );
         append( fibre_sa, suffix( indices, it_index - begin( indices ) ) );
         String< size_t > lcp_values;
         while( it_index < end( indices ) - 1 ){
            lcp_values += lcpLength( suffix( string, *it_index ), suffix( string, *( it_index + 1 ) ) );
            ++it_index;
         }
         append( fibre_lcp, lcp_values );
         fibre_lcp += 0;
      }

      // Cleanup: flatten the Journals to write back to fibres
      flatten( fibre_sa );
      flatten( fibre_lcp );
      flatten( string );
      //Should be ready now!
   }

template< typename TIndex, typename TString >
   void synchronize_SA_LCP( TIndex & index, TString & string, std::vector<Entry>& entries ){

      typename Fibre< TIndex, ESA_SA >::Type &string_sa = indexSA( index );
      String< typename SAValue< TIndex >::Type, Journal< typename SAValue< TIndex >::Type, Alloc<>, Alloc<>, Sloppy > > fibre_sa( string_sa );

      if( !indexSupplied( index, ESA_LCP() ) ){
         indexRequire( index, ESA_LCP() );
      }

      typename Fibre< TIndex, ESA_LCP >::Type &fibre_lcp = indexLCP( index );
      String< size_t > indices;

      for( std::vector<Entry>::iterator it = entries.begin(); it != entries.end(); ++it ){
#ifndef NDEBUG
         it->print();
#endif
         unsigned int last_lcp = 0;
         for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            DEBUG_OUT( "Position: " + fibre_sa[i] );
            if( fibre_sa[i] >= it->position() ){
               if( (int)fibre_sa[i] < (int)it->position() - it->count() ){
                  erase( fibre_sa, i );
                  if( i > 0 ){
                     fibre_lcp[ i - 1 ] = _min( fibre_lcp[i], last_lcp );
                  }
                  erase( fibre_lcp, i );
                  --i;
                  DEBUG_OUT( " erased !" );
               }else{
                  value( fibre_sa, i ) += it->count();
                  DEBUG_OUT( " shifted to " << fibre_sa[i] << " !" );
               }
            }else{
               DEBUG_OUT( " kept !" );
            }
            last_lcp = fibre_lcp[i];
         }
      }

      for( std::vector<Entry>::iterator it = entries.begin(); it != entries.end(); ++it ){
         unsigned int last_lcp = 0;
         for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            if( fibre_sa[i] + _max( fibre_lcp[i], last_lcp ) >= it->position() - 1 && fibre_sa[i] < it->position() ){
               DEBUG_OUT( "Position: " << fibre_sa[i] << " needs recalculation!" );
               indices += fibre_sa[i];
               erase( fibre_sa, i );  // erase werks for journal, and others
               if( i > 0 ){
                  fibre_lcp[ i - 1 ] = _min( fibre_lcp[i], last_lcp );
               }
               erase( fibre_lcp, i );
               --i;
            }
            last_lcp = fibre_lcp[i];
         }

         if( it->insertion() ){
            for( unsigned int i = 0; i < (unsigned int)it->count(); ++i ){
               DEBUG_OUT( "Position: " << ( it->position() + i ) << " with string: " << infix( string, it->position() + i, it->position() + it->count() + i + 4 ) << "...\tneeds adding!");
               indices += it->position() + i;
            }
         }
      }

#ifndef NDEBUG
      std::cout << "> Sorting indices: ";
      for( unsigned int i = 0; i < length( indices ); ++i ){
         std::cout << indices[i] << " ";
      }
      std::cout << " !" << std::endl;
#endif
      suffix_compare_functor< TString > cmp( string );

      std::sort( begin( indices ), end( indices ), cmp );

#ifndef NDEBUG
      std::cout << "Sorted indices( SA and LCP):" << std::endl;
      for( unsigned int i = 0; i < length( indices ); ++i ){
         std::cout << "> " << indices[i] << " " << infix( string, indices[i], 12 ) << "\t" << ( ( i == length( indices ) - 1) ? 0 : lcpLength( suffix( string, indices[i] ), suffix( string, indices[i + 1] ) ) ) << std::endl;
      }
      printSA_LCP( index, string );
      std::cout << "##################################" << std::endl;
#endif

      Iterator< String< size_t> >::Type it_index = begin( indices );

      for( unsigned int j = 0; j < length( fibre_sa ); ++j ){
         DEBUG_OUT( "j: " << j );
         if( suffix( string, fibre_sa[ j ] ) > suffix( string, *it_index ) ){

            DEBUG_OUT( fibre_sa[ j ] << "\t" << suffix( string, fibre_sa[ j ] ) );
            DEBUG_OUT( *it_index << "\t" << suffix( string, *it_index ) << std::endl );

            insert( j, fibre_sa, it_index, 1 );

            resizeSpace( fibre_lcp, 1, j - 1, j - 1 );
            fibre_lcp[j - 1] = lcpLength( suffix( string, fibre_sa[ j - 1 ] ), suffix( string, fibre_sa[ j ] ) );
            fibre_lcp[j] = lcpLength( suffix( string, fibre_sa[ j ] ), suffix( string, fibre_sa[ j + 1 ] ) );
            if( ++it_index == end( indices ) ){
               break;
            }
         }
      }

      size_t last_pos = length( fibre_sa );
      appen( fibre_sa, suffix( indices, it_index - begin( indices ) ) );
      DEBUG_OUT( length( fibre_sa ) );
      Iterator< String< size_t> >::Type it_stop = end( indices ) - 1;
      while( it_index < it_stop ){
         fibre_lcp[ last_pos - 1 ] = lcpLength( suffix( string, fibre_sa[ last_pos - 1 ] ), suffix( string, fibre_sa[ last_pos ] ) );
         fibre_lcp += lcpLength( suffix( string, fibre_sa[ last_pos ] ), suffix( string, fibre_sa[ last_pos + 1 ] ) );
         ++it_index;
      }
      fibre_lcp += 0;

      flatten( fibre_sa );

      std::cout << std::endl << "-------------------------------------------------------" << std::endl;

      for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
         std::cout << fibre_sa[i] << "\t" << suffix( string, fibre_sa[i] ) << "$\t" << fibre_lcp[i] << std::endl;
      }

      std::cout << std::endl << "-------------------------------------------------------" << std::endl;

#ifndef NDEBUG
      std::cout << std::endl << "New Suffix Array and LCP-Table:" << std::endl;
      printSA_LCP( index, string );
      std::cout << std::endl << "-------------------------------------------------------" << std::endl;
#endif
   }

   template< typename TIndex, typename TString >
   void synchronize_ChildTab( TIndex & index, TString & string, std::vector<Entry>& entries ){
      if( !indexSupplied( index, ESA_ChildTab() ) ){
         return;
      }
      entries.begin();
      std::cout << string << "$" << std::endl;
   }

}
