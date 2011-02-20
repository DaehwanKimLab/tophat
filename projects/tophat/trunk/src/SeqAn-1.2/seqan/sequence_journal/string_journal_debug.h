namespace seqan {
    template< typename TIndex, typename TString >
    void printSA( TIndex & index, TString & string ){

        typename Fibre< TIndex, ESA_SA >::Type &fibre_sa = indexSA( index );

        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << fibre_sa[i] << "\t" << suffix( string, fibre_sa[i] ) << "$" << std::endl;
        }

    }

    template< typename TIndex, typename TString >
    void printSA_LCP( TIndex & index, TString & string ){

        typename Fibre< TIndex, ESA_SA >::Type &fibre_sa = indexSA( index );

        if( !indexSupplied( index, ESA_LCP() ) ){
            indexRequire( index, ESA_LCP() );
        }

        typename Fibre< TIndex, ESA_LCP >::Type &fibre_lcp = indexLCP( index );

        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << fibre_sa[i] << "\t" << suffix( string, fibre_sa[i] ) << "$\t" << fibre_lcp[i] << std::endl;
        }

    }

}
