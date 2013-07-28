//  Code created: 2013-07-27 20:01:13


// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP rowCumsums2_source( SEXP matr) ;
}

// definition

//# The C code was posted by Romain Francois at
//# http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-October/001198.html

SEXP rowCumsums2_source( SEXP matr ){
BEGIN_RCPP
     NumericMatrix input( matr ) ;  
          NumericMatrix output  = clone<NumericMatrix>( input ) ;  
       
          int nr = input.nrow(), nc = input.ncol() ;  
          NumericVector tmp( nr );  
          for( int i=0; i<nc; i++){  
              tmp = tmp + input.column(i) ;  
              NumericMatrix::Column target( output, i ) ;  
              std::copy( tmp.begin(), tmp.end(), target.begin() ) ;  
          }  
          return output ;
END_RCPP
}



