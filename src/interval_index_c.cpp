

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
SEXP interval_index_C( SEXP matr, SEXP rn) ;
}

// definition

SEXP interval_index_C( SEXP matr, SEXP rn ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     Rcpp::NumericVector RN(rn) ;  
     	// rn random number for plausible value imputation  
       
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
       
     // create output vectors  
     NumericVector IND (NR) ;  
     IND.fill(0);  
       
     for (int nn=0;nn<NR;++nn){  
      	for (int cc=0 ; cc < NC ; ++cc ){  
     	    if ( MATR(nn,cc) > RN[nn] ){  
     	    	    IND(nn) = cc + 1 ;  
     	    	    break ;   
     	    	           }  
     		}  
     	}  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     return( wrap(IND) );  
     // return List::create(_["maxval"] = MAXVAL , _["maxind"]=MAXIND ) ;     
     
END_RCPP
}



