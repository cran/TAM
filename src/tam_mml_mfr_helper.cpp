

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
SEXP gresp_extend( SEXP gresp_, SEXP xstep_) ;
}

// definition

SEXP gresp_extend( SEXP gresp_, SEXP xstep_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix gresp(gresp_);          
     Rcpp::NumericVector xstep(xstep_);        
       
     int I=gresp.ncol() ;  
     int N=gresp.nrow();  
       
     Rcpp::NumericMatrix gresp2(N,I);  
       
     for (int ii=0;ii<I;ii++){  
     for (int nn=0;nn<N;nn++){  
     	if (! R_IsNA(gresp(nn,ii)) ) {  
     		if ( gresp(nn,ii) == xstep[ii] ){  
     			gresp2(nn,ii) = 1 ;   
     				}  
     				} else {  
     		gresp2(nn,ii) = NA_REAL ;  
     				}  
     		}  
     	}  
     	  
     	  
     return( wrap(gresp2) ) ;  
       
        			  
     //*************************************************      
     // OUTPUT              
                   
      //return Rcpp::List::create(    
     //    _["d1b" ] = d1b ,  
     //    _["d2b" ] = d2b   
     //    ) ;  
END_RCPP
}




// declarations
extern "C" {
SEXP gresp_na_facets( SEXP gresp_, SEXP rnfacets_, SEXP rnx_) ;
}

// definition

SEXP gresp_na_facets( SEXP gresp_, SEXP rnfacets_, SEXP rnx_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix gresp(gresp_);          
     Rcpp::CharacterVector rnfacets(rnfacets_);        
     Rcpp::CharacterVector rnx(rnx_);  
              
     int I=gresp.ncol() ;  
     int N=gresp.nrow();  
       
     Rcpp::NumericMatrix gresp2 = gresp;  
       
     for (int ii=0;ii<I;ii++){  
     for (int nn=0;nn<N;nn++){  
     	if ( rnfacets[nn] != rnx[ii]   ) {  
     		gresp2(nn,ii) = NA_REAL ;  
     				}  
     		}  
     	}  
     	  
     	  
     return( wrap(gresp2) ) ;  
       
        			  
     //*************************************************      
     // OUTPUT              
                   
      //return Rcpp::List::create(    
     //    _["d1b" ] = d1b ,  
     //    _["d2b" ] = d2b   
     //    ) ;  
END_RCPP
}




