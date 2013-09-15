
// includes from the plugin
#include <RcppArmadillo.h>
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
SEXP theta_sq_cpp( SEXP theta_) ;
}

// definition

SEXP theta_sq_cpp( SEXP theta_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix theta(theta_) ;  
       
     int N = theta.nrow() ;  
     int D = theta.ncol() ;   
       
     arma::cube thetasq(N,D,D) ;  
       
     ////////////////////////////////////////  
     // calculation of squared theta matrix  
       
     // int nn = 0 ;  
     for (int nn=0;nn<N;nn++){  
         for (int dd1=0;dd1<D;dd1++){  
         thetasq(nn,dd1, dd1 ) = pow( theta(nn,dd1) , 2 ) ;   
         for (int dd2=dd1+1;dd2<D;dd2++){  
             thetasq(nn,dd1, dd2 ) = theta(nn,dd1) * theta(nn,dd2 ) ;   
             thetasq(nn,dd2 , dd1 ) = thetasq(nn,dd1,dd2) ;  
             }      
         }  
     }  
       
     //// OUTPUT  
     return ( wrap(thetasq) ) ;
END_RCPP
}




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





