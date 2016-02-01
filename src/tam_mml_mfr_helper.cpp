

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



// declarations
extern "C" {
SEXP a_matrix_cumsum( SEXP index_matr_, SEXP mm_, SEXP SG_) ;
}

// definition

SEXP a_matrix_cumsum( SEXP index_matr_, SEXP mm_, SEXP SG_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix index_matr(index_matr_);  
     Rcpp::NumericMatrix mm(mm_);          
     int SG=as<int>(SG_);  
       
     int NR = mm.nrow();  
     int NR1 = NR + SG;  
     int NC = mm.ncol() ;  
       
     Rcpp::NumericMatrix cumsum_mm(NR1,NC);  
       
     double tmp=0;  
     int ss =0;  
     int rr=0;  
       
     for (int cc=0; cc<NC ; cc++){  
     	ss = 2*SG;  
     	rr=0;	  
     	for( int zz=0; zz < NR ; zz++){   
     		if ( index_matr(zz,0) != ss ){  
     			rr ++ ;  
     			tmp = 0 ;  
     			      }  
     		tmp = tmp + mm( index_matr(zz,1) , cc ) ;  
     		cumsum_mm( rr , cc ) = tmp ;  
     		ss = index_matr(zz,0) ;  
     		rr ++ ;  
     				}	  
     		}  
       
     //*************************************************      
     // OUTPUT              
                   
      return Rcpp::List::create(    
         Rcpp::_["index_matr"] = index_matr ,  
         Rcpp::_["SG"] = SG ,  
         Rcpp::_["cumsum_mm"] = cumsum_mm  
         ) ;  
END_RCPP
}



// declarations
extern "C" {
SEXP colsums_gresp( SEXP gresp_) ;
}

// definition

SEXP colsums_gresp( SEXP gresp_ ){
BEGIN_RCPP  
       
     Rcpp::NumericMatrix gresp(gresp_);                
              
     int NR = gresp.nrow();  
     int NC = gresp.ncol();  
       
     Rcpp::NumericVector sumnull(NC);  
       
     for (int cc=0;cc<NC;cc++){  
     sumnull[cc]=1;  
     	for (int nn=0;nn<NR;nn++){  
     	if ( ! R_IsNA( gresp(nn,cc) ) ){  
     	     if ( gresp(nn,cc) > 0 ){  
     		     sumnull[cc] = 0 ;  
     		     break;  
     				}  
     			}  
     		}  
     	}  
       
       
     return( wrap(sumnull) );		  
     		  
     //*************************************************      
     // OUTPUT              
                   
     // return Rcpp::List::create(    
     //    _["gresp"] = gresp  
     //    ) ;  
END_RCPP
}





