


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
SEXP TAM_CALCEXP( SEXP np, SEXP rprobsL, SEXP AL, SEXP indexIPno, SEXP indexIPlist2, SEXP estxsiindex, SEXP CC, SEXP itemwt) ;
}

// definition

SEXP TAM_CALCEXP( SEXP np, SEXP rprobsL, SEXP AL, SEXP indexIPno, SEXP indexIPlist2, SEXP estxsiindex, SEXP CC, SEXP itemwt ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
       
     int NP = as<int>(np);  
     Rcpp::NumericMatrix rprobs(rprobsL);  
     Rcpp::NumericMatrix A(AL);  
     Rcpp::NumericMatrix INDEXIPNO(indexIPno);  
     Rcpp::NumericVector INDEXIPLIST2(indexIPlist2);  
     Rcpp::NumericVector ESTXSIINDEX(estxsiindex);  
     int C = as<int>(CC);  
     Rcpp::NumericMatrix ITEMWT(itemwt);  
     // Rcpp::Numeric C(CC); => Rcpp does not have a class Numeric!  
       
     ////////////////////////////////////////////////////////////  
     // define output vectors  
     NumericVector XBAR (NP) ;  
     NumericVector XBAR2 (NP) ;  
     NumericVector XXF (NP) ;  
       
     ///////////////////////////////////////////////////////////  
     // DEFINE indices  
       
     int TP = rprobs.ncol();  
     int NEXI = ESTXSIINDEX.size() ;  
     int NR = rprobs.nrow();  
     int II = NR / C ;   
       
     /////////////////////////////////////////////////////////  
     // CALCULATIONS  
       
     // loop over xsi item parameters  
     for (int hh=0;hh<NEXI;++hh){  
       
     double pp = ESTXSIINDEX[hh] - 1 ;	  
     double ll = 0 ; // xbar element  
     double yy = 0 ; // xbar2 element  
     double zz = 0 ; // xxf element  
       
     // loop over theta points  
     for (int tt=0;tt<TP;++tt){  
     // loop over items within an item parameter group  
     double GG = INDEXIPNO(pp,1) - INDEXIPNO(pp,0) + 1 ;  
     for (int gg=0;gg<GG;++gg){  
        double ii=INDEXIPLIST2[ INDEXIPNO(pp,0)-1 + gg ]-1 ;	  
        // loop over categories cc = 1, ... , C  
        double vv1 = 0 ;  // xbar counter  
        double vv2 = 0 ;  // xxf counter  
        for (int cc=0;cc<C;++cc){     
     	// xbar calculation  
     	vv1 += A( ii+cc*II , pp) * rprobs( ii + cc*II , tt ) ;  
     	vv2 += A( ii+cc*II , pp) * A( ii+cc*II , pp) * rprobs( ii + cc*II , tt ) ;	  
     		}  
        ll += vv1*ITEMWT(tt,ii) ; // xbar addition  
        yy += vv1*vv1*ITEMWT(tt,ii) ; // xbar2 addition  
        zz += vv2*ITEMWT(tt,ii) ; // xxf addition     
        		}  
     	}	  
     XBAR[pp] = ll ;  
     XBAR2[pp] = yy ;  
     XXF[pp] = zz ;  
     	} // end xsi index  
       
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
       
     return List::create(_["xbar"]=XBAR , _["xbar2"]=XBAR2 , _["xxf"]=XXF  ); 
END_RCPP
}



