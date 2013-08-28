

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


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//*******************************************************************
// normal density with different person means
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// declarations
extern "C" {
SEXP prior_normal_density_C( SEXP theta_, SEXP mu_, SEXP varInverse_, SEXP coeff_) ;
}

// definition

SEXP prior_normal_density_C( SEXP theta_, SEXP mu_, SEXP varInverse_, SEXP coeff_ ){
BEGIN_RCPP
  
       
       
     Rcpp::NumericMatrix theta(theta_);    
     Rcpp::NumericMatrix mu(mu_);          
     Rcpp::NumericMatrix varInverse(varInverse_);  
     Rcpp::NumericVector COEFF(coeff_);  
       
     int nnodes = theta.nrow() ;  
     int ndim = theta.ncol() ;  
     int nstud = mu.nrow() ;  
       
     double coeff = COEFF[0] ;  
       
     Rcpp::NumericMatrix gwt(nstud,nnodes) ;  
       
       
     Rcpp::NumericVector x1(ndim) ;  
       
     //*****************************************  
     // R code to be simplified  
     //      for ( qq in 1:nnodes ) {  
     //        x1 <- - mu + theta[rep(qq,nstud),]  #x has dimension nstud#  
     //        x <- matrix( rowSums( (x1%*%varInverse) * x1 ) , ncol= 1)  
     //        gwt[,qq] <- coeff*exp(-0.5*x) 		  
     //					}  
       
     // int nn=0 ;  
     // int qq=0 ;  
       
     for (int nn=0;nn<nstud;nn++){  
     for (int qq=0;qq<nnodes;qq++){  
       
         for (int dd=0;dd<ndim;dd++){  
             x1[dd] = theta(qq,dd) - mu(nn,dd) ;  
                 }  
       
         for (int dd1=0;dd1<ndim;dd1++){  // beg dd1  
             gwt(nn,qq) += x1[dd1]*x1[dd1] * varInverse(dd1,dd1) ;  
             for (int dd2=(dd1+1);dd2<ndim;dd2++){ // beg dd2  
                     gwt(nn,qq) += 2*x1[dd1]*x1[dd2] * varInverse(dd1,dd2) ;  
             }   // end dd2  
         }  // end dd1  
       
         gwt(nn,qq) = coeff * exp( -0.5*gwt(nn,qq) ) ;  
         } // end qq  
     } // end nn  
               
     //// OUTPUT  
     return ( wrap(gwt) ) ;  
         
END_RCPP
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//*******************************************************************
// normal density applied for all subjects
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


// declarations
extern "C" {
SEXP prior_normal_densityALL_C( SEXP theta_, SEXP mu_, SEXP varInverse_, SEXP coeff_) ;
}

// definition

SEXP prior_normal_densityALL_C( SEXP theta_, SEXP mu_, SEXP varInverse_, SEXP coeff_ ){
BEGIN_RCPP
  
       
       
     Rcpp::NumericMatrix theta(theta_);    
     Rcpp::NumericMatrix mu(mu_);          
     Rcpp::NumericMatrix varInverse(varInverse_);  
     Rcpp::NumericVector COEFF(coeff_);  
       
     int nnodes = theta.nrow() ;  
     int ndim = theta.ncol() ;  
     // int nstud = mu.nrow() ;  
       
     double coeff = COEFF[0] ;  
       
     Rcpp::NumericVector gwt(nnodes) ;  
       
       
     Rcpp::NumericVector x1(ndim) ;  
       
     //*****************************************  
     // R code to be simplified  
     //      for ( qq in 1:nnodes ) {  
     //        x1 <- - mu + theta[rep(qq,nstud),]  #x has dimension nstud#  
     //        x <- matrix( rowSums( (x1%*%varInverse) * x1 ) , ncol= 1)  
     //        gwt[,qq] <- coeff*exp(-0.5*x) 		  
     //					}  
       
     int nn=0 ;  
     // int qq=0 ;  
       
     for (int qq=0;qq<nnodes;qq++){  
       
         for (int dd=0;dd<ndim;dd++){  
             x1[dd] = theta(qq,dd) - mu(nn,dd) ;  
                 }  
       
         for (int dd1=0;dd1<ndim;dd1++){  // beg dd1  
             gwt[qq] += x1[dd1]*x1[dd1] * varInverse(dd1,dd1) ;  
             for (int dd2=(dd1+1);dd2<ndim;dd2++){ // beg dd2  
                     gwt[qq] += 2*x1[dd1]*x1[dd2] * varInverse(dd1,dd2) ;  
             }   // end dd2  
         }  // end dd1  
       
         gwt[qq] = coeff * exp( -0.5*gwt[qq] ) ;  
         } // end qq  
               
     //// OUTPUT  
     return ( wrap(gwt) ) ;  
         
END_RCPP
}






