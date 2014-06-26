

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
SEXP mml3_calc_Fdes( SEXP Fdes_, SEXP dimFdes_) ;
}

// definition

SEXP mml3_calc_Fdes( SEXP Fdes_, SEXP dimFdes_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericVector XDES(Fdes_);          
     Rcpp::NumericVector dimXdes(dimFdes_);        
       
     // $dimXdes  
     // [1]  6  4 21 19  
       
     int I= dimXdes[0] ;  
     int maxK= dimXdes[1] ;  
     int TP= dimXdes[2] ;  
     int Nlam = dimXdes[3] ;  
     int RR = XDES.size() ;  
       
     Rcpp::NumericMatrix XdesM(RR,5) ;  
       
     int rr = 0 ;  
       
     int ind = 0 ;  
       
     for (int ii=0; ii<I;ii++){  
     for (int kk=0; kk <maxK ; kk++){  
     for ( int tt=0; tt<TP; tt++ ){   
     for ( int ll=0; ll<Nlam ; ll++ ){  
       
        rr=ii+I*kk+I*maxK*tt+I*maxK*TP*ll ;  
       
     // Rcpp::Rcout << "rr = " << rr << " XDES[rr] = " << XDES[rr] << std::endl;  
       
     if ( XDES[rr] != 0 ){  
          XdesM(ind,0) = ii ;	  
          XdesM(ind,1) = kk ;  
          XdesM(ind,2) = tt ;  
          XdesM(ind,3) = ll ;  
          XdesM(ind,4) = XDES[rr] ;  
          ind = ind + 1 ;  
          		}  
          	}  
     }  
     }  
     }  
       
       
     //*************************************************      
     // OUTPUT              
                   
     return Rcpp::List::create(    
         _["NFdesM"] = ind  ,  
         _["FdesM" ] = XdesM  
         ) ;  
END_RCPP
}



// declarations
extern "C" {
SEXP mml3_slca_deriv( SEXP FdesM_, SEXP dimFdes_, SEXP gammaslope_, SEXP probs_, 
	SEXP nik_, SEXP Nik_, SEXP guess_, SEXP probs0_) ;
}

// definition

SEXP mml3_slca_deriv( SEXP FdesM_, SEXP dimFdes_, SEXP gammaslope_, SEXP probs_, 
	SEXP nik_, SEXP Nik_, SEXP guess_, SEXP probs0_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix XdesM(FdesM_);          
     Rcpp::NumericVector dimXdes(dimFdes_);        
     Rcpp::NumericVector Xlambda(gammaslope_);  
     Rcpp::NumericVector probs(probs_) ;  
     Rcpp::NumericVector nik(nik_) ;  
     Rcpp::NumericVector Nik(Nik_) ;  
     Rcpp::NumericVector guess(guess_) ;  
     Rcpp::NumericVector probs0(probs0_) ;  
       
     // $dimXdes  
     // [1]  6  4 21 19  
       
     int I= dimXdes[0] ;  
     int maxK= dimXdes[1] ;  
     int TP= dimXdes[2] ;  
     int Nlam = dimXdes[3] ;  
       
     int RR = XdesM.nrow() ;  
       
     // int PP = I * maxK * TP ;  
       
     Rcpp::NumericVector d1b(Nlam);  
     Rcpp::NumericVector d2b(Nlam);  
       
     int ii=0;  
     int hh=0;  
     int tt=0;  
     int ll=0;  
       
     //  # XdesM     [ii,kk,tt,ll, value ]   
     //        # probs  num [1:I, 1:maxK, 1:TP]  
     //        # n.ik  num [1:I, 1:maxK, 1:TP]  
     //        # N.ik  num [1:I,1:TP]  
     //        # Xdes  num [1:I, 1:maxK, 1:TP, 1:Nlam]           
       
       
     ///*********************  
     // First derivative  
       
     //  for (hh in 1:maxK){  
     //      t1 <- sum( Xdes[ , hh , , ll] * ( n.ik[ , hh , ] - probs[,hh,] * N.ik ) )  
     //       d1.b[ll] <- d1.b[ll] + t1  
     //			}  
     for (int rr = 0 ; rr <RR ; rr++){  
     	// # XdesM     [ii,kk,tt,ll, value ]   
     	ii = XdesM(rr,0);  // item ii  
     	hh = XdesM(rr,1);  // category hh  
     	tt = XdesM(rr,2);  // theta grid point tt  
     	ll = XdesM(rr,3);  // gamma parameter ll    
     	  
     	//*** no guessing parameter  
     	if ( guess[ii] == 0 ){	  
     	  d1b[ll] += XdesM(rr,4) * ( nik[ii+I*hh+I*maxK*tt] -  
     	      probs[ii+I*hh+I*maxK*tt] * Nik[ ii+I*tt ] ) ;  
     	  		  } // end if guess[ii] = 0  
     	//*** with guessing parameter	  		    
     	if ( guess[ii] > 0 ){  
     	   d1b[ll] += XdesM(rr,4) * probs0[ii+I*hh+I*maxK*tt] / probs[ii+I*hh+I*maxK*tt] *   
     	      ( nik[ii+I*hh+I*maxK*tt] - probs[ii+I*hh+I*maxK*tt] * Nik[ ii+I*tt ] ) ;	  		    
     		}  // end guess[ii] > 0  
     		  
     		}  
       
       
     ///*********************  
     // Second derivative  
       
     int NS = I*TP*Nlam ;  
     Rcpp::NumericVector tmp1(NS) ;  
     // tmp1 <- 0  
     // for (hh in 1:maxK ){  
     //   tmp1 <- tmp1 + probs[,hh,] * Xdes[,hh,,ll]  
     //		}			  
       
     // parameter ll; item i ; class tt  
     int vv=0;  
     for (int rr=0;rr<RR;rr++){  
        ii = XdesM(rr,0);  
        hh = XdesM(rr,1);  
        tt = XdesM(rr,2);  
        ll = XdesM(rr,3);  
        vv = ii + I*tt + I*TP*ll ;  
        tmp1[vv] += XdesM(rr,4) * probs[ii+I*hh+I*maxK*tt] ;  
        			}  
       
     // for (hh in 1:maxK){  
     //  t2 <- sum( Xdes[ , hh , , ll] * N.ik * probs[,hh,] *  
     //   ( Xdes[, hh , , ll ] - tmp1 ) )  
     // d2.b[ll] <- d2.b[ll] + t2  
     //		}  
       
     // int rr = 0 ;  
     for (int rr=0;rr<RR;rr++){  
        ii = XdesM(rr,0);  
        hh = XdesM(rr,1);  
        tt = XdesM(rr,2);  
        ll = XdesM(rr,3);  
        vv = ii + I*tt + I*TP*ll ;  
     //  t2 <- sum( Xdes[ , hh , , ll] * N.ik * probs[,hh,] *  
     //   ( Xdes[, hh , , ll ] - tmp1 ) )     
       if ( guess[ii] == 0 ){	  
            d2b[ll] += XdesM(rr,4) * Nik[ii + I*tt ] * probs[ ii + I*hh + I*maxK*tt ] *   
       	              ( XdesM(rr,4) - tmp1[vv] ) ;  
       		}  // end if guess[ii] = 0  
       		  
       if ( guess[ii] > 0 ){	  
            d2b[ll] += pow(XdesM(rr,4),2.0) * probs0[ ii + I*hh + I*maxK*tt ] *  
                    probs0[ ii + I*0 + I*maxK*tt ] * ( guess[ii] * nik[ii+I*hh+I*maxK*tt] /  
                    	     pow( probs[ ii + I*hh + I*maxK*tt ] , 2.0) - Nik[ii + I*tt ] );  
       		}  // end if guess[ii] > 0  		  
       		  
       		  
        			}  
       
        			  
     //*************************************************      
     // OUTPUT              
                   
      return Rcpp::List::create(    
         _["d1b" ] = d1b ,  
         _["d2b" ] = d2b   
         ) ;  
END_RCPP
}



// declarations
extern "C" {
SEXP mml3pl_tam_calcexp( SEXP np, SEXP rprobsL, SEXP AL, SEXP indexIPno, 
	SEXP indexIPlist2, SEXP estxsiindex, SEXP CC, SEXP itemwt, SEXP rprobsL0, 
	SEXP guess, SEXP nik_, SEXP ni_) ;
}

// definition

SEXP mml3pl_tam_calcexp( SEXP np, SEXP rprobsL, SEXP AL, SEXP indexIPno, 
	SEXP indexIPlist2, SEXP estxsiindex, SEXP CC, SEXP itemwt, SEXP rprobsL0, 
	SEXP guess, SEXP nik_, SEXP ni_ ){
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
          Rcpp::NumericMatrix rprobs0(rprobsL0);      
          Rcpp::NumericVector GUESS(guess) ;  
          Rcpp::NumericVector nik(nik_) ;  
          Rcpp::NumericVector ni(ni_) ;  
            
          // Rcpp::Numeric C(CC); => Rcpp does not have a class Numeric!    
              
          ////////////////////////////////////////////////////////////    
          // define output vectors    
          NumericVector XBAR (NP) ;  
          NumericVector iscore (NP) ;  
          NumericVector XBAR2 (NP) ;    
          NumericVector XXF (NP) ;    
              
          ///////////////////////////////////////////////////////////    
          // DEFINE indices    
              
          int TP = rprobs.ncol();    
          int NEXI = ESTXSIINDEX.size() ;    
          int NR = rprobs.nrow();    
          int II = NR / C ;     
            
          // II ... number of items  
          // TP ... number of theta points  
          // CC ... number of categories  
            
          /////////////////////////////////////////////////////////    
          // CALCULATIONS    
              
          // loop over xsi item parameters    
     for (int hh=0;hh<NEXI;++hh){    
              
       double pp = ESTXSIINDEX[hh] - 1 ;	    
       double ll = 0 ; // xbar element    
       double yy = 0 ; // xbar2 element    
       double zz = 0 ; // xxf element    
       double ww = 0 ; // iscore element  
         
       // loop over theta points    
       for (int tt=0;tt<TP;++tt){    
       // loop over items within an item parameter group    
       double GG = INDEXIPNO(pp,1) - INDEXIPNO(pp,0) + 1 ;    
         for (int gg=0;gg<GG;++gg){    
          double ii=INDEXIPLIST2[ INDEXIPNO(pp,0)-1 + gg ]-1 ;	    
          // loop over categories cc = 1, ... , C    
          double vv1 = 0 ;  // xbar counter    
          double vv2 = 0 ;  // xxf counter  
          double vv3 = 0 ;  
          for (int cc=0;cc<C;++cc){  
             if ( GUESS[ii] == 0 ){	  
          	     // xbar calculation    
               vv1 += A( ii+cc*II , pp ) * rprobs( ii + cc*II , tt )  ;  
               vv2 += A( ii+cc*II , pp) * A( ii+cc*II , pp) * rprobs( ii + cc*II , tt ) ;  
               vv3 += A( ii+cc*II , pp) * nik[ ii +cc*II + tt*II*C ] ;  
          	                  }   // end guess[ii] == 0  
             if ( ( GUESS[ii] > 0 ) & ( cc == 1 ) ){	  
          	     // xbar calculation    
               vv1 += A( ii+cc*II , pp ) * rprobs0( ii + cc*II , tt )  ;  
     //          vv2 += A( ii+cc*II , pp) * A( ii+cc*II , pp) * rprobs( ii + cc*II , tt ) ;  
               vv2 += A( ii+cc*II , pp) * A( ii+cc*II , pp) * rprobs0( ii + cc*II , tt ) *            
                         rprobs0( ii , tt ) * ( GUESS[ii] / pow(rprobs( ii + cc*II , tt ),2) *  
                		nik[ ii +cc*II + tt*II*C ] - ni[ii+tt*II] ) ;	  
               vv3 += A( ii+cc*II , pp) * nik[ ii +cc*II + tt*II*C ] *   
               		rprobs0( ii + cc*II , tt ) / rprobs( ii + cc*II , tt ) ;  
          	                  }   // end guess[ii] > 0     	                   	                       	                    
          	                }  // end categories cc  
         if ( GUESS[ii] == 0 ){  
             ll += vv1*ni[ii+tt*II] ;  
             yy += vv1*vv1*ni[ii+tt*II] ; // xbar2 addition    
             zz += -vv2*ni[ii+tt*II] ; // xxf addition  
             ww += vv3 ;  
             		} // end guess ii  
         if ( GUESS[ii] > 0 ){  
             ll += vv1*ni[ii+tt*II] ;  
             zz += vv2 ;  
             ww += vv3 ;  
             		} // end guess ii        		  
           		}   // end item parameter group gg    
          }    // end theta grid tt  
        XBAR[pp] = ll ;    
        iscore[pp] = ww ;  
        XBAR2[pp] = yy ;    
        XXF[pp] = zz ;    
        } // end xsi index    
              
          ///////////////////////////////////////////////////////    
          ///////////// O U T P U T   ///////////////////////////    
              
          return List::create(  
          	     	_["xbar"]=XBAR , _["xbar2"]=XBAR2 , _["xxf"]=XXF ,  
          	     	_["iscore"] = iscore  
          	     	);   
       
       
     
END_RCPP
}






