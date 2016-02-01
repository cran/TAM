//  Code created: 2014-03-19 09:11:51


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
SEXP tam_wle_Bs( SEXP rprobsWLEL, SEXP respind, SEXP BL, SEXP BBL, SEXP BBBL, SEXP ndim, SEXP nitems, SEXP maxK, SEXP nstud) ;
}

// definition

SEXP tam_wle_Bs( SEXP rprobsWLEL, SEXP respind, SEXP BL, SEXP BBL, SEXP BBBL, SEXP ndim, SEXP nitems, SEXP maxK, SEXP nstud){
BEGIN_RCPP
	    /////////////////////////////////////    
	    // INPUT    
		  Rcpp::NumericMatrix RPROBS(rprobsWLEL);
	    Rcpp::NumericMatrix RESPIND(respind);
	    Rcpp::NumericMatrix CBL(BL);
	    Rcpp::NumericMatrix CBB(BBL);
	    Rcpp::NumericMatrix CBBB(BBBL);
	    
	    ///////////////////////////////////////////////////////////
	    // INPUT indices
	    int cndim = as<int>(ndim);
	    int cnitems = as<int>(nitems);
	    int cmaxK = as<int>(maxK);
	    int cnstud = as<int>(nstud);
	    int citstud = cnitems*cnstud;
	    
	    ////////////////////////////////////////////////////////////
	    // define output vectors
	    NumericMatrix B_bari (citstud, cndim);
	    NumericMatrix BB_bari (citstud, cndim*cndim);
	    NumericMatrix BBB_bari (citstud, cndim);
	    
	    NumericMatrix B_Sq (citstud, cndim*cndim);
	    NumericMatrix B2_B (citstud, cndim);
	    NumericMatrix B_Cube (citstud, cndim);
	    
	    /////////////////////////////////////////////////////////
	    // CALCULATIONS
	    
	    for(int ii=0; ii<cnitems; ii++){// item loop
	    	for(int jj=0; jj<cnstud; jj++){// student loop
	    		for(int dd1=0; dd1<cndim; dd1++){// dimension loop 1
	    		
	    			B_bari( cnstud*ii+jj , dd1 ) = 0;
	    			BBB_bari( cnstud*ii+jj , dd1 ) = 0;
	    			for(int cc=0; cc<cmaxK; cc++){// category loop
	    				B_bari( cnstud*ii+jj , dd1 ) += CBL( cnitems*cc+ii , dd1 )*RPROBS( cnitems*cc+ii , jj )*RESPIND( jj , ii );
	    				BBB_bari( cnstud*ii+jj , dd1 ) += CBBB( cnitems*cc+ii , dd1 )*RPROBS( cnitems*cc+ii , jj )*RESPIND( jj , ii );
	    			}
	    			
	    			B2_B( cnstud*ii+jj , dd1 ) = 0;
	    			B_Cube( cnstud*ii+jj , dd1 ) = 0;
	    			for(int dd2=0; dd2<cndim; dd2++){// category loop
	    				BB_bari(cnstud*ii+jj , cndim*dd1+dd2 )=0;
	    				for(int cc=0; cc<cmaxK; cc++){// category loop
	    					BB_bari( cnstud*ii+jj , cndim*dd2+dd1 ) += CBB( cnitems*cc+ii , cndim*dd2+dd1 )*RPROBS( cnitems*cc+ii , jj )*RESPIND( jj , ii );
	    				}
	    				
	    				B_Sq( cnstud*ii+jj , cndim*dd2+dd1 ) = B_bari( cnstud*ii+jj , dd1 )*B_bari( cnstud*ii+jj , dd2 );
	    				B2_B( cnstud*ii+jj , dd1 ) += BB_bari( cnstud*ii+jj , cndim*dd2+dd1 )*B_bari( cnstud*ii+jj , dd2 );
	    				B_Cube( cnstud*ii+jj , dd1 ) += B_Sq( cnstud*ii+jj , cndim*dd2+dd1 )*B_bari( cnstud*ii+jj , dd2 );
	    				
	    			}
	    		}
	    	}
	    }
	    
	    ///////////////////////////////////////////////////////
	    ///////////// O U T P U T   ///////////////////////////
	    
	    return Rcpp::List::create(
                    Rcpp::_["B_bari"]=B_bari, 
                    Rcpp::_["BB_bari"]=BB_bari, 
                    Rcpp::_["BBB_bari"]=BBB_bari,
                    Rcpp::_["B_Sq"]=B_Sq, 
                    Rcpp::_["B2_B"]=B2_B, 
                    Rcpp::_["B_Cube"]=B_Cube );
END_RCPP
}

// declarations
extern "C" {
	SEXP tam_wle_errinv( SEXP errl, SEXP ndim, SEXP nstud);
}

// definition
SEXP tam_wle_errinv( SEXP errl, SEXP ndim, SEXP nstud ){
BEGIN_RCPP
	/////////////////////////////////////
	// INPUT
	Rcpp::NumericMatrix myERR(errl);
	
	///////////////////////////////////////////////////////////
	// INPUT indices
	int cndim = as<int>(ndim);
	int cnstud = as<int>(nstud);
	
	////////////////////////////////////////////////////////////
	// define output vectors
	arma::mat ERR_j = arma::zeros(cndim, cndim);
	arma::mat ERR_j_inv;
	NumericMatrix ERR_inv (cnstud, cndim*cndim);
	
	/////////////////////////////////////////////////////////
	// CALCULATIONS
	for(int jj=0; jj<cnstud; jj++){// item loop
		for(int dd1=0; dd1<cndim; dd1++){// dimension loop 1
			for(int dd2=0; dd2<cndim; dd2++){// dimension loop 2
				ERR_j(dd1,dd2) = myERR(jj, dd1+dd2*cndim);
			}
		}
		
		ERR_j_inv = arma::mat(inv(ERR_j));
		
		for(int dd1=0; dd1<cndim; dd1++){// dimension loop 1
			for(int dd2=0; dd2<cndim; dd2++){// dimension loop 2
				ERR_inv(jj, dd1+dd2*cndim) = ERR_j_inv(dd1, dd2);
			}
		}
	}
	
	///////////////////////////////////////////////////////
	///////////// O U T P U T   ///////////////////////////
	return ( wrap(ERR_inv) );
END_RCPP
}
