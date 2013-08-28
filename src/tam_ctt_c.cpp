

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
SEXP tam_ctt_C( SEXP tdat, SEXP wle, SEXP maxK, SEXP est_wle) ;
}

// definition

SEXP tam_ctt_C( SEXP tdat, SEXP wle, SEXP maxK, SEXP est_wle ){
BEGIN_RCPP
  
       
       
     Rcpp::CharacterMatrix TDAT(tdat);          
     Rcpp::NumericVector WLE (wle);  
     Rcpp::NumericVector MAXK (maxK);  
     int EST_WLE=as<int>(est_wle) ;  
       
     int N = TDAT.ncol() ;  
     int I = TDAT.nrow() ;  
     int K = MAXK[0];  
       
     Rcpp::NumericMatrix des( I*K , 9 ) ;  
     Rcpp::CharacterVector desV( I*K ) ;  
       
     int start_cc = 0 ;  
     int cc1 ;  
       
     for ( int ii=0; ii<I ; ii++){  
       
     // int ii = 0;  // select item ii  
       
     double wles1 = 0 ;  
     double wles2 = 0 ;  
     double wles3 = 0 ;                  
     double wles4 = 0 ;   
     double freq = 0 ;                  
     double freq2 = 0 ;   
     double mw = 0 ;  
       
     Rcpp::CharacterVector TDAT_ii = TDAT.row(ii) ;  
     Rcpp::CharacterVector uii = Rcpp::unique( TDAT_ii ) ;  
       
     // Rcpp::IntegerVector uii = table( TDAT_ii ) ;  
     int NC_ii = uii.size() ;  
       
     //Rcpp::CharacterVector categii = uii.attr("names") ;  
     Rcpp::CharacterVector categii = uii ;  
       
     for ( int cc=0; cc < NC_ii ; cc++ ){  
            cc1 = cc + start_cc ;      
            des( cc1  , 0 ) = ii+1 ;          // item number  
            desV[cc1 ] = categii[cc] ;    // category label   
            wles1 = 0 ;  
            wles2 = 0 ;  
            wles3 = 0 ;  
            wles4=0;      
            freq = 0 ;                  
            freq2 = 0 ;                  
       
         for (int nn=0; nn < N ; nn++){                                      
             std::string cx = Rcpp::as<std::string>(TDAT_ii[nn]);  
             std::string cy = Rcpp::as<std::string>(categii[cc]);        
             if ( cx != "NA" ){      
                 if ( cx == cy){  
                    if ( EST_WLE==1){                   
                     wles1 = wles1 + WLE[nn] ;      
                     wles2 = wles2 + WLE[nn]*WLE[nn] ;                  
                             }  
                     freq ++ ;  
                         } else {   
                    if ( EST_WLE==1){                                           
                     wles3 = wles3 + WLE[nn]  ;  
                     wles4 = wles4 + WLE[nn]*WLE[nn] ;                          
                    }  
                     freq2 ++ ;  
                         }  
                     }  
                 }  // end case nn              
         des(cc1,3) = freq ;    // frequency of students at category cc  
         des(cc1,2) = freq2 ;   // frequency of students not at category cc       
         des(cc1,4) = wles1 ;    // score sum of students at category cc    
         des(cc1,5) = wles3 ;    // score sum of students not at category cc                      
         des(cc1,6) = wles2 ;    // sum of squares WLE at category cc       
         // calculate N  
         des(cc1,1) = des(cc1,2) + des(cc1,3) ;  
         if ( EST_WLE==1){                               
         // calculate WLE mean total              
         mw = ( wles1 + wles3 ) / des( cc1,1) ;  
         // calculate SD total              
         des(cc1,8) = sqrt( ( wles2 + wles4 - des(cc1,1)*pow( mw , 2 )  )/ ( des(cc1,1) - 1 )  ) ;              
         // calculate WLE means              
         des(cc1,4) = des(cc1,4) / des(cc1,3) ;  // M at category cc1  
         des(cc1,5) = des(cc1,5) / des(cc1,2) ;  // M not at category cc1  
         // calculate SD of WLE          
         des(cc1,6) = sqrt( ( des(cc1,6) - des(cc1,3)*pow( des(cc1,4) , 2 )  ) / ( des(cc1,3) - 1 ) ) ;     
         // calculate point-biserial correlation              
         des(cc1,7) =  ( des(cc1,4) - des(cc1,5) )/des(cc1,8) * sqrt( des(cc1,2)*des(cc1,3) /   
                         ( des(cc1,1) * ( des(cc1,1)-1 ) ) ) ;   
         }          
                     }   // end category cc  
         start_cc = start_cc + NC_ii ;  
                 }  
       
     ///////////////////////////////////////////  
     // OUTPUT:  
     return Rcpp::List::create(  
            _["des"] = des ,  _["desV"] = desV   
            		) ;  
       
       
               
               
     //  Rcpp::Rcout << "cx " << cx  << std::endl ;                                       
     //   Rcpp::Rcout << "cy " << cy  << std::endl ;                                       
             
END_RCPP
}



