calc_exp_TK <-
function(rprobs , A , itemwt , p , ip , 
                        nitems , resp.ind){
  #... TK: multiple category option -> na.rm = TRUE    
  xbari <- sapply(ip, function(i) colSums(A[i,,p] * rprobs[i,,] , na.rm = TRUE))
  #... TK: multiple category option -> na.rm = TRUE
  xxfi <- sapply(ip, function(i) colSums(A[i,,p]^2 * rprobs[i,,] , na.rm = TRUE))
#*** ARb 2013-04-30: Write these loops more efficient in C.
# Note that outside of the loop the operations are
#        xbar[p] <- sum(res$xbar)
#        xxf[p] <- sum(res$xxf)
#        xbar2[p] <- sum( res$xbar2 )
# Therefore, one can simply use the sum score of xbar, xxf
# and xbar2 inside the C function.
  xbar <- rowSums(xbari * itemwt[,ip]) 
  xxf <- rowSums(xxfi * itemwt[,ip])
  xbar2 <- rowSums(xbari^2 * itemwt[,ip])
  
  return(list("xbar" = xbar, "xxf" = xxf, "xbar2" = xbar2))
}
# This function is not used anymore.
###########################################################

###########################################################
# faster version of calc_exp_TP3
calc_exp_TK3 <- function( rprobs , A , np , est.xsi.index , itemwt ,
	indexIP.no , indexIP.list2 ){
	CC <- dim(rprobs)[2]
	TP <- dim(rprobs)[3]
	NXSI <- dim(A)[3]
	NI <- dim(A)[1]
	# restructure rprobs and AL
	rprobsL <- matrix( rprobs , NI*CC , TP )
	AL <- matrix( A , NI*CC , NXSI )	
	AL[ is.na(AL) ] <- 0
	rprobsL[ is.na(rprobsL) ] <- 0
	# call Rcpp code
#	res <- tam.calcexp( np , rprobsL , AL ,	indexIP.no , indexIP.list2 , 
#			est.xsi.index , CC , itemwt)
	res <- .Call( "TAM_CALCEXP" , np , rprobsL , AL ,	indexIP.no , 
			indexIP.list2 , est.xsi.index , CC , itemwt , PACKAGE="TAM" )
	return(res)
		}
		
##------		
###########################################################
calc_posterior.v2 <-
  function(rprobs , gwt , resp , nitems , 
           resp.ind.list , normalization = TRUE , 
           thetasamp.density = NULL , snodes = 0 , resp.ind=NULL){   

# a0 <- Sys.time()		   
    if ( snodes == 0 ){ 
      fx <- gwt  
    } else {
      # calculate individual 'sampling weight'
      swt <- fx <- gwt / outer( rep(1,nrow(gwt)) , thetasamp.density )
	  # This is essentially equal to one.
#    	swt <-fx <- gwt
    } 
    nstud <- nrow(fx)
    # using c Code here
    storage.mode(resp) <- "integer"
# cat("vor calcfx") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
    fx <- .Call("calcfx", fx, rprobs, resp.ind.list, resp)
# cat("nach calcfx") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
    # numerical integration
    if ( snodes == 0 ){ 
      rfx <- rowSums(fx)
      if (normalization ){
        hwt <- fx / rfx } else {   hwt <- fx }
    }
    # Monte Carlo integration
    if ( snodes > 0 ){ 
#      rfx <- rowMeans(fx)
		rfx <- rowSums(fx)		
      if (normalization ){
		 hwt <- fx / rfx 	
			} else { hwt <- fx }
    }
    res <-  list("hwt" = hwt , "rfx" = rfx )
    if ( snodes > 0 ){ 
		
		res[["swt" ]] <- fx
				}
# cat(" in  posterior rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	    
    return(res)
  }
#####################################################################

#####################################################################
# calc_prob
# Calculation of probabilities
calc_prob.v5 <-
  function(iIndex, A, AXsi, B, xsi, theta, 
           nnodes, maxK, recalc=TRUE){
    
    if(recalc){
      AXsi.tmp <- array( tensor( A[iIndex,,, drop = FALSE], xsi, 3, 1 ) , 
                         dim = c( length(iIndex) , maxK , nnodes ) )
      AXsi[iIndex,] = AXsi.tmp[,,1]
    } else {
      AXsi.tmp <- array( AXsi, dim = c( length(iIndex) , maxK , nnodes ) )
    }
    
    Btheta <- array(0, dim = c(length(iIndex) , maxK , nnodes) )
    for( dd in 1:ncol(theta) ) 
      Btheta <- Btheta + array(B[iIndex,,dd ,drop = FALSE] %o% theta[,dd] , dim = dim(Btheta))
    
    rprobs <- ( rr <- exp(Btheta+AXsi.tmp) )/aperm( array( rep( colSums( aperm( rr , c(2,1,3) ) ,
					dims=1 , na.rm = TRUE) ,                                                           maxK ), dim=dim(rr)[c(1,3,2)] ) , c(1,3,2) )
    return(list("rprobs" = rprobs, "AXsi" = AXsi))
  }

#############################################################
#############################################################
stud_prior.v2 <-
  function(theta , Y , beta , variance , nstud , 
           nnodes , ndim , YSD ){
    if(ndim == 1) {
      ##################################
      # SINGLE DIMENSION
      gwt <- matrix(dnorm(rep(theta, each = nstud), mean= Y%*%beta, sd = sqrt(variance)),
                    nrow = nstud)
    } else {
      ###################################
      # MULTIPLE DIMENSIONS
      mu <- Y%*%beta     #mean vector for each student: dimensions nstud by ndim 
      eps <- 10^(-7)
      #	variance[ diag(variance) ] <- diag(variance) + eps
      diag(variance) <- diag(variance) + eps
      varInverse <- solve(variance)
      coeff <- 1/sqrt( (2*pi)^ndim * det(variance) ) 
      gwt <- matrix( 0 , nrow=nstud , ncol=nnodes )  
      ###*****
	  ### ARb 2013-08-27
	  ### different calculations depend if there are
	  ### person-specific predictors or not  
# a0 <- Sys.time()
##      for ( qq in 1:nnodes ) {
##        x1 <- - mu + theta[rep(qq,nstud),]  #x has dimension nstud#
##        x <- matrix( rowSums( (x1%*%varInverse) * x1 ) , ncol= 1)
##        gwt[,qq] <- coeff*exp(-0.5*x) 		
##					}
#  cat(" * prior Ysd") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    if ( YSD ){
		gwt <- prior.normal.density.R( theta_=theta , mu_=mu , 
			     varInverse_=varInverse ,  coeff_=coeff) 
			   }	
#  cat(" * prior nnodes2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
    if ( ! YSD){
		gwt <- prior.normal.densityALL.R( theta_=theta , mu_=mu , 
			     varInverse_=varInverse ,  coeff_=coeff) 
		gwt <- matrix( gwt , nrow=nstud , ncol=nnodes , byrow=TRUE )
				}
#  cat(" * prior nnodes3") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
    }
    return(gwt)
  }
#####################################################################  
  
#*************************
# auxiliary functions for calculation of prior functions  
prior.normal.density.R <- 
function( theta_ , mu_ , varInverse_ , coeff_){
	.Call("prior_normal_density_C",  theta_ , mu_ , 
			varInverse_ , coeff_ , PACKAGE = "TAM")
} 

prior.normal.densityALL.R <- 
function( theta_ , mu_ , varInverse_ , coeff_){
	.Call("prior_normal_densityALL_C",  theta_ , mu_ , 
		varInverse_  , coeff_ , PACKAGE = "TAM")
} 
  
###########################################################################  
tam.jml.xsi <-
  function ( resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
             ItemScore, theta, xsi, Msteps, pweightsM,
             est.xsi.index
  ){
    
    #Update item parameters
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    r <- matrix(0,nrow=nitems,ncol=maxK) 
    rr <- array(0,dim=c(nitems,maxK,maxK))
    AA <- array (0, dim=c(nitems,maxK,maxK))
    
    maxChangeP <- 0
    errorP <- rep(0, max(est.xsi.index))
    convergeAllP <- FALSE
    p_loop <- est.xsi.index
    convergeP <- rep(FALSE,max(est.xsi.index))
    
    # begin loop 
    iterP <- 1
    old_increment <- rep(5,max(p_loop))
    cat(" Item parameter estimation |")
    while (!convergeAllP & ( iterP <= Msteps ) ) {  
      res.p <- calc_prob.v5( iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , TRUE )        	
      rprobs <- res.p[["rprobs"]]               
      
      #compute probability weights, summed over students, so that there is no cycling through students for parameter estimation (p loop)
      for (k1 in 1:maxK) {
        r[,k1] <- colSums(t(rprobs[,k1,]) * resp.ind * pweightsM, na.rm=TRUE)
        for (k2 in 1:maxK) {
          rr[,k1,k2] <- colSums(t(rprobs[,k1,]*rprobs[,k2,]) * resp.ind * pweightsM, na.rm=TRUE)
        }
      }
      
      for (p in p_loop ) {
        A_bari <- rowSums(A[,,p] * r, na.rm=TRUE)
        AA_bari <- rowSums(A[,,p]^2 * r, na.rm=TRUE)
        for (k1 in 1:maxK) {
          for (k2 in 1:maxK) {
            AA[,k1,k2] <- A[,k1,p]*A[,k2,p] * rr[,k1,k2]
          }
        }
        A_Sq <- apply(AA, 1, sum, na.rm=TRUE)
        
        
        expected <- sum (A_bari, na.rm=TRUE) # sum over items 
        err <- sum (AA_bari - A_Sq, na.rm=TRUE)   #sum over the items  
        err_inv <- abs (1/err)
        scores <- ItemScore[p] - expected 
        increment <-  err_inv*scores
        if (maxChangeP < abs(increment)) {
          maxChangeP <- abs(increment)
        }
        
        while (abs(increment) > abs(old_increment[p])){
          increment <- increment / 2.0  #damping the increment
        }
        
        xsi[p] <- xsi[p] + increment
        old_increment[p] <- increment
        errorP[p] <- sqrt(err_inv)
        if ( max(abs(increment)) < convM ) {
          convergeP[p] <- TRUE
        }
        #      cat( paste( "Iteration in parameter estimation ", iterP, " parameter ", p, "\n" )  ) 
        flush.console()
      }  # parameter p 
      iterP <- iterP + 1 
      p_loop <- est.xsi.index[convergeP[est.xsi.index]==FALSE]
      convergeAllP <- (sum(convergeP[est.xsi.index]) == length(est.xsi.index))  
      cat("-")  
    } # end of all parameters convergence
    
    res <- list( "xsi" = xsi , "errorP" = errorP, "maxChangeP" = maxChangeP)
    return (res)  
  }

######################################################################  
######################################################################
  
tam.jml.xsi2 <-
  function ( resp , resp.ind, A, A.0 , B, nstud, nitems, maxK, convM, 
             ItemScore, theta, xsi, Msteps, pweightsM,
             est.xsi.index , rp3 , rp3.sel , rp3.pweightsM 
  ){
    
    
    #Update item parameters
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    r <- matrix(0,nrow=nitems,ncol=maxK) 
    rr <- array(0,dim=c(nitems,maxK,maxK))
    AA <- array (0, dim=c(nitems,maxK,maxK))
    
    maxChangeP <- 0
    errorP <- rep(0, max(est.xsi.index))
    convergeAllP <- FALSE
    p_loop <- est.xsi.index
    PP <- length(p_loop)
    convergeP <- rep(FALSE,max(est.xsi.index))
    old_xsi <- xsi
    # begin loop 
    iterP <- 1
    old_increment <- rep(5,max(p_loop))
    cat(" Item parameter estimation |")
    while (!convergeAllP & ( iterP <= Msteps ) ) {
      
      res.p <- calc_prob.v5( iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta[ rp3.sel$caseid ,,drop=FALSE ] , 
                             nrow(rp3.sel) , maxK , TRUE )      		
      rprobs <- res.p[["rprobs"]]               
      #compute probability weights, summed over students, so that there is no cycling through students for parameter estimation (p loop)
      for (k1 in 1:maxK) {
        r[,k1] <- colSums(t(rprobs[,k1,]) * resp.ind[ rp3.sel$caseid , ] * rp3.pweightsM, na.rm=TRUE)
        for (k2 in 1:maxK) {
          rr[,k1,k2] <- colSums(t(rprobs[,k1,]*rprobs[,k2,]) * resp.ind[ rp3.sel$caseid , ] * rp3.pweightsM, na.rm=TRUE)
        }
      }
      
      A_Sq <- AA_bari <- A_bari <- matrix( 0 , PP , nitems )
      
      for (kk in 1:maxK){ 
        A_bari <- A_bari + t( A.0[ , kk , ] * r[ , kk ] )
        AA_bari <- AA_bari + t( A.0[ , kk , ]^2 * r[ , kk ] )		
      }
      for (kk1 in 1:maxK){ 
        for (kk2 in 1:maxK){ 
          A_Sq <- A_Sq + t( A.0[,kk1,] * A.0[,kk2,] * rr[ , kk1 , kk2 ] )	
        }
      }
      # A							[ nitems , maxK , length(xsi) ]	
      # A_Sq, AA_bari, A_bari		[ length(xsi) , nitems ]
      # r							[ nitems , maxK ]
      # rr						[ nitems , maxK , maxK ]	
      #    for (p in p_loop ) {  # begin p loop
      #      A_bari[p,] <- rowSums(A[,,p] * r, na.rm=TRUE)
      #      AA_bari[p,] <- rowSums(A[,,p]^2 * r, na.rm=TRUE)
      #      for (k1 in 1:maxK) {
      # for (k2 in 1:maxK) {
      # AA[,k1,k2] <- A[,k1,p]*A[,k2,p] * rr[,k1,k2]
      # }
      # }
      # A_Sq[p,] <- apply(AA, 1, sum, na.rm=TRUE)
      # }	# end parameter p loop
      
      expected <- rowSums (A_bari, na.rm=TRUE) # sum over items
      err <- rowSums(AA_bari - A_Sq, na.rm=TRUE)   #sum over the items
      
      err_inv <- abs (1/( abs(err) + 10^(-10) ))
      scores <- ItemScore * ( ! convergeP ) - expected
      increment <-  err_inv*scores
      ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
      increment <- ifelse( abs( increment) > abs(old_increment)  , 
                           increment/(2*ci) , 
                           increment )
      xsi <- xsi + increment
      old_increment <- increment      
      errorP <- sqrt(err_inv)
      convergeP[ abs(increment) < convM ] <- TRUE
      flush.console()
      iterP <- iterP + 1 
      p_loop <- est.xsi.index[convergeP[est.xsi.index]==FALSE]
      convergeAllP <- (sum(convergeP[est.xsi.index]) == length(est.xsi.index))  
      cat("-")  
    } # end of all parameters convergence
    
    res <- list( "xsi" = xsi , "errorP" = errorP, 
                 "maxChangeP" = max(abs( xsi - old_xsi ) ) )
    return (res)  
  }
