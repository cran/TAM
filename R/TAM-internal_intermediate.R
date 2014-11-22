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
           thetasamp.density = NULL , snodes = 0 , resp.ind=NULL,
		   avoid.zerosum=FALSE){   

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
	
	if (avoid.zerosum ){	
		fxs <- rowSums( fx )
		m1 <- min( fxs[ fxs > 0 ] ) / 1E3 / ncol(fx )
		fx[ fxs == 0 , ] <- m1
		# fx <- fx + m1
		}
	
# cat("nach calcfx") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
    # numerical integration
    if ( snodes == 0 ){ 
#      rfx <- rowSums(fx)
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
# print( rowSums(hwt))	
#print( cbind( rfx , rowSums(hwt ))[1:100,] )
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

#cat("-----")	
#print(AXsi.tmp[,,1])
    
    Btheta <- array(0, dim = c(length(iIndex) , maxK , nnodes) )
    for( dd in 1:ncol(theta) ) 
      Btheta <- Btheta + array(B[iIndex,,dd ,drop = FALSE] %o% theta[,dd] , dim = dim(Btheta))
    
    rprobs <- ( rr <- exp(Btheta+AXsi.tmp) )/aperm( array( rep( colSums( aperm( rr , c(2,1,3) ) ,
					dims=1 , na.rm = TRUE) ,    maxK ), dim=dim(rr)[c(1,3,2)] ) , c(1,3,2) )
#cat("****")	
#rprobs[ is.na(rprobs) ] <- 0
#print(round(rprobs,3))
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
	  # eps <- 1E-3
	  
	    
	  
	# @@@ ARb 2014-10-19: Stabilization of the covariance matrix
		svd_var <- svd(variance)
		d0 <- d <- svd_var$d
		eps2 <- .05
		ind <- which( d < eps2)
		if (length(ind)>0){
			d[ ind ] <- eps2
			d <- d / sum(d) * sum(d0)
			variance <- svd_var$u %*% diag(d) %*% t( svd_var$v  )
							}


      #	variance[ diag(variance) ] <- diag(variance) + eps
#      diag(variance) <- diag(variance) + eps
      varInverse <- solve(variance)
		detvar <- det(variance)

      coeff <- 1/sqrt( (2*pi)^ndim * detvar ) 
	  # coeff <- 1 
      gwt <- matrix( 0 , nrow=nstud , ncol=nnodes )  
	  

	  
      ###*****
	  ### ARb 2013-08-27
	  ### different calculations depend if there are
	  ### person-specific predictors or not  
# a0 <- Sys.time()
# manuell <- TRUE
#if ( manuell ){
#      for ( qq in 1:nnodes ) {
#        x1 <- - mu + theta[rep(qq,nstud),]  #x has dimension nstud#
#        x <- matrix( rowSums( (x1%*%varInverse) * x1 ) , ncol= 1)
#        gwt[,qq] <- coeff*exp(-0.5*x) 		
#					}
# gwt0 <- gwt
#						}						
						
#  cat(" * prior Ysd") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    if ( YSD ){
		gwt <- prior.normal.density.R( theta_=theta , mu_=mu , 
			     varInverse_=varInverse ,  coeff_=coeff) 
			   }	
#  cat(" * prior nnodes2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						
    if ( ! YSD){
		gwt <- prior.normal.densityALL.R( theta_=theta , mu_=mu , 
			     varInverse_=varInverse ,  coeff_=coeff) 
# gwt <- mvtnorm::dmvnorm( theta , mean=mu[1,1:2] , sigma=variance )	  
		# gtw <- gwt / sum(gwt) 
		gwt <- matrix( gwt , nrow=nstud , ncol=nnodes , byrow=TRUE )
				}
#	 gwt <- gwt / matrix( rowSums(gwt) , nrow=nstud , ncol=nnodes , byrow=TRUE )			
				
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
