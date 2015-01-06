


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
