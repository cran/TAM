#########################################################
# calculation of B matrix
.mml.3pl.computeB <- function( E , gammaslope ){
       dimE <- dim(E)
	   I <- dimE[1]
	   maxK <- dimE[2]
	   D <- dimE[3]
	   B <- array( 0 , dim= c(I , maxK , D ) )
	   for (dd in 1:D){
       for (cc in 1:maxK){
			B[ , cc,dd] <- E[,cc,dd,] %*% gammaslope
						}
					}
		return(B)
			}
#############################################################

#####################################################################
# calc_prob
# Calculation of probabilities
.mml.3pl.calc_prob.v5 <-
  function(iIndex, A, AXsi, B, xsi, theta, 
           nnodes, maxK, recalc=TRUE , guess){
		   
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
					dims=1 , na.rm = TRUE) ,    maxK ), dim=dim(rr)[c(1,3,2)] ) , c(1,3,2) )
					
					
	# include guessing	
	rprobs0 <- rprobs
	ind <- which(guess > 0 )
	if ( length(ind) > 0 ){
		rprobs[ ind , 2 , ] <- guess[ind] + ( 1-guess[ind] ) * rprobs0[ind,2,]	
		# include guessing here
		rprobs[ ind , 1 , ] <- 1 - rprobs[ ind , 2 , ]
							}
    return(list("rprobs" = rprobs, "AXsi" = AXsi , "rprobs0"=rprobs0 ))
  }



###########################################################################
# reduced skillspace estimation
.mml.3pl.skillspace <- function( Ngroup, pi.k , 
			delta.designmatrix , G , delta , delta.fixed ,			
			eps=10^(-10) ){		
	Z <- delta.designmatrix	
	delta0 <- delta
	ND <- nrow(delta)
	covdelta <- list(1:G)
	for (gg in 1:G){
		ntheta1 <- Ngroup[gg] * pi.k[,gg]
		ntheta1 <- ntheta1 / sum(ntheta1 )		
		lntheta <- log(ntheta1+eps)
		mod <- lm( lntheta ~ 0 + Z , weights = ntheta1 )
		covbeta <- vcov(mod)		
		beta <- coef(mod)
		if ( ! is.null( delta.fixed ) ){
		# delta.fixed: 1st column: parameter index
		#              2nd column: group index
		#              3rd column: parameter value 
		    ind.gg <- which( delta.fixed[ ,2] == gg )
			if ( length(ind.gg) > 0 ){
				beta[ delta.fixed[ind.gg,1] ] <- delta.fixed[ind.gg,3]
									}
									}
		pi.k[,gg] <- exp( Z %*% beta ) / Ngroup[gg]
		pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
		delta[,gg] <- beta
		covdelta[[gg]] <- covbeta
					}
	res <- list( "pi.k"=pi.k , "delta"=delta , 
			"covdelta" = covdelta )			
			}
##########################################################################
# calculation of expected counts
.mml.3pl.expected.counts <- function( datindw , nitems , maxK , ntheta , hwt){
			# calculate expected counts
			n.ik <- array( 0 , dim=c(nitems , maxK , ntheta ) )
			N.ik <- array( 0 , dim=c( nitems , ntheta ) )
			for (kk in 1:maxK){   # kk <- 1
				dkk <- datindw[[kk]]
				g1 <- crossprod( dkk , hwt ) 
				n.ik[,kk,] <- g1
				N.ik <- N.ik + g1
						}
			res <- list("n.ik"=n.ik , "N.ik" = N.ik )
			return(res)
			}
#####################################################################


#######################################################################
# moments of distribution
.mml.3pl.distributionmoments <- function( D , G , pi.k , theta.k ){
	mean.trait <- sd.trait <- skewness.trait <- matrix( 0 , nrow=D , ncol=G )		
	for (dd in 1:D){
		for (gg in 1:G){
		mean.trait[dd,gg] <- sum( theta.k[,dd] * pi.k[ , gg ] )
		sd.trait[dd,gg] <- sqrt( sum( theta.k[,dd]^2 * pi.k[ , gg ] ) - mean.trait[dd,gg]^2 )
		skewness.trait[dd,gg] <- sum( ( theta.k[,dd] - mean.trait[dd,gg] )^3 * pi.k[ , gg ] ) /
						sd.trait[dd,gg]^3
					}
			}
	rownames(skewness.trait) <- rownames(sd.trait) <- 
				rownames(mean.trait) <- colnames(theta.k)
	colnames(skewness.trait) <- colnames(sd.trait) <- 
				colnames(mean.trait) <- paste0("Group",1:G)
	#*****
	# correlation matrices
	correlation.trait <- as.list(1:G)
	names(correlation.trait) <- colnames(mean.trait)
	for (gg in 1:G){
		# gg <- 1
		mean.gg <- rep(0,D)
		Sigma.gg <- diag(0,D)	
		for (dd in 1:D){
			mean.gg[dd] <- sum( pi.k[,gg] * theta.k[,dd] )
				}
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
		#		dd1 <- 1 ; 	dd2 <- 1
				Sigma.gg[dd1,dd2] <- sum( pi.k[,gg] * (theta.k[,dd1] - mean.gg[dd1] )*(theta.k[,dd2] - mean.gg[dd2] ) ) 
#				Sigma.gg[dd1,dd2] <- Sigma.gg[dd1,dd2] - mean.gg[dd1] * mean.gg[dd2]
				Sigma.gg[dd2,dd1] <- Sigma.gg[dd1,dd2]
									}
						}
		rownames(Sigma.gg) <- colnames(Sigma.gg) <- rownames(mean.trait)
		correlation.trait[[gg]] <- cov2cor(Sigma.gg + diag(10^(-5),D) )
					}	
	res <- list( "mean.trait"=mean.trait , "sd.trait" = sd.trait , 
				"skewness.trait" = skewness.trait , "correlation.trait"=correlation.trait)
    return(res)				
		}