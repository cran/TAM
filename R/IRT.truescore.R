
###########################################################################
# IRT.truescore
IRT.truescore <- function( object , iIndex = NULL , theta=NULL){
		irf <- IRT.irfprob(object)
		theta0 <- attr( irf , "theta" )[,1]
		dim_irf <- dim(irf)
		if ( is.null(iIndex) ){
			iIndex <- seq( 1 , dim_irf[1] )
						}
		irf <- irf[ iIndex ,, ]
		K <- dim_irf[2] - 1
		irf[ is.na(irf)] <- 0
		
		#*******
		# compute expected scores
		irf_K <- 0+0*irf
		for (kk in 1:K){
			irf_K[,kk+1,] <- kk
						}
		TP <- dim_irf[3]
		I <- dim(irf)[1]	
		vec <- rep(0,TP)
		for (ii in 1:I){
			# ii <- 1
			vec <- vec +  colSums( irf[ii,,] * irf_K[ii,,]  )
						}
		dfr <- data.frame( "theta" = theta0 , "truescore" = vec)
		#*******
		if ( ! is.null(theta) ){
			theta <- theta[ ( theta >= min(theta0) ) & ( theta <= max(theta0) ) ]
			TP1 <- length(theta)
			v1 <- 1
			v2 <- TP
			for (tt in 1:TP){
				v1 <- ifelse( theta > theta0[tt] , tt , v1 )	
				v2 <- ifelse( theta < theta0[TP-tt+1] , TP-tt+1 , v2 )
							}
			vec <- dfr$truescore[ v1 ] +  ( dfr$truescore[ v2 ] - dfr$truescore[ v1 ] ) *
						( theta - dfr$theta[ v1 ] ) / ( dfr$theta[ v2 ] - dfr$theta[ v1 ] )
			dfr <- data.frame( "theta" = theta , "truescore" = vec)						
							}
		return(dfr)						
			}
###########################################################################			