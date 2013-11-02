##################################
# Information criteria
.TAM.ic <- function( nstud , deviance , xsi , xsi.fixed ,
	beta , beta.fixed , ndim , variance.fixed , G , irtmodel ,
	B_orig=NULL , B.fixed , E , est.variance , resp ,
		est.slopegroups=NULL ){
	# 2PL estimation
	# c("2PL","GPCM","GPCM.design","2PL.groups") )	
  #***Model parameters
  ic <- data.frame("n" = nstud , "deviance" = deviance )
  dev <- deviance
	# xsi parameters
	ic$Nparsxsi <- length(xsi)
	if ( ! is.null( xsi.fixed) ){ 
			ic$Nparsxsi <- ic$Nparsxsi - nrow(xsi.fixed ) }
	# B slopes
	ic$NparsB <- 0
	if ( irtmodel == "2PL" ){
		ic$NparsB <- sum( B_orig != 0 )
						    }
	if ( irtmodel == "GPCM" ){
		ic$NparsB <- ncol(resp)
						    }
	if ( irtmodel == "GPCM.design" ){
		ic$NparsB <- ncol(E)
						    }							
	if ( irtmodel == "2PL.groups" ){
#		ic$NparsB <- sum( B_orig != 0 )
		ic$NparsB <- length( unique( est.slopegroups ) )
		# This is not yet correct for multiple dimensions and multiple
		# categories
						    }								
	if ( ! is.null( B.fixed ) ){
		ic$NparsB <- ic$NparsB - nrow(B.fixed )
					}
	
	# beta regression parameters
	ic$Nparsbeta <- dim(beta)[1] * dim(beta)[2]
	if ( ! is.null( beta.fixed) ){ 
			ic$Nparsbeta <- ic$Nparsbeta - nrow(beta.fixed ) }
	# variance/covariance matrix
	ic$Nparscov <- ndim + ndim*(ndim-1)/2
	if ( ! est.variance ){ ic$Nparscov <- ic$Nparscov - ndim }
	if ( ! is.null( variance.fixed) ){ 
			ic$Nparscov <- max(0 , ic$Nparscov - nrow(variance.fixed ) )
									 }
	# total number of parameters
	ic$Npars <- ic$np <- ic$Nparsxsi + ic$NparsB + ic$Nparsbeta + ic$Nparscov
    	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
		# adjusted BIC 
		ic$aBIC <- dev + ( log( ( ic$n -2 ) / 24 ) )*ic$np
        # CAIC (consistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )	  
		
	return(ic)
	}
##################################################
# create table of item parameters
.TAM.itempartable <- function( resp , maxK , AXsi , B , ndim ,
			resp.ind , rprobs,n.ik,pi.k){
				
#	item1 <- data.frame( "item" = colnames(resp) )
	item1 <- data.frame( "item" = dimnames(B)[[1]] )
	item1$N <- colSums(resp.ind )
	item1$M <- colSums( resp.ind * resp , na.rm=TRUE) / colSums( resp.ind )
	#****
	# Item fit
	# probs ... [ classes , items , categories ]
	probs <- aperm( rprobs , perm=c(3,1,2))
	pi.k <- matrix( pi.k , ncol=1 )
#	res <- .tam.itemfit.rmsea( n.ik , pi.k , probs )
	#####
	# Exploratory analyses show that item fit rmsea
	# does not seem to be sensitive
#	item1$rmsea <- res
	for (kk in 1:(maxK-1)){ # kk <- 1
		item1[ , paste0("AXsi_.Cat" , kk) ] <- - AXsi[,kk+1]
						}
	for (kk in 1:(maxK-1)){ # kk <- 1
		for (dd in 1:ndim){
			item1[ , paste0("B.Cat" , kk,".Dim",dd) ] <- B[,kk+1,dd]
							}
					}				
    item1 <- item1[ item1$N > 0 , ]	
	item1 <- item1[ order( paste( item1$item)) , ]		
	rownames(item1) <- NULL
	return(item1)
		}
#######################################################
# calculate counts
.tam.calc.counts <- function( resp, theta , resp.ind , 
	group , maxK , pweights , hwt ){
	TP <- nrow(theta)
	I <- ncol(resp)
	if ( is.null( group )){ group <- rep(1 , nrow(resp)) }
	G <- length( unique( group ))
	# n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	n.ik <- array( 0 , dim=c(TP,I,maxK , G ))
	for (gg in 1:G){	# gg <- 1
	ind.gg <- which( group == gg ) 		
		for (kk in 1:(maxK)){   #		kk <- 1	# category 0 ( -> 1 )
    		dkk2 <- ( resp[ ind.gg , ]  == (kk-1) ) * resp.ind[ ind.gg ] * 
					pweights[ind.gg]
			# t( t(A) * B ) = t(B) * A = crossprod(B,A)
#			n.ik[,,kk,gg] <- t( t(dkk2) %*% hwt[ind.gg,] )
			n.ik[,,kk,gg] <- crossprod( hwt[ind.gg,] , dkk2 )
						}						
					}
	# calculate pi.k
	pi.k <- matrix( pweights , nrow=nrow(resp) , ncol= ncol(hwt) )
	pi.k <- colSums( pi.k * hwt ) / colSums( pi.k )
	pi.k <- pi.k / sum( pi.k )
	res <- list( "n.ik" = n.ik , "pi.k" = pi.k)
	return(res)
	}
#####################################
