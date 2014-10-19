##################################
# Information criteria
.mml.3pl.TAM.ic <- function( nstud , deviance , xsi , xsi.fixed ,
	beta , beta.fixed , ndim , variance.fixed , G , irtmodel ,
	B_orig=NULL , B.fixed , E , est.variance , resp ,
	est.slopegroups=NULL , skillspace , delta , est.guess , fulldesign ,
	est.some.slopes , gammaslope , gammaslope.fixed, gammaslope.constr.V ,
	gammaslope.constr.Npars ){
  #***Model parameters
  ic <- data.frame("n" = nstud , "deviance" = deviance )
  dev <- deviance
	# xsi parameters
	ic$Nparsxsi <- length(xsi)
	if ( ! is.null( xsi.fixed) ){ 
			ic$Nparsxsi <- ic$Nparsxsi - nrow(xsi.fixed ) }
	# B slopes
	ic$NparsB <- 0
    if ( est.some.slopes ){
		ic$NparsB <- length(gammaslope)
		if ( ! is.null(gammaslope.constr.V) ){
		   ic$NparsB <- ic$NparsB - ncol(gammaslope.constr.V)		
							}
						}
	if ( ! is.null( gammaslope.fixed ) ){
		ic$NparsB <- ic$NparsB - nrow(gammaslope.fixed )
					}
	ic$NparsB <- ic$NparsB - gammaslope.constr.Npars
	
	
	# beta regression parameters
	ic$Nparsbeta <- 0
	ic$Nparscov <- 0
	if ( skillspace=="normal"){
		ic$Nparsbeta <- dim(beta)[1] * dim(beta)[2]
		if ( ! is.null( beta.fixed) ){ 
				ic$Nparsbeta <- ic$Nparsbeta - nrow(beta.fixed ) }
		# variance/covariance matrix
		ic$Nparscov <- ndim + ndim*(ndim-1)/2
		if ( ! est.variance ){ ic$Nparscov <- ic$Nparscov - ndim }
		if ( ! is.null( variance.fixed) ){ 
				ic$Nparscov <- max(0 , ic$Nparscov - nrow(variance.fixed ) )
										 }
								}												 
	ic$Nguess <- length( setdiff( unique(est.guess) , 0 ) )									 
	if ( skillspace != "normal" ){								 
		ic$Ndelta <- prod( dim(delta) )
        if ( fulldesign ){ 
             ic$Ndelta <- ic$Ndelta - ncol(delta)
							}
							} else { ic$Ndelta <- 0 }
	# total number of parameters
	ic$Npars <- ic$np <- ic$Nparsxsi + ic$NparsB + ic$Nparsbeta + ic$Nparscov + 
		          ic$Nguess + ic$Ndelta	
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
.mml.3pl.TAM.itempartable <- function( resp , maxK , AXsi , B , ndim ,
			resp.ind , rprobs,n.ik,pi.k , guess , est.guess ,
			order.items=FALSE){
				
	item1 <- data.frame( "item" = colnames(resp) )
#	item1 <- data.frame( "item" = dimnames(B)[[1]] )
	item1$N <- colSums(resp.ind )
	item1$M <- colSums( resp.ind * resp , na.rm=TRUE) / colSums( resp.ind )
	#****
	# Item fit
	# probs ... [ classes , items , categories ]
	probs <- aperm( rprobs , perm=c(3,1,2))
	pi.k <- matrix( pi.k , ncol=1 )
    if ( is.null( est.guess) ){ est.guess <- 0 }
	item1$est.guess <- est.guess
	item1$guess <- guess
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
	if ( order.items ){
		item1 <- item1[ order( paste( item1$item)) , ]		
					}
	rownames(item1) <- NULL
	return(item1)
		}
#######################################################
