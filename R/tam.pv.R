tam.pv <-
function( tamobj , nplausible = 10 , 
			ntheta = 2000 , 
			normal.approx = FALSE , 
            samp.regr = FALSE , np.adj = 8 ){
    #####################################################
    # INPUT:
    # tamobj ... result from tam analysis
    # nplausible ... number of plausible values
    # ntheta ... number of simulated theta values
    # samp.regr ... sample regression coefficients?
	#        when sampling regression coefficients,
	#        plausible values are used for recalculating
	#        regression coefficients
	#        (sampling of regression coefficients only
	#          works in the unidimensional case)
	# normal.approx ... use normal distribution as an 
	#					approximation of the posterior
    ####################################################
	# 2012-07-28:  person weights included in sampling of regression coefficients
# a0 <- Sys.time()	
    type <- "nonparm"		# there is no type='normal' up to now implemented
    B <- tamobj$B
    A <- tamobj$A
    Y <- tamobj$Y
	YSD <- tamobj$YSD
    nitems <- tamobj$nitems
    xsi <- ( tamobj$xsi )[,1]
    beta <- tamobj$beta
    variance <- tamobj$variance
    nstud <- tamobj$nstud
	nthetal <- rep( 1 , ntheta )
	nnodes <- ntheta
    AXsi <- tamobj$AXsi
	ndim <- tamobj$ndim
	pweights <- tamobj$pweights
	maxK <- tamobj$maxK
    # simulate theta 
	if ( ndim == 1 ){
		MEAP <- mean( tamobj$person$EAP )
		SDEAP <- sqrt( var( tamobj$person$EAP ) + mean( tamobj$person$SD.EAP^2 ) )
				}
	if ( ndim > 1 ){
		tp1 <- tamobj$person
		ind <- grep("EAP\\.Dim" , colnames(tp1) )
		ind <- ind[ seq( 1 , length(ind) , 2 ) ]
		dat1 <- tp1[ ,  ind ]
		mu1 <- as.vector( colMeans( dat1 ) )
		var1 <- apply( dat1 , 2 , var ) / tamobj$EAP.rel
		Sigma1 <- cov2cor(variance)
		Sigma1 <- np.adj * diag( sqrt( var1) ) %*% Sigma1 %*% diag( sqrt( var1 ))
					}
    # create pv matrix (uni- and multidimensional case)
    pv <- matrix( 0 , nrow=nstud , ncol= nplausible*ndim)     
    NPV <- nplausible
    pp <- 1
#	if (samp.regr){
		cat("|")
		cat( paste( rep("*" , nplausible ) , collapse="") )
		cat("|\n|") ; flush.console()
#				}
# cat("start routine") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
	###################################################
	# routine for drawing plausible values
	while ( pp <= NPV ){
	
	   #***************************
       # 1-dimensional PV imputation	   
	   if (ndim == 1){
			# unidimensional theta simulation
			if ( ! normal.approx){			
#				theta <- matrix( rnorm( ntheta , mean = MEAP , sd = SDEAP )  , ncol= 1)	
				theta <- matrix( rnorm( ntheta , mean = MEAP , sd = np.adj*SDEAP )  , ncol= 1)	
				theta <- theta[ order( theta[,1] ) , , drop=FALSE]
								} else {
				 theta <- matrix( SDEAP * seq( - 5 , 5 , len=ntheta ) + MEAP , ncol=1 )		 
							}
				} else {
		#*****************************
        # multidimensional PV imputation		
		  # adapt for multidimensional case here!!
	      theta <- mvrnorm( ntheta , mu = mu1 , Sigma = Sigma1 )
					}
# cat("start prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							
 	    res <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
 	                         nnodes=nnodes, maxK=maxK )
		rprobs <- res[["rprobs"]]
		AXsi <- res[["AXsi"]]
# cat("calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
		# calculate student's prior distribution    	
		gwt <- stud_prior.v2( theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                          nnodes=nnodes , ndim=ndim , YSD=YSD)
# cat("stud prior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  
		# posterior distribution
		hwt <- calc_posterior.v2( rprobs=rprobs , gwt=gwt , resp=tamobj$resp , nitems=nitems , 
		                          resp.ind.list=tamobj$resp.ind.list , normalization=TRUE , 
		                          thetasamp.density=NULL , snodes=0 )$hwt
        hwt1 <- hwt			
# cat("posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  		
#		hwt1 <- rowcumsums( hwt1)
# cat("rowcumsums orig") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  		
		hwt1 <- rowCumsums.TAM(hwt1) # include this function in later versions!!
# cat("rowcumsums TAM") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  		
		if (  samp.regr ){
			#*****
			# no normal approximation
			if ( ! normal.approx){			
				rn1 <- runif( nstud )
				ind <- rowSums( hwt1 < outer( rn1 , nthetal ) ) +1
#				ind <- rowSums( hwt1 < outer( rn1 , nthetal ) )
#				ind <- ifelse( ind > ntheta , ntheta , ind )
#				ind <- ifelse( ind == 0 , 1 , ind )
				pv[,pp] <- theta1 <- theta[ind]				
									}
			#*****
			# normal approximation in unidimensional case
			if ( normal.approx){
				thetaM <- matrix( theta[,1] , nstud , ntheta , byrow=TRUE )
				EAP <- rowSums( thetaM * hwt )
				SDPost <- sqrt( rowSums( thetaM^2 * hwt ) - EAP^2 )
				pv[,pp] <- theta1 <- rnorm( nstud , mean = EAP , sd = SDPost )		
					}

				pp <- pp + 1
				#!!!				
				# adapt for multidimensional case here!		
				# up to now only a warning is included
				if ( ndim > 1 ){ 
					stop("PV imputation with sampled regression coefficients only possible for one dimension" ) 
								}	
				# sampling of regression coefficients
				modlm <- lm( theta1  ~  -1 + as.matrix(Y) , weights = pweights)
				beta <- modlm$coef			# sampled regression coefficients
				smod <- summary(modlm)
				variance <- smod$sigma^2	# sampled residual variance
				vcovmod <- vcov( modlm )		
				# sampling of a new regression coefficient from the posterior distribution				
				beta <- mvrnorm( 1 , mu = beta , Sigma = vcovmod )
				beta <- matrix( beta , ncol=1 )
				# no sampling of residual variance from the posterior
								}                  
			#********************								
			# no sampling of regression cofficients
     		if ( ! samp.regr ){
			for ( pp in 1:nplausible ){
			 #****
			 # no normal approximation
			 if ( ( ! normal.approx ) | ndim > 1 ){
				rn1 <- runif( nstud )
				ind <- interval_index( hwt1 , rn1 )	
				# Correction MWu 2012-09-18
                pv[ , (pp-1)*(ndim) + 1:ndim ] <- theta[ind , ]			
									}
              #****									 
			  # normal approximation
			  if ( normal.approx & ( ndim == 1) ){
				thetaM <- matrix( theta[,1] , nstud , ntheta , byrow=TRUE )
				EAP <- rowSums( thetaM * hwt )
				SDPost <- sqrt( rowSums( thetaM^2 * hwt ) - EAP^2 )
				pv[,pp] <- rnorm( nstud , mean = EAP , sd = SDPost )					  
									}
               if (pp!=nplausible){ cat("-") ; flush.console() }
#				pv[,pp] <- theta[ind,]
#!!!
# adapt for multidimensional case!
								}
			NPV <- nplausible / 2
						}		
			cat("-" ) ; flush.console()
					}
			cat("|\n")
# cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  			
	# label the pv matrix	
	colnames(pv) <- paste("PV" , rep(1:nplausible,each=ndim) , 
					".Dim" , rep(1:ndim,nplausible) , sep="")   
    pv <- data.frame( "pid" =tamobj$pid , pv )					
    res <- list( "pv" = pv , "hwt" = hwt , "hwt1" = hwt1 ,
                "theta" = theta  )
    return(res)
    }

	
##################################################################
##################################################################