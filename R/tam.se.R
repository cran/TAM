tam.se <-
  function( tamobj , ...){
    if(class(tamobj) == "tam.mml"){
      res <- tam.mml.se( tamobj, ...)
    }
    if(class(tamobj) == "tam.jml"){
      # res <- tam.jml.se( tamobj, ...)
    }
    
    return( res )
    }

tam.mml.se <-
function( tamobj , numdiff.parm = .001){
    h <- numdiff.parm
    B <- tamobj$B
    A <- tamobj$A
    Y <- tamobj$Y
    nitems <- tamobj$nitems
    xsi <- ( tamobj$xsi )[,1]
    beta <- tamobj$beta
    variance <- tamobj$variance
    nstud <- tamobj$nstud
    AXsi <- tamobj$AXsi
    resp <- tamobj$resp
	ndim <- tamobj$ndim	
    theta <- tamobj$theta
	maxK <- tamobj$maxK
	thetawidth <- diff(theta[,1] )
	thetawidth <- ( ( thetawidth[ thetawidth > 0 ])[1] )^ndim
    hwt <- tamobj$hwt
    resp <- tamobj$resp
    resp.ind <- tamobj$resp.ind
    resp.ind.list <- tamobj$resp.ind.list 
    nnodes <- tamobj$nnodes	
    ntheta <- length(theta)
	# multiplication parameters for numerical differentiation
	mult <- c(0,1,-1)
	##############################################
	# Item parameters
	##############################################	
    ip <- length(xsi)
    # create result object for item parameters
    se.xsi <- rep( 0 , ip )
    cat("Item parameters\n|")
	cat(paste( rep("*" , ip) , collapse=""))
	cat("|\n|")    
    # compute likelihood
    # prior distribution for each student (normal density)
	res0a <- stud_prior.v2( theta=theta , Y=Y , beta=beta , variance=variance , 
                          nstud=nstud , nnodes=nnodes , ndim=ndim )
				  						  
	ll <- matrix( 0 , nrow=nstud , ncol=3 )	
    for (pp in 1:ip){
        # pp <- 10  # parameter pp
		if (pp == 1){ vec <- 1:3 } else { vec <- 2:3 }
        for (mm in vec){
            xsi0 <- xsi	
            xsi0[ pp ] <- xsi0[ pp] + mult[mm] * numdiff.parm
            ll0 <- 0
            # calculate probabilities
			res0 <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                            xsi=xsi0 , theta=theta , nnodes=nnodes , maxK=maxK )
            rprobs <- res0[["rprobs"]]
			# posterior distribution
		   # calculate student's likelihood
		   res0b <- calc_posterior.v2(rprobs=rprobs , gwt=res0a , resp=resp , 
                                  nitems=nitems , resp.ind.list=resp.ind.list , 
                                  normalization = FALSE , thetasamp.density = NULL , 
                                  snodes = 0 )$hwt
		   # calculate individual log likelihood																						
			ll[,mm] <- log( rowSums( res0b * thetawidth ) )
				}
			se.xsi[pp] <- sqrt( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
            cat("-") ; 
			flush.console()
			}
	cat("|\n")	
	##############################################
	# Regression parameters
	##############################################	
    # create result object for item parameters
    se.beta <- 0*beta
	nreg <- nrow(beta)
    cat("Regression parameters\n|")
	cat(paste( rep("*" , nreg*ndim) , collapse=""))
	cat("|\n|")    
	# compute response probabilities
    for (pp in 1:nreg){
	for (dd in 1:ndim){	
#		ll <- matrix( 0 , nrow=nstud , ncol=3 )
		for (mm in 2:3){
			 beta0 <- beta 
			 beta0[ pp ,dd] <- beta0[pp,dd] + mult[mm] * numdiff.parm
			# compute likelihood
			# prior distribution for each student (normal density)
			res0a <- stud_prior.v2( theta=theta , Y=Y , beta=beta0 , variance=variance , 
			                        nstud=nstud , nnodes=nnodes , ndim=ndim )
			# calculate probabilities
            res0 <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                  xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK )
            rprobs <- res0[["rprobs"]]
			# posterior distribution
            res0b <- calc_posterior.v2( rprobs=rprobs , gwt=res0a , resp=resp , 
                                        nitems=nitems , resp.ind.list=resp.ind.list ,
                                        normalization=FALSE , thetasamp.density=NULL , 
                                        snodes=0)$hwt
			ll[,mm] <- log( rowSums( res0b *thetawidth ) )
						}
			se.beta[pp,dd] <- sqrt( - 1 / sum( tamobj$pweights * ( ll[,2] + ll[,3] - 2*ll[,1] )/h^2 ) )
			cat("-") ; 
			flush.console()
				} }
	#-----------------------------------------------------------
	# handle fixed parameters
	if ( ! is.null( tamobj$xsi.fixed ) ){
		se.xsi[ tamobj$xsi.fixed[,1] ] <- 0
								}
	if ( ! is.null( tamobj$beta.fixed ) ){
		se.beta[ tamobj$beta.fixed[,1:2] ] <- 0
								}	
	
	#-----------------------------------------------------------
	cat("|\n")	
    xsi <- data.frame( tamobj$item[,1:2] , 
					"est" = xsi , "se" = se.xsi )		
	beta <- data.frame( "beta" = beta , "se" = se.beta )
	colnames(beta) <- c( paste("est.Dim" , 1:ndim , sep="")	, paste("se.Dim" , 1:ndim , sep="")	)
    flush.console()
    res <- list( "xsi" = xsi , "beta" = beta )
    return(res)
        }
