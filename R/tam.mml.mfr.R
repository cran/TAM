tam.mml.mfr <-
  function( resp , Y=NULL , group = NULL ,  irtmodel ="1PL" ,
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.inits = NULL , 
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = FALSE , formulaA=~item+item:step, constraint="cases",
            A=NULL , B=NULL , B.fixed = NULL , 
            Q=NULL , facets=NULL, est.slopegroups=NULL , E = NULL , 
            pweights = NULL , control = list() 
            # control can be specified by the user 
  ){
    
    #------------------------------------
    # INPUT:
    # Y ... matrix of regression predictors
    # group ... group indicators
    # formulaY ... a formula object for creating regression indicators
    # dataY ... corresponding data frame for formula object of Y
    # pid ... person ID
    # xsi.fixed  matrix L * 2 where L<=np
    #		1st column: xsi parameter label
    #		2nd column: fixed item parameter value
    # beta.fixed ... L*3 matrix
    # 	1st column: beta parameter integer predictor in Y
    #   2nd column: dimension
    #   3rd column: which parameter should be fixed?
    # B.fixed ... fixed B entries
    #			first three columns are B entries, fourth entry is the value to be fixed
    # est.variance ... should variance be estimated? Relevant for 2PL model
    #  2PL ... default FALSE
    # control:
    #		control = list( nodes = seq(-6,6,len=15) , 
    #		          			convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 30 ,            
    #                   maxiter = 1000 , progress = TRUE) 
    # progress ... if TRUE, then display progress
    #-------------------------------------
    
    s1 <- Sys.time()
    # display
    disp <- "....................................................\n"
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
	
    # attach control elements
    e1 <- environment()
    con <- list( nodes = seq(-6,6,len=21) , snodes = 0 , QMC=TRUE , 
                 convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 4 ,            
                 maxiter = 1000 , max.increment = 1 , 
                 min.variance = .001 , progress = TRUE , ridge=0,seed= NULL ,
				 xsi.start0=TRUE , increment.factor=1)  	
    con[ names(control) ] <- control  
    Lcon <- length(con)
    con1a <- con1 <- con ; 
    names(con1) <- NULL
    for (cc in 1:Lcon ){
      assign( names(con)[cc] , con1[[cc]] , envir = e1 ) 
    }
  if ( !is.null(con$seed)){ set.seed( con$seed )	 }
  
  if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; flush.console()
				}    
  if ( ! is.null(group) ){ 
    con1a$QMC <- QMC <- FALSE
	con1a$snodes <- snodes <- 0
					}
				
    resp <- as.matrix(resp)
    nullY <- is.null(Y)

# cat("read data" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    # create design matrices
    design <- designMatrices.mfr(resp, formulaA=formulaA, facets=facets,  
                                 constraint=constraint, ndim=ndim,
                                 Q=Q, A=A, B=B , progress=progress)
    A <- design$A$A.3d.0
    cA <- design$A$A.flat.0
    B <- design$B$B.3d.0
    Q <- design$Q$Q.flat.0
    X <- design$X$X
    X.red <- design$X$X.noStep
    gresp <- design$gresp$gresp
    gresp.noStep <- design$gresp$gresp.noStep
	xsi.constr <- design$xsi.constr
    if ( is.null( pid ) ){ pid <- 1:(nrow(gresp) ) }

	
  if (progress){ 
      cat("    * Created Design Matrices   (", 
			paste(Sys.time()) , ")\n") ; flush.console()	  
				}    
	
	
# cat("design matrix" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
	
	#***
	# preprocess data if multiple person IDs do exist
	tp <- max( table( pid ))
	if ( tp > 1){
		persons <- sort( unique( pid ) )
		NP <- length( persons )
		person.ids <- sapply( persons , FUN = function( pp){ which( pid == pp ) } )
		PP <- matrix( NA , nrow=NP , ncol=tp)
		for (pos in 1:tp){
			#pos <- 1
			PP[,pos] <- unlist( lapply( person.ids , FUN = function( vv){ vv[pos] } ) )
					}

		gresp0 <- matrix( NA , nrow=NP , ncol= ncol(gresp) )
		colnames(gresp0) <- colnames(gresp)
		gresp0.noStep <- matrix( NA , nrow=NP , ncol= ncol(gresp.noStep) )
		colnames(gresp0.noStep) <- colnames(gresp.noStep)
		for (pos in 1:tp){
			# pos <- 2
			ind.pos <- which( ! is.na( PP[,pos]  ) )
			PP.pos <- PP[ind.pos,pos]
			
			g1 <- gresp[ PP.pos , ]
			g0 <- gresp0[ ind.pos , ]
			g0[ ! is.na(g1) ] <- g1[ ! is.na(g1) ]
			gresp0[ ind.pos , ] <- g0
			
			g1 <- gresp.noStep[ PP.pos , ]
			g0 <- gresp0.noStep[ ind.pos , ]
			g0[ ! is.na(g1) ] <- g1[ ! is.na(g1) ]
			gresp0.noStep[ ind.pos , ] <- g0
					}
		gresp0 -> gresp
		gresp0.noStep -> gresp.noStep
		pid <- persons	
		  if (progress){ 
			  cat("    * Arranged Response Data with Multiple Person Rows   (", 
					paste(Sys.time()) , ")\n") ; flush.console()	  
						}  		
				}
	###################################################

	
    nitems <- nrow( X.red )
    nstud <- nrow(gresp)        # number of students
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    
	if (progress){ 
	  cat("    * Response Data:" , nstud , "Persons and " , 
			ncol(gresp.noStep) , "Generalised Items \n" )  ;
	  flush.console()	  
					}  	
	
    #!! check dim of person ID pid
    if ( is.null(pid) ){ pid <- seq(1,nstud) }
   
    
    # normalize person weights to sum up to nstud
    pweights <- nstud * pweights / sum(pweights)
    # a matrix version of person weights
    pweightsM <- outer( pweights , rep(1,nitems) )
    
    # calculate ndim if only B or Q are supplied
    if ( ! is.null(B) ){ ndim <- dim(B)[3] } 
    if ( ! is.null(Q) ){ ndim <- dim(Q)[2] }
    
    betaConv <- FALSE         #flag of regression coefficient convergence
    varConv <- FALSE          #flag of variance convergence
    nnodes <- length(nodes)^ndim
    if ( snodes > 0 ){ nnodes <- snodes }
    
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
    
    ################################
    # number of parameters
    np <- dim(A)[[3]]
    
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      xsi <- xsi.inits 
    } else { xsi <- rep(0,np)   } 
    if ( ! is.null( xsi.fixed ) ){
      xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2]
      est.xsi.index <- setdiff( 1:np , xsi.fixed[,1] )
    } else { est.xsi.index <- 1:np }
	est.xsi.index -> est.xsi.index0
   
    # variance inits  
    # initialise conditional variance 
    if ( !is.null( variance.inits ) ){
      variance <- variance.inits
    } else variance <- diag( ndim ) 
    if ( !is.null(variance.fixed) ){
      variance[ variance.fixed[,1:2 ,drop=FALSE] ] <- variance.fixed[,3]
      variance[ variance.fixed[,c(2,1) ,drop=FALSE] ] <- variance.fixed[,3]	
    }
    
    # group indicators for variance matrix
    if ( ! is.null(group) ){ 
      groups <- sort(unique(group))
      G <- length(groups)
      # user must label groups from 1, ... , G
      if ( length( setdiff( 1:G , groups)  ) > 0 ){
        stop("Label groups from 1, ...,G\n")
      }
      var.indices <- rep(1,G)
      for (gg in 1:G){
        var.indices[gg] <- which( group == gg )[1]				
      }
    } else { G <- 1 }
    
    
    # beta inits
    # (Y'Y)
    if ( ! is.null( formulaY ) ){
      formulaY <- as.formula( formulaY )
      Y <- model.matrix( formulaY , dataY )[,-1]   # remove intercept
    }
    
    if ( ! is.null(Y) ){ 
      Y <- as.matrix(Y)
      nreg <- ncol(Y)
      if ( is.null( colnames(Y) ) ){
        colnames(Y) <- paste("Y" , 1:nreg , sep="")
      }
      if ( ! nullY ){ 		
        Y <- cbind(1,Y)          #add a "1" column for the Intercept
        colnames(Y)[1] <- "Intercept"
      }
    } else 
    { 
      Y <- matrix( 1 , nrow=nstud , ncol=1 ) 
      nreg <- 0
    }
    if ( G > 1 & nullY ){	
      Y <- matrix( 0 , nstud , G )
      colnames(Y) <- paste("group" , 1:G , sep="")
      for (gg in 1:G){ Y[,gg] <- 1*(group==gg) }
      nreg <- G - 1
    }
    
    
    
    W <- t(Y * pweights) %*% Y
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    
    
    #initialise regressors
    if ( is.null(beta.fixed) & (  is.null(xsi.fixed) ) ){
      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
      if ( ndim > 1){ 
        for ( dd in 2:ndim){
          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
        }}}  
    
    if ( ! is.null( beta.inits ) ){ 
      beta <- beta.inits 
    } else {
      beta <- matrix(0, nrow = nreg+1 , ncol = ndim)        
    }


    #***
	#*** ARb 2013-03-26   o Change setup of calculation

    # define response indicator matrix for missings
#    resp.ind <- 1 - is.na(resp)
#    resp.col.ind <- as.numeric(X.red[,ifelse("item" %in% colnames(X.red), "item", 1)]) 
    resp.ind.list <- list( 1:nitems )
    gresp.ind <- 1 - is.na( gresp ) 
    gresp.noStep.ind <- 1 - is.na( gresp.noStep )	 
#	nomiss <- sum( is.na(gresp.noStep) == 0 )  	#*** included nomiss in M step regression
	resp.ind <- gresp.noStep.ind
	nomiss <- mean( gresp.noStep.ind ) == 1
    for (i in 1:nitems){ 
#      resp.ind.list[[i]] <- which( resp.ind[,resp.col.ind[i]] == 1)  
       resp.ind.list[[i]] <- which( gresp.noStep.ind[,i] == 1)  
    }
#    resp[ is.na(resp) ] <- 0 	# set all missings to zero
	 gresp0.noStep <- gresp.noStep
     gresp[ is.na( gresp) ] <- 0
	 gresp.noStep[ is.na( gresp.noStep) ] <- 0
#    gresp.noStep.ind <- resp.ind[ ,resp.col.ind]
#    gresp.ind <- resp.ind[ ,as.numeric(X[,ifelse("item" %in% colnames(X), "item", 1)])]
    
    #@@ ARb:Include Starting values for xsi??
    #       xsi <- - qnorm( colMeans( resp ) )
    AXsi <- matrix(0,nrow=nitems,ncol=maxK )  #A times xsi
    
    # Create an index linking items and parameters
    indexIP <- colSums(aperm(A, c(2,1,3)) != 0, na.rm = TRUE)
    # define list of elements for item parameters
    indexIP.list <- list( 1:np )
    for ( kk in 1:np ){
      indexIP.list[[kk]] <- which( indexIP[,kk] > 0 )
    }
    lipl <- cumsum( sapply( indexIP.list , FUN = function(ll){ length(ll) } ) )
	indexIP.list2 <- unlist(indexIP.list)
	indexIP.no <- as.matrix( cbind( c(1 , lipl[-length(lipl)]+1 ) , lipl ) )
	
	
	#******** addition ARb 2013-05-03
	
	
    # These sufficient statistics must be changed
    # to make it more general
    # First extension:  pweights and dependent on A; needs to be further extended (e.g., different number of categories)
    # Second extension: multiple category option       -> resp \in 0:maxKi (see method definition calc_posterior_TK)
    #                                                  -> length(ItemScore) = np (see diff computation in M Step)
    #                   multiple category option Bugfix
    #                                                  -> dim(cResp) = (nstud, nitems*maxK)
    #                                                  -> adapt dim(A) to dim(cResp) for sufficient statistic (cf. print.designMatrices)

# cat("before item score" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1		
	#*********************** 
#	  col.index <- rep( 1:nitems , each = maxK )
#	  cResp <- gresp.noStep[ , col.index  ]*gresp.noStep.ind[ , col.index ]
#	  cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
	  ##@@@## missings should be taken into account
#	  cResp <- cResp * gresp.noStep.ind[ , col.index ] 
	#***********************	
    #...TK - 24.08.2012:  cA returned from designMatrices
    # sufficient statistics
#    ItemScore <- (cResp %*% cA) %t*% pweights
	col.index <- rep( 1:nitems , each = maxK )
	cResp <- (gresp.noStep+1)*gresp.noStep.ind
    cResp <- cResp[ , col.index  ]
    cResp <- 1 * ( cResp == matrix( rep(1:(maxK), nitems) , nrow(cResp) , 
		    		ncol(cResp) , byrow=TRUE ) )
    cA <- t( matrix( aperm( A , c(2,1,3) ) , nrow = dim(A)[3] , byrow = TRUE ) )
    cA[is.na(cA)] <- 0					
  if ( sd(pweights) > 0 ){ 
	ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
			} else { 
	ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
				   }
				   
# cat("calc ItemScore" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    if (progress){ 
        cat("    * Calculated Sufficient Statistics   (", 
			paste(Sys.time()) , ")\n") ; flush.console()	  
				}   			
    # starting values for xsi
    gresp.ind.tmp <- gresp.noStep.ind[ , col.index  ]
#    gresp.ind.tmp[,- grep(paste0("-step",(maxK-1)),colnames(gresp))] <- 0
    ItemMax <- (gresp.ind.tmp %*% cA) %t*% pweights
	ItemMax <- as.vector( t( colSums( gresp.ind.tmp * pweights ) ) %*% cA )    

	
	# cat("calc ItemMax" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1

	
#    xsi[est.xsi.index] <- - log(abs(ItemScore[est.xsi.index]/(ItemMax[est.xsi.index]-ItemScore[est.xsi.index])))  
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
				(ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )
	# starting values of zero
	if( xsi.start0 ){ xsi <- 0*xsi }

	
	#log of odds ratio of raw scores  
	xsi[ is.na(xsi) ] <- 0
    if ( ! is.null(xsi.inits) ){   xsi <- xsi.inits  }
    if ( ! is.null(xsi.fixed) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }
    
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
    # nodes
    if ( snodes == 0 ){ 
      theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )
      #we need this to compute sumsig2 for the variance
      theta2 <- matrix(theta.sq(theta), nrow=nrow(theta),ncol=ncol(theta)^2)            
      # grid width for calculating the deviance
      thetawidth <- diff(theta[,1] )
      thetawidth <- ( ( thetawidth[ thetawidth > 0 ])[1] )^ndim 
      thetasamp.density <- NULL
    } else {
    # sampled theta values
	if (QMC){
#		library(sfsmisc)					
		r1 <- QUnif (n=snodes, min = 0, max = 1, n.min = 1, p=ndim, leap = 409)						
		theta0.samp <- qnorm( r1 )
			} else {
			theta0.samp <- matrix( mvrnorm( snodes , mu = rep(0,ndim) , 
						Sigma = diag(1,ndim ) )	,
                        nrow= snodes , ncol=ndim )			
				}
    thetawidth <- NULL
    }
    
    
    deviance <- 0  
    deviance.history <- matrix( 0 , nrow=maxiter , ncol = 2)
    colnames(deviance.history) <- c("iter" , "deviance")
    deviance.history[,1] <- 1:maxiter
    
    iter <- 0 
    a02 <- a1 <- 999	# item parameter change
    a4 <- 0
	
    	hwt.min <- 0
		rprobs.min <- 0
		AXsi.min <- 0
		B.min <- 0
		deviance.min <- 0
		itemwt.min <- 0

  ##**SE
	se.xsi <- 0*xsi
	se.B <- 0*B
	se.xsi.min <- se.xsi
	se.B.min <- se.B
	
    # display
    disp <- "....................................................\n"
    # define progress bar for M step    
	
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & (iter < maxiter) ) { 

# a0 <- Sys.time()	
      iter <- iter + 1
      if (progress){ 
        cat(disp)	
        cat("Iteration" , iter , "   " , paste( Sys.time() ) )
        cat("\nE Step\n") ; flush.console()
      }
      # calculate nodes for Monte Carlo integration	
      if ( snodes > 0){
        theta <- beta[ rep(1,snodes) , ] +  t ( t(chol(variance)) %*% t(theta0.samp) )
        # calculate density for all nodes
        thetasamp.density <- dmvnorm( theta , mean = as.vector(beta[1,]) , sigma = variance )
        # recalculate theta^2
        theta2 <- matrix( theta.sq(theta) , nrow=nrow(theta) , ncol=ncol(theta)^2 )   
      }			
      olddeviance <- deviance
      # calculation of probabilities
      res <- calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                          nnodes=nnodes , maxK=maxK , recalc=TRUE )	
# cat("calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  
						  
      rprobs <- res[["rprobs"]]
      AXsi <- res[["AXsi"]]
      # calculate student's prior distribution
      gwt <- stud_prior.v2(theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                           nnodes=nnodes , ndim=ndim)
#print( head(gwt))						   

      # calculate student's likelihood
      res.hwt <- calc_posterior.v2(rprobs=rprobs , gwt=gwt , resp=gresp.noStep , nitems=nitems , 
                                   resp.ind.list= resp.ind.list , normalization=TRUE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
						resp.ind=resp.ind	)	
      hwt <- res.hwt[["hwt"]]  
# print( head(hwt))

	  
# cat("calc posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		  
      if (progress){ cat("M Step Intercepts   |"); flush.console() }
      # collect old values for convergence indication
      oldxsi <- xsi
      oldbeta <- beta
      oldvariance <- variance 
      # M step: estimation of beta and variance
      resr <- mstep.regression( resp=gresp.noStep , hwt=hwt , resp.ind=gresp.noStep.ind , pweights=pweights , 
                                pweightsM=pweightsM , Y=Y , theta=theta , theta2=theta2 , YYinv=YYinv , 
                                ndim=ndim , nstud=nstud , beta.fixed=beta.fixed , variance=variance , 
                                Variance.fixed=variance.fixed , group=group ,  G=G , snodes = snodes ,
								nomiss=nomiss)
      
      beta <- resr$beta
#cat("m step regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        
      variance <- resr$variance	
      if( ndim == 1 ){  # prevent negative variance
        variance[ variance < min.variance ] <- min.variance 
      }
      itemwt <- resr$itemwt
      
      # constraint cases (the design matrix A has no constraint on items)
      if ( max(abs(beta-oldbeta)) < conv){    
        betaConv <- TRUE       # not include the constant as it is constrained
      }
      if (G == 1){
        diag(variance) <- diag(variance) + 10^(-10)
        #		if( det(variance) < 10^(-20) ){
        #		  stop("\n variance is close to singular or zero. Estimation cannot proceed")
        #@@ARb: I would not prefer to stop the program but adding a small
        #       constant in the diagonal.
        #				} 
      }
    
      if (max(abs(variance-oldvariance)) < conv) varConv <- TRUE     
      ######################################
      #M-step item parameters
      converge <- FALSE
      Miter <- 1
      old_increment <- rep( max.increment , np )
	  est.xsi.index <- est.xsi.index0
      while (!converge & ( Miter <= Msteps ) ) {	
        # Only compute probabilities for items contributing to param p
        if (Miter > 1){ 
          res.p <- calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK)					
          rprobs <- res.p[["rprobs"]]            
        }
		#	indexIP.list2 <- unlist(indexIP.list)
		# ==> Vector with item indices for parameter estimation
		#	indexIP.no <- cbind( c(1 , lipl[-length(lipl)]+1 ) , lipl )
		# ==> Vector with start and end indices for item parameter estimation

		res <- calc_exp_TK3( rprobs , A , np , est.xsi.index , itemwt ,
			indexIP.no , indexIP.list2 )
        xbar <- res$xbar
		xbar2 <- res$xbar2
		xxf <- res$xxf

			        
        # Compute the difference between sufficient statistic and expectation
        diff <- as.vector(ItemScore) - xbar
        #Compute the Newton-Raphson derivative for the equation to be solved
        deriv <- xbar2 - xxf 			
        increment <- diff*abs(1/( deriv + 10^(-20) ) )
        if ( !is.null( xsi.fixed) ){ increment[ xsi.fixed[,1] ] <- 0 } 
        #!!!	  necesessary to include statement to control increment?
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
        increment <- ifelse( abs( increment) > abs(old_increment)  , 
                             increment/(2*ci) , 
                             increment )
        # Code from Margaret
        #    	for (p in est.xsi.index ) {
        #			 while (abs(increment[p]) > abs(old_increment[p])){
        #				increment[p] <- increment[p] / 2.0  #damping the increment
        #						}
        #						}
        
        old_increment <- increment
        xsi <- xsi+increment   # update parameter p
#		est.xsi.index <- which( abs(increment) > convM )		
        if ( max(abs(increment)) < convM ) { 
				converge <- TRUE 
						}
        Miter <- Miter + 1						

		##**SE
		 se.xsi <- sqrt( 1 / abs(deriv) )
		 if ( ! is.null( xsi.fixed) ){ se.xsi[ xsi.fixed[,1] ] <- 0 } 
		##**        
        # progress bar
        if (progress){ 
      #    cat( paste( rep("-" , sum( mpr == p ) ) , collapse="" ) )
		  cat("-")
          flush.console()
			}
      } # end of all parameters loop
      
	#***
	# decrease increments in every iteration
	if( increment.factor > 1){max.increment <-  1 / increment.factor^iter }	  
	  

# cat("m step item parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        
      # calculate deviance
      if ( snodes == 0 ){ 
        deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )
      } else {
#        deviance <- - 2 * sum( pweights * log( res.hwt$rfx ) )
      deviance <- - 2 * sum( pweights * log( rowMeans( res.hwt$swt ) ) )
      }
      deviance.history[iter,2] <- deviance
      a01 <- abs( ( deviance - olddeviance ) / deviance  )
      a02 <- abs( ( deviance - olddeviance )  )	
      
	if( deviance - olddeviance > 0 ){ 
		xsi.min.deviance <- xsi.min.deviance 
		beta.min.deviance <- beta.min.deviance
		variance.min.deviance <- variance.min.deviance
		hwt.min <- hwt.min
		rprobs.min <- rprobs.min
		AXsi.min <- AXsi.min
		B.min <- B.min
		deviance.min <- deviance.min
		itemwt.min <- itemwt.min
		se.xsi.min <- se.xsi.min
		se.B.min <- se.B.min
					}   else { 
		xsi.min.deviance <- xsi 
		beta.min.deviance <- beta
		variance.min.deviance <- variance	
        hwt.min <- hwt	
		AXsi.min <- AXsi	
		B.min <- B
		deviance.min <- deviance
		itemwt.min <- itemwt
		se.xsi.min <- se.xsi
		se.B.min <- se.B
				}
      
      a1 <- max( abs( xsi - oldxsi ))	
      a2 <- max( abs( beta - oldbeta ))	
      a3 <- max( abs( variance - oldvariance ))
      if (progress){ 
        cat( paste( "\n  Deviance =" , round( deviance , 4 ) ))
      devch <- -( deviance - olddeviance )
	  cat( " | Deviance change:", round( devch  , 4 ) )
	  if ( devch < 0 & iter > 1 ){ cat ("   Deviance increases!") }
        cat( "\n  Maximum intercept parameter change:" , round( a1 , 6 ) )
        cat( "\n  Maximum regression parameter change:" , round( a2 , 6 ) )  
        if ( G == 1 ){ 
          cat( "\n  Variance: " , round( variance[ ! lower.tri(variance)] , 4 ) , " | Maximum change:" , round( a3 , 6 ) )  
        } else {
          cat( "\n  Variance: " , round( variance[var.indices] , 4 ) ,
               " | Maximum change:" , round( a3 , 6 ) )  		
        }					
        cat( "\n  beta ",round(beta,4)  )
        cat( "\n" )
        flush.console()
      }
    } # end of EM loop
    #############################################################
	#############################################################
  		xsi.min.deviance -> xsi 
		beta.min.deviance -> beta
		variance.min.deviance -> variance	
        hwt.min -> hwt	
		AXsi.min -> AXsi	
		B.min -> B
		deviance.min -> deviance
		itemwt.min -> itemwt
		se.xsi.min -> se.xsi	
		se.B.min -> se.B
  #******
  	#***
	resp <- gresp0.noStep
	resp.ind <- gresp.noStep.ind
	
  ##**SE  
  # standard errors of AXsi parameters
  # check for missing entries in A
	se.AXsi <- 0*AXsi
	A1 <- A
	A1[ is.na(A) ] <- 0
	se.xsiD <- diag( se.xsi^2 )
	for (kk in 1:maxK){  # kk <- 1
	se.AXsi[,kk] <- sqrt( diag( A1[,kk,] %*% se.xsiD %*% t( A1[,kk,]) ) )
			}

	##*** Information criteria
	ic <- .TAM.ic( nstud , deviance , xsi , xsi.fixed ,
		beta , beta.fixed , ndim , variance.fixed , G ,
		irtmodel , B_orig=NULL , B.fixed , E , est.variance =TRUE ,
		resp )
		
	#***
	# calculate counts
	res <- .tam.calc.counts( resp, theta , resp.ind , 
		group , maxK , pweights , hwt )
	n.ik <- res$n.ik
	pi.k <- res$pi.k 

	#****
	# collect item parameters
	item1 <- .TAM.itempartable( resp , maxK , AXsi , B , ndim ,
				 resp.ind , rprobs,n.ik,pi.k)
		
  
    #####################################################
    # post ... posterior distribution	
    # create a data frame person	
    person <- data.frame( "pid"=pid , "case" = 1:nstud , "pweight" = pweights )


#    person$score <- rowSums( resp * resp.ind )
    resp2 <- resp
	resp2[ is.na(resp2) ] <- 0
    person$score <- rowSums( resp * resp.ind , na.rm=TRUE)

    # use maxKi here; from "design object"
    nstudl <- rep(1,nstud)

    person$max <- rowSums( outer( nstudl , apply( resp2 ,2 , max , na.rm=T) ) * resp.ind )
		
    # calculate EAP
    # EAPs are only computed in the unidimensional case for now,
    # but can be easily adapted to the multidimensional case
    if ( snodes == 0 ){ 
      hwtE <- hwt 
    } else { 	
#		hwtE <- hwt / snodes 
		hwtE <- hwt
			}
    if ( ndim == 1 ){
      person$EAP <- rowSums( hwtE * outer( nstudl , theta[,1] ) )
      person$SD.EAP <- sqrt( rowSums( hwtE * outer( nstudl , theta[,1]^2 ) ) - person$EAP^2)
      #***
      # calculate EAP reliability
      # EAP variance
      EAP.variance <- weighted.mean( person$EAP^2 , pweights ) - ( weighted.mean( person$EAP , pweights ) )^2
      EAP.error <- weighted.mean( person$SD.EAP^2 , pweights )
      EAP.rel <- EAP.variance / ( EAP.variance + EAP.error )	
    } else { 
      EAP.rel <- rep(0,ndim)
      names(EAP.rel) <- paste("Dim",1:ndim , sep="")
      for ( dd in 1:ndim ){
        #	dd <- 1  # dimension
        person$EAP <- rowSums( hwtE * outer( nstudl , theta[,dd] ) )
        person$SD.EAP <- sqrt(rowSums( hwtE * outer( nstudl , theta[,dd]^2 ) ) - person$EAP^2)	
        #***
        # calculate EAP reliability
        # EAP variance
        EAP.variance <- weighted.mean( person$EAP^2 , pweights ) - ( weighted.mean( person$EAP , pweights ) )^2
        EAP.error <- weighted.mean( person$SD.EAP^2 , pweights )
        EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )	
        colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP.Dim" , dd , sep="")
        colnames(person)[ which( colnames(person) == "SD.EAP" ) ] <- paste("SD.EAP.Dim" , dd , sep="")				
      }
      person <- data.frame( "pid" = pid , person )
    }
    ############################################################
    s2 <- Sys.time()
	
   item <- data.frame( "xsi.index" = 1:np , 
                        "xsi.label" = dimnames(A)[[3]] , 
                        "est" = xsi )
    if (progress){
      cat(disp)
      cat("Item Parameters\n")
      item2 <- item
      item2[,"est"] <- round( item2[,"est"] , 4 )
      print(item2)
      cat("...................................\n")
      cat("Regression Coefficients\n")
      print( beta , 4  )
      cat("\nVariance:\n" ) # , round( varianceM , 4 ))
      if (G==1 ){ 
        varianceM <- matrix( variance , nrow=ndim , ncol=ndim ) 
        print( varianceM , 4 )	
      } else { 
        print( variance[ var.indices] , 4 )	}
      if ( ndim > 1){ 
        cat("\nCorrelation Matrix:\n" ) # , round( varianceM , 4 ))	
        print( cov2cor(varianceM) , 4 )	
      }
      cat("\n\nEAP Reliability:\n")
      print( round (EAP.rel,3) )
      cat("\n-----------------------------")
      devmin <- which.min( deviance.history[,2] )
      if ( devmin < iter ){
        cat(paste("\n\nMinimal deviance at iteration " , devmin , 
                  " with deviance " , round(deviance.history[ devmin , 2 ],3) , sep="") , "\n")
        cat("The corresponding estimates are\n")
        cat("  xsi.min.deviance\n  beta.min.deviance \n  variance.min.deviance\n\n")
      }
      cat( "\nStart: " , paste(s1))
      cat( "\nEnd: " , paste(s2),"\n")
      print(s2-s1)
      cat( "\n" )
    }
    
	####################################
	# collect xsi parameters
	
	xsiFacet <- as.data.frame( (xsi.constr$xsi.table)[,1:2]	)
 	obji <- data.frame( "parameter" = dimnames(A)[[3]] , 
			"xsi"=xsi , "se.xsi"=se.xsi ) 		
	rownames(obji) <- paste(obji$parameter)
	rownames(xsiFacet) <- paste( xsi.constr$xsi.table[,1] )
	xsi1 <- merge( x = xsiFacet , y= obji , by="parameter" , all=TRUE )
	A1 <- xsi.constr$xsi.constraint %*% xsi
	xsi1[ match( rownames(xsi.constr$xsi.constraint) , xsi1$parameter) , 
				"xsi" ] <- A1
	xsi1 <- xsi1[ match( xsiFacet$parameter , xsi1$parameter) , ]
    xsi.facets <- xsi1
	rownames(xsi.facets) <- NULL
	xsi <- obji[,-1]
	rownames(xsi) <- dimnames(A)[[3]]
	
  # Output list
  deviance.history <- deviance.history[ 1:iter , ]
  res <- list( "xsi" = xsi , "xsi.facets" = xsi.facets , 
			   "beta" = beta , "variance" = variance ,
			   "item" = item1 , 
               "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
               "post" = hwt ,  "rprobs" = rprobs , "itemweight" = itemwt ,
               "theta" = theta , 
			   "n.ik" = n.ik , "pi.k" = pi.k ,
			   "Y" = Y , "resp" = resp , 
               "resp.ind" = resp.ind , "group" = group , 
			   "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
               "formulaY" = formulaY , "dataY" = dataY , 
               "pweights" = pweights , 
               "time" = c(s1,s2,s2-s1) , "A" = A , "B" = B  ,
			   "se.B" = se.B , 
               "nitems" = nitems , "maxK" = maxK , "AXsi" = AXsi ,
			   "se.AXsi" = se.AXsi , 
               "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
               "hwt" = hwt , "ndim" = ndim ,
               "xsi.fixed" = xsi.fixed , "beta.fixed" = beta.fixed ,
			   "formulaA"=formulaA , "facets"=facets ,
               "variance.fixed" = variance.fixed ,
               "nnodes" = nnodes , "deviance" = deviance ,
			   "ic" = ic , 
               "deviance.history" = deviance.history ,
               "control" = con1a , "irtmodel" = irtmodel ,
			   "iter" = iter
#			   "xsi.min.deviance" = xsi.min.deviance ,
#			   "beta.min.deviance" = beta.min.deviance , 
# "variance.min.deviance" = variance.min.deviance 
						)
  class(res) <- "tam.mml"
  return(res)
  }
