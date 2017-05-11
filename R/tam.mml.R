tam.mml <- function( resp , Y=NULL , group = NULL ,  irtmodel ="1PL" ,
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.inits = NULL , 
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = TRUE , constraint="cases" , 
            A=NULL , B=NULL , B.fixed = NULL , 
            Q=NULL , est.slopegroups=NULL , E = NULL , 
            pweights = NULL , 
			userfct.variance = NULL , variance.Npars = NULL , 
			item.elim = TRUE , verbose = TRUE , 
			control = list() 
            # control can be specified by the user 
  ){
    
    s1 <- Sys.time()
	CALL <- match.call()
		
    # display
    disp <- "....................................................\n"    
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    R <- trim_increment <- NULL
	fac.oldxsi <- acceleration <- NULL
	
	#**** handle verbose argument
	args_CALL <- as.list( sys.call() )
	control$progress <- tam_args_CALL_search( args_CALL=args_CALL , variable="verbose" , 
								default_value = TRUE )				
	#*******
    
    if ( is.null(A)){ printxsi <- FALSE  } else { printxsi <- TRUE }
    
	#--- attach control elements
    e1 <- environment()
	tam_fct <- "tam.mml"
	
	res <- tam_mml_control_list_define(control=control, envir=e1, tam_fct=tam_fct)
	con <- res$con
	con1a <- res$con1a

	resp <- as.matrix(resp)
	resp0 <- resp <- add.colnames.resp(resp)
	
	
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; flush.console()
    }  
    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }
    # define design matrix in case of PCM2
#    if (( irtmodel=="PCM2" ) & (is.null(Q)) & ( is.null(A)) ){ 
#      A <- .A.PCM2( resp ) 
#    }  
	#**** constraints
	if ( constraint == "items" ){
		irtmodel <- "PCM2"
						}
						
						
    if (( irtmodel=="PCM2" ) & ( is.null(A)) ){ 
      A <- .A.PCM2( resp , constraint=constraint , Q=Q  ) 	    
    }  
    
    
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
    
    nullY <- is.null(Y)
    
    nitems <- ncol(resp)       # number of items
    if (is.null(colnames(resp))){
      colnames(resp) <- paste0( "I" , 100+1:nitems )
    }
    nstud <- nrow(resp)        # number of students
	#*****
	nstud1 <- sum(1*( rowSums( 1 - is.na(resp) ) > 0 ))
	
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons and " , 
          nitems , "Items \n" )  ;
      flush.console()	  
    }  	  
    
    #!! check dim of person ID pid
    if ( is.null(pid) ){ pid <- seq(1,nstud) }else{ pid <- unname(c(unlist(pid))) }
    
    # print( colSums( is.na(resp)) )
    
    # normalize person weights to sum up to nstud
	pweights0 <- pweights	
	pweights <- nstud * pweights / sum(pweights)		
	#@@--

				
	# a matrix version of person weights
    pweightsM <- outer( pweights , rep(1,nitems) )
    
    # calculate ndim if only B or Q are supplied
    if ( ! is.null(B) ){ ndim <- dim(B)[3] } 
    if ( ! is.null(Q) ){ ndim <- dim(Q)[2] }
    
    betaConv <- FALSE         #flag of regression coefficient convergence
    varConv <- FALSE          #flag of variance convergence
    nnodes <- length(nodes)^ndim
    if ( snodes > 0 ){ nnodes <- snodes }

	#--- print information about nodes
	res <- tam_mml_progress_proc_nodes( progress=progress, snodes=snodes, nnodes=nnodes, 
					skillspace="normal", QMC=QMC)  	
    
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
    
    # create design matrices
    modeltype <- "PCM"
    if (irtmodel=="RSM"){  modeltype <- "RSM" }
	#****
	# ARb 2015-12-08
	maxKi <- NULL
	if ( ! (item.elim ) ){
		maxKi <- rep( maxK - 1 , ncol(resp) )		
				}
    #*** 
    design <- designMatrices( modeltype= modeltype , maxKi=maxKi , resp=resp , 
                              A=A , B=B , Q=Q , R=R, ndim=ndim  , constraint=constraint )
    A <- design$A
    B <- design$B
    cA <- design$flatA
    cA[is.na(cA)] <- 0
    if (progress){ 
      cat("    * Created Design Matrices   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }    
    
    design <- NULL				
   
    #   #---2PL---
    #   B_orig <- B  #keep a record of generated B before estimating it in 2PL model 
    #   #---end 2PL---
    
    ################################
    # number of parameters
    np <- dim(A)[[3]]
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      #    xsi <- xsi.inits 
      xsi <- rep(0,np)
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]	
    } else { xsi <- rep(0,np)   } 
    
    
    if ( ! is.null( xsi.fixed ) ){
      xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2]
      est.xsi.index <- setdiff( 1:np , xsi.fixed[,1] )
				} else { 
		est.xsi.index <- 1:np 
						}
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
      group <- match( group , groups )
      # user must label groups from 1, ... , G
      #    if ( length( setdiff( 1:G , groups)  ) > 0 ){
      #      stop("Label groups from 1, ...,G\n")
      #				}							
      var.indices <- rep(1,G)
      for (gg in 1:G){
        var.indices[gg] <- which( group == gg )[1]				
      }
    } else { 
      G <- 1 
      groups <- NULL
    }  
    # beta inits
    # (Y'Y)
    if ( ! is.null( formulaY ) ){
      formulaY <- stats::as.formula( formulaY )
      Y <- stats::model.matrix( formulaY , dataY )[,-1]   # remove intercept
      nullY <- FALSE
    }
    
    
    #  if ( ! is.null(Y) ){ 
    if (! nullY){
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
      #		colnames(Y) <- paste("group" , 1:G , sep="")
      colnames(Y) <- paste("group" , groups , sep="")
      for (gg in 1:G){ Y[,gg] <- 1*(group==gg) }
      nreg <- G - 1
    }
    
    W <- t(Y * pweights) %*% Y
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
    if ( is.null(beta.fixed) & (  is.null(xsi.fixed) ) & ( constraint=="cases") ){
      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
      if (  ndim > 1){ 
        for ( dd in 2:ndim){
          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
        }}}
    
    #****
    # ARb 2013-08-20: Handling of no beta constraints	
    # ARb 2013-08-24: correction	
    if( ! is.matrix(beta.fixed) ){
      if ( ! is.null(beta.fixed) ){
        if ( ! beta.fixed   ){ beta.fixed <- NULL }
      }
    }
    #****
    
    beta <- matrix(0, nrow = nreg+1 , ncol = ndim)  
    if ( ! is.null( beta.inits ) ){ 
      beta[ beta.inits[,1:2] ] <- beta.inits[,3]
    }	
    #	if ( ! is.null( beta.inits ) ){ 
    #		beta <- beta.inits 
    #  } else {
    #		beta <- matrix(0, nrow = nreg+1 , ncol = ndim)        
    #			}
    
    
    # define response indicator matrix for missings
    resp.ind <- 1 - is.na(resp)
    nomiss <- sum( is.na(resp) ) == 0
    #*** included nomiss in M step regression
    resp.ind.list <- list( 1:nitems )
    for (i in 1:nitems){ resp.ind.list[[i]] <- which( resp.ind[,i] == 1)  }
    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    #stop("here")  
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
    
    
    # These sufficient statistics must be changed
    # to make it more general
    # First extension:  pweights and dependent on A; needs to be further extended (e.g., different number of categories)
    # Second extension: multiple category option       -> resp \in 0:maxKi (see method definition calc_posterior_TK)
    #                                                  -> length(ItemScore) = np (see diff computation in M Step)
    #                   multiple category option Bugfix
    #                                                  -> dim(cResp) = (nstud, nitems*maxK)
    #                                                  -> adapt dim(A) to dim(cResp) for 
    #							sufficient statistic (cf. print.designMatrices)
    #**************************************************	
    # ARb 2013-04-29
    # more efficient calculation of sufficient statistics
    col.index <- rep( 1:nitems , each = maxK )
    #  cResp <- (resp[ , col.index  ]+1) *resp.ind[ , col.index ]
    cResp <- (resp +1) *resp.ind
    
    cResp <- cResp[ , col.index  ]
    #   cResp <- 1 * t( t(cResp) == rep(1:(maxK), nitems) )
    cResp <- 1 * ( cResp == matrix( rep(1:(maxK), nitems) , nrow(cResp) , 
                                    ncol(cResp) , byrow=TRUE ) )
    if ( sd(pweights) > 0 ){ 
      ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
    } else { 
      ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
    }
    if (progress){ 
      cat("    * Calculated Sufficient Statistics   (", 
          paste(Sys.time()) , ")\n") ; flush.console()	  
    }    				   
    
    #**************************************************
    # starting values for xsi
    maxAi <-  - (apply(-(A) , 3 , rowMaxs , na.rm=TRUE))  
    personMaxA <- resp.ind %*% maxAi
    #ItemMax <- personMaxA %t*% pweights  
	ItemMax <- crossprod( personMaxA , pweights ) 
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
                                      (ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )
    # starting values of zero
    if( xsi.start0 ){ xsi <- 0*xsi }
    
    #log of odds ratio of raw scores  
    if ( ! is.null(xsi.inits) ){    
		xsi[ xsi.inits[,1] ] <- xsi.inits[,2]
    }
    if ( ! is.null( xsi.fixed ) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }
    
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
# cat("b200"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1  									  
    
	#--- create grid of nodes for numeric or stochastic integration
	res <- tam_mml_create_nodes( snodes=snodes, nodes=nodes, ndim=ndim, theta=theta, QMC=QMC ) 
	theta <- res$theta
	theta2 <- res$theta2
	thetawidth <- res$thetawidth
	theta0.samp <- res$theta0.samp
	thetasamp.density <- res$thetasamp.density
        
    deviance <- 0  
	deviance.history <- tam_deviance_history_init(maxiter=maxiter)
    
    iter <- 0 
    a02 <- a1 <- 999	# item parameter change
    a4 <- 0
    
    ##**SE
    se.xsi <- 0*xsi
    se.B <- 0*B
    
    YSD <- max( apply( Y , 2 , sd ) )
    if (YSD > 1E-15 ){ YSD <- TRUE } else { YSD <- FALSE }
    
	#*****
	#@@@@ speed gains, further auxiliary objects, 2015-06-26
	Avector <- as.vector(A)
	Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }	
	if (  is.null(beta.fixed) ){ unidim_simplify <- FALSE }
	
	# unidim_simplify <- FALSE
	#@@@@
		
	#--- warning multiple group estimation
	res <- tam_mml_warning_message_multiple_group_models( ndim=ndim, G=G)
			
	#--- acceleration
	res <- tam_acceleration_inits(acceleration=acceleration, G=G, xsi=xsi, 
				variance=variance)	
	xsi_acceleration <- res$xsi_acceleration
	variance_acceleration <- res$variance_acceleration						
			
    # define progress bar for M step
    mpr <- round( seq( 1 , np , len = 10 ) )    
    hwt.min <- 0
    rprobs.min <- 0
    AXsi.min <- 0
    B.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0
    se.xsi.min <- se.xsi
    se.B.min <- se.B
    
	
	#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & (iter < maxiter) ) { 
      	  	  
		iter <- iter + 1
		#--- progress
		res <- tam_mml_progress_em0(progress=progress, iter=iter, disp=disp)		
		#--- calculate nodes for Monte Carlo integration	
		if ( snodes > 0){
			res <- tam_mml_update_stochastic_nodes( theta0.samp=theta0.samp, variance=variance, 
						snodes=snodes, beta=beta, theta=theta ) 
			theta <- res$theta
			theta2 <- res$theta2
			thetasamp.density <- res$thetasamp.density
		}			
		olddeviance <- deviance
# a0 <- Sys.time()	

		#--- calculation of probabilities
		res <- tam_mml_calc_prob( iIndex=1:nitems, A=A, AXsi=AXsi, B=B, xsi=xsi, theta=theta, 
					nnodes=nnodes, maxK=maxK, recalc=TRUE ) 
		rprobs <- res$rprobs
		AXsi <- res$AXsi	  
#cat("calc_prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
      
		#--- calculate student's prior distribution
		gwt <- tam_stud_prior( theta=theta, Y=Y, beta=beta, variance=variance, nstud=nstud, 
					nnodes=nnodes, ndim=ndim, YSD=YSD, unidim_simplify=unidim_simplify, 
					snodes=snodes ) 
# cat("stud_prior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 

		#--- calculate student's likelihood
		res.hwt <- tam_calc_posterior( rprobs=rprobs, gwt=gwt, resp=resp, nitems=nitems, 
						resp.ind.list=resp.ind.list, normalization=TRUE, 
						thetasamp.density=thetasamp.density, snodes=snodes, resp.ind=resp.ind, 
						avoid.zerosum=TRUE ) 
		hwt <- res.hwt$hwt
# cat("calc_posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 
      
# cat("before mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						           

		#--- M step: estimation of beta and variance
		resr <- tam_mml_mstep_regression( resp=resp, hwt=hwt, resp.ind=resp.ind, 
					pweights=pweights, pweightsM=pweightsM, Y=Y, theta=theta, theta2=theta2, 
					YYinv=YYinv, ndim=ndim, nstud=nstud, beta.fixed=beta.fixed, variance=variance, 
					Variance.fixed=variance.fixed, group=group, G=G, snodes=snodes, nomiss=nomiss, 
					thetasamp.density=thetasamp.density, iter=iter, min.variance=min.variance, 
					userfct.variance=userfct.variance, variance_acceleration=variance_acceleration, 
					est.variance=est.variance, beta=beta ) 
# cat("mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						       
		beta <- resr$beta     
		variance <- resr$variance	
		itemwt <- resr$itemwt	  
	    variance_acceleration <- resr$variance_acceleration
		variance_change <- resr$variance_change
		beta_change <- resr$beta_change		
		if ( beta_change < conv){ betaConv <- TRUE }
		if ( variance_change < conv){ varConv <- TRUE }
      
		#--- M-step item intercepts
		res <- tam_mml_mstep_intercept( A=A, xsi=xsi, AXsi=AXsi, B=B, theta=theta, 
					nnodes=nnodes, maxK=maxK, Msteps=Msteps, rprobs=rprobs, np=np, 
					est.xsi.index0=est.xsi.index0, itemwt=itemwt, indexIP.no=indexIP.no, 
					indexIP.list2=indexIP.list2, Avector=Avector, max.increment=max.increment, 
					xsi.fixed=xsi.fixed, fac.oldxsi=fac.oldxsi, ItemScore=ItemScore, 
					convM=convM, progress=progress, nitems=nitems, iter=iter, 
					increment.factor=increment.factor, xsi_acceleration=xsi_acceleration,
					trim_increment=trim_increment ) 
		xsi <- res$xsi
		se.xsi <- res$se.xsi
		max.increment <- res$max.increment
		xsi_acceleration <- res$xsi_acceleration
		xsi_change <- res$xsi_change	
# cat("mstep item parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 	
		
		#--- compute deviance
		res <- tam_mml_compute_deviance( loglike_num=res.hwt$rfx, loglike_sto=res.hwt$rfx, 
					snodes=snodes, thetawidth=thetawidth, pweights=pweights, deviance=deviance, 
					deviance.history=deviance.history, iter=iter ) 			  
		deviance <- res$deviance
		deviance.history <- res$deviance.history
		a01 <- rel_deviance_change <- res$rel_deviance_change
		a02 <- deviance_change <- res$deviance_change
	  	  
		if (con$dev_crit == "relative" ){ a02 <- a01 }
      
		if( deviance < deviance.min ){ 	
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
		a1 <- xsi_change
		a2 <- beta_change
		a3 <- variance_change
		devch <- -( deviance - olddeviance )	  
		
		#** print progress
		res <- tam_mml_progress_em( progress=progress, deviance=deviance, 
					deviance_change=deviance_change, iter=iter, 
					rel_deviance_change=rel_deviance_change, xsi_change=xsi_change, 
					beta_change=beta_change, variance_change=variance_change, B_change=0,
					devch=devch ) 
      
# cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
      
    } # end of EM loop
    #******************************************************
	
	# a0 <- Sys.time()	
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
    # a0 <- Sys.time()  
    # cat("restructure output") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
    
	#******
	# generate input for fixed parameters
	xsi.fixed.estimated <- generate.xsi.fixed.estimated( xsi=xsi , A=A )
	B.fixed.estimated <- generate.B.fixed.estimated( B=B )
	
	#**** standard errors AXsi
	se.AXsi <- tam_mml_se_AXsi( AXsi=AXsi, A=A, se.xsi=se.xsi, maxK=maxK ) 
    
    ##*** Information criteria
    ic <- tam_mml_ic( nstud=nstud1, deviance=deviance, xsi=xsi, xsi.fixed=xsi.fixed, 
				beta=beta, beta.fixed=beta.fixed, ndim=ndim, variance.fixed=variance.fixed, 
				G=G, irtmodel=irtmodel, B_orig=NULL, B.fixed=B.fixed, E=E, 
				est.variance=TRUE, resp=resp, est.slopegroups=NULL, 
				variance.Npars=variance.Npars, group=group ) 
#    cat("TAM.ic") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			

    #***
    # calculate counts
	res <- tam_calc_counts( resp=resp, theta=theta, resp.ind=resp.ind, group=group, 
				maxK=maxK, pweights=pweights, hwt=hwt )
    n.ik <- res$n.ik
    pi.k <- res$pi.k 	
#    cat("calculate counts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				

    #****
    # collect item parameters
	item1 <- tam_itempartable( resp=resp, maxK=maxK, AXsi=AXsi, B=B, ndim=ndim, 
					resp.ind=resp.ind, rprobs=rprobs, n.ik=n.ik, pi.k=pi.k ) 									
#     cat("tam itempartable") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1								 

	#**** collect all person statistics
	res <- tam_mml_person_posterior( pid=pid, nstud=nstud, pweights=pweights, 
				resp=resp, resp.ind=resp.ind, snodes=snodes, 
				hwtE=hwt, hwt=hwt, ndim=ndim, theta=theta ) 
	person <- res$person
	EAP.rel <- res$EAP.rel	
#    cat("person parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				  

    s2 <- Sys.time()
    if ( is.null( dimnames(A)[[3]] ) ){  
      dimnames(A)[[3]] <- paste0("Xsi" , 1:dim(A)[3] )
    }
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
        print( stats::cov2cor(varianceM) , 4 )	
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
    
    # collect xsi parameters
    obji <- data.frame( "xsi"=xsi , "se.xsi"=se.xsi ) 
    rownames(obji) <- dimnames(A)[[3]]	
    xsi <- obji
    
	#**** calculate individual likelihood
    res.hwt <- tam_calc_posterior( rprobs=rprobs, gwt=1+0*gwt, resp=resp, nitems=nitems, 
					resp.ind.list=resp.ind.list, normalization=FALSE, 
					thetasamp.density=thetasamp.density, snodes=snodes, resp.ind=resp.ind ) 
    res.like <- res.hwt$hwt
		
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "xsi" = xsi ,
                 "beta" = beta , "variance" = variance ,
                 "item" = item1 , 
                 "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
                 "post" = hwt ,  "rprobs" = rprobs , "itemweight" = itemwt ,
                 "theta" = theta , 
                 "n.ik" = n.ik , "pi.k" = pi.k ,
                 "Y" = Y , "resp" = resp0 , 
                 "resp.ind" = resp.ind , "group" = group , 
                 "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
                 "groups" = if ( is.null(group)){1} else { groups } , 			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights0 , 
                 "time" = c(s1,s2,s2-s1) , "A" = A , "B" = B  ,
                 "se.B" = se.B , 
                 "nitems" = nitems , "maxK" = maxK , "AXsi" = AXsi ,
                 "AXsi_" = - AXsi ,
                 "se.AXsi" = se.AXsi , 
                 "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
                 "hwt" = hwt ,  "like" = res.like , 
				 "ndim" = ndim ,
                 "xsi.fixed" = xsi.fixed , 
				 "xsi.fixed.estimated" = xsi.fixed.estimated , 
				 "B.fixed.estimated" = B.fixed.estimated ,
				 "beta.fixed" = beta.fixed , "Q" = Q  ,
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "irtmodel" = irtmodel ,
                 "iter" = iter ,
                 "printxsi" = printxsi ,
                 "YSD"=YSD , CALL =CALL
                 #			   "design"=design
                 #			   "xsi.min.deviance" = xsi.min.deviance ,
                 #			   "beta.min.deviance" = beta.min.deviance , 
                 # "variance.min.deviance" = variance.min.deviance 
    )
    class(res) <- "tam.mml"
    return(res)
  }


# tam.mml.output <- function(){
# 	}
