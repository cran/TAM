tam.mml.mfr <-
  function( resp , Y=NULL , group = NULL ,  irtmodel ="1PL" ,
            formulaY = NULL , dataY = NULL , 
            ndim = 1 , pid = NULL ,
            xsi.fixed=NULL ,  xsi.setnull = NULL , 
			xsi.inits = NULL , 			
            beta.fixed = NULL , beta.inits = NULL , 
            variance.fixed = NULL , variance.inits = NULL , 
            est.variance = TRUE , formulaA=~item+item:step, constraint="cases",
            A=NULL , B=NULL , B.fixed = NULL , 
            Q=NULL , facets=NULL, est.slopegroups=NULL , E = NULL , 
            pweights = NULL , verbose = TRUE , control = list() ,
            delete.red.items=TRUE
            # control can be specified by the user 
  ){
    
    
	CALL <- match.call()
    a0 <- Sys.time()    
    s1 <- Sys.time()
    # display
    disp <- "....................................................\n"
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    resp_orig <- resp
    B00 <- B 
	B <- trim_increment <- NULL
	fac.oldxsi <- acceleration <- NULL	
	
	#**** handle verbose argument
	args_CALL <- as.list( sys.call() )
	control$progress <- tam_args_CALL_search( args_CALL=args_CALL , variable="verbose" , 
								default_value = TRUE )				
	#*******	
  	#--- attach control elements
    e1 <- environment()
	tam_fct <- "tam.mml.mfr"	
	res <- tam_mml_control_list_define(control=control, envir=e1, tam_fct=tam_fct)
	con <- res$con
	con1a <- res$con1a

	# userfct.variance is not allowed in tam.mml.mfr
	userfct.variance <- NULL
	
    #***
    fac.oldxsi <- max( 0 , min( c( fac.oldxsi , .95 ) ) )
    
    if ( constraint=="items" ){ beta.fixed <- FALSE }
    
    pid0 <- pid <- unname(c(unlist(pid)))
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; utils::flush.console()
    }    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }
    
	resp <- as.matrix(resp)
	resp <- add.colnames.resp(resp)		
    itemnames <- colnames(resp)
	
    nullY <- is.null(Y)
    
    if ( ! is.null(facets) ){ facets <- as.data.frame(facets) }
    
# cat("read data" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    # create design matrices
    if(ncol(resp)>1){
      maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    } else {
		if(ncol(resp)==1){
		  item.ind <- grep("Item", names(facets), ignore.case=TRUE)
		  if(!is.null(item.ind)){
				if ( length(item.ind) == 0 ){
					item.ind <- NULL 
								}
							}		   
		  if(!is.null(item.ind)){ 
			maxKi <- stats::aggregate( resp , facets[,item.ind,drop=FALSE] , 
								max, na.rm=TRUE )[,2]
		  }else{
			maxKi <- stats::aggregate( resp , facets[,1,drop=FALSE] , 
								max, na.rm=TRUE )[,2]
			
		  }
		}		
	}
	
	#*****************
	# handle formula and facets
    resp00 <- resp

	res <- tam_mml_mfr_dataprep( formulaA=formulaA, xsi.setnull=xsi.setnull, B=B, 
				Q=Q, resp=resp, pid=pid, facets=facets, beta.fixed=beta.fixed ) 
	formulaA <- res$formula_update
	xsi.setnull <- res$xsi.setnull	
	beta.fixed <- res$beta.fixed
	facets <- res$facets
	PSF <- res$PSF
	pid <- res$pid

	diffKi <- FALSE
	var_ki <- stats::var( maxKi )
	if ( is.na(var_ki) ){ var_ki <- 0 }

	
	if ( var_ki > 1E-3 ){ 
	     diffKi <- TRUE
		 design <- designMatrices.mfr2(resp, formulaA=formulaA, facets=facets,  
                                 constraint=constraint, ndim=ndim,
                                 Q=Q, A=A, B=B , progress=progress)
		 xsi.elim <- design$xsi.elim
		if ( ! is.null(xsi.elim) ){	
			 if ( nrow(xsi.elim) > 0 ){
				 xsi.elim2 <- cbind( xsi.elim[,2] , 99 )		 
				 xsi.fixed <- rbind( xsi.fixed , xsi.elim2 )
										}
									}
			# set first beta coefficient to zero
			if ( is.null( beta.fixed ) ){
				dimB <- dim(design$B$B.3d.0	)	
			    beta.fixed <- cbind( 1 , 1:dimB[3] , 0)
							}
					} else {				
         design <- designMatrices.mfr(resp, formulaA=formulaA, facets=facets,  
                                 constraint=constraint, ndim=ndim,
                                 Q=Q, A=A, B=B , progress=progress)	
								 
							}	 
																										
													
    A <- design$A$A.3d.0	
    cA <- design$A$A.flat.0	
    B <- design$B$B.3d.0
    Q <- design$Q$Q.flat.0
    X <- design$X$X
    X.red <- design$X$X.noStep
    gresp <- design$gresp$gresp
    gresp.noStep <- design$gresp$gresp.noStep
    xsi.constr <- design$xsi.constr

	#****************************
	items00 <- colnames(resp00)
	I00 <- length(items00)
	D <- dim(B00)[3]
	if ( ! is.null(B00) ){		
			rownames_A <- dimnames(A)[[1]]
			for (ii in 1:I00){
				#		ii <- 1
				ind <- grep( paste0( items00[ii] , "-" ) , rownames_A )
				if ( length(ind) > 0 ){
					I2 <- length(ind)
					for (vv in 1:I2){
						B[ ind[vv] , , 1:D] <- B00[ii , , 1:D ,drop=FALSE] * ( B[ ind[vv] , , 1:D] != 0  )
								}						
							}		
						}									
					}   # end is.null()
			#*****************

	
    if ( is.null( pid ) ){ pid <- 1:(nrow(gresp) ) }
    
	
#Revalpr(" head(gresp.noStep)")	
	
    design <- NULL
    
    if (progress){ 
      cat("    * Created Design Matrices   (", 
          paste(Sys.time()) , ")\n") ; utils::flush.console()	  
    }    
    
#    cat(" --- design matrix ready" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    
    #***
    # preprocess data if multiple person IDs do exist
    tp <- max( table( pid ))

    if ( tp > 1){
      persons <- sort( unique( pid ) )
      NP <- length( persons )
      # cat("*** multiple persons start" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1		
      #****
      # ARb 2013-08-23: added simplify=TRUE
      person.ids <- sapply( persons , FUN = function( pp){ which( pid == pp ) } ,
                            simplify=FALSE)
      #print(person.ids[[5]] )					
      #cat("*** multiple persons sapply function" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1				
      PP <- matrix( NA , nrow=NP , ncol=tp)
      for (pos in 1:tp){
        #pos <- 1
        PP[,pos] <- unlist( lapply( person.ids , FUN = function( vv){ vv[pos] } ) )
      }

     
#     cat("*** multiple persons lapply function" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
      gresp0 <- matrix( NA , nrow=NP , ncol= ncol(gresp) )
      colnames(gresp0) <- colnames(gresp)
      gresp0.noStep <- matrix( NA , nrow=NP , ncol= ncol(gresp.noStep) )
      colnames(gresp0.noStep) <- colnames(gresp.noStep)
      grespNA <- ( ! is.na( gresp ) )
      grespnoStepNA <- ( ! is.na( gresp.noStep ) )		
      #cat("*** multiple persons start pos" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1			
	  
	  #***
	  # check multiple rows
	  m1 <- rowsum( 1-is.na(gresp.noStep) , pid )
	  	  
# Revalpr("table( rowsum( m1 > 1 )[,1] )")	  
	  
	  h1 <- sum(m1>1)
	  if (h1>0){
		cat("* Combinations of person identifiers and facets are not unique.\n")
		cat("* Use an extended 'formulaA' to include all \n")
		cat("  relevant facets and the argument 'xsi.setnull'.\n")
		cat("  See the help page of 'tam.mml' (?tam.mml) Example 10a.\n") 
		stop()			
				}
	  	  
      for (pos in 1:tp){
        ind.pos <- which( ! is.na( PP[,pos]  ) )
        PP.pos <- PP[ind.pos,pos]
        g1 <- gresp[ PP.pos , ]
        g0 <- gresp0[ ind.pos , ]
        #			ig1 <- ( ! is.na(g1) )
        ig1 <- grespNA[ PP.pos , ]
        # * this check is time-consuming! release it to rcpp
        g0[ ig1 ] <- g1[ ig1 ]
        gresp0[ ind.pos , ] <- g0
        g1 <- gresp.noStep[ PP.pos , ]
        g0 <- gresp0.noStep[ ind.pos , ]
        #			ig1 <- ( ! is.na(g1) )
        ig1 <- grespnoStepNA[ PP.pos , ]
        g0[ ig1 ] <- g1[ ig1 ]
        gresp0.noStep[ ind.pos , ] <- g0

	
      }
      #cat("*** multiple persons loop over pos" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1				
      gresp0 -> gresp
      gresp0.noStep -> gresp.noStep
      pid <- persons	
      if (progress){ 
        cat("    * Arranged Response Data with Multiple Person Rows   (", 
            paste(Sys.time()) , ")\n") ; utils::flush.console()	  
      }  		
    }
    ###################################################
    
	# set some xsi effects to null
	if ( ! is.null(xsi.setnull) ){	
		xsi.labels <- dimnames(A)[[3]]
		xsi0 <- NULL	
		N1 <- length(xsi.setnull)
		for (nn in 1:N1){
			ind.nn <- grep( xsi.setnull[nn] , xsi.labels )	
			l1 <- cbind( ind.nn , 0 )
			xsi0 <- rbind( xsi0 , l1 )	
			colnames(xsi0) <- NULL
						}
		xsi.fixed <- rbind( xsi.fixed , xsi0 )
		i2 <- duplicated(xsi.fixed[,1])
		if ( sum(i2) > 0 ){
			xsi.fixed <- xsi.fixed[ - i2  , ]
								}							
						}

#cat("process data in case of multiple persons" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1


    nitems <- nrow( X.red )
    nstud <- nrow(gresp)        # number of students
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons and " , 
          ncol(gresp.noStep) , "Generalized Items (" , paste(Sys.time()) ,")\n" )  ;
      utils::flush.console()	  
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

	#--- print information about nodes
	res <- tam_mml_progress_proc_nodes( progress=progress, snodes=snodes, nnodes=nnodes, 
					skillspace="normal", QMC=QMC)  		
	
    #********* 
    # maximum no. of categories per item. Assuming dichotomous
    maxK <- max( resp , na.rm=TRUE ) + 1 
    
    ################################
    # number of parameters
    np <- dim(A)[[3]]
    
    # xsi inits
    if ( ! is.null(xsi.inits) ){
      #      xsi <- xsi.inits 
      xsi <- rep(0,np) 
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]				
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
    
    #    if ( ! is.null(Y) ){ 
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
      colnames(Y) <- paste("group" , groups , sep="")
      for (gg in 1:G){ Y[,gg] <- 1*(group==gg) }
      nreg <- G - 1
    }
    # redefine Y in case of multiple persons
    if (tp>1){
      if ( nrow(gresp) != nrow(Y) ){
        Ypid <- rowsum( Y , pid0 )
        Y <- Ypid / Ypid[,1]
      }
    }
       
    #    W <- t(Y * pweights) %*% Y
    W <- crossprod(Y * pweights ,  Y )
    
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
    if ( is.null(beta.fixed) & (  is.null(xsi.fixed) ) ){
      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
      if ( ndim > 1){ 
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
    #*****
    
    beta <- matrix(0, nrow = nreg+1 , ncol = ndim)  
    if ( ! is.null( beta.inits ) ){ 
      beta[ beta.inits[,1:2] ] <- beta.inits[,3]
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
    #***
    miss.items <- rep(0,nitems)
    for (i in 1:nitems){ 
      #      resp.ind.list[[i]] <- which( resp.ind[,resp.col.ind[i]] == 1)  
      resp.ind.list[[i]] <- which( gresp.noStep.ind[,i] == 1)  
      miss.items[i] <- i * ( length(resp.ind.list[[i]]) == 0 )
    }
    #    resp[ is.na(resp) ] <- 0 	# set all missings to zero
    gresp0.noStep <- gresp.noStep
    gresp[ is.na( gresp) ] <- 0
    gresp.noStep[ is.na( gresp.noStep) ] <- 0
    #    gresp.noStep.ind <- resp.ind[ ,resp.col.ind]
    #    gresp.ind <- resp.ind[ ,as.numeric(X[,ifelse("item" %in% colnames(X), "item", 1)])]
        
    #*****
    # ARb 2013-09-09: deletion of items
    miss.items <- miss.items[ miss.items > 0 ]	
    if ( length(miss.items) == 0 ){ delete.red.items <- FALSE }
    if (delete.red.items){					
      miss.itemsK <- NULL
      for (kk in 1:maxK ){
        miss.itemsK <- c( miss.itemsK , ( miss.items - 1 )* maxK + kk )
      }
      
      miss.itemsK <- sort(miss.itemsK)
      gresp <- gresp[ , - miss.itemsK ]
      gresp.noStep <- gresp.noStep[ , - miss.items ]
      gresp.noStep.ind <- gresp.noStep.ind[ , - miss.items ]		
      A <- A[ - miss.items , , , drop=FALSE]
      B <- B[ - miss.items , , ,drop=FALSE]	
      resp.ind.list <- resp.ind.list[ - miss.items ]
      resp.ind <- resp.ind[ , - miss.items ]
      nitems <- ncol(gresp.noStep)
      pweightsM <- outer( pweights , rep(1,nitems) )
      
      if (progress){ 
        cat("    * Reduced Response Data:" , nstud , "Persons and " , 
            ncol(gresp.noStep) , "Generalized Items (" , paste(Sys.time()) ,")\n" )  ;
        utils::flush.console()	  
      } 		
      
    }
    
    #****
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
    if ( stats::sd(pweights) > 0 ){ 
      ItemScore <- as.vector( t( colSums( cResp * pweights ) ) %*% cA )
    } else { 
      ItemScore <- as.vector( t( colSums( cResp) ) %*% cA )			
    }
    
    # cat("calc ItemScore" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
       
    if (progress){ 
      cat("    * Calculated Sufficient Statistics   (", 
          paste(Sys.time()) , ")\n") ; utils::flush.console()	  
    }   			
    # starting values for xsi
    gresp.ind.tmp <- gresp.noStep.ind[ , col.index  ]
    #    gresp.ind.tmp[,- grep(paste0("-step",(maxK-1)),colnames(gresp))] <- 0
    # ItemMax <- (gresp.ind.tmp %*% cA) %t*% pweights
	ItemMax <- crossprod(gresp.ind.tmp %*% cA , pweights )
    ItemMax <- as.vector( t( colSums( gresp.ind.tmp * pweights ) ) %*% cA )    
    
    
    # cat("calc ItemMax" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
    
    
    #    xsi[est.xsi.index] <- - 
    #  log(abs(ItemScore[est.xsi.index]/(ItemMax[est.xsi.index]-ItemScore[est.xsi.index])))  
    xsi[est.xsi.index] <- - log(abs(( ItemScore[est.xsi.index]+.5)/
                                      (ItemMax[est.xsi.index]-ItemScore[est.xsi.index]+.5) ) )							  
    # starting values of zero
    if( xsi.start0 == 1){ 
			xsi <- 0*xsi 
					}
    if( xsi.start0 == 2){ 
		ind1 <- which( dimnames(A)[[3]] %in% colnames(resp) )
		ind2 <- which( dimnames(A)[[3]] %in% paste0( "step" ,1:9) )
		ind3 <- setdiff( seq(1,length(xsi) ) , union(ind1,ind2) )
		xsi[ind3] <- 0
					}
					
    
    
    #log of odds ratio of raw scores  
    xsi[ is.na(xsi) ] <- 0
    if ( ! is.null(xsi.inits) ){  
      #			xsi <- xsi.inits  
      xsi[ xsi.inits[,1] ] <- xsi.inits[,2]			
    }
    if ( ! is.null(xsi.fixed) ){   xsi[ xsi.fixed[,1] ] <- xsi.fixed[,2] }
    
    xsi.min.deviance <- xsi
    beta.min.deviance <- beta
    variance.min.deviance <- variance
    
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
    
    hwt.min <- 0
    rprobs.min <- 0
    AXsi.min <- 0
    B.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0
    
	#*****
	#@@@@ 2015-06-26
	Avector <- as.vector(A)
	Avector[ is.na(Avector) ] <- 0
	unidim_simplify <- TRUE
    YSD <- max( apply( Y , 2 , stats::sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }		
	if (G > 1){ unidim_simplify <- FALSE }
	if ( YSD){ unidim_simplify <- FALSE }	
	if (  is.null(beta.fixed) ){ unidim_simplify <- FALSE }
	#@@@@	
	
	#--- acceleration
	res <- tam_acceleration_inits(acceleration=acceleration, G=G, xsi=xsi, 
				variance=variance)	
	xsi_acceleration <- res$xsi_acceleration
	variance_acceleration <- res$variance_acceleration					
	
	#--- warning multiple group estimation
	res <- tam_mml_warning_message_multiple_group_models( ndim=ndim, G=G)
	
    ##**SE
    se.xsi <- 0*xsi
    se.B <- 0*B
    se.xsi.min <- se.xsi
    se.B.min <- se.B
    

    
    devch <- 0
    
    # display
    disp <- "....................................................\n"
    # define progress bar for M step        
 # cat("rest  " ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1								 

   
    ##############################################################   
    ##############################################################   
    ##############################################################   
    #Start EM loop here
    while ( ( (!betaConv | !varConv)  | ((a1 > conv) | (a4 > conv) | (a02 > convD)) )  & 
              (iter < maxiter) ) { 

		# a0 <- Sys.time()	
		iter <- iter + 1
		#--- progress
		res <- tam_mml_progress_em0(progress=progress, iter=iter, disp=disp)
		# calculate nodes for Monte Carlo integration	
		if ( snodes > 0){
			res <- tam_mml_update_stochastic_nodes( theta0.samp=theta0.samp, variance=variance, 
						snodes=snodes, beta=beta, theta=theta ) 
			theta <- res$theta
			theta2 <- res$theta2
			thetasamp.density <- res$thetasamp.density  
		}			
		olddeviance <- deviance
		#--- calculation of probabilities
		res <- tam_mml_calc_prob(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                          nnodes=nnodes , maxK=maxK , recalc=TRUE )	
      # cat("calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							  
		rprobs <- res$rprobs
		AXsi <- res$AXsi
		
		#--- calculate student's prior distribution
		gwt <- tam_stud_prior( theta=theta, Y=Y, beta=beta, variance=variance, nstud=nstud, 
					nnodes=nnodes, ndim=ndim, YSD=YSD, unidim_simplify=unidim_simplify, 
					snodes=snodes ) 					   
      
		#--- calculate student's likelihood
		res.hwt <- tam_calc_posterior( rprobs=rprobs, gwt=gwt, resp=gresp.noStep, nitems=nitems, 
						resp.ind.list=resp.ind.list, normalization=TRUE, 
						thetasamp.density=thetasamp.density, snodes=snodes, resp.ind=resp.ind, 
						avoid.zerosum=TRUE ) 
		hwt <- res.hwt$hwt

# cat("calc posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		  

		#--- M step: estimation of beta and variance
		resr <- tam_mml_mstep_regression( resp=gresp.noStep, hwt=hwt, 
					resp.ind=gresp.noStep.ind, pweights=pweights, pweightsM=pweightsM, 
					Y=Y, theta=theta, theta2=theta2, YYinv=YYinv, ndim=ndim, nstud=nstud, 
					beta.fixed=beta.fixed, variance=variance, Variance.fixed=variance.fixed, 
					group=group, G=G, snodes=snodes, nomiss=nomiss, iter=iter, 
					min.variance=min.variance, userfct.variance=userfct.variance, 
					variance_acceleration=variance_acceleration, est.variance=est.variance, 
					beta=beta ) 
																									
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
	    devch <- - ( deviance - olddeviance )
		
		#--- print progress
		res <- tam_mml_progress_em( progress=progress, deviance=deviance, deviance_change=deviance_change,
					iter=iter, rel_deviance_change=rel_deviance_change, xsi_change=xsi_change, 
					beta_change=beta_change, variance_change=variance_change, B_change=0,
					devch=devch )       
					
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
    
	#****
	# look for non-estimable xsi parameters
#    xsi[ xsi == 99 ] <- NA	

	#******
	# generate input for fixed parameters
	xsi.fixed.estimated <- generate.xsi.fixed.estimated( xsi , A )
	B.fixed.estimated <- generate.B.fixed.estimated(B)
	
	#**** standard errors AXsi
	se.AXsi <- tam_mml_se_AXsi( AXsi=AXsi, A=A, se.xsi=se.xsi, maxK=maxK )     
    
    ##*** Information criteria
    ic <- tam_mml_ic( nstud=nstud, deviance=deviance, xsi=xsi, xsi.fixed=xsi.fixed, 
				beta=beta, beta.fixed=beta.fixed, ndim=ndim, 
				variance.fixed=variance.fixed, G=G, irtmodel=irtmodel, B_orig=NULL, 
				B.fixed=B.fixed, E=E, est.variance=TRUE, resp=resp, 
				est.slopegroups=NULL, variance.Npars=NULL, group=group ) 
    
    #***
    # calculate counts
	res <- tam_calc_counts( resp=gresp.noStep, theta=theta, resp.ind=gresp.noStep.ind, 
				group=group, maxK=maxK, pweights=pweights, hwt=hwt ) 				
    n.ik <- res$n.ik
    pi.k <- res$pi.k 
        
    #****
    # collect item parameters    
    item1 <- tam_itempartable( resp=gresp.noStep , maxK=maxK , AXsi=AXsi, B=B, 
					ndim=ndim, resp.ind=gresp.noStep.ind, 
					rprobs=rprobs, n.ik=n.ik, pi.k=pi.k, order=TRUE )    
  
	
	#**** collect all person statistics
	res <- tam_mml_person_posterior( pid=pid, nstud=nstud, pweights=pweights, 
				resp=gresp.noStep, resp.ind=gresp.noStep.ind, snodes=snodes, 
				hwtE=hwt, hwt=hwt, ndim=ndim, theta=theta ) 
	person <- res$person
	EAP.rel <- res$EAP.rel			
	
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
    
    ####################################
    # collect xsi parameters
    xsiFacet <- as.data.frame( (xsi.constr$xsi.table)[,1:2]	)
    obji <- data.frame( "parameter" = dimnames(A)[[3]] , 
                        "xsi"=xsi , "se.xsi"=se.xsi ) 		
    rownames(obji) <- paste(obji$parameter)
    rownames(xsiFacet) <- paste( xsi.constr$xsi.table[,1] )

    xsi1 <- merge( x = xsiFacet , y= obji , by="parameter" , all=TRUE )
    A1 <- xsi.constr$xsi.constraint %*% xsi
	
	incl <- match( rownames(xsi.constr$xsi.constraint) , xsi1$parameter)
    xsi1[ incl , "xsi" ] <- A1

    xsi1 <- xsi1[ match( xsiFacet$parameter , xsi1$parameter) , ]
    xsi.facets <- xsi1
    rownames(xsi.facets) <- NULL
    i1 <- grep( "Intercept" , xsi.facets$parameter)

    if ( length(i1) > 0 ){
      xsi.facets <-  xsi.facets[ - i1 , ] 
    }

	XX <- xsi.constr$xsi.constraint
	incl <- match( rownames(XX) , xsi.facets$parameter)
	vcov_xsi <- diag( obji$se.xsi^2 )
	se2 <- sqrt( diag( XX %*% vcov_xsi %*% t(XX) ))
	xsi.facets[ incl , "se.xsi"] <- se2
	
	#@@@@@@@@@@@@@@@@@@@@@@@@ control xsi.facets
	if( xsi.constr$intercept_included ){	
		ind <- which( paste(xsi.facets$facet) == "item" )
		n1 <- length(ind)
			if ( n1 > 0 ){
				itemc <- itemnames[n1]
				itemo <- paste0("item" , n1 )		
				g1 <- which( paste(xsi.facets$parameter) == itemo )
				if ( length(g1) > 0 ){
					xsi.facets$parameter[g1] <- itemc
									}

				g1 <- grep( paste0(itemo , ":") , paste(xsi.facets$parameter)  )
				if ( length(g1) > 0 ){
					xsi.facets$parameter <- gsub( paste0(itemo , ":") , paste0(itemc , ":")  ,
										paste(xsi.facets$parameter) )
									}

				g1 <- grep( paste0(itemo , "-") , dimnames(A)[[1]]  )
				if ( length(g1) > 0 ){
					dimnames(A)[[1]] <- gsub( paste0(itemo , "-") , paste0(itemc , "-")  ,
										dimnames(A)[[1]] )
									}

									
						}
				}	
	#@@@@@@@@@@@@@@@@@@@@@@@@@
	
    xsi <- obji[,-1]
    rownames(xsi) <- dimnames(A)[[3]]
    
    if(delete.red.items) resp <- resp[,-miss.items]
    colnames(resp) <- dimnames(A)[[1]]
    
	 res.hwt <- tam_calc_posterior(rprobs=rprobs , gwt=1+0*gwt , resp=resp , nitems=nitems , 
                                   resp.ind.list=resp.ind.list , normalization=FALSE , 
                                   thetasamp.density=thetasamp.density , snodes=snodes ,
                                   resp.ind=resp.ind )	
      res.like <- res.hwt[["hwt"]] 		

	
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
                 "groups" = if ( is.null(group)){1} else { groups } , 			   			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights , 
                 "time" = c(s1,s2,s2-s1) , "A" = A , "B" = B  ,
                 "se.B" = se.B , 
                 "nitems" = nitems , "maxK" = maxK , "AXsi" = AXsi ,
                 "AXsi_" = - AXsi ,			   
                 "se.AXsi" = se.AXsi , 
                 "nstud" = nstud , "resp.ind.list" = resp.ind.list ,
                 "hwt" = hwt , "like" = res.like , "ndim" = ndim ,
                 "xsi.fixed" = xsi.fixed , 
				 "xsi.fixed.estimated" = xsi.fixed.estimated , 
				 "B.fixed.estimated" = B.fixed.estimated , 
				 "beta.fixed" = beta.fixed , "Q"=Q,
                 "formulaA"=formulaA , "facets"=facets ,
				 "xsi.constr" = xsi.constr , 
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a , "irtmodel" = irtmodel ,
                 "iter" = iter , "resp_orig" = resp_orig ,
                 "printxsi"=TRUE , "YSD"=YSD , "PSF" = PSF ,
				 CALL = CALL 
                 #			   "design"=design
                 #			   "xsi.min.deviance" = xsi.min.deviance ,
                 #			   "beta.min.deviance" = beta.min.deviance , 
                 # "variance.min.deviance" = variance.min.deviance 
    )
    class(res) <- "tam.mml"
    return(res)
  }
