
###################################################################
# latent regression
tam.latreg <- function( like , theta=NULL , Y=NULL , group=NULL , 
				formulaY = NULL , dataY = NULL , 
				beta.fixed = NULL , beta.inits = NULL , 
				variance.fixed = NULL , variance.inits = NULL , 
				est.variance = TRUE , pweights = NULL , pid=NULL , 
				userfct.variance = NULL , variance.Npars = NULL , 
				verbose = TRUE , control = list() 
  ){
       
    s1 <- Sys.time()
	CALL <- match.call()
    # display
    disp <- "....................................................\n"    
    increment.factor <- progress <- nodes <- snodes <- ridge <- xsi.start0 <- QMC <- NULL
    maxiter <- conv <- convD <- min.variance <- max.increment <- Msteps <- convM <- NULL 
    pweightsM <- R <- NULL
	
	#**** handle verbose argument
	args_CALL <- as.list( sys.call() )
	control$progress <- tam_args_CALL_search( args_CALL=args_CALL , variable="verbose" , 
								default_value = TRUE )				
	#*******
  	#--- attach control elements
    e1 <- environment()
	tam_fct <- "tam.latreg"	
	res <- tam_mml_control_list_define(control=control, envir=e1, tam_fct=tam_fct)
	con <- res$con
	con1a <- res$con1a
	
	if ( is.null(theta) ){
	   theta <- attr( like , "theta" )
	}
						
    nodes <- theta 
	ndim <- ncol(theta)		 
	
    if (progress){ 
      cat(disp)	
      cat("Processing Data     ", paste(Sys.time()) , "\n") ; 
	  utils::flush.console()
    }  
    
    if ( ! is.null(group) ){ 
      con1a$QMC <- QMC <- FALSE
      con1a$snodes <- snodes <- 0
    }

    
    
    if ( !is.null(con$seed)){ set.seed( con$seed )	 }
    nullY <- is.null(Y)
    
    nstud <- nrow(like)        # number of students
    if ( is.null( pweights) ){
      pweights <- rep(1,nstud) # weights of response pattern
    }
    if (progress){ 
      cat("    * Response Data:" , nstud , "Persons \n" )  ;
      utils::flush.console()	  
    }  	  
    
    #!! check dim of person ID pid
    if ( is.null(pid) ){ 
			pid <- seq(1,nstud) 
				} else { 
			pid <- unname(c(unlist(pid))) 
						}
       
    # normalize person weights to sum up to nstud
    pweights <- nstud * pweights / sum(pweights)
        
    betaConv <- FALSE         #flag of regression coefficient convergence
    varConv <- FALSE          #flag of variance convergence
    # nnodes <- length(nodes)^ndim
	nnodes <- nrow(nodes)

	
	if ( snodes > 0 ){ nnodes <- snodes }
    
	#--- print information about nodes
	res <- tam_mml_progress_proc_nodes( progress=progress, snodes=snodes, nnodes=nnodes, 
					skillspace="normal", QMC=QMC)  	
	
    #*********
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
    
    # W <- t(Y * pweights) %*% Y
	W <- crossprod(Y * pweights ,  Y )
    if (ridge > 0){ diag(W) <- diag(W) + ridge }
    YYinv <- solve( W )
    
    #initialise regressors
#    if ( is.null(beta.fixed)  ){
#      beta.fixed <- matrix( c(1,1,0) , nrow= 1) 
#      if (  ndim > 1){ 
#        for ( dd in 2:ndim){
#          beta.fixed <- rbind( beta.fixed , c( 1 , dd , 0 ) )
#        }}}
    
    #****
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
    
    YSD <- max( apply( Y , 2 , stats::sd ) )
    if (YSD > 10^(-15) ){ YSD <- TRUE } else { YSD <- FALSE }
    
    # define progress bar for M step
#    mpr <- round( seq( 1 , np , len = 10 ) )
    
	#--- warning multiple group estimation
	res <- tam_mml_warning_message_multiple_group_models( ndim=ndim, G=G)	
	
    hwt.min <- 0
    deviance.min <- 1E100
    itemwt.min <- 0

	nomiss <- TRUE
	Variance.fixed <- variance.fixed
	res.hwt <- list()
    
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
 
      # cat("calc_prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
      
		# calculate student's prior distribution
		gwt <- tam_stud_prior(theta=theta , Y=Y , beta=beta , variance=variance , nstud=nstud , 
                           nnodes=nnodes , ndim=ndim,YSD=YSD, unidim_simplify=FALSE)
	    hwt <- like * gwt
		res.hwt$rfx <- rowSums(hwt)
		hwt <- hwt / rowSums(hwt)	 
      # cat("calc_posterior") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						 
            
		#-- M step: estimation of beta and variance
		if (progress){ 
	       cat("M Step Latent Regression")
		   utils::flush.console() 
		}
		oldbeta <- beta
		oldvariance <- variance 				
		resr <- tam_mml_mstep_regression( resp=NULL, hwt=hwt, resp.ind=NULL, 
					pweights=pweights, pweightsM=pweightsM, Y=Y, theta=theta, theta2=theta2, 
					YYinv=YYinv, ndim=ndim, nstud=nstud, beta.fixed=beta.fixed, variance=variance, 
					Variance.fixed=variance.fixed, group=group, G=G, snodes=snodes, nomiss=nomiss, 
					thetasamp.density=thetasamp.density, iter=iter, min.variance=min.variance, 
					userfct.variance=userfct.variance, variance_acceleration=NULL, 
					est.variance=est.variance, beta=beta, latreg_use=TRUE ) 	  
      # cat("mstep regression") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1						       

		beta <- resr$beta     
		variance <- resr$variance	
		variance_change <- resr$variance_change
		beta_change <- resr$beta_change
		
		if ( beta_change < conv){ betaConv <- TRUE }
		if ( variance_change < conv){ varConv <- TRUE }		  
      	  
		#--- compute deviance
		res <- tam_mml_compute_deviance( loglike_num=res.hwt$rfx, 
					loglike_sto=rowMeans(res.hwt$swt), snodes=snodes, 
					thetawidth=thetawidth, pweights=pweights, deviance=deviance, 
					deviance.history=deviance.history, iter=iter ) 			  
		deviance <- res$deviance
		deviance.history <- res$deviance.history
		a01 <- rel_deviance_change <- res$rel_deviance_change
		a02 <- deviance_change <- res$deviance_change	  
	        
      if( deviance < deviance.min ){ 	 
        beta.min.deviance <- beta
        variance.min.deviance <- variance	
        hwt.min <- hwt	
        deviance.min <- deviance
      }
      
		a1 <- 0
		a2 <- beta_change
		a3 <- variance_change		  
		devch <- - ( deviance - olddeviance )
		
		#** print progress
		res <- tam_mml_progress_em( progress=progress, deviance=deviance, 
					deviance_change=deviance_change, iter=iter, 
					rel_deviance_change=rel_deviance_change, xsi_change=0, 
					beta_change=beta_change, variance_change=variance_change, B_change=0, 
					is_latreg=TRUE, devch=devch ) 
        
      # cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
      
      
    } # end of EM loop
    #******************************************************

    beta.min.deviance -> beta
    variance.min.deviance -> variance	
    hwt.min -> hwt	
    deviance.min -> deviance

    ##*** information criteria
	ic <- tam_latreg_ic( nstud=nstud, deviance=deviance, beta=beta, 
				beta.fixed=beta.fixed, ndim=ndim, variance.fixed=variance.fixed, 
				G=G, est.variance=est.variance, variance.Npars=NULL, group=group ) 
				
	#**** collect all person statistics
	res <- tam_mml_person_posterior( pid=pid, nstud=nstud, pweights=pweights, 
				resp=NULL, resp.ind=NULL, snodes=snodes, 
				hwtE=hwt, hwt=hwt, ndim=ndim, theta=theta ) 
	person <- res$person
	EAP.rel <- res$EAP.rel		
	
    #cat("person parameters") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				  
    ############################################################
    s2 <- Sys.time()
	
	

    if (progress){
      cat(disp)
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
    
    # Output list
    deviance.history <- deviance.history[ 1:iter , ]
    res <- list( "beta" = beta , "variance" = variance ,
                 "person" = person , pid = pid , "EAP.rel" = EAP.rel , 
                 "post" = hwt , "theta" = theta , 
                 "Y" = Y ,  "group" = group , 
                 "G" = if ( is.null(group)){1} else { length(unique( group ) )} , 
                 "groups" = if ( is.null(group)){1} else { groups } , 			   
                 "formulaY" = formulaY , "dataY" = dataY , 
                 "pweights" = pweights , 
                 "time" = c(s1,s2,s2-s1) , 
                 "nstud" = nstud , 
				 "hwt" = hwt ,  "like" = like , 
				 "ndim" = ndim ,
                 "beta.fixed" = beta.fixed , 
                 "variance.fixed" = variance.fixed ,
                 "nnodes" = nnodes , "deviance" = deviance ,
                 "ic" = ic , 
                 "deviance.history" = deviance.history ,
                 "control" = con1a ,    "iter" = iter ,
                 "YSD"=YSD , CALL = CALL 

    )
    class(res) <- "tam.latreg"
    return(res)
  }


# tam.mml.output <- function(){
# 	}
