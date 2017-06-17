
tam_mml_mstep_intercept <- function( A, xsi, AXsi, B, theta , nnodes , maxK,
		Msteps, rprobs, np , est.xsi.index0, itemwt, indexIP.no , indexIP.list2 , 
		Avector, max.increment, xsi.fixed, fac.oldxsi, ItemScore, convM,
		progress, nitems, iter, increment.factor, xsi_acceleration, 
		trim_increment = "cut" , prior_list_xsi=NULL, eps = 1E-20)
{	  
    converge <- FALSE
    Miter <- 1
		
    old_increment <- rep( max.increment , np )
    est.xsi.index <- est.xsi.index0
	oldxsi <- old_xsi <- xsi
	increments_msteps <- rep(NA,Msteps)
    if (progress){ 
		cat("M Step Intercepts   |")
		utils::flush.console() 
	}	
	#----------------------------------------------------
	# begin algorithm M-step
    while ( ! converge & ( Miter <= Msteps ) ) {	       
        if (Miter > 1){ 
			res.p <- tam_mml_calc_prob( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK)					
			rprobs <- res.p[["rprobs"]]            
        }

		res <- tam_calc_exp( rprobs=rprobs, A=A, np=np, est.xsi.index=est.xsi.index, 
					itemwt=itemwt, indexIP.no=indexIP.no, indexIP.list2=indexIP.list2, 
					Avector=Avector ) 						 
        xbar <- res$xbar
        xbar2 <- res$xbar2
        xxf <- res$xxf	        
		
        # Compute the difference between sufficient statistic and expectation
        diff <- as.vector(ItemScore) - xbar
        #Compute the Newton-Raphson derivative for the equation to be solved
        deriv <- xbar2 - xxf
		
		#-- include prior distributions
		res <- tam_evaluate_prior( prior_list = prior_list_xsi , parameter = xsi )
		d1 <- res$d1
		d2 <- res$d2
		logprior_xsi <- res$d0
		
		diff <- diff + d1
		deriv <- abs(deriv) + abs( d2 )		
		#-- define increments
        increment <- diff*abs( 1/( deriv + eps ) )
        if ( ! is.null( xsi.fixed) ){ 
			increment[ xsi.fixed[,1] ] <- 0 
		} 
		#--- trim increments
		increment <- tam_trim_increment( increment=increment, max.increment=old_increment, 
							trim_increment=trim_increment)		
		old_increment <- increment        
        ##**SE
        se.xsi <- sqrt( 1 / abs(deriv) )
        if ( ! is.null( xsi.fixed) ){ 
			se.xsi[ xsi.fixed[,1] ] <- 0 
		} 
        ##**	        
        xsi <- xsi+increment

		max_change <- max(abs(increment))				
		increments_msteps[Miter] <- max_change		
        if ( max_change < convM ){ converge <- TRUE }

        Miter <- Miter + 1						        
        # stabilizing the algorithm | ARb 2013-09-10
        if (fac.oldxsi > 0 ){
			xsi <-  (1-fac.oldxsi) * xsi + fac.oldxsi *oldxsi
        }			   
        # progress bar
        if (progress){ 
			cat("-")
			utils::flush.console()
        }
    } # end of all parameters loop
	#--------------------------------------
    #*** decrease increments in every iteration
    if( increment.factor > 1){
		max.increment <-  1 / increment.factor^iter 
	}   
	#*** acceleration of xsi parameter
	if ( xsi_acceleration$acceleration != "none" ){		
		xsi_acceleration <- tam_accelerate_parameters( xsi_acceleration=xsi_acceleration , 
								xsi=xsi , iter=iter , itermin=3)
		xsi <- xsi_acceleration$parm
	}	
	#*** maximum xsi parameter change
	xsi_change <- max( abs( xsi - old_xsi) )
	#------------------
	# OUTPUT
	res <- list(xsi=xsi, max.increment = max.increment, se.xsi = se.xsi, Miter=Miter,
				xsi_acceleration=xsi_acceleration, xsi_change = xsi_change,
				Miter=Miter, increments_msteps=increments_msteps, logprior_xsi=logprior_xsi )
	return(res)	
}

