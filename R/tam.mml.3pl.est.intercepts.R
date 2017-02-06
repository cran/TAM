

##########################################################################
# estimation of item intercepts
tam.mml.3pl.est.intercepts <- function( max.increment , np , est.xsi.index0 , 
		Msteps , nitems , A , AXsi , B , xsi , guess , theta , nnodes , maxK ,
		progress , itemwt , indexIP.no , indexIP.list2 ,
		ItemScore , fac.oldxsi , rprobs , xsi.fixed , convM , rprobs0 ,
		n.ik , N.ik  , xsi.prior , indexIP.list)
{
      converge <- FALSE
      Miter <- 1
	  eps <- 1e-10
	  oldxsi <- xsi
      old_increment <- rep( max.increment , np )
      est.xsi.index <- est.xsi.index0
	  
	  #***************************************************
      while ( !converge & ( Miter <= Msteps ) ) {	
        #      xbar2 <- xxf <- xbar <- rep(0,np)		
		# numerical differentiation parameter
		h <- 1E-4
		
		#----- switch with respect to existence of guessing parameter
		guess_exists <- max( guess ) > eps
		
		if ( ! guess_exists ){
			# Only compute probabilities for items contributing to param p
			if (Miter > 1){ 
			  res.p <- tam_mml_3pl_calc_prob( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
									 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK , 
									 guess=guess )					
			  rprobs <- res.p[["rprobs"]]
			  rprobs0 <- res.p$rprobs0		  
			}

			res <- .mml.3pl.calc_exp_TK3( rprobs , A , np , est.xsi.index , itemwt ,
								 indexIP.no , indexIP.list2 , rprobs0 , guess ,
								 n.ik , N.ik )
			xbar <- res$xbar
			xbar2 <- res$xbar2
			xxf <- res$xxf
			ItemScore <- res$iscore					
			# Compute the difference between sufficient statistic and expectation
			diff <- as.vector(ItemScore) - xbar
			#Compute the Newton-Raphson derivative for the equation to be solved
			deriv <- xbar2 - xxf 	
		}
		
		if ( guess_exists){
			NX <- length(xsi)		
			ll0 <- rep( NA , NX )
			ll1m <- ll1p <- NA*ll0
			for (xx in 1:NX){
				# xx <- 1		
				iIndex <- 1:nitems				
				# iIndex <- indexIP.list[[xx]]				
				ll0[xx] <- tam_mml_3pl_calc_total_ll( iIndex= iIndex , A=A , B=B , xsi=xsi , theta=theta ,
								nnodes=nnodes , guess=guess , n.ik=n.ik , eps=eps )	
				ll1p[xx] <- tam_mml_3pl_calc_total_ll( iIndex= iIndex , A=A , B=B , 
								xsi= tam_mml_3pl_vec_add_increment( vec=xsi , h=h , index=xx ) , 
								theta=theta, nnodes=nnodes , guess=guess , n.ik=n.ik , eps=eps )
				ll1m[xx] <- tam_mml_3pl_calc_total_ll( iIndex= iIndex , A=A , B=B , 
								xsi= tam_mml_3pl_vec_add_increment( vec=xsi , h=-h , index=xx ) , 
								theta=theta, nnodes=nnodes , guess=guess , n.ik=n.ik , eps=eps )
			}
			res <- tam_mml_3pl_difference_quotient( d0=ll0 , d0p=ll1p , d0m=ll1m , h=h)
			diff <- res$d1
			deriv <- res$d2
		}
		
		#***********************
		# xsi prior
		  if ( ! is.null(xsi.prior) ){
			  d0  <- log( stats::dnorm( xsi , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d0p  <- log( stats::dnorm( xsi + h , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d0m  <- log( stats::dnorm( xsi - h , mean=xsi.prior[,1] , 
						       sd=xsi.prior[,2] ) + eps)
			  d1 <- ( d0p - d0 ) / h
			  d2 <- ( ( d0p - d0 ) - ( d0 - d0m ) ) / h^2		
              diff <- diff + d1
              deriv <- deriv + d2			  
								}			
		#************************
		
        increment <- diff*abs(1/( deriv + 10^(-20) ) )
        if ( !is.null( xsi.fixed) ){ increment[ xsi.fixed[,1] ] <- 0 } 
        #!!!	  necessary to include statement to control increment?
        ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
        increment <- ifelse( abs( increment) > abs(old_increment)  , 
                             increment/(2*ci) , 
                             increment )
        
        old_increment <- increment
        
        
        ##**SE
        se.xsi <- sqrt( 1 / abs(deriv) )
        if ( ! is.null( xsi.fixed) ){ se.xsi[ xsi.fixed[,1] ] <- 0 } 
        ##**

        xsi <- xsi+increment   # update parameter p
        #	  est.xsi.index <- which( abs(increment) > convM )
        if ( max(abs(increment)) < convM ) { converge <- TRUE }
        Miter <- Miter + 1						
        
        # stabilizing the algorithm | ARb 2013-09-10
        if (fac.oldxsi > 0 ){
          xsi <-  (1-fac.oldxsi) * xsi + fac.oldxsi *oldxsi
        }	  	  
        
        # progress bar
        if (progress){ 
          #        cat( paste( rep("-" , sum( mpr == p ) ) , collapse="" ) )
          cat("-") ; utils::flush.console()
        }
      } # end of all parameters loop
	  
	  res <- list( "xsi" = xsi , "se.xsi" = se.xsi )
	  return(res)
		}
#######################################################################

.mml.3pl.est.intercepts <- tam.mml.3pl.est.intercepts

