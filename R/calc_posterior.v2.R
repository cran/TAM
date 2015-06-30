


###########################################################
calc_posterior.v2 <-
  function(rprobs , gwt , resp , nitems , 
           resp.ind.list , normalization = TRUE , 
           thetasamp.density = NULL , snodes = 0 , resp.ind=NULL,
		   avoid.zerosum=FALSE , logprobs=FALSE ){   

# a0 <- Sys.time()		   


    if ( snodes == 0 ){ 
      fx <- gwt  
    } else {
      # calculate individual 'sampling weight'
      swt <- fx <- gwt / outer( rep(1,nrow(gwt)) , thetasamp.density )
	  # This is essentially equal to one.
#    	swt <-fx <- gwt
    } 
    nstud <- nrow(fx)
    # using c Code here
    storage.mode(resp) <- "integer"
	
	fx0 <- fx


#***
# eps <- .001
# rprobs <- rprobs + eps
#***
	
	
# cat("vor calcfx") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
    fx <- .Call("calcfx", fx, rprobs, resp.ind.list, resp)
# cat("nach calcfx") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
#Revalpr("head(fx0)")	
	# logprobs
#	eps <- 1E-10
#	rprobs <- log(rprobs + eps )
#	fx1 <- .Call("calcfx_logprobs", fx=fx0, rprobs, resp.ind.list, resp)
#	fx <- exp(fx1)
# cat("nach calcfx_logprobs") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
	
	if (avoid.zerosum ){	
		fxs <- rowSums( fx )
		#m1 <- min( fxs[ fxs > 0 ] ) / 1E3 / ncol(fx )
		m1 <- max( min( fxs[ fxs > 0 ] , na.rm=TRUE) , 1E-200 ) / 1E3 / ncol(fx )	
		fx[ (fxs == 0) , ] <- m1
		fx[ is.na(fxs) , ] <- m1 
		# fx <- fx + m1
		}
		
# cat("nach calcfx (2)") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
    # numerical integration
    if ( snodes == 0 ){ 
#      rfx <- rowSums(fx)
      rfx <- rowSums(fx)
      if (normalization ){
        hwt <- fx / rfx } else {   hwt <- fx }
    }
    # Monte Carlo integration
    if ( snodes > 0 ){ 
#      rfx <- rowMeans(fx)
		rfx <- rowSums(fx)		
      if (normalization ){
		 hwt <- fx / rfx 	
			} else { hwt <- fx }
    }

    res <-  list("hwt" = hwt , "rfx" = rfx )
    if ( snodes > 0 ){ 
		
		res[["swt" ]] <- fx
				}
# cat(" in  posterior rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	    
    return(res)
  }
#####################################################################
