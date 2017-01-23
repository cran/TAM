

######################################################################
  
tam.jml.xsi2 <-
  function ( resp , resp.ind, A, A.0 , B, nstud, nitems, maxK, convM, 
             ItemScore, theta, xsi, Msteps, pweightsM,
             est.xsi.index , rp3 , rp3.sel , rp3.pweightsM 
  ){
    
    
    #Update item parameters
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    r <- matrix(0,nrow=nitems,ncol=maxK) 
    rr <- array(0,dim=c(nitems,maxK,maxK))
    AA <- array (0, dim=c(nitems,maxK,maxK))
    
    maxChangeP <- 0
    errorP <- rep(0, max(est.xsi.index))
    convergeAllP <- FALSE
    p_loop <- est.xsi.index

    PP <- length(p_loop)
    # convergeP <- rep(FALSE,max(est.xsi.index))
    old_xsi <- xsi
	PP1 <- length(xsi)
	convergeP <- rep(FALSE,PP1)
	nonest.xsi.index <- setdiff( seq(1,PP1) , est.xsi.index )
	convergeP[ nonest.xsi.index ] <- TRUE
	
    # begin loop 
    iterP <- 1
    # old_increment <- rep(5,max(p_loop))
	old_increment <- rep(5,PP1)
    cat(" Item parameter estimation |")
    while (!convergeAllP & ( iterP <= Msteps ) ) {
   
# a0 <- Sys.time()
		  res.p <- calc_prob.v5( iIndex = 1:nitems , A , AXsi , 
								 B , xsi , theta[ rp3.sel$caseid ,,drop=FALSE ] , 
								 nrow(rp3.sel) , maxK , TRUE )      		
		  rprobs <- res.p[["rprobs"]]       

# cat("one iteration xsi -- calc prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 		  
		  #compute probability weights, summed over students, so that there is no cycling through students for parameter estimation (p loop)
		  M1 <- resp.ind[ rp3.sel$caseid , ] * rp3.pweightsM 
		  t_rprobs <- base::aperm( rprobs , dim=c(3,2,1) )
		  for (k1 in 1:maxK) {
		    M1_k1 <- t_rprobs[,k1,] * M1 
			r[,k1] <- base::colSums( M1_k1 , na.rm=TRUE)

			for (k2 in 1:maxK) {
			  rr[,k1,k2] <- base::colSums( M1_k1 * t_rprobs[,k2,] , na.rm=TRUE)
			}
		  }		  
		  

# cat("one iteration xsi -- matrix 1 (V2)") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 		  
		  
		  A_Sq <- AA_bari <- A_bari <- matrix( 0 , PP1 , nitems )
		  t_r <- base::t( r )
          t_rr <- base::aperm( rr , dim=c(3,2,1) ) 		  
		  t_A.0 <- base::aperm( A.0 , dim=c(3,2,1) )
		  for (kk in 1:maxK){ 
			A0_r_kk <- A.0[ , kk , ] * r[ , kk ]
			A_bari <- A_bari + t( A0_r_kk )
			# A_bari <- A_bari + t_A.0[ , kk , ] * t_r[ kk ,  ]
			# AA_bari <- AA_bari + t( A.0[ , kk , ]^2 * r[ , kk ] )		
			AA_bari <- AA_bari + t( A.0[ , kk , ] * A0_r_kk )		
		  }
		  for (kk1 in 1:maxK){ 
			# A0_rr_kk1 <- A.0[,kk1,] * rr[, 
			for (kk2 in 1:maxK){ 
			  fac <- 1
			  # A_Sq <- A_Sq + fac * t( A.0[,kk1,] * A.0[,kk2,] * rr[ , kk1 , kk2 ] )	
			  
			  if ( kk1 < kk2 ){ fac <- 2 }
			  if (kk1 <= kk2 ){
				# A_Sq <- A_Sq + fac * t( A.0[,kk1,] * A.0[,kk2,] * rr[ , kk1 , kk2 ] )	
				A_Sq <- A_Sq + fac * t( A.0[,kk1,] * A.0[,kk2,] * rr[ , kk1 , kk2 ] )	
				}
			}
		  }
# cat("one iteration xsi -- matrix 2 (V2)") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
 		  
			  # A							[ nitems , maxK , length(xsi) ]	
			  # A_Sq, AA_bari, A_bari		[ length(xsi) , nitems ]
			  # r							[ nitems , maxK ]
			  # rr						[ nitems , maxK , maxK ]	
			  #    for (p in p_loop ) {  # begin p loop
			  #      A_bari[p,] <- rowSums(A[,,p] * r, na.rm=TRUE)
			  #      AA_bari[p,] <- rowSums(A[,,p]^2 * r, na.rm=TRUE)
			  #      for (k1 in 1:maxK) {
			  # for (k2 in 1:maxK) {
			  # AA[,k1,k2] <- A[,k1,p]*A[,k2,p] * rr[,k1,k2]
			  # }
			  # }
			  # A_Sq[p,] <- apply(AA, 1, sum, na.rm=TRUE)
			  # }	# end parameter p loop
		  
		  expected <- rowSums (A_bari, na.rm=TRUE) # sum over items
		  err <- rowSums(AA_bari - A_Sq, na.rm=TRUE)   #sum over the items

# cat("one iteration xsi -- matrix 3") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 		  
		  
		  err_inv <- abs (1/( abs(err) + 10^(-10) ))		  
		  scores <- ItemScore * ( ! convergeP ) - expected
		  increment <-  err_inv*scores
		  ci <- ceiling( abs(increment) / ( abs( old_increment) + 10^(-10) ) )
		  increment <- ifelse( abs( increment) > abs(old_increment)  , 
							   increment/(2*ci) , 
							   increment )
		  increment[ nonest.xsi.index ] <- 0
		  xsi <- xsi + increment
		  old_increment <- increment      
		  # xsi[ nonest.xsi.index ] <- old_xsi[ nonest.xsi.index ]
		  
		  errorP <- sqrt(err_inv)
		  convergeP[ abs(increment) < convM ] <- TRUE
		  utils::flush.console()
		  iterP <- iterP + 1 
		  p_loop <- est.xsi.index[convergeP[est.xsi.index]==FALSE]
		  convergeAllP <- (sum(convergeP[est.xsi.index]) == length(est.xsi.index))  	  
		  cat("-")  
# cat("one iteration xsi -- update") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1 		  
    } # end of all parameters convergence

    res <- list( "xsi" = xsi , "errorP" = errorP, 
                 "maxChangeP" = max(abs( xsi - old_xsi ) ) )
    return (res)  
  }
