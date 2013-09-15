tam.wle <- function( tamobj, ... ){
  if(class(tamobj) == "tam.mml"){
    res <- tam.mml.wle( tamobj, ...)
  }
  if(class(tamobj) == "tam.jml"){
    res <- tam.jml.WLE( tamobj, ...)
  }
  
  return( res )
}

tam.mml.wle <-
function( tamobj, WLE=TRUE , adj=.3 , Msteps=20 , 
				convM = .0001 ){
  #########################################################
  # INPUT:
  # tamobj ... result from tam analysis
  # (WLE = TRUE) will produce WLE. Otherwise it will be MLE
  # 
  #########################################################
#  adj <- 0.3
#  Msteps <- 20
#  convM <- .0001
  B <- tamobj$B
  A <- tamobj$A
  nitems <- tamobj$nitems
  xsi <- ( tamobj$xsi )[,1]
  nstud <- tamobj$nstud
  AXsi <- tamobj$AXsi
  ndim <- tamobj$ndim
  maxK <- tamobj$maxK
  resp <- tamobj$resp
  resp[ is.na(resp) ] <- 0  
  resp.ind <- tamobj$resp.ind  
  col.index <- rep( 1:nitems , each = maxK )
  cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
  cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
  cB <- t( matrix( aperm( B , c(2,1,3) ) , nrow = dim(B)[3] , byrow = TRUE ) )
  cB[is.na(cB)] <- 0
  
  #Compute person sufficient statistics (total score on each dimension)
  PersonScores <- cResp %*% cB
  
  #Compute possible maximum score for each item on each dimension
  maxBi <- apply(B , 3 , rowMaxs , na.rm = TRUE)
  
  #Compute possible maximum score for each person on each dimension
  PersonMax <- resp.ind %*% maxBi
  PersonMax[ PersonMax == 0 ] <- 2 * adj
  
  #Adjust perfect scores for each person on each dimension
  PersonScores[PersonScores==PersonMax] <- PersonScores[PersonScores==PersonMax] - adj
  
  #Adjust zero scores for each person on each dimension
  PersonScores[PersonScores==0] <- PersonScores[PersonScores==0] + adj
  
  #Calculate Axsi. Only need to do this once.
  for (i in 1:nitems) {
    for (k in 1:maxK){
      AXsi[i,k] <- ( A[i,k,] %*% xsi )
    }
  }
 
  #Initialise theta (WLE) values for all students
  theta <- log(PersonScores/(PersonMax-PersonScores)) #log of odds ratio of raw score
  
  ######################################
  #Compute WLE
  #similar to the M step in the tam function, but each student's theta vector is now one node.
  converge <- FALSE
  Miter <- 0
  BB <- array (0, dim=c(nitems,maxK,ndim,ndim))
  BBB <- array (0, dim=c(nitems,maxK,ndim)) 
  for (i in 1:nitems) {
    for (k in 1:maxK) {
      BB[i,k,,] <- B[i,k,] %*% t(B[i,k,])
      BBB[i,k,] <- BB[i,k,,] %*% B[i,k,]
    }
  }
  increment <- array(0, dim=c(nstud,ndim))
  old_increment <- 3 + increment
  
  
  # Begin iterations
  while (!converge & ( Miter <= Msteps ) ) {  
    resWLE <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                           B , xsi , theta , nstud, maxK , recalc=FALSE )      	
    rprobsWLE <- resWLE[["rprobs"]] 
    B_bari <- array(0,dim=c(nstud, nitems, ndim))
    BB_bari <- array(0, dim=c(nstud, nitems, ndim, ndim))
    BBB_bari <- array(0,dim=c(nstud, nitems, ndim))
    for (d1 in 1:ndim) {
      B_bari[,,d1] <- sapply(1:nitems, function(i) colSums(B[i,,d1] * rprobsWLE[i,,] , na.rm = TRUE)) * resp.ind
      for (d2 in 1:ndim) {
        BB_bari[,,d1,d2] <- sapply(1:nitems, function(i) colSums(BB[i,,d1,d2] * rprobsWLE[i,,] , na.rm = TRUE)) *resp.ind
      }
      BBB_bari[,,d1] <- sapply(1:nitems, function(i) colSums(BBB[i,,d1] * rprobsWLE[i,,] , na.rm = TRUE)) *resp.ind  
    }
 
    B_Sq <- array(0,dim=c(nstud, nitems, ndim, ndim))
    B2_B <- array(0,dim=c(nstud, nitems, ndim))
    B_Cube <- array(0,dim=c(nstud, nitems, ndim))
    for (d1 in 1:ndim) {      
      B2_B[,,d1] <- 0
      B_Cube[,,d1] <- 0
      for (d2 in 1:ndim) {
        B_Sq[,,d1,d2] <- B_bari[,,d1]*B_bari[,,d2]
        B2_B[,,d1] <- B2_B[,,d1] + BB_bari[,,d1,d2]*B_bari[,,d2]   
        B_Cube[,,d1] <- B_Cube[,,d1] + B_Sq[,,d1,d2]*B_bari[,,d2]
      }
    }
    expected <- colSums(aperm(B_bari,c(2,1,3)))
    err <- colSums(aperm(BB_bari,c(2,1,3,4))) - colSums(aperm(B_Sq, c(2,1,3,4)))  #sum over the items
    if (ndim == 1) {
#      err_inv <- apply(err,1,function(x) 1/x )
      err_inv <- 1 / err
    } else {
      #err_inv <- aperm(apply(err,1,solve),c(2,1))   
      err_inv <- aperm(apply(err,1,function(ee){
						ee1 <- ee		
						diag(ee1) <- diag(ee1) + 10^(-15)
						solve(ee1)	
							}
					),c(2,1))            
    }
    err_inv <- array(abs(err_inv),dim=c(nstud,ndim,ndim))
    warm <- -3*B2_B + 2*B_Cube + BBB_bari
    warmadd <- colSums(aperm(warm,c(2,1,3)))  #sum over the items
    scores <- PersonScores - expected
    if (WLE) {
      warmaddon <- array(0,dim=c(nstud,ndim))
      for (d1 in 1:ndim) {
        warmaddon[,d1] <- 0
        for (d2 in 1:ndim) {
          warmaddon[,d1] <- warmaddon[,d1] + err_inv[,d1,d2]*warmadd[,d2]
        }
      }
      scores <- scores + warmaddon/2.0      
    }
    increment <- array(0, dim=c(nstud,ndim))
    for (d1 in 1:ndim) {
      increment[,d1] <- 0
      for (d2 in 1:ndim) {
        increment[,d1] <- increment[,d1] + err_inv[,d1,d2]*scores[,d2]
      }
    }
	# dampening the increment
    for ( d1 in 1:ndim){ 
#	   increment[,d1] <- ifelse( abs(increment[,d1]) > 3 , sign( increment[,d1] )*3 , increment[,d1] )
	  ci <- ceiling( abs(increment[,d1]) / ( abs( old_increment[,d1]) + 10^(-10) ) )	   
      increment[,d1] <- ifelse( abs( increment[,d1]) > abs(old_increment[,d1])  , 
								increment[,d1]/(2*ci) , 
								increment[,d1] )	   
	  old_increment[,d1] <- increment[,d1] 
	   #***
	   # avoid NaNs in increment
	   increment[,d1] <- ifelse( is.na(increment[,d1] ) , 0 , increment[,d1] )
		# increment[abs(increment)>3] <- sign(increment[abs(increment)>3])*3	
				}
    theta <- theta + increment
    if ( max(abs(increment)) < convM ) {
      converge <- TRUE
			}
    Miter <- Miter + 1 
    cat( paste( "Iteration in WLE/MLE estimation ", Miter, 
			"  | Maximal change " , round( max(abs(increment)) , 4) , "\n" )  ) 
    flush.console()
  }  # end of Newton-Raphson
  
  #standard errors of theta estimates
  if (ndim == 1) {
    error <- apply(err_inv,1,sqrt) 
  } else {    
    error <- aperm(apply(sqrt(err_inv),1,diag), c(2,1))
  }
  
  # The output contains 
  #   Person Scores on the test, by dimension
  #   Person possible maximum score, by dimension (Each person could take 
  #    different items, so possible maximum could vary)
  #   WLE or MLE estimate, by dimension
  #   Standard errors of WLE/MLE estimates, by dimension
  
  if ( ndim> 1){
	colnames(error) <- paste0("error.Dim" , substring( 100+1:ndim , 2) )
			}
  res <- data.frame( "pid" = tamobj$pid , 
				"N.items" = rowSums(resp.ind) , 
				"PersonScores" = PersonScores, 
				"PersonMax" = PersonMax, "theta" = theta , error )
					
  if (ndim==1){ colnames(res)[4:5] <- c("PersonMax" , "theta") }
  if (ndim>1){  
	colnames(res)[ 1:ndim + 2] <- paste0("PersonScores.Dim" , substring( 100+1:ndim , 2) )	
	ind <- grep( "theta" , colnames(res) )	
	colnames(res)[ind] <- 	paste0("theta.Dim" , substring( 100+1:ndim , 2) )	
		}
	####################
	# correct personMax set theta and standard error to missing		
	# if there are no observations on one dimension
	ind1 <- grep("PersonMax" , colnames(res))
    check1 <- ( res[ , ind1 , drop=FALSE] == 2*adj )
	ind2 <- grep("theta" , colnames(res))
	D <- length(ind1)
    for (ii in 1:D){
		res[ check1[,ii] , ind2[ii] ] <- NA
					}
	ind2 <- grep("error" , colnames(res))
	    for (ii in 1:D){
		res[ check1[,ii] , ind2[ii] ] <- NA
					}
  	#***
	# WLE reliability
	if ( ndim==1 ){
		ind <- which( res$N.items > 0 )
		v1 <- var( theta[ind] , na.rm=TRUE )	
		v2 <- mean( error[ind]^2 , na.rm=TRUE)
		# WLE_Rel = ( v1 - v2 ) / v1 = 1 - v2 / v1
		WLE.rel <- 1 - v2 / v1
	  cat("----\nWLE Reliability =" , round(WLE.rel,3) ,"\n" )
	  res$WLE.rel <- rep( WLE.rel , nrow(res) )
				}
	if ( ndim>1 ){
		cat("\n-------\n")
		for (dd in 1:ndim){
#	dd <- 1
		v1 <- var( res[,paste0("theta.Dim" , substring( 100+1:ndim , 2))[dd] ] , na.rm=TRUE)
		v2 <- mean( res[,paste0("error.Dim" , substring( 100+1:ndim , 2))[dd] ]^2 , na.rm=TRUE)
#		v2 <- mean( error^2 )
		res[ ,paste0("WLE.rel.Dim" , substring( 100+ dd , 2)) ]	<- h1 <- 1 - v2 / v1
	  cat(paste0("WLE Reliability (Dimension" , dd , ") = " , round(h1,3) ) , "\n" )
#	  res$WLE.rel <- rep( WLE.rel , nrow(res) )
					}
				}

				
#  res <- list( "PersonScores" = PersonScores, "PersonMax" = PersonMax, "theta" = theta , "error" =  error )
  return(res)
}

################################################################
################################################################
################################################################
tam.jml.WLE <-
  function ( tamobj, resp , resp.ind, A, B, nstud, nitems, maxK, convM, 
             PersonScores, theta, xsi, Msteps, WLE=FALSE
  ){
    
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    B1 <- B[,,1]
    BB <- array (0, dim=c(nitems,maxK))
    BBB <- array (0, dim=c(nitems,maxK)) 
    B_bari <- array(0,dim=c(nstud, nitems))
    BB_bari <- array(0, dim=c(nstud, nitems))
    BBB_bari <- array(0,dim=c(nstud, nitems))
    B_Sq <- array(0,dim=c(nstud, nitems))
    B2_B <- array(0,dim=c(nstud, nitems))
    B_Cube <- array(0,dim=nstud)
    
    #Calculate Axsi. Only need to do this once for ability estimates.
    for (i in 1:nitems) {
      for (k in 1:maxK){
        AXsi[i,k] <- ( A[i,k,] %*% xsi )
      }
    }
    cat("\n MLE/WLE estimation        |")
   
    #Compute WLE
    #similar to the M step in the tam function, but each student's theta is now one node.
    convergeWLE <- FALSE
    iterWLE <- 0
    BB <- B1^2
    BBB <- BB * B1  
    # BB	[ nitems , maxK ]
    # BBB	[ nitems , maxK ]
    maxChangeWLE <- 0
    thetaOld <- theta
    while (!convergeWLE & ( iterWLE <= Msteps ) ) {  
      resWLE <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                             B , xsi , theta , nstud, maxK , recalc=FALSE )      	
      rprobsWLE <- resWLE[["rprobs"]] 
      rprobsWLE[ is.na( rprobsWLE ) ] <- 0
      
      
      B_bari <- B1[,1] * rprobsWLE[ , 1 , ]
      BB_bari <- BB[,1] * rprobsWLE[ , 1 , ]
      for (kk in 2:maxK){ 
        B_bari <- B_bari + B1[,kk]*rprobsWLE[,kk,] 
        BB_bari <- BB_bari + BB[,kk] * rprobsWLE[ , kk , ]			
      }
      B_bari <- t(B_bari) * resp.ind 
      
      #    B_bari.OLD <- sapply(1:nitems, function(i) colSums(B1[i,] * rprobsWLE[i,,] , na.rm = TRUE)) * resp.ind
      # B1		[ nitems , maxK ]
      # rprobsWLE [ nitems , maxK , nstud ]
      # resp.ind  [ nstud , nitems ]	
      # colSums(B1[i,] * rprobsWLE[i,,] , na.rm = TRUE))
      #	-> colSums( [ maxK , nstud ] ) = [nstud]
      # B_bari	[ nstud , nitems ]
      
      BB_bari <- t(BB_bari ) * resp.ind
      B_Sq <- B_bari^2
      expected <- rowSums(B_bari, na.rm=TRUE)
      err <- rowSums(BB_bari, na.rm=TRUE) - rowSums(B_Sq, na.rm=TRUE)  #sum over the items
      err_inv <- abs(1/err)
      scores <- PersonScores - expected    
      
      if (WLE) {
        BBB_bari <- BBB[,1] * rprobsWLE[ , 1 , ]
        for (kk in 2:maxK){ 
          BBB_bari <- BBB_bari + BBB[,kk] * rprobsWLE[ , kk , ]			
        }	  
        BBB_bari <- t(BBB_bari ) * resp.ind				
        B2_B <- BB_bari*B_bari   
        B_Cube <- B_Sq*B_bari
        warm <- -3*B2_B + 2*B_Cube + BBB_bari
        warmadd <- rowSums(warm, na.rm=TRUE)                 #sum over the items
        warmaddon <- err_inv*warmadd
        scores <- scores + warmaddon/2.0      
      }     
      increment <-  err_inv*scores
      
      if (maxChangeWLE < max(abs(increment))) {
        maxChangeWLE <- max(abs(increment))
      }
      
      increment[abs(increment)>3] <- sign(increment[abs(increment)>3])*3
      
      theta <- theta + increment
      if ( max(abs(increment)) < convM ) {
        convergeWLE <- TRUE
      }
      iterWLE <- iterWLE + 1 
      cat( "-"  ) 
      flush.console()
    }  # end of Newton-Raphson   
    cat("\n")
    meanChangeWLE <- mean(theta - thetaOld)
    #standard errors of theta estimates
    errorWLE <- sqrt(err_inv)
	
    res <- list( "theta" = theta , "errorWLE" = errorWLE, "meanChangeWLE" = meanChangeWLE)
    return (res)
  }
