tam.fit <- function( tamobj, ... ){
  if(class(tamobj) == "tam.mml"){
    res <- tam.mml.fit( tamobj, ...)
  }
  if(class(tamobj) == "tam.jml"){
    res <- tam.jml.fit( tamobj, ...)
  }
  
  return( res )
}

tam.mml.fit <-
function( tamobj, FitMatrix=NULL , progress = FALSE ){
  #####################################################
  # INPUT:
  # tamobj ... result from tam analysis
  # FitMatrix is the fit design matrix. If it's NULL, then we will use the A matrix.
  # MW: We need to check whether this works when there are missing responses.
  # progress ... fit progress
  ####################################################
  #   
  resp <- tamobj$resp
  rprobs <- tamobj$rprobs
  
  indexIP.list <- tamobj$indexIP.list
  hwt <- tamobj$hwt
  resp.ind <- tamobj$resp.ind
  nnodes <- nrow(tamobj$theta)
  pweights <- tamobj$pweights
  nstud <- tamobj$nstud
  nitems <- tamobj$nitems
  maxK <- tamobj$maxK
  
  if ( is.null(FitMatrix) ) {
    FitMatrix <- tamobj$A 
  } 
  col.index <- rep( 1:nitems , each = maxK )
  cResp <- resp[ , col.index  ]*resp.ind[ , col.index ]
  cResp <- 1 * t( t(cResp) == rep(0:(maxK-1), nitems) )
  cF <- t( matrix( aperm( FitMatrix , c(2,1,3) ) , nrow = dim(FitMatrix)[3] , byrow = TRUE ) )
  cF[is.na(cF)] <- 0
  # sufficient statistics by person by parameter, for Fit design matrix.
  ParamScore <- (cResp %*% cF)
  
  np <- dim(FitMatrix)[3]  #The parameter dimension
  indexIP <- colSums( aperm( FitMatrix, c(2,1,3) ) != 0, na.rm = TRUE )
  
  # define list of elements for item parameters
  indexIP.list <- list( 1:np )
  for ( kk in 1:np ){ 
    indexIP.list[[kk]] <- which( indexIP[,kk] > 0 )
  }
  
  Outfit <- rep(0,np)
  Infit <- rep(0,np)
  Outfit_t <- rep(0,np)
  Infit_t <- rep(0,np)
  if (progress){
	cat(paste( "|" , paste(rep("*" , np ), collapse="") , "|\n|" ,sep="") )
				}
  for (p in 1:np) {
	if (progress){ cat("-") ; flush.console() }
    ip <- indexIP.list[[p]]
    xbari <- sapply( ip, function(i) colSums(FitMatrix[i,,p] * rprobs[i,,] , na.rm = TRUE ))
    #... TK: multiple category option -> na.rm = TRUE
    xxfi <- sapply( ip, function(i) colSums(FitMatrix[i,,p]^2 * rprobs[i,,] , na.rm = TRUE ))
    vari <- xxfi - xbari^2
	  xxxfi <- sapply( ip, function(i) colSums((FitMatrix[i,,p])^3 * rprobs[i,,] , na.rm = TRUE ))
	  xxxxfi <- sapply( ip, function(i) colSums(FitMatrix[i,,p]^4 * rprobs[i,,] , na.rm = TRUE ))
	  C4i <- xxxxfi - 4*xbari*xxxfi + 6*(xbari^2)*xxfi - 3*(xbari^4)
    Vz2i <- C4i - vari^2
    Uz2i <- C4i/(vari^2) - 1
	
    xbar <- resp.ind[,ip]%*%t( xbari )
    var1 <- resp.ind[,ip]%*%t( vari )
	  Vz2 <- resp.ind[,ip]%*%t( Vz2i )
	  Uz2 <- resp.ind[,ip]%*%t( Uz2i )

    Ax <- matrix(rep(ParamScore[,p],nnodes),nrow=nstud, ncol=nnodes)
	N <- nrow(hwt)
#	c_hwt<- aperm(apply (hwt,1,cumsum),c(2,1))
	c_hwt <- rowCumsums.TAM(hwt)

	
	rn1 <- runif(N)
	nthetal <- rep(1,ncol(c_hwt))
	
#    j <- apply(c_hwt,1, function (x) {findInterval(runif(1),x)})
#    j <- sapply(1:N, function (ii) {findInterval(rn1[ii],c_hwt[ii,])})
	j <- rowSums( c_hwt < outer( rn1 , nthetal ) )
	j <- rowSums( c_hwt < rn1)
	j[ j == 0 ] <- 1
	NW <- ncol(c_hwt)
    j <- j + 1
    j[ j > NW ] <- NW	
    s <- cbind(seq(1:N),j)
    wt_numer <- ( Ax[s] - xbar[s] )^2
    wt_denom <- var1[s]
    z2 <- wt_numer/wt_denom
    varz2 <- Uz2[s]
    wt_var <- Vz2[s]

  #Outfit MNSQ (unweighted fit)
  #calculate number of students per item parameter	
    nstud.ip <- sum( rowMeans( resp.ind[ , ip , drop=FALSE],na.rm = TRUE ), na.rm=TRUE )
#	  z2[z2 > 10*sd(z2)] <- 10*sd(z2)  #Trim extreme values
	  Outfit[p] <- sum( z2*pweights, na.rm = TRUE  ) / nstud.ip
	
  #Infit MNSQ (weighted fit)
    Infit[p] <- sum( wt_numer*pweights,na.rm = TRUE )/sum(wt_denom*pweights,na.rm = TRUE  )
	
  #Infit t
	  vf <- sum(wt_var*pweights,na.rm = TRUE )/(sum(wt_denom*pweights,na.rm = TRUE)^2 ) 
	  Infit_t[p] <- (Infit[p]^(1/3)-1) * 3/sqrt(vf) + sqrt(vf)/3   
    
  #Outfit t
	  vf2 <- sum(varz2*pweights,na.rm = TRUE )/(nstud.ip^2)
    Outfit_t[p] <- (Outfit[p]^(1/3)-1) * 3/sqrt(vf2) + sqrt(vf2)/3
  }
  if (progress){ cat("|\n") ; flush.console() }
  res <- data.frame("Outfit" = Outfit , 
				"Outfit_t" = Outfit_t, "Infit" = Infit , 
				"Infit_t" = Infit_t,   row.names = dimnames(FitMatrix)[[3]])
#data.frame( "Outfit" = round(Outfit,2) , "Outfit_t" = round(Outfit_t,1), "Infit" = round(Infit,2), Infit_t = round(Infit_t,1) )    
  return ( res )
}


tam.jml.fit <-
  function ( tamobj, resp , resp.ind, A, B, nstud, nitems, maxK, 
             ItemScore, theta, xsi, Msteps, pweightsM,
             est.xsi.index
  ){
    AXsi <- matrix(0, nrow=nitems, ncol=maxK) 
    BB <- array (0, dim=c(nitems,maxK))
    B_Sq <- array(0,dim=c(nstud, nitems))
    C4 <- array(0,dim=c(nitems,nstud))
    for (k in 1:maxK){ 
      AXsi[,k] <- A[,k,] %*% xsi
    } 
    B.0 <- B
    B.0[ is.na(B.0) ] <- 0					
    B1 <- B.0[,,1]					
    BB <- B1^2
    theta.unique <- unique( theta[,1] )
    NU <- length(theta.unique)
    B_bari <- array(0,dim=c(NU, nitems))
    BB_bari <- array(0, dim=c(NU, nitems))  
    res <- calc_prob.v5(iIndex = 1:nitems , A , AXsi , 
                        B , xsi , theta= matrix( theta.unique , nrow=NU , ncol=1) , 
                        NU, maxK , recalc=FALSE )        
    rprobs <- res[["rprobs"]] 
    #  rprobs <- rprobs[ , , match( theta[,1] , theta.unique)  ]
    rprobs[ is.na( rprobs) ] <- 0 
    # rprobs [ nitems , maxK , nstud ]
    # B1	[ nitems , maxK ]
    # BB	[ nitems , maxK ]
    # B_bari, BB_bari 	[nstud , nitems ]
    
    for (kk in 1:maxK){ 
      B_bari <- B_bari + t( B1[,kk]*rprobs[,kk,] )
      BB_bari <- BB_bari + t( BB[,kk] * rprobs[ , kk , ] )
    }
    ind.theta <- match( theta[,1] , theta.unique)				 
    rprobs <- rprobs[ , ,  ind.theta ]				 
    B_bari <- B_bari[ ind.theta , ] 
    BB_bari <- BB_bari[ ind.theta , ] 				 
    B_bari <- B_bari * resp.ind 		
    BB_bari <- BB_bari  * resp.ind
    #  B_bari <- sapply(1:nitems, function(i) colSums(B1[i,] * rprobs[i,,] , na.rm = TRUE)) * resp.ind
    #  BB_bari <- sapply(1:nitems, function(i) colSums(BB[i,] * rprobs[i,,] , na.rm = TRUE))  * resp.ind
    B_var <- BB_bari - (B_bari^2)
    z_sq <- (resp - B_bari)^2/B_var
    zsd <- sd(as.numeric(z_sq),na.rm=TRUE)
    z_sq[z_sq > 10*zsd] <-  10*zsd   #Trim extreme values
    B_bariM <- aperm(outer(B_bari,rep(1,maxK)),c(2,3,1))
    B1M <- outer(B1,rep(1,nstud))
    # B1M	[ nitems , maxK , nstud ]
    # B_bariM [ nitems , maxK , nstud ]
    #  C4 <- sapply(1:nitems, function(i) colSums((B1M[i,,]-B_bariM[i,,])^4 * rprobs[i,,] , na.rm = TRUE)) * resp.ind  
    # rprobs[i,,]  [ maxK , nstud ]
    # C4	[ nstud , nitems ]
    # B1M	[ nitems , maxK , nstud ]
    for (kk in 1:maxK){
      C4 <- C4 + ( B1M[,kk,] - B_bariM[,kk,] )^4 * rprobs[,kk,]
    }
    
    C4 <- t(C4) * resp.ind			 
    #  outfitPerson <- apply(z_sq, 1, mean, na.rm = TRUE)
    outfitPerson <- rowMeans( z_sq , na.rm=TRUE )
    #  outfitItem <- apply(z_sq * pweightsM, 2, mean, na.rm = TRUE)
    outfitItem <- colMeans(z_sq * pweightsM, na.rm = TRUE)
    
    infitPerson <- rowSums((resp - B_bari)^2, na.rm=TRUE)/rowSums(B_var, na.rm=TRUE)
    infitItem <- colSums((resp - B_bari)^2 * pweightsM, na.rm=TRUE)/colSums(B_var * pweightsM, na.rm=TRUE) 
    
    z_sq.ind <- !is.na(z_sq)
    #  var_outfit <- rowSums(C4/(B_var^2), na.rm=TRUE)/(nitems^2) - 1/nitems
    var_outfit <- rowSums(C4/(B_var^2), na.rm=TRUE)/(rowSums(z_sq.ind)^2) - 1/rowSums(z_sq.ind)
    outfitPerson_t <- (outfitPerson^(1/3) - 1) * (3/sqrt(var_outfit)) + sqrt(var_outfit)/3
    
    #  var_outfit <- colSums(C4/(B_var^2), na.rm=TRUE)/(nstud^2) - 1/nstud
    var_outfit <- colSums(C4/(B_var^2), na.rm=TRUE)/(colSums(z_sq.ind)^2) - 1/colSums(z_sq.ind)
    outfitItem_t <- (outfitItem^(1/3) - 1) * (3/sqrt(var_outfit)) + sqrt(var_outfit)/3
    
    var_infit <- rowSums(C4-B_var^2, na.rm=TRUE)/((rowSums(B_var, na.rm=TRUE))^2)
    infitPerson_t <- (infitPerson^(1/3) - 1) * (3/sqrt(var_infit)) + sqrt(var_infit)/3  
    
    var_infit <- colSums(C4-B_var^2, na.rm=TRUE)/((colSums(B_var, na.rm=TRUE))^2)
    infitItem_t <- (infitItem^(1/3) - 1) * (3/sqrt(var_infit)) + sqrt(var_infit)/3
    
    res <- list( "outfitPerson" = outfitPerson , "outfitItem" = outfitItem,
                 "infitPerson" = infitPerson , "infitItem" = infitItem,
                 "outfitPerson_t" = outfitPerson_t , "outfitItem_t" = outfitItem_t,
                 "infitPerson_t" = infitPerson_t , "infitItem_t" = infitItem_t)
    return (res)  
  }
