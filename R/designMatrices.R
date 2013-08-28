designMatrices <-
  function( modeltype = c( "PCM" , "RSM" ) , 
            maxKi = NULL , resp = resp , ndim = 1 ,
            A = NULL , B = NULL , Q = NULL , R = NULL , ... ){
    
    modeltype <- match.arg(modeltype)
#a0 <- Sys.time();
    I <- ncol(resp)
    A.draft <- A	# if ! is.null(A), it is necessary
    if( is.null(maxKi) ){
      if( !is.null(resp) ){
        resp[is.na(resp)] <- 0
        maxKi <- apply( resp , 2 , max , na.rm=TRUE )
      } else 
        #... TK: 24.07.2012 -- check     
        if( !is.null(A) ){   
          np <- ncol(A)
          maxKi <- -colSums(A)
          maxKi <- maxKi[ - (which( (maxKi - 1) > 0) + maxKi[ which( (maxKi - 1) > 0) ]-1) ]
        } else return( warning("Not enough information to generate design matrices") )
    }
	# stop processing if there are items with a maximum score of 0
	i11 <- names(maxKi)[ maxKi == 0 ]
    if ( length(i11) > 0 ){
		stop( cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) )
					}
	
    nI <- length(maxKi)
    maxK <- max(maxKi)
    item <- rep( 1:nI , maxKi+1 )
#cat("g100"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1      
    if(modeltype %in%  c("PCM","RSM") ){
      
      cat <- unlist( lapply ( maxKi, seq , from=0 ) )
      np <- sum(maxKi)
      repnP <- cat[ cat != 0 ]
      revCat <- unlist( lapply ( maxKi, seq , to=1 ) )
      
      # Q Matrix
      if( is.null(Q) ){
        if( ndim > 1 ) warning("random q matrix")
        Q.draft <- matrix( 0 , nrow = nI , ncol = ndim )
        Q.draft[ cbind( 1:nI , sample(1:ndim, nI, replace=TRUE) ) ] <- 1 
      }else{
        Q.draft <- Q
      } 
      ndim <- dim(Q.draft)[2]
#cat("g150"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1            
      # B Matrix
      if( is.null(B) ){      
        B.draft <- array( 0 , dim = c(nI, maxK+1 , ndim) , 
                          dimnames = list(
                            #							paste( "Item", ifelse( 1:nI<10 , "0" , "" ), 1:nI , sep = "" ) , 
                            colnames(resp)  , 
                            paste( "Cat", 0:maxK , sep = "" ) , 
                            paste( "Dim" , ifelse( (1:ndim) < 10 , "0" , "" ), 1:ndim , sep = "")))
        
        for(dd in ndim){
          ind <- cbind( rep( 1:nI , maxKi+1 ) ,
                        cat + 1 , 
                        rep( dd , sum(maxKi+1) ) )
          B.draft[ind] <- cat*rep(Q.draft[,dd], maxKi+1)
        }
        
      } else { 
        B.draft <- B 
      }
      
      if ( ! is.null(Q) ){
        for (dd in 1:dim(Q)[2] ){
          for (zz in 1:dim(B.draft)[2] ){ 
            B.draft[ , zz , dd ] <- (zz-1)*Q[,dd]
          } 
        }
      }
#cat("g200"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1         
      # A Matrix
      if( is.null(A) || length( dim(A) ) < 3 ){
        A.draft <- array(, dim = c( nI, maxK+1 , np ) ,
                         dimnames = list( paste( "Item", ifelse( 1:nI<10 , "0" , "" ), 1:nI , sep = "" ) , 
                                          paste( "Category", 0:maxK , sep = "" ) , 
                                          paste( "Xsi", ifelse( (1:np) < 10 , "0" , "" ), 1:np , sep = "" ) ))
        
        ind.0 <- cbind( rep( item , np ) ,
                        rep( cat + 1 , np ) ,
                        rep( 1 : np , each = nI+np ) )
        
        i_pars <- cbind( pars_start <- rep(c(0, cumsum(maxKi))+1, c(maxKi, 0)),
                         pars_end <- pars_start+repnP-1 )
        
        ind.1 <- cbind( "item" = rep( item[ cat != 0 ] , repnP ), 
                        "category" = rep( repnP+1 , repnP) ,
                        "xsi" = unlist( apply( i_pars , 1 , 
								function(i_par) seq(from = i_par[1], to = i_par[2]) ) )
					)
        
        A.draft[ ind.0 ] <- 0
        A.draft[ ind.1 ] <- -1
        
        
        # item labels for Partial Credit Model
        l0 <- unlist(sapply( maxKi , FUN = function(cc){ seq( 1, cc)  } ))
        l1 <- paste( rep( colnames(resp) , maxKi) , "_Cat" , l0 , sep="" )
        dimnames(A.draft)[[3]] <- l1
        
        # item labels for Rasch model
        if (maxK == 1 ){
          dimnames(A.draft)[[3]] <- colnames(resp)
        }
        
      }    
      
    }

#cat("g300"); a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1       
    if(modeltype == "RSM"){
      if( is.null(A) ) 
#        return( warning("Not enough information to generate design matrices") )
#          nP <- sum(maxKi) + length(rater)
	Nxsi <- I + maxK - 1
	Kitem <- maxKi+1
	A <- array( 0 , dim=c( I , maxK+1 , Nxsi ) )
	vv <- 1
	for (ii in 1:I){
		A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
		if ( Kitem[ii] <= maxK ){
			A[ ii , ( Kitem[ii] + 1 ):(maxK+1) ,  ] <- NA                
							}
		vv <- vv+1
					}
    for (ii in 1:I){					
	  if ( Kitem[ii] > 2 ){
		for (kk in 1:(Kitem[ii] - 2) ){
			A[ ii , 1 + ( kk:(Kitem[ii]-2) ) , I+kk ] <- - 1
						}
					}
				}
				
	dimnames(A)[[1]] <- colnames(resp)
	vars <- colnames(resp)
	vars <- c(vars , paste0( "Cat"  , 1:(maxK-1) ) )
	dimnames(A)[[3]] <- vars			
	A.draft <- A


	
    }
    
    flatA <- t( matrix( aperm( A.draft , c(2,1,3) ) , nrow = dim(A.draft)[3] , byrow = TRUE ) )
    colnames(flatA) <- dimnames(A.draft)[[3]]
    
    flatB <- t( matrix( aperm( B.draft , c(2,1,3) ) , nrow = dim(B.draft)[3] , byrow = TRUE ) )
    colnames(flatB) <- dimnames(B.draft)[[3]]
    rownames(flatB) <- rownames(flatA) <- t(outer(dimnames(B.draft)[[1]] , 
					dimnames(B.draft)[[2]] , paste , sep ="."))
    
    out <- list( "item" = item , "maxKi" = maxKi , "cat" = cat , 
                 "A" = A.draft , "flatA" = flatA , "B" = B.draft , 
                 "flatB" = flatB , "Q" = Q , "R" = R )
    class(out) <- "designMatrices"
    return(out)
  }
  
#############################################################
print.designMatrices <-
  function( X , ... ){
    x <- X
    BB <- x$flatB
    colnames(BB) <- paste("B_", colnames(BB), sep ="")
    out <- cbind( x$flatA, BB )
    
    NAs <- apply( x$flatA , 1 , function(fA) all(is.na(fA)) )
    out <- out[!NAs, ]
    
    print(out)
    invisible( out )
  }


rownames.design <- function(X){
  Y <- apply(X, 2, as.numeric )
  Y <- sapply(1:ncol(Y), function(vv) 
    paste( colnames(Y)[vv], add.lead(Y[,vv], ceiling(log( max(as.numeric(Y[,vv])), 10)) ), sep ="" )
  )
  
  rownames(X) <- apply(Y, 1, paste, collapse = "-")
  return(X)
} 

rownames.design2 <- function(X){
  Y <- apply(X, 2, as.numeric )
  Y <- sapply(1:ncol(Y), function(vv) 
    # paste( colnames(Y)[vv], add.lead(Y[,vv], ceiling(log( max(as.numeric(Y[,vv])), 10)) ), sep ="" )
	paste( colnames(Y)[vv], add.lead(Y[,vv], 1) , sep ="" )
  )
  
  rownames(X) <- apply(Y, 1, paste, collapse = "-")
  return(X)
} 


###########################################################
.A.matrix <-
  function( resp, formulaA = ~ item + item*step, facets = NULL,  
            constraint = c("cases", "items") ){
    
	### redefine facets matrix
	facets0 <- facets
	NF <- length(facets)
	facet.list <- as.list( 1:NF )
	names(facet.list) <- colnames(facets)
	if (NF==0){ facet.list <- NULL }
	if (NF>0){
	for (ff in 1:NF){
#		ff <- 2
		uff <- sort( unique( facets[,ff] ) )
		facets[,ff] <- match( facets[,ff] , uff )
		facet.list[[ff]] <- data.frame(
			"facet.label" = paste0( colnames(facets)[ff] , uff ) , 
			"facet.index" = paste0( colnames(facets)[ff] , seq(1,length(uff) ) ) )
						}
					}
    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
    
    maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    maxK <- max( maxKi )
    nI <- ncol( resp )
    
	# stop processing if there are items with a maximum score of 0
	i11 <- names(maxKi)[ maxKi == 0 ]
    if ( length(i11) > 0 ){
		stop( cat( "Items with maximum score of 0:" , paste(i11 , collapse=" " ) ) )
					}
	
    tf <- terms( formulaA )	
    fvars <- as.vector( attr(tf,"variables"), mode = "character" )[-1]
#cat("fvars 212") ; print(fvars)	
    otherFacets <- setdiff( fvars, c("item", "step") )
    contr.list <- as.list( rep( "contr.sum", length(fvars) ) )
    names( contr.list ) <- fvars
	#****
	# ARb: 2013-03-27
	# no contrasts for items
	nitems <- ncol(resp)
#	contr.list[["item"]] <- diag(1,nitems)
#******	
    ### prepare data-Object for model.matrix()
    expand.list <- 
      as.vector( c( list( if( "item" %in% fvars ) factor(1:nI),
                          if( "step" %in% fvars ) factor(1:maxK) ),
                    if( length( otherFacets ) == 1){
                      list( factor( 1:max(facets[, otherFacets]) ) )
#  					list( factor( unique( facets[ , otherFacets] ) ) )
                    }else if( length( otherFacets ) > 1 ){
                      apply( as.matrix( facets[, otherFacets] ), 2, 
							function(ff) as.factor(1:max(ff)) 
									)
                    }                     
      ) )
	  	  
    expand.list <- expand.list[ !unlist( lapply(expand.list, is.null) ) ]
    names( expand.list ) <- fvars
	for (vv in seq(1 , length(expand.list) ) ){
		expand.list[[vv]] <- paste( expand.list[[vv]] ) 
					}
 					
    X <- rownames.design2( expand.grid(expand.list) )
    ### constraints and formulaA
    if( constraint == "cases" ) formulaA <- update.formula(formulaA, ~0+.)
	NX <- ncol(X)
	for (ff in 1:NX){
		uff <- length( unique(X[,ff] ) )
		if (uff==1){ cat(paste0("          - facet " ,
					colnames(X)[ff] , " does only have one level!" ) , "\n") }
				}


    mm <- - model.matrix(formulaA, X, contrasts = contr.list)
#    mm <- - model.matrix(formulaA, X )
    if( constraint == "items" ) mm <- mm[,-1]
	
	############################################################
	###*** ARb 2013-03-28
	### generate all interactions	
	xsi.constr <- .generate.interactions(X , facets , formulaA , mm )									
	###############################################################
		
    ### Postprocessing
    # model.matrix _ case: step in fvars
    if( "step" %in% fvars ){
      if( ncol( attr(tf, "factors") ) == 1 ){
        return( warning("Can't proceed the estimation: 
                        Factor of order 1 other than step must be specified.") )
      } 
      if( all( attr(tf, "factors")["step",] != 1 ) ){
        return( warning("Can't proceed the estimation: 
                        Lower-order term is missing.") )
      } 
      
      A <- NULL
      
      stepgroups <- unique( gsub( "-step([[:digit:]])*", "-step([[:digit:]])*", rownames(X) ) )
      X.out <- data.frame(as.matrix(X), stringsAsFactors = FALSE)

     
      for( sg in stepgroups ){
        # sg <- stepgroups[1]
        mm.sg.temp <- rbind( 0, apply( mm[ grep(sg, rownames(mm)) ,], 2, cumsum ) )
        rownames(mm.sg.temp)[1] <- gsub("-step([[:digit:]])*", "-step0", sg, fixed = TRUE)
        A <- rbind(A, mm.sg.temp)
        
        x.sg.temp <- X.out[grep(sg, rownames(X.out))[1], ]
        x.sg.temp[,"step"] <- 0
        rownames(x.sg.temp) <- gsub("-step([[:digit:]])*", "-step0", sg, fixed = TRUE)
        X.out <- rbind(X.out, x.sg.temp)
      } 
      
    } else 
      # model.matrix _ case: step not in fvars
    {
      
      rownames(mm) <- paste( rownames(X) , "-step1", sep = "")
      A <- mm
      
      for( kk in setdiff(0:maxK, 1) ){
        mm.k.temp <- mm*kk
        rownames(mm.k.temp) <- paste( rownames(X) , "-step", kk , sep ="")
        A <- rbind(A, mm.k.temp)
      }
      
      X.out <- expand.grid( c( expand.list, list("step"=factor(0:maxK)) ) )
      X.out <- rownames.design2( data.frame(as.matrix(X.out), stringsAsFactors = FALSE) )
      
    }# end step in fvars

	# facet design
	facet.design <- list( "facets" = facets , "facets.orig" = facets0 , 
			"facet.list" = facet.list[otherFacets])
	A <- A[ ! duplicated( rownames(A) ) , ]
    A <- A[order(rownames(A)), ,drop = FALSE]      
    X.out <- X.out[order(rownames(X.out)), ,drop = FALSE]
	
    return(list( "A"=A, "X"=X.out, "otherFacets"=otherFacets , "xsi.constr"=xsi.constr ,
			"facet.design" = facet.design ) )
  }


  
 
####################################################
# create ConQuest parametrization for 
# partial credit model
.A.PCM2 <- function( resp ){
	Kitem <- apply( resp , 2 , max , na.rm=T ) + 1
	maxK <- max(Kitem)
	I <- ncol(resp)
	Nxsi <- sum(Kitem) - I
	A <- array( 0 , dim=c( I , maxK , Nxsi ) )
	vv <- 1
	for (ii in 1:I){
		A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
		if ( Kitem[ii] < maxK ){
			A[ ii , ( Kitem[ii] + 1 ):maxK , ] <- NA                
							}
		vv <- vv+1
					}
	for (ii in 1:I){
	  if ( Kitem[ii] > 2 ){
		for (kk in 2:(Kitem[ii] - 1) ){
			A[ ii , kk:(Kitem[ii]-1) , vv ] <- - 1
			vv <- vv + 1
						}
					}
				}
	dimnames(A)[[1]] <- colnames(resp)
	vars <- colnames(resp)
	vars <- c(vars , unlist( sapply( (1:I)[Kitem>2] , FUN = function(ii){
		paste0( colnames(resp)[ii] , "_step" , 1:(Kitem[ii] - 2) ) } )	) )
	dimnames(A)[[3]] <- vars			
	return(A) 
	}
#############################################################



####################################################
# create ConQuest parametrization for 
# partial credit model
.A.PCM3 <- function( resp ){
	Kitem <- apply( resp , 2 , max , na.rm=T ) + 1
	maxK <- max(Kitem)
	I <- ncol(resp)	
	Nxsi <- I + sum( Kitem > 2 )
	A <- array( 0 , dim=c( I , maxK , Nxsi ) )
	vv <- 1
	for (ii in 1:I){
		A[ ii , 2:Kitem[ii] , vv ] <- - ( 2:Kitem[ii] - 1 )
		if ( Kitem[ii] < maxK ){
			A[ ii , ( Kitem[ii] + 1 ):maxK , ] <- NA                
							}
		vv <- vv+1
					}
	for (ii in 1:I){
	  if ( Kitem[ii] > 2 ){
			Kii <- Kitem[ii]-1
			A[ ii , 1:(Kii+1) , vv ] <- ( 0:Kii ) * ( Kii - ( 0:Kii) )						
			vv <- vv + 1
						}
				}
	dimnames(A)[[1]] <- colnames(resp)
	vars <- colnames(resp)
#	vars <- c(vars , unlist( sapply( 1:I , FUN = function(ii){
#		paste0( colnames(resp)[ii] , "_step" , 1:(Kitem[ii] - 2) ) } )	) )
	vars1 <- paste0( vars[ Kitem > 2 ] , "_disp" )
	vars <- c( vars , vars1 )
	dimnames(A)[[3]] <- vars			
	return(A) 
	}
#############################################################
