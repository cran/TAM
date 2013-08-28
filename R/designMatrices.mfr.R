
  
#########################################################################
designMatrices.mfr <-
  function( resp, formulaA = ~ item + item:step, facets = NULL,  
            constraint = c("cases", "items"), ndim = 1,
            Q=NULL, A=NULL, B=NULL , progress=FALSE ){
 z0 <- Sys.time()
    ### Basic Information and Initializations
    constraint <- match.arg(constraint)
	## restructure formulaA
	t1 <- attr( terms( formulaA ) , "term.labels" )
	t2 <- intersect( c("item" , "step" , "item:step") , t1 )
#	if ( length( t2 ) == 0 ){ t2 <- "item" 
#		cat("   ~ Included 'item' in model formula\n")
#			}
# print("s20")
	
	formulaA <- paste(  paste( c(t2 , setdiff(t1 , t2 ) ) , collapse= " + " ) )
	formulaA <- as.formula( paste( " ~ " , formulaA ) )	
	
	#********************************
	# change formate in facets
	FF <- ncol(facets)
	NFF <- nrow(facets)
	if (progress){ 
		cat( "        o Check facets\n") ; flush.console();
			}
	
	for (ff in 1:FF){			
		# ff <- 1
		cff <- nchar(facets[,ff] )
		Mff <- max(cff)
		sff <- paste( rep("_" , Mff ) , collapse="" )
		if( min(cff) < Mff ){
		    facets.ff0 <- facets[,ff]
#			facets[,ff] <- paste0( facets[,ff] , substring( sff , 1 , Mff - cff ) )
			facets[,ff] <- paste0( "_" , facets[,ff] , substring( sff , 1 , Mff - cff ) )
			if (progress){
    			u1 <- unique( setdiff( paste(facets[,ff]) , paste( facets.ff0 ) ) )			
				p1 <- paste0( "          * Changed levels of facet ", colnames(facets)[ff], ":" )
				p1 <- paste( p1 , paste( paste0("'",u1,"'") , collapse= " " ) )
				cat(p1, "\n")
						}
				}
			}

			
	#********************************	
#    resp[ is.na(resp) ] <- 0
    maxKi <- apply( resp , 2 , max , na.rm=TRUE )
    maxK <- max( maxKi )
    I <- nI <- ncol( resp )
    item <- rep( 1:nI , maxKi+1 )
    if ( is.null( colnames(resp) ) ){
		colnames(resp) <- paste0( "item" , 1:nI )
				}
					
    # A Matrix
    if( is.null(A) ){
      AX <- .A.matrix( resp, formulaA = formulaA, facets = facets, constraint = constraint )
	  
      A <- AX$A; X <- AX$X; otherFacets <- AX$otherFacets
	  xsi.constr <- AX$xsi.constr
	  facet.design <- AX$facet.design
	  facet.list <- facet.design$facet.list
	  facets <- facet.design$facets
      X.noStep <- unique(X[,- grep("step", colnames(X)), drop = FALSE ])     
      rownames(X.noStep) <- gsub("-step([[:digit:]])*", "", rownames(X.noStep))
		} 		
	if (progress){ 
		cat( "        o Created A Matrix\n") ; flush.console();
			}
    # Q Matrix
    if( is.null(Q) ){
      if( ndim > 1 ) warning("random q matrix")
      Q <- matrix( 0 , nrow = nI , ncol = ndim )
      Q[ cbind( 1:nI , sample(1:ndim, nI, replace=T) ) ] <- 1 
    } 
    Q <- Q[ as.numeric(X[,ifelse("item" %in% colnames(X), "item", 1)]) ,,drop=FALSE]
    dimnames(Q) <- list( rownames(X), 
                         paste( "Dim" , add.lead( 1:ndim, ceiling(log(ndim, 10)) ), sep ="") 
    )
    # ndim
    ndim <- dim(Q)[2]
    
    # B Matrix
    if( is.null(B) ){ 
      B <- Q * apply(X, c(1,2), as.numeric)[,"step"]
    }
	if (progress){ 
		cat( "        o Created B Matrix\n") ; flush.console();
			}	
    # gresp
    ind.resp.cols <- as.numeric(X[, ifelse("item" %in% colnames(X), "item", 1) ])

	
    gresp <- resp[,ind.resp.cols]	
#    gresp <- 1* ( gresp == t(array( X[, "step"], dim = dim(t(gresp)) )) )	
	gresp <- 1*(gresp==matrix( as.numeric(X[,"step"]) , nrow(gresp) , ncol(gresp) , byrow=TRUE ))
# cat("calc gresp (version 2)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    

    ind.resp.cols <- as.numeric(X.noStep[, ifelse("item" %in% colnames(X.noStep), "item", 1) ])
    gresp.noStep <- resp[,ind.resp.cols]
# cat("calc gresp.noStep" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1    
    if( length(otherFacets) > 0 ){
      rnFacets <- rownames( rownames.design2( as.matrix(facets[,otherFacets]) ))
      rnX <-      rownames( rownames.design2( as.matrix(X[,otherFacets]) ))
      rnX.noStep <-      rownames( rownames.design2( as.matrix(X.noStep[,otherFacets]) ))   
# cat("rownames.design2  X" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1		  
#*** ARb 2013-03-26:
#*** Set all entries in gresp and gresp.noStep to missing
#*** if they are not observed.
#      gresp <- gresp * (1* outer(rnFacets, rnX, "=="))
	   gresp[ outer(rnFacets, rnX, "!=") ] <- NA
#      gresp.noStep <- gresp.noStep * (1* outer(rnFacets, rnX.noStep, "=="))
       gresp.noStep[ outer(rnFacets, rnX.noStep, "!=") ] <- NA
# cat("gresp NA " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
   
    }        
    colnames(gresp) <- rownames(X)
    X$empty <- 1* (colSums( gresp, na.rm=TRUE ) == 0)
    colnames(gresp.noStep) <- rownames(X.noStep)
    X.noStep$empty <- 1* (colSums( gresp.noStep, na.rm=TRUE ) == 0)
    ### output
    ind <- X[,"empty"] == 1
    nStep <- maxK+1
    nGenit <- nrow(X)
# cat("gresp" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1
    
    .generate.3d <- function(x){      
      return( aperm( array( as.matrix(x), dim= c( nStep, nGenit/nStep, ncol(x) ), 
                            dimnames= list( paste("_step",0:maxK, sep= ""), 
                                            unique(gsub("-step([[:digit:]])*", "", rownames(x))),
                                            colnames(x) ) )
                     , c(2,1,3) 
						) )
					}
	# generate B
    .generateB.3d <- function(x){  
		x2 <- array( 0 , c(nStep , nGenit/nStep , ncol(x) ) )
        dimnames(x2) <- list( paste("_step",0:maxK, sep= ""), 
                                            unique(gsub("-step([[:digit:]])*", "", rownames(x))),
                                            colnames(x) )
		for (dd in seq(1 , ncol(x) ) ){
		  for (ss in 0:(nStep-1)){
		    str.ss <- paste0("-step",ss )
			iss <- grep(  str.ss , rownames(x) , fixed=TRUE )
			str.ss2 <- gsub( str.ss , "" , rownames(x)[iss] )
			x2[ss+1,str.ss2,dd] <- as.vector(x[ iss , dd])
#			x2[ ss,,dd] <- matrix( x[,dd] , nrow=nStep , ncol= dim(x2)[2] , byrow=TRUE)
							}
						}
		x2 <- aperm( x2 , c(2,1,3) )
		return(x2)
					}
	#***
	# debugging ind manually
	ind <- FALSE * ind
	#*************
	# rename items
	itemren <- data.frame( "item" = colnames(resp) , "itemren" = paste0( "item" , 1:nI ) )
# cat(".....\nbefore rename" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
# print("g100")	
	A <- .rename.items( matr=A , itemren )
# cat(".rename.items (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1		
    xsi.table <- xsi.constr$xsi.table
#	A <- .rename.items3( matr=A , facet.list , I )	
#cat(".rename.items3 (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
	A <- .rename.items3a( matr=A , facet.list , I , cols=TRUE , xsi.table )	
# cat(".rename.items3a (A)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
	B <- .rename.items( matr=B , itemren )	
# cat(".rename.items (B)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1			
    dimnames(B)[[1]] <- dimnames(A)[[1]]	
#	B <- .rename.items3( matr=B , facet.list )	
	gresp <- t( .rename.items( matr=t(gresp) , itemren , cols=FALSE)	)
# cat(".rename.items (gresp)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1				
#	gresp <- t( .rename.items3( matr=t(gresp) , facet.list , cols=FALSE)	)	
	dimnames(gresp)[[2]] <- dimnames(A)[[1]]	
	gresp.noStep <- t( .rename.items( matr=t(gresp.noStep) , itemren , cols=FALSE)	)	
# cat("h2" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
#	gresp.noStep <- t( .rename.items3( matr=t(gresp.noStep) , facet.list , I , cols=FALSE)	)	
	gresp.noStep <- t( .rename.items3a( matr=t(gresp.noStep) , facet.list , I , cols=FALSE ,
				xsi.table )	)		
 #cat(".rename.items (gresp.noStep)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
	Q <- .rename.items( matr=Q , itemren , cols=FALSE)
	dimnames(Q)[[1]] <- dimnames(A)[[1]]
# cat(".rename.items (Q)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1						
# 	Q <- .rename.items3( matr=Q , facet.list , cols=FALSE)
	X <- .rename.items( matr=X , itemren , cols=FALSE)
	dimnames(X)[[1]] <- dimnames(A)[[1]]
#	X <- .rename.items3( matr=X , facet.list , cols=FALSE)	
 # cat(".rename.items (Q,X)" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1					
	X.noStep <- .rename.items( matr=X.noStep , itemren , cols=FALSE)
#	X.noStep <- .rename.items3( matr=X.noStep , facet.list , cols=FALSE)	
	#***
	G1 <- xsi.constr$xsi.table 	
	G1$parameter <- .rename.items2( paste( G1$parameter) , itemren) 	
# cat(".rename.items2 (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1							
#	G1$parameter <- .rename.items2a( paste( G1$parameter) , facet.list , I) 	
#cat(".rename.items2a (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1						
	G1$parameter <- .rename.items2b( paste( G1$parameter) , facet.list , I , xsi.table ) 	
	xsi.constr$xsi.table <- G1	
# cat(".rename.items2b (G1$parameter)  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1							
	#***
	G1 <- xsi.constr$xsi.constraints
	rownames(G1) <- .rename.items2( rownames(G1) , itemren) 	
	colnames(G1) <- .rename.items2( colnames(G1) , itemren) 
# cat(".rename.items2 (colnames(G1))  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1							
	colnames(G1) <- .rename.items2b( colnames(G1) , facet.list , I , xsi.table , sel1=1) 
	rownames(G1) <- .rename.items2b( rownames(G1) , facet.list , I , xsi.table , sel1=2) 		
 #cat(".rename.items2a (colnames(G1))  " ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1								
	G1 -> xsi.constr$xsi.constraints
# cat("rename items" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1

	if (progress){ 
		cat( "        o Relabeled Variable Names\n") ; flush.console();
			}
	
    # A
    A.flat.0 <- A.flat <- A; A.flat.0[ind,] <- 0
    A.3d <- .generateB.3d( A.flat )
    A.flat <- A.flat[!ind,]
    A.3d.0 <- .generateB.3d( A.flat.0 )
    # B                  
    B.flat.0 <- B.flat <- B; B.flat.0[ind,] <- 0
    B.3d <- .generateB.3d( B.flat )
    B.flat <- B.flat[!ind,]
    B.3d.0 <- .generateB.3d( B.flat.0 )
   
    # Q                  
    Q.flat.0 <- Q.flat <- Q; Q.flat.0[ind,] <- 0
    Q.3d <- .generateB.3d( Q.flat )
    Q.flat <- Q.flat[!ind,]
    Q.3d.0 <- .generateB.3d( Q.flat.0 )     
# cat("output mfr" ) ; z1 <- Sys.time() ; print(z1-z0) ; z0 <- z1	
    # out
    out <- list( "gresp" = list("gresp"=gresp, "gresp.noStep"=gresp.noStep), 
                 "A" = list("A.flat"=A.flat, "A.flat.0"=A.flat.0, 
                            "A.3d"=A.3d, "A.3d.0"=A.3d.0), 
                 "B" = list("B.flat"=B.flat, "B.flat.0"=B.flat.0, 
                            "B.3d"=B.3d, "B.3d.0"=B.3d.0), 
                 "Q" = list("Q.flat"=Q.flat, "Q.flat.0"=Q.flat.0,
                            "Q.3d"=Q.3d, "Q.3d.0"=Q.3d.0), 
                 "X" = list("X"=X, "X.noStep"=X.noStep) ,
				 "xsi.constr" = xsi.constr 
    )
    class(out) <- "designMatrices.mfr"
    return(out)
  }

  

.generate.interactions <- function(X , facets , formulaA , mm ){	
	d1 <- d0 <- X	
	h1 <- sapply( colnames(d1) , FUN = function(vv){
				length(grep( vv , paste(formulaA) )) } )
	h1 <- colnames(d1)[ h1 == 0 ]
	d0 <- d0[ , ! ( colnames(d1) %in% h1 ) , drop=FALSE]
	M2 <- model.matrix( object= formulaA , data= d1 , 
			contrasts.arg = lapply( d0 , contrasts, contrasts=FALSE) )
	h2 <- colnames(M2)
	h1 <- colnames(mm)
	# extract facets
	xsi.table <- data.frame( "parameter" = h2 )
	xsi.split <- sapply( xsi.table$parameter , FUN = function(ll){ 
		l1 <- as.vector( unlist( strsplit( paste(ll) , split=":" ) ) )
		v1 <- l1
		for (ii in 1:length(l1) ){
			for (cc in colnames(X) ){ 
				kk <- grep( cc , l1[ii] )
				if (length(kk)>0){ v1[ii] <- cc }
								}
						} 
		v1 <- paste0( v1 , collapse=":" )
		return(v1)		
					} )
	xsi.table$facet <- unlist(xsi.split)
	xsi.table$facet.order <- sapply( xsi.table$parameter , FUN = function(ll){ 
		length( as.vector( unlist( strsplit( paste(ll) , split=":" ) ) ) ) } )
	xsi.table$constraint <- 1 - 1*(xsi.table$parameter %in% h1)
	xsi.table$facet.index <- match( xsi.table$facet , unique( xsi.table$facet ) )
#	xsi.table$orig.index <- seq(1,nrow(xsi.table))
#	xsi.table[ order( paste( xsi.table$facet.index+100 , xsi.table$parameter ) ) , ]

	facets.unique <- unique( xsi.table$facet )
	b1 <- xsi.table[ xsi.table$constraint == 1 , "parameter" ]
	c1 <- xsi.table[ xsi.table$constraint == 0 , "parameter" ]	
    xsi.constraints <- matrix( NA  , nrow=length(b1) , ncol=length(c1) )
	rownames(xsi.constraints) <- paste(b1)
	colnames(xsi.constraints) <- paste(c1)
# b1 <- b1[3]
	############################
	# loop over terms
    for (bb in b1 ){
		#bb <- b1[9]
		v1 <- 0
		mult <- 1
		xsi.table.bb <- xsi.table[ xsi.table$parameter == bb , ]
		x0 <- x1 <- xsi.table[ xsi.table$facet %in% xsi.table.bb$facet , ]
		if ( xsi.table.bb$facet.order==1){ 
			xsi.constraints[paste(bb),] <- 0		
			xsi.constraints[ paste(bb) , paste( x1[ x1$constraint == 0 , "parameter" ] ) ] <- -1
										}
		if ( xsi.table.bb$facet=="item:step"){ 
			v1 <- 1
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					ll[1] == s2[1] } ) )
			x1 <- x1[ g1 , ]	
			mult <- 1
# cat("......",bb,"......\n")
# print(x1)			
			varsc <- paste( x1[ x1$constraint == 0 , "parameter" ] )
			if ( length(varsc) == 0){
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						ll[2] == s2[2] } ) )
				x1 <- x0[ g1 , ]	
				varsc <- paste(x1[ x1$constraint == 0 , "parameter" ])
				mult <- 1
				if ( length(varsc) == 0){
					varsc <- x1[ , "parameter" ]
					varsc <- setdiff( varsc , paste(bb) )
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					mult <- -1
									}				
						}
			
			xsi.constraints[ paste(bb) , varsc ] <- -1*mult
										}
		##################
		### order 2
		if ( xsi.table.bb$facet.order==2 & v1 ==0 ){ 							
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )	
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					ll[2] == s2[2] } ) )
			x1 <- x1[ g1 , ]	
			varsc <- x1[ x1$constraint == 0 , "parameter" ]
			if ( length(varsc) == 0){
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						ll[1] == s2[1] } ) )
				x1 <- x0[ g1 , ]	
				varsc <- x1[ x1$constraint == 0 , "parameter" ]	
				if ( length(varsc) == 0){
					varsc <- x1[ , "parameter" ]
					varsc <- setdiff( varsc , paste(bb) )
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					mult <- -1
									}
								}							
			xsi.constraints[ paste(bb) , paste( varsc ) ] <- -1	* mult									
										}
		#########################
		### order 3
		if ( xsi.table.bb$facet.order==3 & v1 ==0 ){ 			
			mult <- 1	
			xsi.constraints[paste(bb),] <- 0
			s2 <- unlist( strsplit( paste(xsi.table.bb$parameter) , split=":" ) )	
			s20 <- strsplit( paste(x1$parameter) , split=":" )		
			g1 <- unlist( lapply( s20 , FUN = function(ll){
					( ll[2] == s2[2] ) & (ll[3] == s2[3] )  } ) )
			x1 <- x1[ g1 , ]
			varsc <- x1[ x1$constraint == 0 , "parameter" ]
			if ( length(varsc) == 0 ){			
				g1 <- unlist( lapply( s20 , FUN = function(ll){
						( ll[1] == s2[1] ) & (ll[3]==s2[3]) } ) )
				x1 <- x0[ g1 , ]	
				varsc <- x1[ x1$constraint == 0 , "parameter" ]		
				
				if ( length(varsc) == 0 ){			
					varsc <- x1[ , "parameter" ]			
					varsc <- setdiff( varsc , paste(bb) )	
					h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)
					varsc <- names(h1)[  h1 != 0 ]
					varsc <- na.omit( varsc)
					mult <- -1			
						if ( length(varsc) == 0 ){
							g1 <- unlist( lapply( s20 , FUN = function(ll){
									( ll[1] == s2[1] ) & (ll[2]==s2[2]) } ) )
							x1 <- x0[ g1 , ]	
							varsc <- x1[ x1$constraint == 0 , "parameter" ]		
							mult <- 1	
							
							if ( length(varsc) == 0 ){ 	
#							varsc <- setdiff( varsc , paste(bb) )		
								varsc <- x1[ , "parameter" ]
   							    varsc <- setdiff( varsc , paste(bb) )	
						        h1 <- colSums( xsi.constraints[ varsc , , drop=FALSE]	)					
								varsc <- names(h1)[  h1 != 0 ]
								varsc <- na.omit( varsc)
								mult <- -1					
													}
						}
										}																
							}		
						
			if ( length(varsc) > 0 ){
				xsi.constraints[ paste(bb) , paste( varsc ) ] <- -1	* mult	
					} else {
				xsi.constraints[ paste(bb) , ] <- NA	
						}
	
		
							}
						}
		xsi.constraints[ rowSums( abs(xsi.constraints) ) == 0 , ] <- NA
		res <- list( "xsi.constraints" = xsi.constraints , "xsi.table" = xsi.table )
#print(res) ;  stop("here")		
		
		return(res)
}						
######################
# rename item names
.rename.items <- function( matr , itemren , cols=TRUE ){
	rM <- rownames(matr)
	cM <- colnames(matr)
	I <- nrow(itemren)
	for ( ii in 1:I){
		rM <- gsub( paste0( itemren[ii,2] , "-") , paste0( itemren[ii,1] , "-") , rM )
		if (cols){
			cM <- gsub( paste0( itemren[ii,2] , ":") , paste0( itemren[ii,1] , ":") , cM )	
			cM[ cM == itemren[ii,2] ] <- paste(itemren[ii,1])
				}
					}
	rM -> rownames(matr) 
	if ( cols){ cM -> colnames(matr) }
	return(matr)
		}
#############################################################		
.rename.items2 <- function( vec , itemren ){
	cM <- vec
	I <- nrow(itemren)
	for ( ii in 1:I){
			cM <- gsub( paste0( itemren[ii,2] , ":") , paste0( itemren[ii,1] , ":") , cM )	
			cM[ cM == itemren[ii,2] ] <- paste(itemren[ii,1])
				}
	return(cM)
		}
#############################################################
.rename.items3 <- function( matr , facet.list , I , cols=TRUE  ){
### check for equalities in rM and cM in all entries!!!!
	rM <- rownames(matr)
	rMsplit <- strsplit( rM , split="-" )	
	RR <- length(rMsplit)
	FF <- length(facet.list)
	for (rr in 1:RR){
		rr1 <- rMsplit[[rr]]
		if (FF > 0 ){
		for (ff in 1:FF){ # begin ff
#		for (ff in seq(1,FF,1) ){
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)
				for (ii in 1:I ){ 
					rr1[ rr1 == itemren[ii,2] ] <- paste(itemren[ii,1])
					rMsplit[[rr]] <- rr1
								}
							} # end ff
							} # end if FF
					}
	rM <- unlist( lapply( rMsplit , FUN = function(ll){ paste( ll , collapse="-") } )	)
	rownames(matr) <- rM	
	#****************************************
	if ( cols){
		cM <- colnames(matr)
		cMsplit <- strsplit( cM , split=":" )	
		RR <- length(cMsplit)
		FF <- length(facet.list)
		for (rr in 1:RR){
			rr1 <- cMsplit[[rr]]
			if (FF>0){
			for (ff in 1:FF){ # begin ff
				itemren <- facet.list[[ff]]
				I <- nrow(itemren)
					for (ii in 1:I ){ 
						rr1[ rr1 == itemren[ii,2] ] <- paste(itemren[ii,1])
						cMsplit[[rr]] <- rr1
									}
								} # end ff
							}
						}
		cM <- unlist( lapply( cMsplit , FUN = function(ll){ paste( ll , collapse=":") } )	)
		colnames(matr) <- cM	
			}
	return(matr)
		}
#############################################################
.rename.items2a <- function( vec , facet.list , I ){
### check for equalities!!!
	cM <- vec
	FF <- length(facet.list)
	rM <- cM
    if ( ! is.null(rM) ){ 
		rMsplit <- strsplit( rM , split=":" )	
		RR <- length(rMsplit)
		FF <- length(facet.list)
		for (rr in 1:RR){
			rr1 <- rMsplit[[rr]]
			if (FF>0){
			for (ff in 1:FF){
				itemren <- facet.list[[ff]]
				I <- nrow(itemren)
					for (ii in 1:I ){ 
						rr1[ rr1 == itemren[ii,2] ] <- paste(itemren[ii,1])
						rMsplit[[rr]] <- rr1
									}
								}
							}
						}	
		rM <- unlist( lapply( rMsplit , FUN = function(ll){ paste( ll , collapse=":") } )	)
		cM <- rM
			}
	return(cM)
		}
		
		
#############################################################		
#############################################################
.rename.items3a <- function( matr , facet.list , I , cols=TRUE ,
			xsi.table ){
### check for equalities in rM and cM in all entries!!!!
	rM <- rownames(matr)
	rMsplit <- strsplit( rM , split="-" )
	RR <- length(rMsplit)
	rMM <- matrix( unlist(rMsplit) , nrow=RR , byrow=TRUE)
	rMM.ncol <- ncol(rMM)
	FF <- length(facet.list)	

	if (FF>0){ 
	for (ff in 1:FF){ # ff <- 1
		itemren <- facet.list[[ff]]
		I <- nrow(itemren)		
		for (ii in 1:I){ # ii <- 1
		if (rMM.ncol>1){
			for (kk in 2:rMM.ncol){# kk <- 3
			rMM[ rMM[,kk] == itemren[ii,2] , kk ] <- paste(itemren[ii,1])
							}
						}
						}
					}
				}   # end if FF>0	
	rM <- rMM[,1]
    if (rMM.ncol>1){
		for (rr in 2:rMM.ncol){ rM <- paste( rM , rMM[,rr] , sep="-" ) 	}
				}
	#****************************************
	if ( cols){
		rM <- colnames(matr)
		rMsplit <- unlist( strsplit( rM , split=":" ) )
	    xsi.table <- xsi.table[xsi.table$constraint==0,]
		XT <- nrow(xsi.table)
		F0 <- max(xsi.table$facet.order)
		index <- sapply( 1:XT , FUN = function(xx){
				m1 <- cbind( xx , 1:xsi.table[xx,"facet.order"] )
				matrix( t(m1) , ncol=1 , byrow=FALSE)
								} )					
		index <- matrix( unlist(index) , ncol=2 , byrow=T)
		rMMsub <- matrix("" , nrow=XT, ncol=F0) 
		rMMsub[ index ] <- rMsplit
		FF <- length(facet.list)	
		if (FF>0){ 
		for (ff in 1:FF){ # ff <- 1
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)		
			for (ii in 1:I){ # ii <- 1
			for (kk in 1:F0){# kk <- 3
				rMMsub[ rMMsub[,kk] == itemren[ii,2] , kk ] <- paste(itemren[ii,1])
								}
							}
						}
					}	# end if FF>0
		rM <- unlist( sapply( 1:XT , FUN = function(kk){
			paste( rMMsub[ kk , seq(1 , xsi.table$facet.order[kk] ) ] , collapse=":" ) } )
						)
		colnames(matr) <- rM					
			}
	return(matr)
		}
################################################################################		
################################################################################
.rename.items2b <- function( vec , facet.list , I , xsi.table , sel1=0 ){
### check for equalities!!!
#	cM <- vec	
	rM <- vec
	if ( ! is.null(rM)){
		rMsplit <- unlist( strsplit( rM , split=":" ) )
		if (sel1==1){ xsi.table <- xsi.table[xsi.table$constraint==0,] }
		if (sel1==2){ xsi.table <- xsi.table[xsi.table$constraint==1,] }		
		XT <- nrow(xsi.table)
		F0 <- max(xsi.table$facet.order)
		index <- sapply( 1:XT , FUN = function(xx){
				m1 <- cbind( xx , 1:xsi.table[xx,"facet.order"] )
				matrix( t(m1) , ncol=1 , byrow=FALSE)
								} )					
		index <- matrix( unlist(index) , ncol=2 , byrow=T)
		rMMsub <- matrix("" , nrow=XT, ncol=F0) 
		rMMsub[ index ] <- rMsplit
		FF <- length(facet.list)	
		if (FF>0){ 
		for (ff in 1:FF){ # ff <- 1
			itemren <- facet.list[[ff]]
			I <- nrow(itemren)		
			for (ii in 1:I){ # ii <- 1
			for (kk in 1:F0){# kk <- 3
				rMMsub[ rMMsub[,kk] == itemren[ii,2] , kk ] <- paste(itemren[ii,1])
								}
							}
						}
					}	
		rM <- unlist( sapply( 1:XT , FUN = function(kk){
			paste( rMMsub[ kk , seq(1 , xsi.table$facet.order[kk] ) ] , collapse=":" ) } )
						)
				}
	
	return(rM)
		}
		
		
