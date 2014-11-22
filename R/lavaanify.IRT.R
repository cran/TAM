
########################################################################
# lavaanify extension
lavaanify.IRT <- function( lavmodel , items=NULL , data = NULL){
#a0 <- Sys.time()
    res <- lavaanify.sirt.v1( lavmodel = lavmodel )		
# cat("sirt.v1 first") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1	
	ind <- grep( "__" , lavmodel )
	
	if ( length(ind) > 0 ){
	    if ( is.null(items) ){
              items <- colnames(data)
							}	
		res <- lavpartable.grep.underbrace( lavpartable=res$lavpartable , items )		
#cat("grep under") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
		res <- remove.duplicated.variances.lavsyn(res)				
#cat("remove dupl") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
		res <- lavaanify.sirt.v1( lavmodel = res)
#cat("sirt.v1 second") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
		lavpar <- res$lavpartable
		lavsyn <- res$lavaan.syntax
		lavpar0 <- lavpar
		ind1 <- which( paste(lavpar$op) == "?="  )
		ind2 <- which( !( paste(lavpar$rhs) %in% c("g1","s1")  ) )
		ind <- intersect( ind1 , ind2 )
		if ( length(ind) > 0 ){
		  lavpar0 <- lavpar0[ - ind , ]
						}
		res$lavpartable <- lavpar0
		lavsyn1 <- unlist( strsplit( lavsyn , "\n") )
		cn <- items
		v2 <- paste0( cn , "?=1*" , cn )
		v2 <- c(v2,paste0( cn , "?=0*" , cn ))
		lavsyn1 <- lavsyn1[ ! (lavsyn1 %in% v2 ) ]
		lavsyn1 <- paste0( lavsyn1 , collapse="\n" )
		res$lavaan.syntax <- lavsyn1
# cat("rest ind") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							
				}	
				
	# eliminate some entries of "?=" from parameter table
	# e.g. I1 ?= 1*I1
	lavpar0 <- res$lavpartable
	ind <- which( ( lavpar0$lhs == lavpar0$rhs ) & ( lavpar0$op == "?=" ) )
	if ( length(ind) > 0 ){	lavpar0 <- lavpar0[ -ind, ] }
	res$lavpartable <- lavpar0
# cat("all out") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1					
	return(res)
			}

###################################################
lavaanify.sirt.v1 <- function( lavmodel ){
# z0 <- Sys.time()
	syn <- lavmodel
	syn <- strsplit( syn , " " )[[1]]
	syn <- syn[ syn != "" ]
	syn <- gsub( ";" , "\n" , syn )
	#*****	
	syn <- split_syn_string( syn , "\\n" )	
	syn[ syn == "\\n" ] <- "\n"
	#***
	dfr1 <- data.frame( "index" = 1:length(syn) , "syntax"=syn )
	# look for specific strings and breaks
	dfr1$eqind <- 0
	N1 <- nrow(dfr1)
	vv <- 1
	for (ii in 1:N1){
	#    ii <- 1
		dfr1[ ii , "eqind" ] <- vv
		if ( length( grep( "\n" , dfr1$syntax[ii] ) ) > 0 ){ vv <- vv + 1 }
					}										
	syn0 <- lavmodel

	
	#***************************************************************
	# handling of guessing and slipping parameters						
	dfr1$guess_slip <- 0
	ind <- grep( "\\?=" , dfr1$syntax , perl=FALSE)
	if ( length(ind) > 0 ){
		dfr1$guess_slip[ ind ] <- 1
							}
	eqgroups <- dfr1$eqind[ which( dfr1$guess_slip == 1 ) ]
	dfr1$guess_slip[ dfr1$eqind %in% eqgroups ] <- 1
	
	# create "normal" lavaan syntax
	dfr2 <- dfr1[ dfr1$guess_slip == 0 , ]
	lavmodel1 <- paste0( dfr2$syntax , collapse="")
	lavpartable1 <- lavaan::lavaanify( as.character(lavmodel1 ) , warn = FALSE , debug=FALSE )
    # lavpartable1 <- lavaanify_in_sirt( as.character(lavmodel1 ) , warn = FALSE , debug=FALSE )	
	res1 <- change.grep.lavpartable( lavpartable1 )
# cat("**** change grep") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		
	# create a new model syntax here!!!
    if ( res1$changed ){
	   syn0 <- lavpartable2lavsyntax( res1$lavpartable )	   
# cat("**** lavpartable2lavsyntax") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1			   
	   lavpartable1 <- lavaan::lavaanify( as.character( syn0 ) , warn = FALSE , debug=FALSE )
# cat("**** lavaanify changed sirt.v1") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		   
						}
	
	
	# stop()
	# create lavaan parameter table for guessing/slipping parameters
	dfr2 <- dfr1[ dfr1$guess_slip == 1 , ]
	ug <- unique( dfr2$eqind)
	vecstr <- c("\\+" , "\\\n" , "\\?=" , "\\*" )

	for (uu in ug){
#		uu <- ug[1]
		syn.temp <- paste0( dfr2$syntax[ dfr2$eqind == uu ] , collapse="")
		syn.temp <- split_syn_string_vec( syn=syn.temp , vecstr = vecstr )
        syn.temp[ syn.temp == "g1" ] <- "t1"
        syn.temp[ syn.temp == "s1" ] <- "t2"
		syn.temp[ syn.temp == "\\?=" ] <- "|"
        syn.temp[ syn.temp == "\\+" ] <- "+"
		syn.temp[ syn.temp == "\\*" ] <- "*"
		syn.temp[ syn.temp == "\\\n" ] <- "\n"
        syn.temp <- paste0( syn.temp , collapse="")
		# h1 <- lavaanify_in_sirt( syn.temp)
		h1 <- lavaan::lavaanify( syn.temp)	
        h1$op <- "?="
        h1[ h1$rhs == "t1" ,"rhs"] <- "g1"		
        h1[ h1$rhs == "t2" ,"rhs"] <- "s1"	
		h0 <- h1
		h1$free <- h1$free + max(lavpartable1$free)
		h1$free[ h0$free == 0 ] <- 0
		h1$eq.id <- h1$eq.id + max(lavpartable1$eq.id)
		h1$eq.id[ h0$eq.id == 0 ] <- 0
		h1$unco <- h1$unco + max(lavpartable1$unco)
		h1$unco[ h0$unco == 0 ] <- 0
		h1$id <- h1$id + max(lavpartable1$id)
		lavpartable1 <- rbind( lavpartable1 , h1 )
				}										
#cat("**** guess / slip") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
	res <- list("lavpartable" = lavpartable1 , "lavaan.syntax"=syn0 )			
	return(res)	
			}
##################################################################			
			
#######################################################################
# create a lavaan syntax from lavaan parameter table
# This syntax only works for single groups.
lavpartable2lavsyntax <- function( lavpartable ){
    LL <- nrow(lavpartable )    
    syn0 <- paste0( lavpartable$lhs , lavpartable$op )
    lavpartable$prefix <- ""
    lavpartable$prefix <- ifelse( paste(lavpartable$label ) != "" , 
                paste0(lavpartable$label , "*"  ) , lavpartable$prefix ) 
    lavpartable$prefix <- ifelse( ( lavpartable$free == 0 ) & ( paste(lavpartable$label ) == "" ) , 
                paste0(lavpartable$ustart , "*"  ) , lavpartable$prefix )              
    syn0 <- paste0( syn0 , lavpartable$prefix , lavpartable$rhs  )
    syn0 <- paste0( syn0 , collapse="\n")
    return(syn0)
            }			
#######################################################################
			
################################################################
# change lavaan parameter label
change.grep.lavpartable <- function( lavpartable ){
    labels1 <- lavpartable$label
    labels1 <- sort( unique( labels1 ) )
    labels1 <- grep( "__" , labels1 , value=TRUE )
    LL <- length(labels1)
    if (LL > 0 ){
        for (ll in 1:LL){
            # ll <- 1
            label.group <- labels1[ll]
            lav.changed <- extend.label.group( label.group )
            ind.ll <- which( lavpartable$label %in% label.group )
			if ( length(ind.ll) == length(lav.changed) ){ 
				lavpartable[ ind.ll , "label" ] <- lav.changed 
									}
                        }
                }
	res <- list("lavpartable" = lavpartable , "changed" = LL > 0 )				
    return(res)
            }			

#######################################################
# extend labels arranged in a group
extend.label.group <- function( label.group ){
    str1 <- unlist( strsplit( label.group , split="__" ) )
    M1 <- min( nchar( str1) )
    ii <- 1
    for (mm in 1:M1){
        if ( substring( str1[1] , mm , mm) == substring( str1[2] , mm , mm) ){
            ii <- ii + 1 
                } else {
            break 
                    }
            }
    
    l1 <- as.numeric( substring( str1 , ii )    )
    str2 <- seq( l1[1] , l1[2] )
    str2 <- paste0( substring( str1 , 1 , ii - 1 )[1] , str2 )
    return( str2 )
        }
#######################################################
			
			
# remove duplicated variances
remove.duplicated.variances.lavsyn <- function( res0) {
	res0 <- gsub( " " , "" , res0 )
	res1 <- strsplit( res0 , split="\n")[[1]]
	res1 <- data.frame( "syn" = res1 , "sel" = 0 )
	res1$variance <- 0
	ind <- grep( "~~" , res1$syn )
	if ( length(ind) > 0 ){
		res1[ ind,"variance"] <- 1
		res1$variance.obs <- ""
		l1 <- res1[ ind, "syn" ]
		l1 <- strsplit( paste(l1) , split="~~")
		l2 <- lapply( l1 , FUN = function(ll){ ll[1] } )
		res1$variance.obs[ind] <- unlist(l2)
		l3 <- duplicated( res1[ind , "variance.obs"] )
		if ( sum(l3) > 0 ){
			res1 <- res1[ - ind[ l3 ] , ]
							}
					}
	# recreate  lavaan syntax
	lav2 <- paste( res1$syn , collapse="\n") 			
	return(lav2)
			}
			
#################################################################
# grep for "__" operator, meaning I01__I10
lavpartable.grep.underbrace <- function( lavpartable , items ){
    lav2 <- lavpartable
	LL <- nrow(lav2)
	syn <- NULL
    for (ll in 1:LL){
        # ll <- 13
        lav2.ll <- lav2[ll,]
		##*** ARb 2014-09-21: bug fix for "~1" operator
		if (lav2.ll$op == "~1" ){
			lav2.ll$rhs <- "1" 
			lav2.ll$op <- "~"
						}			
        ind.ll <- grep( "__" , c( lav2.ll$lhs  , lav2.ll$rhs  ) )
        v11 <- v1 <- lav2.ll$lhs            
            v10 <- strsplit( v1 , "__" )[[1]]  
            if ( length(v10) > 1 ){        
                v11 <- items[ seq( which( items == v10[1] ) , which( items == v10[2] )  ) ]            
                                    }        
        syn0 <- paste0( v11 , " " , lav2.ll$op , " ")
        g1 <- ifelse( lav2.ll$label != "" ,  paste0( lav2.ll$label , "*" ) , "" )
        g1 <- paste0( g1 , "" ,  ifelse( lav2.ll$free == 0 ,  paste0( lav2.ll$ustart , "*" ) , "" ) )
        
        syn0 <- paste0( syn0 , g1 , "" )            
        v11 <- v1 <- lav2.ll$rhs            
            v10 <- strsplit( v1 , "__" )[[1]]
			
            if ( length(v10) > 1 ){        
                v11 <- items[ seq( which( items == v10[1] ) , which( items == v10[2] )  ) ]            
                                    }
        syn0 <- paste0( syn0 , paste0( v11 , "\n")    )
        syn0 <- paste0( syn0 , collapse = "" )		
        syn <- paste0( syn , syn0 , collapse = "")   
	
                }   

    return(syn)
        }
########################################################################

#########################################################################		
# split_syn_string vectorized input
split_syn_string_vec <- function( syn, vecstr ){
	for (vv in vecstr){
    syn <- split_syn_string( syn , vv  )    
            }
	return(syn)
		}
########################################################################

		
		

########################################################################
# cleans syntax in a vector from strings vv
split_syn_string <- function( syn , vv ){
	syn <- as.list(syn )
	syn.vv <- grep( vv , syn )
	LL <- length(syn.vv)
	if (LL>0){
		for (ii in 1:LL){
			ll <- syn.vv[ii]
 			syn.ll <- syn[[ll]]		
			syn[[ll]] <- split_conc( syn.ll , vv )
							}
					}
	syn <- unlist(syn)
	return(syn)
			}
########################################################################

########################################################################
# splits a string syn.ll and concatanates it with string vv
split_conc <- function( syn.ll , vv ){
	g1 <- strsplit( syn.ll , vv , perl=FALSE )[[1]] 
	Lg1 <- length(g1)
	vec <- NULL
	if (Lg1 == 1 ){ vec <- c( g1 , vv ) }
	if (Lg1 > 1 ){
		vec <- rep("" , Lg1 + (Lg1-1) )
		vec[ seq( 1 , 2*Lg1 , 2 ) ] <- g1
		vec[ seq( 2 , 2*Lg1 - 1 , 2 ) ] <- vv	
		Ls1 <- nchar(syn.ll)
		if ( substring( syn.ll , Ls1 , Ls1 ) == vv ){
		       vec <- c( vec , vv )
							}
			}	
	return(vec)
			}
########################################################################					