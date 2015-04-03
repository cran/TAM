
############################################################
mfr.dataprep <- function( formulaA , xsi.setnull , B , Q ,
		resp, pid, facets , beta.fixed ){
	
	tA <- terms( formulaA )
	tlab <- attr(tA , "factors")
	
	# redefine formula
	stlab <- apply( tlab , 1 , sum )		
	ind <- which( stlab == 0 )
	nullfacets <- names(stlab)[ind]

	#********************
	# data restructuring for non-identifiable combinations
	facets_labs <- setdiff( rownames(tlab) , c("item" , "step") )
	
	# create combination of pid and facets
	FF <- length(facets_labs)	
	combi <- pid
	if ( FF>0){
		for (ff in 1:FF){
			combi <- paste0( combi , "-" , facets[ , facets_labs[ff] ] )
						}
					}
	dups_combi <- any( duplicated( combi ) )
    PSF <- FALSE	
    if ( dups_combi ){	
		NC <- max( table( table( combi) ) )		
		l1 <- sapply( unique(combi) , FUN = function(cc){
				  N1 <- which( combi == cc ) 
				  m1 <- t( cbind( N1 , seq(1,length(N1) ) ) )
				  return(m1)
							} )			
		l1 <- matrix( unlist(l1) , ncol=2 , byrow=TRUE)
#		facets$psf <- paste0("PF",10^(round(log(NC,10)+1 )) + l1[,2] )
		N1 <- nchar( paste(max(table(combi))))	
		facets$psf <- paste0("PF", 10^N1 + l1[,2] )			
		nullfacets <- c( nullfacets , "psf" )
        PSF <- TRUE
        cat("   -- Created a pseudo facet 'psf' (with zero effects)\n")
		cat("   -- because of non-unique person-facet combinations.\n") 		
					}

	# new formula
	formula_update <- paste( c( attr( tA , "term.labels") , nullfacets ) , collapse=" + ")
	formula_update <- as.formula( paste0( "~ " , formula_update ) )
	xsi.setnull <- unique( c( xsi.setnull , nullfacets ) )
	
	if ( length(xsi.setnull)==0 ){
			xsi.setnull <- NULL
						}						

	#********************
	# dimensions for beta fixed
	D <- 1
	if ( ! is.null(B) ){
		D <- dim(B)[[3]]
					}
	if ( ! is.null(Q) ){
		D <- dim(Q)[[2]]
					}
	if ( is.null(beta.fixed) ){				
		beta.fixed <- cbind( 1 , 1:D , 0 )
							}

	res <- list( "formula_update" = formula_update , 
				"xsi.setnull" = xsi.setnull ,
				"beta.fixed" = beta.fixed ,
				"facets"=facets , "PSF" = PSF )
	return(res)
		}
############################################################		