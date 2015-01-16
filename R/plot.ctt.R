
#########################################################################
# plot function for empirical item characteristic curves
plotctt <- function( resp , theta , Ncuts = NULL , ask =FALSE , ... ){
	if ( ! is.null(Ncuts) ){
		stepw <- 1/Ncuts
		cuts <- quantile( theta , seq( stepw , 1-stepw , length=Ncuts-1 ) , na.rm=TRUE)
		eps <- .001
		cuts <- round( c( min(theta, na.rm=TRUE) - eps  , cuts , max( theta , na.rm=TRUE ) + eps ) , 2 )	
		resp <- resp[ ! is.na( theta) , , drop=FALSE ]
		theta <- theta[ !is.na(theta) ]
		theta.cuts <- cut( theta , cuts )
				} else { theta.cuts <- theta }
	I <- ncol(resp)
	par( mfrow=c(1,1))
	for (ii in 1:I){
		# ii <- 25		
		y <- resp[,ii]
		unique.y <- sort( unique(y) )
		unique.y <- unique.y[ unique.y != "NA" ]
		L1 <- length(unique.y)
		item <- colnames(resp)[ii]		
		dfr <- NULL
		for (ll in 1:L1){
			# ll <- 1
			dfr.ll <- aggregate( 1 * ( y == unique.y[ll] ) , list( theta.cuts ) , mean , na.rm=TRUE )
			dfr.ll <- data.frame("cat" = unique.y[ll] , dfr.ll )
			colnames(dfr.ll)[-1] <- c("theta.cut" , "prob")
			dfr <- rbind( dfr , dfr.ll )
						}				
		main <- paste0('Trace lines for item ', item)
		vkey <- paste0("Cat " , unique.y)
		print( 
			lattice::xyplot( prob ~ theta.cut , data=dfr , group=cat , ylim=c(-.1 , 1.1) , type="o" ,  
				auto.key= TRUE  ,       
				par.settings = list(superpose.symbol = list(pch = 1:L1))	,			
				ylab = expression(P(theta)), xlab = expression(theta) , main=main , lty=1:L1 , pch=1:L1 ,
				...
							)  
					)
		par( ask=ask )
# stop()		
			}
		}
###################################################################################