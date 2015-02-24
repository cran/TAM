


###################################################
# summary for tam.fit
summary.msq.itemfit <- function( object , ... ){
    object <- object$itemfit
	ind <- grep( "fitgroup" , colnames(object) )
	obji <- object
	for ( vv in seq(ind+1,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}
	return(obji)
		}
###################################################