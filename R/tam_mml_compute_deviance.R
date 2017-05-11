
tam_mml_compute_deviance <- function( loglike_num , loglike_sto , snodes,
	thetawidth , pweights, deviance=NA, deviance.history=NULL , iter=NULL )
{
	olddeviance <- deviance
    # calculate deviance
    if ( snodes == 0 ){ 
        deviance <- - 2 * sum( pweights * log( loglike_num * thetawidth ) )
    } else {
		deviance <- - 2 * sum( pweights * log( loglike_sto ) )
	}
	#----- deviance change
	rel_deviance_change <- abs( ( deviance - olddeviance ) / deviance  )
    deviance_change <- abs( ( deviance - olddeviance )  )	
	#----- deviance_history
	if (!is.null(deviance.history)){
		deviance.history[iter,2] <- deviance		
	}
	#----- OUTPUT
	res <- list( deviance = deviance, deviance_change=deviance_change ,
				rel_deviance_change=rel_deviance_change, 
				deviance.history=deviance.history)
	return(res)
}	
	