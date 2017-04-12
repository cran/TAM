
tam_mml_3pl_deviance <- function( hwt0 , rfx , res.hwt , pweights , snodes)
{
	rfx <- NULL
    #---- calculate deviance
    if ( snodes == 0 ){ 
		rfx <- rowSums( hwt0 )
		deviance <- - 2 * sum( pweights * log( rfx ) )
    } else {
		deviance <- - 2 * sum( pweights * log( res.hwt$rfx   ) )				
    }
	res <- list( deviance = deviance , rfx = rfx )
	return(res)
}