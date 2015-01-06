


###########################################################
# object of class tam (tam.mml)
IRT.irfprob.tam <- function( object , ... ){
    ll <- object$rprobs
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.irfprob.tam.mml <- IRT.irfprob.tam
# IRT.irfprob.tam.mfr <- IRT.irfprob.tam		
###########################################################

###########################################################
# object of class tam (tam.mml)
IRT.irfprob.tam.mml.3pl <- function( object , ... ){
    ll <- object$rprobs
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	res <- list( "delta" = object$delta , 
	             "delta.designmatrix" = object$delta.designmatrix )
	attr(ll,"skillspace") <- res	
	attr(ll,"G") <- 1
    return(ll)
        }
###########################################################

###########################################################
# objects of class tamaan
IRT.irfprob.tamaan <- function( object , ... ){
	if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
			res0 <- IRT.irfprob.tam( object , ... )			
			}
	if (object$tamaanify$method == "tam.mml.3pl"){
			res0 <- IRT.irfprob.tam.mml.3pl( object , ... )			
							}
	return(res0)	
			}
###################################################################			