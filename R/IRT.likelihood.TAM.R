

###########################################################
# likelihood
# object of class tam (and tam.mml)
IRT.likelihood.tam <- function( object , ... ){
    ll <- object$like
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.likelihood.tam.mml <- IRT.likelihood.tam 		
IRT.likelihood.tam.mml.3pl <- IRT.likelihood.tam.mml 
###########################################################


###########################################################
# posterior
# object of class tam (and tam.mml)
IRT.posterior.tam <- function( object , ... ){
    ll <- object$hwt
    attr(ll,"theta") <- object$theta
	attr(ll,"prob.theta") <- object$pi.k
	attr(ll,"G") <- 1
    return(ll)
        }
IRT.posterior.tam.mml <- IRT.posterior.tam 		
IRT.posterior.tam.mml.3pl <- IRT.posterior.tam.mml 
###########################################################