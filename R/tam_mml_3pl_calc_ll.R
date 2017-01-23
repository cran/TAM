

tam_mml_3pl_calc_ll <- function( n.ik , probs, eps )
{
	maxK <- base::dim(n.ik)[2]
	probs <- probs + eps
	l1 <- 0
	for (kk in 1:maxK){
		l1 <- l1 + base::sum( n.ik[,kk,,drop=FALSE] * base::log( probs[,kk,,drop=FALSE] ) )
	}
	base::return(l1)
}