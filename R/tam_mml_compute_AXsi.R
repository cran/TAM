
tam_mml_compute_AXsi <- function( A , xsi )
{
	dimA <- base::dim(A)	
	AXsi <- matrix( NA , nrow=dimA[1] , ncol=dimA[2] )
	for (kk in 1:dimA[2]){
		AXsi[,kk] <- A[,kk,] %*% xsi	
	}
	base::return(AXsi)
}