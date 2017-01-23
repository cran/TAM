
tam_mml_3pl_difference_quotient <-
	function( d0 , d0p , d0m , h)
{
	d1 <- ( d0p - d0 ) / h
	d2 <- ( ( d0p - d0 ) - ( d0 - d0m ) ) / h^2	
	res <- base::list( d1 = d1 , d2 = d2 )
	base::return(res)
}	
	