
#####################################################
# computes reliability for one-dimensional WLEs
WLErel <- function( theta , error , w = rep(1,length(theta) )){
    v1 <- weighted_var( x = theta , w = w  )
    v2 <- weighted_mean( error^2 , w=w  )
    # WLE_Rel = ( v1 - v2 ) / v1 = 1 - v2 / v1
    rel <- 1 - v2 / v1
	return(rel)
}
#######################################################	
EAPrel <- function( theta , error , w = rep(1,length(theta) )){
    v1 <- weighted_var( x = theta , w = w  )
    v2 <- weighted_mean( error^2 , w=w )	
    # v1 / (v1+v2) = 1 - v2 / ( v1 + v2 )
    rel <- v1 / ( v1 + v2 )
	return(rel)
}
#######################################################	
				
