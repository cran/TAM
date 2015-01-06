

#######################################################################
# create a lavaan syntax from lavaan parameter table
# This syntax only works for single groups.
lavpartable2lavsyntax <- function( lavpartable ){
    LL <- nrow(lavpartable )    
    syn0 <- paste0( lavpartable$lhs , lavpartable$op )
    lavpartable$prefix <- ""
    lavpartable$prefix <- ifelse( paste(lavpartable$label ) != "" , 
                paste0(lavpartable$label , "*"  ) , lavpartable$prefix ) 
    lavpartable$prefix <- ifelse( ( lavpartable$free == 0 ) & ( paste(lavpartable$label ) == "" ) , 
                paste0(lavpartable$ustart , "*"  ) , lavpartable$prefix )              
    syn0 <- paste0( syn0 , lavpartable$prefix , lavpartable$rhs  )
    syn0 <- paste0( syn0 , collapse="\n")
    return(syn0)
            }			
#######################################################################
	