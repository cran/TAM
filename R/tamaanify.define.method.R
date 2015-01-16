
###################################################
# define R function (method)
tamaanify.define.method <- function(res , tam.method ){

	itemtypes <- paste0( res$items$itemtype )
	l1 <- strsplit( itemtypes , split="," , fixed=TRUE )
	itemtypes <- unlist( lapply( l1 , FUN = function(ll){ ll[1] } ) )

	m1 <- mean(  itemtypes %in% c("Rasch" , "PCM" ) )
	
	if ( m1 == 1 ){  
			res$method <- "tam.mml" 
					}
	if ( ! is.null(tam.method) ){
		res$method <- tam.method 
					}				
	al <- res$ANALYSIS.list$type
	if ( al %in% c("LCA", "LOCLCA","MIXTURE" , "OLCA") ){
		res$method <- "tam.mml.3pl"
						}																	
	return(res)
			}
###################################################			

