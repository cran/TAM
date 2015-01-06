
###################################################
# define R function (method)
tamaanify.define.method <- function(res , tam.method ){

	m1 <- mean( paste0( res$items$itemtype ) %in% c("Rasch" , "PCM" ) )
	
	if ( m1 == 1 ){  res$method <- "tam.mml" }
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

