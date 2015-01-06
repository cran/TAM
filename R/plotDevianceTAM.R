###############################################################################
plotDevianceTAM   <- function ( tam.obj , omitUntil = 1, reverse = TRUE ) {
        stopifnot(class(tam.obj) %in% c("tam.mml","tam.mml.2pl","tam.mml.mfr","tam.mml.3pl","tamaan") )
        if(omitUntil>0)  {
				devChange <- diff(tam.obj$deviance.history[-c(1:omitUntil),2])
						} else { 
					devChange <- diff(tam.obj$deviance.history[,2]) 
							}
        if(reverse)      {devChange <- -1 *  devChange }
        devChange <- data.frame ( nr = 1:length(devChange), devChange)
		xm        <- ceiling( max(devChange[,1])/10 )*10      
					### maximum auf x Achse
		xt        <- NULL;	for ( i in c( 1:30 ) ) xt <- c ( xt , (xm/10) %% i == 0 )
		xt        <- max ( which ( xt ) )                         
					### ticks auf x achse so setzen dass schön 10er
		cex       <- 0.85 - ( length(devChange[,1]) / 1000 )       
					### Punktgröße setzen (cex). initial: .85
		if ( cex < 0.40 ) cex <- 0.40                             
					### mit jeden 100 Iteration 0.1 runter, aber mind. 0.4
		plot ( devChange[,1] , devChange[,2] , type = "o" , 
				main = "Deviance Change Plot", xlab = "Iteration" , 
				xlim = c(min(devChange[,1]) ,max(devChange[,1])) ,  xaxp = c(0,xm,xt) , 
				ylab = "Deviance Change" , pch = 20 , cex = cex , lwd = 0.75 )
		abline ( a=0 , b=0 )                                       ### Linie bei 0
		dcr       <- devChange[devChange[,2]<0,]                   ### Punkte unter 0 rot
        points( dcr[,1] , dcr[,2] , pch=20, cex = cex , col="red") 
			}
###############################################################################