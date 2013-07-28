###########################################################
# plotting tam expected scores curves
#..........................................................
plot.tam <- function(x, low=-3, high=3, ngroups=6, wle=NULL, export=TRUE, 
                       observed=TRUE, overlay=FALSE , ask=FALSE , ... ) {
  tamobj <- x
  ndim <- tamobj$ndim
  tammodel <- "mml"
  if(is.null(ndim)) {
    ndim <- 1
    tammodel <- "jml"
  }
  if (ndim > 1 ) {
    stop ("Expected scores curves are only available for uni-dimensional models")
  }
  
  nitems <- tamobj$nitems
  nnodes <- 100
  theta <- matrix(seq(low, high, length=nnodes), nrow=nnodes, ncol=ndim)
  iIndex <- 1:nitems
  A <- tamobj$A
  B <- tamobj$B
  if (tammodel == "mml") {
    xsi <- tamobj$xsi$xsi
  }  else {
    xsi <- tamobj$xsi
  }
  maxK <- tamobj$maxK
  resp <- tamobj$resp
  resp.ind <- tamobj$resp.ind
  resp[resp.ind==0] <- NA
  AXsi <- matrix(0,nrow=nitems,ncol=maxK )
  res <- calc_prob.v5(iIndex=1:nitems , A=A , AXsi=AXsi , B=B , xsi=xsi , theta=theta , 
                      nnodes=nnodes , maxK=maxK , recalc=TRUE )  
  rprobs <- res[["rprobs"]]
  AXsi <- res[["AXsi"]]
  cat <- 1:maxK - 1
  expScore <- sapply(1:nitems, function(i) colSums(cat*rprobs[i,,], na.rm=TRUE))
  
  if (is.null(wle)) {
    if (tammodel == "mml") {
      wleobj <- tam.wle(tamobj)
      wle <- wleobj$theta
      } 
    else {
      wle <- tamobj$WLE    # model is jml
    }
  }
  wleSorted <- sort(wle)
  ncases <- length(wleSorted)
  groupnumber <- round(seq(1:ncases) / (ncases/ngroups) + 0.5)
  
  aggr <- aggregate(wleSorted, list(groupnumber), mean)
  theta2 <- aggr$x
  
  d <- data.frame(wle, resp)
  d1 <- d[order(wle),]
  d2 <- d1[-1]
  obScore <- apply(d2,2, function(x) aggregate(x, list(groupnumber), mean, na.rm=TRUE))
  
  for (i in 1:nitems) {
    if (i==1 || !overlay) {
	  ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )
      plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability",
#         main=paste("Expected Scores Curve - Item ", i)
		 main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] )	 ,
			ylim=ylim2 , ...
				)
    } else {
      lines(theta, expScore[,i],type="l", col=i, lwd=3, pch=1) 
    }
    if (observed) {
      lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
      }
	  par(ask=ask)	
  }
  #*****
  # export item plots
  dir.create( "Plots" )
  if (export) {
    for (i in 1:nitems) {
	  itemlab <- colnames(tamobj$resp)[i]
	  ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )	  
      png(filename=paste("Plots\\Item_", itemlab ,".png", sep=""))  
      plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability", 
          main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] ) ,
				ylim=ylim2 , ... )
      if (observed ) {
        lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
        }
#      Sys.sleep(0.01)
    }
    for (dv in dev.list()) if(!(is.null(dv)) && (dv!=2)) dev.off(dv)
  }
}

plot.tam.mml <- plot.tam
plot.tam.jml <- plot.tam