###########################################################
# plotting tam expected scores curves
#..........................................................
plot.tam <- function(x, items=1:x$nitems, type="expected" ,
					low=-3, high=3, ngroups=6, 
                     wle=NULL, export=TRUE, export.type="png", 
                     export.args=list(), observed=TRUE, overlay=FALSE , 
                     ask=FALSE, ...) {
  
  old.opt.dev <- getOption("device")
  old.opt.err <- c(getOption("show.error.messages"))
  old.par.ask <- par("ask")
  
  on.exit(options("device"=old.opt.dev))
  on.exit(options("show.error.messages"=old.opt.err), add=TRUE)
  on.exit(par("ask"=old.par.ask), add=TRUE)
  
  tamobj <- x
  ndim <- tamobj$ndim
  tammodel <- "mml"
  if(is.null(ndim)) {
    ndim <- 1
    tammodel <- "jml"
  }
  if (ndim > 1 ) {
   if ( type=="expected"){
    stop ("Expected scores curves are only available for uni-dimensional models")
						}
                  }
  
  nitems <- tamobj$nitems

  nnodes <- 100
  if (ndim == 1 ){
	theta <- matrix(seq(low, high, length=nnodes), nrow=nnodes, ncol=ndim)
				} else {
#	theta <- tamobj$theta
	nnodes <- 40
    nodes <- seq(low, high, length=nnodes)
	theta <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , ncol = ndim ) ) ) )	
	nnodes <- nrow(theta)	
	B <- tamobj$B
					}

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
  
  if ( type == "expected" ){
	  
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
	  wleSorted <- sort(wle, na.last=FALSE)
	  ncases <- length(wleSorted)
	  groupnumber <- round(seq(1:ncases) / (ncases/ngroups) + 0.5)
	  
	  aggr <- aggregate(wleSorted, list(groupnumber), mean)
	  theta2 <- aggr$x
	  
	  d <- data.frame(wle, resp)
	  d1 <- d[order(wle),]
	  d2 <- d1[-1]
	  obScore <- apply(d2,2, function(x) aggregate(x, list(groupnumber), mean, na.rm=TRUE))
			}
			
  
  #*************************************************
  # begin plot function
  
  for (i in (1:nitems)[items]) {
	#***********************************************************
    #** expected item response curves	
    if ( type=="expected"){
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
	}
	#***********************************************************
    if ( type=="items"){	

		rprobs.ii <- rprobs[i,,]
		rprobs.ii <- rprobs.ii[ rowMeans( is.na(rprobs.ii) ) < 1 , ]		
		K <- nrow(rprobs.ii)		
		if ( ndim == 1 ){ theta0 <- theta }
		dat2 <- NULL
		#************
		if ( ndim > 1 ){
				B.ii <- B[i,,]	
				ind.ii <- which( colSums( B.ii ) > 0 )[1]
				rprobs0.ii <- rprobs.ii
				rprobs0.ii <- aggregate( t(rprobs0.ii) , list( theta[,ind.ii] ) , mean )
				theta0 <- rprobs0.ii[,1,drop=FALSE]
				rprobs.ii <- t( rprobs0.ii[,-1] )						
						}
		#**************
		for (kk in 1:K){
			# kk <- 1
			dat2a <- data.frame( "Theta" = theta0[,1] , "cat" = kk , "P" = rprobs.ii[kk,] )
			dat2 <- rbind(dat2 , dat2a)
						}
		main <- paste("Item", colnames(x$resp)[i] )				
		auto.key <- NULL				
		simple.key <- paste0("Cat" , 1:K -  1)
		auto.key <- simple.key
		dat2$time <- dat2$cat
		dat2$time <- paste0("Cat" , dat2$time )
		
		simple.key <- FALSE
		Kpercol <- K
		# floor(K/Kpercol)+1
		auto.key <- list( # columns = Kpercol , 
						    lines=TRUE , points=FALSE , rows=2)
        h1 <- xyplot(P ~ Theta, dat2, group = time, type = 'l', auto.key = auto.key,
                                main = main, ylim = c(-0.1,1.1), simple.key = simple.key , 
								xlim=c(low,high) , 
								ylab = expression(P(theta)), xlab = expression(theta), ... ) 
		plot(h1)		# plot
							
			}
	#***************		
    par(ask=ask)	
  }             # end item ii
  #*************************************************
  
  #*****
  # export item plots
  if (export) {
    
    if(!file.exists("Plots")) dir.create( "Plots" )
    export.type.dev <- switch(export.type,
                              "ps"="postscript",
                              "emf"=if (.Platform$OS.type == "windows") "win.metafile" else "x11",
                              "wmf"=if (.Platform$OS.type == "windows") "win.metafile" else "x11",
                              export.type)
    export.type.ff <- switch(export.type,
                             "postscript"="ps",
                             "win.metafile"="wmf",
                             "x11"="wmf",
                             export.type)
    
    
    options(show.error.messages = FALSE)
    options("device"=export.type.dev)
    
    for (i in (1:nitems)[items]) {
      
      
      itemlab <- colnames(tamobj$resp)[i]
      dev.err <- try({ 
        do.call("dev.new", 
                args=list("filename"=file.path("Plots", paste("Item_", itemlab, ".", export.type.ff, sep="")), 
                          export.args))
      })
      
      if(!is.null(dev.err)){
        warning( dev.err[1], "  --> No file created."  )
      }else{
        
		#***************************************************
		# expected response functions
		if (type=="expected"){
			ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )        
			plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability", 
				 main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] ) ,
				 ylim=ylim2 , ... )        
			if (observed ) {
			  lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
						}
							}
		if ( type=="items" ){
		    
			rprobs.ii <- rprobs[i,,]
			rprobs.ii <- rprobs.ii[ rowMeans( is.na(rprobs.ii) ) < 1 , ]		
			K <- nrow(rprobs.ii)		
			if ( ndim == 1 ){ theta0 <- theta }
			dat2 <- NULL
			#************
			if ( ndim > 1 ){
					B.ii <- B[i,,]	
					ind.ii <- which( colSums( B.ii ) > 0 )[1]
					rprobs0.ii <- rprobs.ii
					rprobs0.ii <- aggregate( t(rprobs0.ii) , list( theta[,ind.ii] ) , mean )
					theta0 <- rprobs0.ii[,1,drop=FALSE]
					rprobs.ii <- t( rprobs0.ii[,-1] )						
							}
			#**************
			for (kk in 1:K){
				# kk <- 1
				dat2a <- data.frame( "Theta" = theta0[,1] , "cat" = kk , "P" = rprobs.ii[kk,] )
				dat2 <- rbind(dat2 , dat2a)
							}
			main <- paste("Item", colnames(x$resp)[i] )				
			auto.key <- NULL				
			simple.key <- paste0("Cat" , 1:K -  1)
			auto.key <- simple.key
			dat2$time <- dat2$cat
			dat2$time <- paste0("Cat" , dat2$time )
			
			simple.key <- FALSE
			Kpercol <- K
			# floor(K/Kpercol)+1
			auto.key <- list( # columns = Kpercol , 
								lines=TRUE , points=FALSE , rows=2)
			h1 <- xyplot(P ~ Theta, dat2, group = time, type = 'l', auto.key = auto.key,
									main = main, ylim = c(-0.1,1.1), simple.key = simple.key , 
									xlim=c(low,high) , 
									ylab = expression(P(theta)), xlab = expression(theta), ... ) 
			plot(h1)		# plot
								
							
							}
							
        dev.off(dev.cur())
        
      }
      
      #       options("device"=old.opt.dev)
      #       options(show.error.messages = as.character(old.opt.err))        
    }
    
    #*****
    # Print path
    if(is.null(dev.err)){ cat("....................................................\n",
                              "Plots exported in", export.type, "format into folder:\n", 
                              file.path(getwd(), "Plots")) ; flush.console() }
  }
  
}

plot.tam.mml <- plot.tam
plot.tam.jml <- plot.tam
