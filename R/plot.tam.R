###########################################################
# plotting tam expected scores curves
#..........................................................
plot.tam <- function(x, items=1:x$nitems, low=-3, high=3, ngroups=6, 
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
  wleSorted <- sort(wle, na.last=FALSE)
  ncases <- length(wleSorted)
  groupnumber <- round(seq(1:ncases) / (ncases/ngroups) + 0.5)
  
  aggr <- aggregate(wleSorted, list(groupnumber), mean)
  theta2 <- aggr$x
  
  d <- data.frame(wle, resp)
  d1 <- d[order(wle),]
  d2 <- d1[-1]
  obScore <- apply(d2,2, function(x) aggregate(x, list(groupnumber), mean, na.rm=TRUE))
  
  for (i in (1:nitems)[items]) {
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
        
        ylim2 <- c(0,max( tamobj$resp[,i] , na.rm=TRUE ) )
        
        plot(theta, expScore[,i], ,col=12, type="l", lwd=3, las=1, ylab="Score", xlab="Ability", 
             main=paste("Expected Scores Curve - Item ", colnames(tamobj$resp)[i] ) ,
             ylim=ylim2 , ... )
        
        if (observed ) {
          lines(theta2,obScore[[i]]$x, type="o", lwd=2, pch=1)
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
