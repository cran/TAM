## File Name: tampv2datalist.R
## File Version: 9.146


##***  converts a pv object into a list of datasets
tampv2datalist <- function( tam.pv.object, pvnames=NULL, Y=NULL,
            Y.pid="pid", as_mids=FALSE, stringsAsFactors=FALSE )
{
    pv <- tam.pv.object$pv
    ndim <- tam.pv.object$ndim
    nplausible <- tam.pv.object$nplausible
    Y00 <- data.frame( pid=tam.pv.object$pid, pweights=tam.pv.object$pweights,
                stringsAsFactors=stringsAsFactors )

    if ( ! is.null(Y) ){
        Y <- as.data.frame(Y)
        if(  sum( colnames(Y) %in% Y.pid )==0 ){
            Y[, Y.pid] <- pv$pid

        }
    }
    if ( is.null(pvnames) ){
        pvnames <- paste0("PV.Dim", 1:ndim )
    }
    # create list of multiply imputed datasets
    datalist <- list( 1:nplausible )
    for (ii in 1:nplausible){
        pv1 <- pv[, 1 + 1:ndim + (ii-1)*ndim, drop=FALSE]
        colnames(pv1) <- pvnames
        Y0 <- data.frame( Y00, pv1 )
        if ( ! is.null(Y) ){
            Y0 <- merge( x=Y, y=Y0, by.x=Y.pid, by.y="pid", all=TRUE )
        }
        dat <- Y0
        datalist[[ii]] <- dat
    }
    if (as_mids){
        require_namespace_msg(pkg="miceadds")
        datalist <- miceadds::datalist2mids(datalist)
    }

    return(datalist)
}

