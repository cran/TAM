## File Name: IRTLikelihood.cfa.R
## File Version: 9.201

#---- IRTLikelihood for fitted CFA model
IRTLikelihood.cfa <- function( data, cfaobj=NULL,
        theta=NULL, L=NULL, nu=NULL, psi=NULL,
        snodes=NULL, snodes.adj=2, version=1)
{
    if ( ! is.null(cfaobj) ){
        cfaobj <- cfa.extract.itempars(object=cfaobj)
        L <- cfaobj$L
        nu <- cfaobj$nu
        psi <- cfaobj$psi
        obs.vars <- cfaobj$obs.vars
        data <- data[, obs.vars, drop=FALSE ]
    }
    D <- ncol(L)
    #*****
    # create theta if it is not provided
    if ( is.null(theta) ){
        if ( is.null(snodes) ){
            if ( D==1){ snodes <- 100 }
            if ( D==2){ snodes <- 200 }
            if ( D>2){ snodes <- 1000 }
        }
        theta0 <- snodes.adj * seq(-3,3,len=21)
        if (D>2){
            r1 <- tam_import_sfsmisc_QUnif(n=snodes, min=0, max=1, n.min=1, p=D, leap=409)
            theta <- stats::qnorm( r1 )
            for ( dd in 1:D){
                theta[,dd] <- snodes.adj*theta[,dd]
            }
        }
        if (D==1){
            theta <- matrix( theta0, ncol=1 )
        }
        if (D==2){
            theta <- expand.grid( theta0, theta0 )
        }
    }

    #**** R version
    if (version==0){
        TP <- nrow(theta)
        N <- nrow(data)
        colnames(theta) <- colnames(L)
        hwt <- matrix( 1, nrow=N, ncol=TP )
        I <- ncol(data)

        for (ii in 1:I){
            #ii <- 1
            term <- matrix( nu[ii], nrow=N, ncol=TP)
            for (dd in 1:D){
                term <- term + matrix( L[ii,dd] * theta[,dd], nrow=N, ncol=TP, byrow=TRUE )
            }
            h1 <- stats::dnorm( data[,ii], mean=term, sd=sqrt( psi[ii,ii] ) )
            ind <- which( ! is.na( data[,ii] ) )
            hwt[ind,] <- hwt[ind,] * h1[ind,]
        }
    }
    #***** Rcpp version
    if (version==1){
        data <- as.matrix(data)
        nu <- as.vector(nu)
        psi <- as.matrix(psi)
        L <- as.matrix(L)
        theta <- as.matrix(theta)
        hwt <- tam_rcpp_irt_likelihood_cfa( data=data, nu=nu, psi=psi, L=L, theta=theta )
        hwt <- hwt$hwt
    }
    res <- hwt
    attr(res,"theta") <- theta
    attr(res,"prob.theta") <- NA
    attr(res,"G") <- 1
    return(res)
}
