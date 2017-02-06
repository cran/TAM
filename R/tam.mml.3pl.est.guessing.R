####################################################
# estimation guessing parameters
.mml.3pl.est.guessing <- function( guess , Msteps , convM , 
	  nitems , A , AXsi , B, xsi , theta , nnodes , maxK ,
	  n.ik , N.ik , est.guess , old.increment.guess , guess.prior ,
	  progress , max.guess 	  )
{
	  old_increment <- old.increment.guess	  
	  if (progress){ 
			cat("\nM Step Guessing     |")
			utils::flush.console() 
	  }
	  eps <- 1e-10
	  eps10 <- 1E-30
#	  eps1 <- 1e-5
	  guess.logit0 <- guess.logit <- stats::qlogis( guess + eps)	  
	  ind.guess <- which( est.guess != 0 )	  
	  n1ij <- n.ik[,2,]
	  nij <- N.ik
	  Miter <- 1
	  converge <- FALSE
	  guess_old <- guess	 
	  h <- 1E-4	
	  # h <- 1E-3
	  
	  while ( !converge & ( Miter <= Msteps ) ){	

	      res.p <- .mml.3pl.calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK , 
								 guess=guess )					
          rprobs00 <- rprobs <- res.p[["rprobs"]]
          rprobs0 <- res.p$rprobs0
		  pij <- rprobs[ , 2 , ]
		  pij <- pij + eps
	      pij0 <- rprobs0[ , 2 , ]		 
		  guess <- ifelse( ( guess < h ) & ( est.guess != 0 ) , 4*h , guess )
		  # guess.logit <- stats::qlogis( guess + eps)	  
		  #****
		  # derivatives guessing priors
		  if ( ! is.null(guess.prior) ){
#			  d0  <- log( stats::dbeta( stats::plogis(guess.logit) , guess.prior[,1] , guess.prior[,2] ) + eps)
#			  d0p <- log( stats::dbeta( stats::plogis(guess.logit+h), guess.prior[,1] , guess.prior[,2] ) + eps)
#			  d0m <- log( stats::dbeta( stats::plogis(guess.logit-h), guess.prior[,1] , guess.prior[,2] ) + eps)
			  d0  <- log( stats::dbeta( guess , guess.prior[,1] , guess.prior[,2] ) + eps)
			  d0p <- log( stats::dbeta( guess + h , guess.prior[,1] , guess.prior[,2] ) + eps)
			  d0m <- log( stats::dbeta( guess - h, guess.prior[,1] , guess.prior[,2] ) + eps)
			  d1 <- ( d0p - d0 ) / h
			  d2 <- ( ( d0p - d0 ) - ( d0 - d0m ) ) / h^2		  		  		  
		  }

		  
	      rprobs1 <- .mml.3pl.calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK , 
								 guess=guess + h )$rprobs		
	      rprobs2 <- .mml.3pl.calc_prob.v5( iIndex=1:nitems , A=A , AXsi=AXsi , B=B , 
                                 xsi=xsi , theta=theta , nnodes=nnodes, maxK=maxK , 
								 guess=guess - h )$rprobs		
		   #*** calculate log-likelihood
			n0ij <- nij - n1ij
			eps2 <- 1E-5
			l0 <- tam_mml_3pl_calc_ll_est_guessing( n0ij , n1ij , probs=rprobs00, eps=eps2 )			
			l1p <- tam_mml_3pl_calc_ll_est_guessing( n0ij , n1ij , probs=rprobs1, eps=eps2 )						
			l1m <- tam_mml_3pl_calc_ll_est_guessing( n0ij , n1ij , probs=rprobs2, eps=eps2 )						
			res <- tam_mml_3pl_difference_quotient( d0=l0 , d0p=l1p , d0m=l1m , h=h)
			der1 <- res$d1
			der2 <- res$d2
	  
		  #***
		  # first derivative
#		  der1 <- rowSums( ( n1ij - pij * nij ) * guess / pij )
		  # second derivative
		  #der2 <- rowSums( guess*(1-guess) * ( ( guess + pij0 ) / pij^2 * n1ij - nij ) )
#		  der2 <- rowSums( guess*(1-guess) * ( ( guess + pij0 ) / ( pij^2 + eps ) * n1ij - nij ) )
		  if ( ! is.null( guess.prior) ){	  
				der1 <- der1 + d1
				der2 <- der2 + d2 
		}
									
		  # aggregation over group of parameters
		  der1 <- stats::aggregate( der1 , list( est.guess ) , sum )
		  der2 <- stats::aggregate( der2 , list( est.guess ) , sum )	  	  		  
		  der1 <- der1[ der1[,1] != 0  , , drop=FALSE]
		  der2 <- der2[ der2[,1] != 0  , , drop=FALSE]
		  increment0 <- increment <- der1[,2] / ( abs(der2[,2]) + eps )
		  ci <- ceiling( abs(increment) / ( abs( old_increment) + eps ) )
		  increment <- ifelse( abs( increment) > abs(old_increment)  , 
								 increment/(2*ci) , increment )	  

		  increment <- increment[ est.guess ]
	      guess[ ind.guess ] <- guess[ind.guess ] + increment
		  guess <- ifelse( ( guess < h ) & ( est.guess != 0 ) , 4*h , guess )
		  guess <- ifelse( guess > max.guess , max.guess , guess )
		  old_increment <- max(abs(increment))
		  if ( old_increment < convM){ 
			converge <- TRUE 
			}
		  Miter <- Miter + 1
		  if (progress){ cat("-") ; utils::flush.console() }			  
					}
	  #*********************************************************
	  # standard error of logit guessing parameter				
	  # se.guess <- sqrt( 1 / abs(der2[ est.guess , 2] ) )
	  
	  guess.change <- max( abs( guess - guess_old ))
	  se.guess <-  sqrt( 1 / ( abs( der2[ est.guess ,2] ) + eps10 ) )
	  se2 <- 0*guess
	  se2[ ind.guess ] <- se.guess[ est.guess ]
	  # transform standard errors according to delta formula
	  # h = plogis = ( 1 + exp( -x ) )^(-1)
	  # h' = -1 * exp(-x) * ( 1 + exp( -x ) )^(-2)	 = h * ( 1 - h )
	 hast <- guess * ( 1 - guess )
	 se2 <- sqrt( hast^2 ) * se2
	  res <- list( "guess" = guess , "guess.change" = guess.change ,
			se.guess = se2 )
	  return(res)
	  }
################################################################################	  
