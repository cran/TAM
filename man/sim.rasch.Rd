\name{sim.rasch}

\alias{sim.rasch}
\alias{sim.rasch.pweights}
\alias{sim.rasch.missing}

\docType{data}

\title{
Simulated Rasch data
}

\description{
Simulated Rasch data under unidimensional trait distribution
%% ## DON'T FORGET TO SET SEEDS BEFORE SIMULATING ACTUAL SIM-DATA
%% ## (NOT DONE YET)
%% ## Safe simulation script in "inst"-folder
}

\usage{
 data(sim.rasch) 
 data(sim.rasch.pweights)
 data(sim.rasch.missing)
}

\format{
 The format is: \cr \cr
 \code{num [1:2000, 1:40] 1 0 1 1 1 1 1 1 1 1 ...}\cr
 \code{ - attr(*, "dimnames")=List of 2} \cr
 \code{  ..$ : NULL} \cr
 \code{  ..$ : chr [1:40] "I1" "I2" "I3" "I4" ...} \cr
}

\details{
\code{N <- 2000     } \cr
\code{ # simulate predictors} \cr
\code{ Y <- cbind( rnorm( N , sd = 1.5) , rnorm(N , sd = .3 ) )} \cr
\code{ theta <- rnorm( N ) + .4 * Y[,1] + .2 * Y[,2]  # latent regression model} \cr
\code{ # simulate item responses with missing data} \cr
\code{ I <- 40} \cr
\code{ resp[ theta < 0 , c(1,seq(I/2+1 , I)) ] <- NA} \cr
\code{ # define person weights} \cr
\code{ pweights <- c(  rep(3,N/2) , rep( 1, N/2 ) )} \cr
}

\source{
 Simulated data (see Details)
}

%\references{
 %none
%}

\examples{
\dontrun{      
data(sim.rasch)
N <- 2000
Y <- cbind( rnorm( N , sd = 1.5) , rnorm(N , sd = .3 ) )

# Loading Matrix
# B <- array( 0 , dim = c( I , 2 , 1 )  )
# B[1:(nrow(B)), 2, 1] <- 1
B <- designMatrices(resp = sim.rasch)[["B"]]
  
# estimate model
mod1_1 <- tam(resp=sim.rasch , Y=Y)

# standard errors
res1 <- tam.se(mod1_1)

# Compute fit statistics
tam.fit(mod1_1)

# plausible value imputation
# PV imputation has to be adpated for multidimensional case!
pv1 <- tam.pv( mod1_1 , nplausible = 7 , # 7 plausible values 
               samp.regr = TRUE       # sampling of regression coefficients
              )
              
# item parameter constraints
xsi.fixed <- matrix( c( 1, -2,5, -.22,10, 2 ), nrow=3 , ncol=2 , byrow=TRUE)
xsi.fixed
mod1_4 <- tam( resp=sim.rasch , xsi.fixed=xsi.fixed )

# missing value handling
data(sim.rasch.missing)
mod1_2 <- tam(sim.rasch.missing , Y = Y)

# handling of sample (person) weights
data(sim.rasch.pweights)
N <- 1000
pweights <- c(  rep(3,N/2) , rep( 1, N/2 ) )
mod1_3 <- tam( sim.rasch.pweights , control = list(conv = .001) , pweights = pweights )}
}

\keyword{datasets}
