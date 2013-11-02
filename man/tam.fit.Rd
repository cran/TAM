\name{tam.fit}

\alias{tam.fit}
\alias{tam.mml.fit}
%% \alias{tam.jml.fit}

\title{
Item Infit and Outfit Statistic
}

\description{
The item infit and outfit statistic are calculated for
objects of classes \code{tam}, \code{tam.mml} and
\code{tam.jml},  respectively.
}

\usage{
tam.fit(tamobj, ...)

tam.mml.fit(tamobj, FitMatrix=NULL, progress=FALSE)
}

%% tam.jml.fit(tamobj, resp, resp.ind, A, B, nstud, nitems, maxK, ItemScore, 
%%    theta, xsi, Msteps, pweightsM, est.xsi.index)


\arguments{
  \item{tamobj}{
  An object of class \code{tam}, \code{tam.mml} or \code{tam.jml}
}
 \item{FitMatrix}{
 A fit matrix \eqn{F} for a specific hypothesis of fit of the linear function
\eqn{F \xi } (see Simulated Example 3 and Adams & Wu 2007).
  }
 \item{progress}{
    An optional logical indicating whether computation progress should
    be displayed at console.
}
\item{\dots}{Further arguments to be passed}
}

\details{
Item fit is automatically calculated in JML estimation using
\code{\link{tam.jml}}.
}

\value{
In case of \code{tam.mml.fit} a data frame with four columns:
\item{Outfit}{Item outfit statistic}
\item{Outfit_t}{The \eqn{t} value for the outfit statistic}
\item{Infit}{Item infit statistic}
\item{Infit_t}{The \eqn{t} value for the infit statistic}
}

\references{
Adams, R. J., & Wu, M. L. (2007). The mixed-coefficients multinomial logit model.
A generalized form of the Rasch model. In M. von Davier & C. H. Carstensen (Eds.):
\emph{Multivariate and mixture distribution Rasch models: Extensions and applications}
(pp. 55-76). New York: Springer.
}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#############################################################################
# EXAMPLE 1: Dichotomous data sim.rasch
#############################################################################

data(sim.rasch)
# estimate Rasch model
mod1 <- tam.mml(resp=sim.rasch) 
# item fit
fit1 <- tam.fit( mod1 )
##   > print( round( fit1 ,3 ))
##       Outfit Outfit_t Infit Infit_t
##   I1   0.951   -1.549 0.994  -0.135
##   I2   1.033    1.040 1.028   0.755
##   I3   1.012    0.387 1.009   0.273
##   I4   1.041    1.296 1.057   1.614

\dontrun{
#############################################################################
# EXAMPLE 2: Partial credit model data.gpcm
#############################################################################

data( data.gpcm )
# estimate partial credit model in ConQuest parametrization 'item+item*step'
mod2 <- tam.mml( resp=data.gpcm , irtmodel="PCM2" )
summary(mod2)
# estimate item fit
fit2 <- tam.fit(mod2)
## > round( fit2 ,3 )
##                 Outfit Outfit_t Infit Infit_t
##   Comfort        0.968   -0.426 0.967  -0.375
##   Work           1.061    0.867 1.054   0.695
##   Benefit        0.984   -0.196 0.987  -0.165
##   Comfort_step1  0.932   -0.946 0.955  -0.728
##   [...]
##  => The first three rows of the data frame correspond to the fit statistics
#      of first three items Comfort, Work and Benefit.

#############################################################################
# SIMULATED EXAMPLE 3: Fit statistic testing for local independence
#############################################################################

#generate data with local DEPENDENCE and User-defined fit statistics
set.seed(4888)
I <- 40		# 40 items
N <- 1000       # 1000 persons

delta <- seq(-2,2, len=I)
theta <- rnorm(N, 0, 1)
# simulate data
prob <- plogis(outer(theta, delta, "-"))
rand <- matrix(runif(N*I), nrow=N, ncol=I)
resp <- 1*(rand < prob)
colnames(resp) <- paste("I" , 1:I, sep="")
  
#Introduce some local dependence
for (item in c(10, 20, 30)){ 
 #  20% of the item responses for items 11, 21 and 31
 #are made equal to the previous item 
  row <- round(runif(0.2*N)*N + 0.5)     
  resp[row, item+1] <- resp[row, item]
}
  
#run TAM
mod1 <- tam(resp)
  
#User-defined fit design matrix
F <- array(0, dim=c(dim(mod1$A)[1], dim(mod1$A)[2], 6))
F[,,1] <- mod1$A[,,10] + mod1$A[,,11]
F[,,2] <- mod1$A[,,12] + mod1$A[,,13]
F[,,3] <- mod1$A[,,20] + mod1$A[,,21]
F[,,4] <- mod1$A[,,22] + mod1$A[,,23]
F[,,5] <- mod1$A[,,30] + mod1$A[,,31]
F[,,6] <- mod1$A[,,32] + mod1$A[,,33]
fit <- tam.fit(mod1, FitMatrix=F)
fit
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Item fit}
\keyword{Infit}
\keyword{Outfit}
\keyword{Fit statistics}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
