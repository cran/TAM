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
 A fit matrix for a specific hypothesis of fit 
 
 EXPLAIN IT IN DETAIL
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

%\references{
%% ~put references to the literature/web site here ~
%}

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
# EXAMPLE 1: dichotomous data sim.rasch
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
# EXAMPLE 2: partial credit model data.gpcm
#############################################################################

data( data.gpcm )
# estimate partial credit model in ConQuest parametrization 'item+iten*step'
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
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Item fit}
\keyword{Infit}
\keyword{Outfit}
\keyword{Fit statistics}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
