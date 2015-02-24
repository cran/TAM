\name{msq.itemfit}

\alias{msq.itemfit}
\alias{summary.msq.itemfit}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Squared Residual Based Item Fit Statistics (Infit, Outfit)
}


\description{
This function computes computed the outfit and infit statistic for items
or item groups. Contrary to \code{\link{tam.fit}}, this function is 
not based on simulation from individual posterior distributions
but rather on evaluating the individual posterior. 
}


\usage{
msq.itemfit( object , fitindices=NULL )

\method{summary}{msq.itemfit}(object, \dots)
}



\arguments{
  \item{object}{
Object for which the classes \code{IRT.data}, \code{IRT.posterior}
and \code{predict} are defined.
}
\item{fitindices}{Vector with integers defining the item groups for
which the fit should be evaluated.}
\item{\dots}{Further arguments to be passed}
}


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Data frame with outfit and infit statistics
}


%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%% Alexander Robitzsch
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{	
\dontrun{
#############################################################################
# SIMULATED EXAMPLE 1: Simulated data Rasch model
############################################################################# 	

#*** simulate data
library(sirt)
set.seed(9875)
N <- 2000
I <- 20
b <- sample( seq( -2 , 2 , length=I ) )
a <- rep( 1, I )
# create some misfitting items
a[c(1,3)] <- c(.5 , 1.5 )
# simulate data
dat <- sirt::sim.raschtype( rnorm(N) , b=b , fixed.a=a )
#*** estimate Rasch model
mod1 <- tam.mml(resp=dat) 

#--- item fit from "msq.itemfit" function
fit1 <- msq.itemfit(mod1)
summary( fit1 )

#--- item fit using simulation in "tam.fit"
fit0 <- tam.fit( mod1 )
summary(fit0)

# define some item groups for fit assessment
fitindices <- rep( c(1,2) , each=10)
fit2 <- msq.itemfit( mod1 , fitindices )
summary(fit2)
}	
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Item fit}
% \keyword{Individual posterior}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
