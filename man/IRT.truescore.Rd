\name{IRT.truescore}
\alias{IRT.truescore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a \eqn{\theta} Score into a True Score \eqn{\tau ( \theta)}
}
\description{
Converts a \eqn{\theta} score into a true score \eqn{\tau ( \theta)}
}
\usage{
IRT.truescore(object, iIndex = NULL, theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object for which the 
\code{\link[CDM:IRT.irfprob]{CDM::IRT.irfprob}} S3 method is defined 
}
  \item{iIndex}{
Optional vector with item indices
}
  \item{theta}{
Optional vector with \eqn{\theta} values
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame containing \eqn{\theta} values and corresponding
true scores \eqn{\tau( \theta ) }.
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

\seealso{
See also \code{\link[sirt:truescore.irt]{sirt::truescore.irt}} 
for a conversion function for generalized partial credit models.
}

\examples{
#############################################################################
# EXAMPLE 1: True score conversion for a test with polytomous items
#############################################################################

data(data.Students, package="CDM")
dat <- data.Students[ , paste0("mj",1:4) ]

# fit partial credit model
mod1 <- tam.mml( dat ,control=list(maxiter=20) )
summary(mod1)

# true score conversion
tmod1 <- IRT.truescore( mod1 )
round( tmod1 , 4 )
# true score conversion with user-defined theta grid
tmod1b <- IRT.truescore( mod1 , theta=seq( -8,8, len=33 ) )
# plot results
plot( tmod1$theta , tmod1$truescore , type="l" , 
            xlab=expression(theta) , ylab=expression(tau( theta ) ) )
points( tmod1b$theta , tmod1b$truescore , pch=16 , col="brown" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{True scores}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
