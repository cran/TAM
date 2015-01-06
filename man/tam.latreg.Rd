\name{tam.latreg}
\alias{tam.latreg}
\alias{summary.tam.latreg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Latent Regression Model
}

\description{
This function fits a latent regression model \eqn{\bold{\theta} = Y \beta +
\varepsilon}. Only the individual likelihood evaluated at a 
\eqn{\bold{\theta}} grid is needed as the input. Like in
\code{\link{tam.mml}} a multivariate normal distribution is posed
on the residual distribution. Plausible values can be drawn by subsequent
application of \code{\link{tam.pv}} (see Example 1).
}

\usage{
tam.latreg(like, theta, Y = NULL, group = NULL, formulaY = NULL, dataY = NULL, 
   beta.fixed = NULL, beta.inits = NULL, variance.fixed = NULL, 
   variance.inits = NULL, est.variance = TRUE, pweights = NULL, pid = NULL, 
   userfct.variance = NULL, variance.Npars = NULL, control = list())
   
\method{summary}{tam.latreg}(object,file=NULL,\dots)   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{like}{
Individual likelihood. This can be typically extracted from fitted 
item response models by making use of \code{IRT.likelihood}.
}
  \item{theta}{
Used \eqn{\bold{\theta}} grid in the fitted IRT model
}
  \item{Y}{
     A matrix of covariates in latent regression. Note that the
     intercept is automatically included as the first predictor.
}
  \item{group}{
     An optional vector of group identifiers
}
  \item{formulaY}{
     An \R formula for latent regression. Transformations of predictors
     in \eqn{Y} (included in \code{dataY}) can be easily spcified,
     e. g. \code{female*race} or \code{I(age^2)}.
}
  \item{dataY}{
     An optional data frame with possible covariates \eqn{Y} in latent regression.
     This data frame will be used if an \R formula in \code{formulaY}
     is specified.
}
  \item{beta.fixed}{
     A matrix with three columns for fixing regression coefficients.
     1st column: Index of \eqn{Y} value, 2nd column: dimension,
     3rd column: fixed \eqn{\beta} value. \cr
     If no constraints should be inposed on \eqn{\beta}, then
     set \code{beta.fixed=FALSE} (see Example 2, Model \code{2_4}).
}
  \item{beta.inits}{
     A matrix (same format as in \code{beta.fixed}) 
     with initial \eqn{\beta} values
}
  \item{variance.fixed}{
     An optional matrix with three columns for fixing
     entries in covariance matrix:
     1st column: dimension 1, 2nd column: dimension 2,
     3rd column: fixed value
}
  \item{variance.inits}{
     Initial covariance matrix in estimation. All matrix entries have to be
     specified and this matrix is NOT in the same format like
     \code{variance.inits}.
}
  \item{est.variance}{
     Should the covariance matrix be estimated? This argument
     applies to estimated item slopes in \code{tam.mml.2pl}.
     The default is \code{FALSE} which means that latent
     variables (in the first group) are standardized in 2PL estimation.
}
\item{pweights}{
     An optional vector of person weights
}
  \item{pid}{
     An optional vector of person identifiers
}
\item{userfct.variance}{Optional user customized function for variance speicification
(See Simulated Example 17).}
\item{variance.Npars}{Number of estimated parameters of variance matrix
if a \code{userfct.variance} is provided.} 
\item{control}{
List of control parameters, see \code{\link{tam.mml}}.
}
  \item{object}{
Object of class \code{tam.latreg}
}
\item{file}{
A file name in which the summary output will be written 
    }
\item{\dots}{
	Further arguments to be passed
}   
}

%\details{
%}

\value{
Subset of values of \code{\link{tam.mml}}	
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
See also \code{\link{tam.pv}} for plausible value imputation.
}

\examples{
#############################################################################
# EXAMPLE 1: Unidimensional latent regression model with fitted IRT model in
#            sirt package
#############################################################################

library(sirt)
data(data.pisaRead, package="sirt")
dat <- data.pisaRead$data

items <- grep("R4" , colnames(dat) , value=TRUE )    # select test items from data
# define testlets
testlets <- substring( items , 1 , 4 )
itemcluster <- match( testlets , unique(testlets) )
# fit Rasch copula model (only few iterations)
mod <- sirt::rasch.copula2( dat[,items] , itemcluster=itemcluster , mmliter=5)
# extract individual likelihood
like1 <- IRT.likelihood( mod )
# fit latent regression model in TAM
Y <- dat[ , c("migra" , "hisei", "female") ]
mod2 <- tam.latreg( like1 , theta=attr(like1, "theta") , Y = Y , pid=dat$idstud )
summary(mod2)
# plausible value imputation
pv2 <- tam.pv( mod2 ) 
# create list of imputed datasets
Y <- dat[ , c("idstud" , "idschool" , "female" , "hisei" , "migra") ]
pvnames <- c("PVREAD")
datlist <- tampv2datalist( pv2, pvnames = pvnames , Y=Y , Y.pid="idstud")

\dontrun{
#--- fit some models
library(mice)
library(miceadds)
# convert data list into a mice object
mids1 <- miceadds::datalist2mids( datlist )
# perform an ANOVA
mod3a <- with( mids1 , lm(PVREAD ~ hisei*migra) )
summary( pool( mod3a ))
mod3b <- miceadds::mi.anova( mids1 , "PVREAD ~ hisei*migra" )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Latent regression model}
\keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
