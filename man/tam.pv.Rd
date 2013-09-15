\name{tam.pv}
\alias{tam.pv}

\title{
Plausible Value Imputation
}
\description{
Plausible value imputation for objects of the classes \code{tam} and \code{tam.mml}
(Adams & Wu, 2007).
}

\usage{
tam.pv(tamobj, nplausible = 10, ntheta = 2000, normal.approx = FALSE, 
    samp.regr = FALSE , np.adj=8 )
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{tamobj}{
     Object of class \code{tam} or \code{tam.mml}
}
 \item{nplausible}{
     Number of plausible values to be drawn
}
\item{ntheta}{
    Number of ability nodes for plausible value imputation. Note
    that in this function ability nodes are simulated for the
    whole sample, not for every person (contrary to the software ConQuest).
}
  \item{normal.approx}{
     An optional logical indicating whether the individual posterior distributions 
     should be approximated by a normal distribution? 
     The default is \code{FALSE}. In the case \code{normal.approx=TRUE}
     (normal distribution approximation), the number of ability nodes 
     \code{ntheta} can be substantially smaller than 2000, say 200 or 500.
     Note that normal approximation is only implemented for unidimensional
     models.
}
  \item{samp.regr}{
     An optional logical indicating whether regression coefficients 
     should be fixed in the plausible value imputation or
     also sampled from their posterior distribution? 
     The default is \code{FALSE}. Note that sampling of regression
     coefficients is only implemented for unidimensional models.
}
  \item{np.adj}{
     This parameter defines the "spread" of the random theta values
     for drawing plausible values when \code{normal.approx=FALSE}.
     If \eqn{s_{EAP}} denotes the standard deviation of the posterior
     distribution of theta (in the one-dimensional case), then theta
     is simulated from a normal distribution with standard deviation
     \code{np.adj} times \eqn{s_{EAP}}.
        }
\item{\dots}{Further arguments to be passed}	        
}

%\details{
%% For details about the implemented method see Adams and Wu (2007).
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries:
  \item{pv}{
    A data frame containing a person identifier (\code{pid})
    and plausible values denoted by \code{PVxx.Dimyy} which
    is the \code{xx}th plausible value of
    dimension \code{yy}.
        }
  \item{hwt}{Individual posterior distribution evaluated at
    the ability grid \code{theta}
    }
  \item{hwt1}{Cumulated individual posterior distribution}
  \item{theta}{Simulated ability nodes}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
# EXAMPLE 1: Dichotomous unidimensional data sim.rasch
#############################################################################

data(sim.rasch)
resp <- sim.rasch[ 1:500 , 1:15 ]  # select subsample of students and items

# estimate Rasch model
mod <- tam.mml(resp)

# draw 5 plausible values without a normality
# assumption of the posterior and 2000 ability nodes
pv1a <- tam.pv( mod , nplausible=5 , ntheta=2000 )

# draw 5 plausible values with a normality
# assumption of the posterior and 500 ability nodes
pv1b <- tam.pv( mod , nplausible=5 , ntheta=500 , normal.approx=TRUE )

# distribution of first plausible value from imputation pv1
hist(pv1a$pv$PV1.Dim1 )
# boxplot of all plausible values from imputation pv2
boxplot(pv1b$pv[ , 2:6 ] )

\dontrun{
#############################################################################
# EXAMPLE 2: Unidimensional plausible value imputation with
#            background variables; dataset data.pisaRead from sirt package
#############################################################################

data(data.pisaRead, package="sirt")
dat <- data.pisaRead$data
##   > colnames(dat)
##    [1] "idstud"   "idschool" "female"   "hisei"    "migra"    "R432Q01" 
##    [7] "R432Q05"  "R432Q06"  "R456Q01"  "R456Q02"  "R456Q06"  "R460Q01" 
##   [13] "R460Q05"  "R460Q06"  "R466Q02"  "R466Q03"  "R466Q06" 

## Note that reading items have variable names starting with R4

# estimate 2PL model without covariates
items <- grep("R4" , colnames(dat) )	# select test items from data
mod2a <- tam.mml.2pl( resp=dat[,items] )
summary(mod2a)

# fix item parameters for plausible value imputation
   # fix item intercepts by defining xsi.fixed
xsi0 <- mod2a$xsi$xsi
xsi.fixed <- cbind( seq(1,length(xsi0)) , xsi0 )
   # fix item slopes using mod2$B
# matrix of latent regressors female, hisei and migra
Y <- dat[ , c("female" , "hisei" , "migra") ]
mod2b <- tam.mml( resp=dat[,items] , B=mod2a$B , xsi.fixed=xsi.fixed , Y=Y)

# plausible value imputation with normality assumption
# and ignoring uncertainty about regression coefficients
#    -> the default is samp.regr=FALSE
pv2c <- tam.pv( mod2b , nplausible=10 , ntheta=500 )
# sampling of regression coefficients
pv2d <- tam.pv( mod2b , nplausible=10 , ntheta=500 , samp.regr=TRUE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plausible value imputation}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
