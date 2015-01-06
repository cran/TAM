\name{tam.pv}
\alias{tam.pv}

\title{
Plausible Value Imputation
}
\description{
Plausible value imputation for objects of the classes \code{tam} and \code{tam.mml}
(Adams & Wu, 2007). For converting generated plausible values into
a list of multiply imputed datasets see \code{\link{tampv2datalist}}
and the Examples 2 and 3 of this function.
}

\usage{
tam.pv(tamobj, nplausible = 10, ntheta = 2000, normal.approx = FALSE, 
    samp.regr = FALSE, theta.model=FALSE, np.adj=8 )
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
     The normal approximation is implemented for unidimensional and
     multidimensional models.
}
  \item{samp.regr}{
     An optional logical indicating whether regression coefficients 
     should be fixed in the plausible value imputation or
     also sampled from their posterior distribution? 
     The default is \code{FALSE}. Sampled regression coefficients are
     obtained by nonparametric bootstrap.
}
\item{theta.model}{Logical indicating whether the theta grid from the
	\code{tamobj} object should be used for plausible value
	imputation. In case of \code{normal.approx=TRUE}, this should
	be sufficient in many applications.
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
pv2c <- tam.pv( mod2b , nplausible=10 , ntheta=500 , normal.approx=TRUE )
# sampling of regression coefficients
pv2d <- tam.pv( mod2b , nplausible=10 , ntheta=500 , samp.regr=TRUE)
# sampling of regression coefficients, normal approximation using the
# theta grid from the model
pv2e <- tam.pv( mod2b , samp.regr=TRUE , theta.model=TRUE , normal.approx=TRUE)

#--- create list of multiply imputed datasets with plausible values
# define dataset with covariates to be matched
Y <- dat[ , c("idstud" , "idschool" , "female" , "hisei" , "migra") ]
# define plausible value names
pvnames <- c("PVREAD")
# create list of imputed datasets
datlist1 <- tampv2datalist( pv2e , pvnames = pvnames , Y=Y , Y.pid="idstud")
str(datlist1)

# create a matrix of covariates with different set of students than in pv2e
Y1 <- Y[ seq( 1 , 600 , 2 ) , ]
# create list of multiply imputed datasets
datlist2 <- tampv2datalist( pv2e , pvnames = c("PVREAD"), Y=Y1 , Y.pid="idstud")

#--- fit some models in lavaan and semTools
library(lavaan)
library(semTools)

#*** Model 1: Linear regression
lavmodel <- "
   PVREAD ~ migra + hisei
   PVREAD ~~ PVREAD
        "
mod1 <- semTools::lavaan.mi( lavmodel , data = datlist1 , m=0)        
summary(mod1 , standardized=TRUE, rsquare=TRUE)

# apply lavaan for third imputed dataset
mod1a <- lavaan::lavaan( lavmodel , data = datlist1[[3]] )        
summary(mod1a , standardized=TRUE, rsquare=TRUE)

# compare with mod1 by looping over all datasets
mod1b <- lapply( datlist1 , FUN = function(dat0){
    mod1a <- lavaan( lavmodel , data = dat0 )        
    coef( mod1a)
        } )
mod1b        
mod1b <- matrix( unlist( mod1b ) , ncol= length( coef(mod1)) , byrow=TRUE )
mod1b
round( colMeans(mod1b) , 3 )
coef(mod1)   # -> results coincide

#*** Model 2: Path model
lavmodel <- "
   PVREAD ~ migra + hisei
   hisei ~ migra
   PVREAD ~~ PVREAD
   hisei ~~ hisei
        "
mod2 <- semTools::lavaan.mi( lavmodel , data = datlist1 )        
summary(mod2 , standardized=TRUE, rsquare=TRUE)
# fit statistics
inspect( mod2 , what="fit")

#--- using mice
library(mice)
library(miceadds)
# convert datalist into a mids object
mids1 <- miceadds::datalist2mids( datlist1 )
# fit linear regression
mod1c <- with( mids1 , lm( PVREAD ~ migra + hisei ) )
summary( pool(mod1c) )

#############################################################################
# SIMULATED EXAMPLE 3: Multidimensional plausible value imputation
#############################################################################

# (1) simulate some data
set.seed(6778)
library(mvtnorm)
N <- 1000      
Y <- cbind( rnorm( N ) , rnorm(N) )
theta <- rmvnorm( N,mean=c(0,0), sigma=matrix( c(1,.5,.5,1) , 2 , 2 ))
theta[,1] <- theta[,1] + .4 * Y[,1] + .2 * Y[,2]  # latent regression model
theta[,2] <- theta[,2] + .8 * Y[,1] + .5 * Y[,2]  # latent regression model
I <- 20
p1 <- plogis( outer( theta[,1] , seq( -2 , 2 , len=I ) , "-" ) )
resp1 <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
p1 <- plogis( outer( theta[,2] , seq( -2 , 2 , len=I ) , "-" ) )
resp2 <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
resp <- cbind(resp1,resp2)
colnames(resp) <- paste("I" , 1:(2*I), sep="")
  
# (2) define loading Matrix
Q <- array( 0 , dim = c( 2*I , 2 ))
Q[cbind(1:(2*I), c( rep(1,I) , rep(2,I) ))] <- 1

# (3) fit latent regression model
mod <- tam.mml( resp=resp , Y=Y , Q=Q , control=list(maxiter=5) )

# (4) draw plausible values with normal approximation using the orginal theta grid
pv1 <- tam.pv( mod , normal.approx=TRUE , theta.mod = TRUE )

# (5) convert plausible values to list of imputed datasets
Y1 <- data.frame(Y)
colnames(Y1) <- paste0("Y",1:2)
pvnames <- c("PVFA","PVFB")
# create list of imputed datasets
datlist1 <- tampv2datalist( pv1 , pvnames = pvnames , Y=Y1 )
str(datlist1)

# (6) apply statistical models
library(semTools)
# define linear regression
lavmodel <- "
   PVFA ~ Y1 + Y2
   PVFA ~~ PVFA
        "
mod1 <- semTools::lavaan.mi( lavmodel , data = datlist1 )        
summary(mod1 , standardized=TRUE, rsquare=TRUE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plausible value imputation}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
