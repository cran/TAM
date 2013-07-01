\name{tam}

\alias{tam}
\alias{tam.mml}
\alias{tam.mml.2pl}
\alias{tam.mml.mfr}
%%\alias{tam.jml}
%%\alias{tam.jml2}


\title{
Test Analysis Modules: Marginal Maximum Likelihood Estimation
}
\description{
Modules for psychometric test analysis demonstrated 
  with the help of artificial example data. 
  The package includes MML and JML estimation of 
  uni- and multidimensional IRT (Rasch,  2PL, Generalized Partial 
  Credit, Rating Scale, Multi Facets) 
  models, fit statistic computation, standard error estimation, as 
  well as plausible value imputation and weighted likelihood
  estimation of ability. 
  Wrapper and outlining methods (S3, wrapper, plot, print, summary) 
  are still to come.
}
\usage{
tam(resp , irtmodel ="1PL" , formulaA=NULL, ...)

tam.mml(resp, Y = NULL, group = NULL, irtmodel = "1PL", formulaY = NULL,
    dataY = NULL, ndim = 1,  pid = NULL, xsi.fixed = NULL, xsi.inits = NULL, 
    beta.fixed = NULL, beta.inits = NULL, variance.fixed = NULL, 
    variance.inits = NULL, est.variance = FALSE, A = NULL, B = NULL, 
    B.fixed = NULL,  Q = NULL, R = NULL, est.slopegroups = NULL, E = NULL, 
    pweights = NULL, control = list())

tam.mml.2pl(resp, Y = NULL, group = NULL, irtmodel = "2PL", formulaY = NULL, 
    dataY = NULL, ndim = 1, pid = NULL, xsi.fixed = NULL, xsi.inits = NULL,
    beta.fixed = NULL, beta.inits = NULL, variance.fixed = NULL, 
    variance.inits = NULL, est.variance = FALSE, A = NULL, B = NULL, 
    B.fixed = NULL, Q = NULL, R = NULL,  est.slopegroups = NULL, E = NULL, 
    pweights = NULL, control = list())
            
tam.mml.mfr(resp, Y = NULL, group = NULL, irtmodel = "1PL", formulaY = NULL,
    dataY = NULL, ndim = 1, pid = NULL, xsi.fixed = NULL, xsi.inits = NULL, 
    beta.fixed = NULL, beta.inits = NULL, variance.fixed = NULL , 
    variance.inits = NULL , est.variance = FALSE , formulaA=~item+item:step, 
    constraint="cases", A=NULL , B=NULL ,  B.fixed = NULL , Q=NULL , 
    facets=NULL, est.slopegroups=NULL , E = NULL , pweights = NULL , control = list())
}
%% USAGE
%%tam.jml(resp, group = NULL, disattenuate = FALSE, bias = TRUE, 
%%    xsi.fixed = NULL, xsi.inits = NULL,  A = NULL, B = NULL, Q = NULL, 
%%    R = NULL, ndim = 1, pweights = NULL, control = list())

%%tam.jml2(resp, group = NULL, disattenuate = FALSE, bias = TRUE, 
%%    xsi.fixed = NULL, xsi.inits = NULL,  A = NULL, B = NULL, 
%%    Q = NULL, R = NULL, ndim = 1, pweights = NULL, control = list())

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
     Data frame with polytomous item responses \eqn{k=0,...,K}.
     Missing responses must be declared as \code{NA}.
}
  \item{Y}{
     A matrix of covariates in latent regression. Note that the
     intercept is automatically included as the first predictor.
}
  \item{group}{
     An optional vector of group identifiers
}
  \item{irtmodel}{
     For fixed item slopes (in \code{tam.mml}) 
     options include \code{PCM} (partial credit model)
     and \code{RSM} (not yet implemented). \cr
     For estimated item slopes (in \code{tam.mml.2pl})
     options are \code{2PL} (all slopes of item categories
     are estimated), \code{GPCM} (generalized partial credit
     model in which every item gets one and only slope
     parameter per dimension) and \code{2PL.groups}
     (subsets of items get same item slope estimates)
     and a design matrix \code{E} on item slopes in the
     generalized partial credit model (\code{GPCM.design},
     see Examples).
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
  \item{ndim}{
     Number of dimensions (is not needed to determined by the user)
}
  \item{pid}{
     An optional vector of person identifiers
}
  \item{xsi.fixed}{
     A matrix with two columns for fixing \eqn{\xi} parameters. 
     1st column: index of \eqn{\xi} parameter, 2nd column: fixed value 
}
  \item{xsi.inits}{
     A vector with initial value for \eqn{\xi}
}
  \item{beta.fixed}{
     A matrix with two columns for fixing regression coefficients.
     1st column: Index of \eqn{Y} value, 2nd column: dimension,
     3rd column: fixed \eqn{\beta} value
}
  \item{beta.inits}{
     A matrix with initial \eqn{\beta} values
}
  \item{variance.fixed}{
     An optional matrix with three columns for fixing
     entries in covariance matrix:
     1st column: dimension 1, 2nd column: dimension 2,
     3rd column: fixed value
}
  \item{variance.inits}{
     Initial covariance matrix in estimation
}
  \item{est.variance}{
     Should the covariance matrix be estimated? This argument
     applies to estimated item slopes in \code{tam.mml.2pl}.
     The default is \code{FALSE} which means that latent
     variables (in the first group) are standardized in 2PL estimation.
}
\item{constraint}{Set sum constraint for parameter identification 
	for \code{items} or \code{cases} (applies to \code{tam.mml.mfr})
	}
  \item{A}{
     An optional array of dimension \eqn{ I \times (K+1) \times N_\xi}.
     Only \eqn{\xi} parameters are estimated, entries in \eqn{A}
     only correspond to the design.
}
  \item{B}{
     An optional array of dimension \eqn{ I \times (K+1) \times D}.
     In case of \code{tam.mml.2pl} entries of the \eqn{B} matrix can be 
     estimated.
}
  \item{B.fixed}{
  An optional matrix with four columns for fixing \eqn{B} matrix entries in 2PL estimation.
  1st column: item index, 2nd column: category, 3rd column: dimension,
  4th column: fixed value.
}
  \item{Q}{
  An optional \eqn{I \times D} matrix which specifies the loading structure
  of items on dimensions.  
}
  \item{R}{
     Not yet implemented
}
  \item{est.slopegroups}{
  A vector of integers of length \eqn{I} for estimating item slope
  parameters of item groups. This function only applies to
  the generalized partial credit model (\code{irtmodel="2PL.groups"}).
}
  \item{E}{
     An optional design matrix for estimating item slopes
     in the generalized partial credit model (\code{irtmodel="GPCM"})
}
  \item{pweights}{
     An optional vector of person weights
}
\item{formulaA}{Design formula (only applies to \code{tam.mml.mfr}). 
    See Example 8.}
\item{facets}{Data frame with facet entries (only applies to \code{tam.mml.mfr})}
  \item{control}{
  A list of control arguments for the algorithm: \cr
  \code{ list( nodes = seq(-6,6,len=21) , snodes = 0 , QMC=TRUE,}\cr
    \code{        convD = .001 ,conv = .0001 , convM = .0001 , Msteps = 4 ,   }\cr
    \code{        maxiter = 1000 , max.increment = 1 , }\cr
    \code{        min.variance = .001 , progress = TRUE , ridge=0 ,    }\cr
    \code{        seed = NULL , xsi.start0=FALSE)  } \cr
\code{nodes}: the discretized \eqn{\theta} nodes for numerical integration \cr
\code{snodes}: number of simulated \eqn{\theta} nodes for stochastic integration. 
    If \code{QMC=TRUE}, then quasi Monte Carlo integration (Pan & Thompson, 2007)
    is used which is a nonstochastic integration but prevents the very demanding
    numeric integration using Gaussian quadrature. If \code{QMC=FALSE},
    then "ordinary" stochastic integration is used.\cr
    If \code{snodes=0}, numerical integration is used. \cr
\code{xsi.start0}: Should all \eqn{\xi} starting values be initially set to
   zero? This can be especially useful in complex multifaceted designs. \cr 
\code{\dots}: EXPLAIN OTHER ARGUMENTS
}
\item{\dots}{Further arguments to be passed}
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\details{
The multidimensional item response model in \pkg{TAM} is well
described in Adams, Wilson and Wu (1997) or Adams and Wu (2007).

Setup of the model: \cr
The data frame \code{resp} contains \eqn{N} persons (in rows) 
at \eqn{I} items (in columns) with \eqn{K} categories \eqn{k=0,...,K}.
The item reponse model has \eqn{D} dimensions of the \eqn{\theta} ability
vector and can be written as
\deqn{ 
        P( X_{pi} = k | \theta_p )  \propto exp( b_{ik} \theta_p + a_{ik} \xi ) 
        }
%%% is this correct???
Item category thresholds for item \eqn{i} in category \eqn{k} 
are written as a linear combination
\eqn{a_i \xi} where the vector \eqn{\xi} of length \eqn{N_\xi}
contains generalized item parameters and \eqn{A= ( a_{ik} )_{ik}} is 
a three-dimensional design array (specified in \code{A}).

The scoring vector \eqn{b_{ik}} contains the fixed (in \code{tam.mml})
or estimated (in \code{tam.mml.2pl}) scores of item \eqn{i} in category \eqn{k}
on dimension \eqn{d}.

The symbol \eqn{\propto } means that response probabilities are normalized such
that \eqn{ \sum _k P( X_{pi} = k | \theta_p ) = 1 }.

The latent regression model regresses the latent trait \eqn{\theta_p}
on covariates \eqn{Y} which results in
\deqn{ 
    \theta_p = Y \beta + \epsilon_p , \epsilon_p \sim N_D ( 0 , \Sigma )
}{%
  \theta_p = Y \beta + \epsilon_p , \epsilon_p ~ N_D ( 0 , \Sigma )
    }

Where \eqn{\beta} is a \eqn{N_Y} times \eqn{D} matrix of regression
coefficients for \eqn{N_Y} covariates in \eqn{Y}.

Until now, the multiple group model for groups 
\eqn{g=1,...,G} works only in unidimensional item
response models. In this case, variance heterogeneity is allowed
\deqn{ 
    \theta_p = Y \beta + \epsilon_p , \epsilon_p \sim N ( 0 , \sigma_g^2 )
}{%
  \theta_p = Y \beta + \epsilon_p , \epsilon_p ~ N( 0 , \sigma_g^2 )
    }

Wrapper for \code{\link{tam.mml}}, \code{\link{tam.jml}}. 
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:

\item{xsi}{Vector of \eqn{\xi} parameter estimates and their
corresponding standard errors}
\item{xsi.facets}{Data frame of \eqn{\xi} parameters and corresponding
constraints for multifacet models}
\item{beta}{Matrix of \eqn{\beta} regression coefficient estimates}
\item{variance}{Covariance matrix. In case of multiple groups, 
    it is a vector indicating heteroscedastic variances}
\item{item}{Data frame with item parameters}
\item{person}{Matrix with person parameter estimates. 
    \code{EAP} is the mean of the posterior distribution and \code{SD.EAP} 
    the corresponding standard deviation
        }
\item{pid}{Vector of person identifiers}
\item{EAP.rel}{EAP reliability}
\item{post}{Posterior distribution for item response pattern}
\item{rprobs}{A three-dimensional array with estimated response probabilities
        (dimensions are items times categories times theta length)}
\item{itemweight}{Matrix of item weights}
\item{theta}{Theta grid}
\item{n.ik}{Array of expected counts: theta class x item x category x group}
\item{pi.k}{Marginal trait distribution at grid \code{theta}} 
\item{Y}{Matrix of covariates}
\item{resp}{Original data frame}
\item{resp.ind}{Response indicator matrix}
\item{group}{Group identifier}
\item{G}{Number of groups}
\item{formulaY}{Formula for latent regression}
\item{dataY}{Data frame for latent regression}
\item{pweights}{Person weights}
\item{time}{Computation time}
\item{A}{Design matrix \eqn{A} for \eqn{\xi} parameters}
\item{B}{Fixed or estimated loading matrix}
\item{se.B}{Standard errors of \eqn{B} parameters}
\item{nitems}{Number of items}
\item{maxK}{Maximum number of categories}
\item{AXsi}{Estimated item intercepts \eqn{ a_{ik} \xi}}
\item{se.AXsi}{Standard errors of \eqn{a_{ik} \xi} parameters}
\item{nstud}{Number of students}
\item{resp.ind.list}{List of response indicator vectors}
\item{hwt}{Individual posterior distribution}
\item{ndim}{Number of dimensions}
\item{xsi.fixed}{Fixed \eqn{\xi} parameters}
\item{beta.fixed}{Fixed \eqn{\beta} parameters}
\item{formulaA}{Design formula (only applies to \code{tam.mml.mfr})}
\item{facets}{Data frame with facet entries (only applies to \code{tam.mml.mfr})}
\item{variance.fixed}{Fixed covariance matrix}
\item{nnodes}{Number of theta nodes}
\item{deviance}{Final deviance}
\item{ic}{Vector with information criteria}
\item{deviance.history}{Deviance history in iterations}
\item{control}{List of control arguments}
\item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Adams, R. J., Wilson, M. & Wu, M. (1997).
Multilevel item response models: An approach to errors in
variables regression. \emph{Journal of Educational and Behavioral 
Statistics}, \bold{22}, 47-76.

Adams, R. J., & Wu, M. L. (2007). The mixed-coefficients multinomial logit model.
A generalized form of the Rasch model. In M. von Davier & C. H. Carstensen (Eds.):
\emph{Multivariate and mixture distribution Rasch models: Extensions and applications}
(pp. 55-76). New York: Springer.

Pan, J. & Thompson, R. (2007). Quasi-Monte Carlo estimation in
generalized linear mixed models. \emph{Computational Statistics &
Data Analysis}, \bold{51}, 5765-5775.
}

%\author{
%\pkg{TAM} authors
%}

\note{
For more than three dimensions, quasi-Monte Carlo or stochastic integration
is recommended. Choose in \code{control} a value of \code{snodes} (say, larger than 1000).
See Pan and Thompson (2007) for more details.
}

\section{Warning}{
Be cautious when using the number of estimated parameters in \code{ic}.
Extensive checking of errors has not be done yet.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{tam.mml}}, \code{\link{tam.jml}}.

Standard errors are estimated by a rather crude (but quick) approximation.
Use \code{\link{tam.se}} for improved standard errors.

For model comparisons see \code{\link{anova.tam}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#####################################
# Example 1: dichotomous data
# sim.rasch: 2000 persons, 40 items
data(sim.rasch)

#***
# Rasch model (MML estimation)
mod1 <- tam.mml(resp=sim.rasch) 
# extract item parameters
mod1$item    # item difficulties

# WLE estimation
wle1 <- tam.wle( mod1 )
# item fit
fit1 <- tam.fit(mod1)
# plausible value imputation
pv1 <- tam.pv(mod1)
# standard errors
se1 <- tam.se( mod1 )
# summary
summary(mod1)

#***
# 2PL model
mod2 <- tam.mml.2pl(resp=sim.rasch,irtmodel="2PL") 

# extract item parameters
mod2$xsi    # item difficulties
mod2$B      # item slopes

\dontrun{
#****
# constrained 2PL estimation
# estimate item parameters in different slope groups
# items 1-10, 21-30 group 1
# items 11-20 group 2 and items 31-40 group 3
est.slope <- rep(1,40)
est.slope[ 11:20 ] <- 2
est.slope[ 31:40 ] <- 3
mod3 <- tam.mml.2pl( resp=sim.rasch , irtmodel="2PL.groups" , 
    est.slopegroups = est.slope )
mod3$B
summary(mod3)

###############################
# Example 2: Unidimensional calibration with latent regressors
#   simulated data

# (1) simulate data
set.seed(6778)
N <- 2000     # number of persons
# latent regressors Y
Y <- cbind( rnorm( N , sd = 1.5) , rnorm(N , sd = .3 ) )
# simulate theta
theta <- rnorm( N ) + .4 * Y[,1] + .2 * Y[,2]  # latent regression model
# number of items
I <- 40
p1 <- plogis( outer( theta , seq( -2 , 2 , len=I ) , "-" ) )
# simulate response matrix
resp <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
colnames(resp) <- paste("I" , 1:I, sep="")
  
# (2) estimate model
mod2_1 <- tam(resp=resp , Y=Y)
summary(mod2_1)

#################################
# Example 3: Multiple group estimation

# (1) simulate data
set.seed(6778)
N <- 3000
theta <- c( rnorm(N/2,mean=0,sd = 1.5) , rnorm(N/2,mean=.5,sd = 1)  )
I <- 20
p1 <- plogis( outer( theta , seq( -2 , 2 , len=I ) , "-" ) )
resp <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
colnames(resp) <- paste("I" , 1:I, sep="")
group <- rep(1:2 , each=N/2 )

# (2) estimate model
mod3_1 <- tam.mml( resp ,  group = group )
summary(mod3_1)

###########################################################
# Example 4: Multidimensional estimation
# with two dimensional theta's - simulate some bivariate data, 
# and regressors
# 40 items: first 20 items load on dimension 1, 
#           second 20 items load on dimension 2

# (1) simulate some data
set.seed(6778)
library(mvtnorm)
N <- 1000      
Y <- cbind( rnorm( N ) , rnorm(N) )
theta <- rmvnorm( N,mean=c(0,0), sigma=matrix( c(1,.5,.5,1) , 2 , 2 ))
theta[,1] <- rnorm( N ) + .4 * Y[,1] + .2 * Y[,2]  # latent regression model
theta[,2] <- rnorm( N ) + .8 * Y[,1] + .5 * Y[,2]  # latent regression model
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

# (3) estimate models 

# Rasch model: without regressors
mod4_1 <- tam.mml( resp=resp , Q=Q )

# Rasch model: set covariance between dimensions to zero
variance_fixed <- cbind( 1 , 2 , 0 )
mod4_2 <- tam.mml( resp=resp , Q=Q , variance.fixed = variance_fixed )
summary(mod4_2)

# 2PL model
mod4_3 <- tam.mml.2pl( resp=resp , Q=Q , irtmodel="2PL" )

# Rasch model with 2000 quasi monte carlo nodes and a maximum of 15 iterations
# -> nodes are useful for more than 3 or 4 dimensions
mod4_4 <- tam.mml( resp=resp , Q=Q , 
        control=list(snodes=2000,maxiter=15) )

# Rasch model with 2000 stochastic nodes
mod4_5 <- tam.mml( resp=resp , Q=Q , 
        control=list(snodes=2000,QMC=FALSE,maxiter=15) )

# estimate two dimensional Rasch model with regressors
mod4_6 <- tam.mml( resp=resp , Y=Y , Q=Q )

#########################################################################
# Example 5: Two dimensional estimation with within item dimensionality

# (1) simulate data
set.seed(4762)
N <- 2000 # 2000 persons
Y <- rnorm( N )
theta <- rmvnorm( N,mean=c(0,0), sigma=matrix( c(1,.5,.5,1) , 2 , 2 ))
I <- 10
# 10 items load on the first dimension
p1 <- plogis( outer( theta[,1] , seq( -2 , 2 , len=I ) , "-" ) )
resp1 <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
# 10 items load on the second dimension
p1 <- plogis( outer( theta[,2] , seq( -2 , 2 , len=I ) , "-" ) )
resp2 <- 1 * ( p1 > matrix( runif( N*I ) , nrow=N , ncol=I ) )
# 20 items load on both dimensions
p1 <- plogis( outer( 0.5*theta[,1] + 1.5*theta[,2] , seq( -2 , 2 , len=2*I ) , "-" ) )
resp3 <- 1 * ( p1 > matrix( runif( N*2*I ) , nrow=N , ncol=2*I ) )
#Combine the two sets of items into one response matrix
resp <- cbind(resp1, resp2, resp3 )
colnames(resp) <- paste("I" , 1:(4*I), sep="")

# (2) define loading matrix 
Q <- cbind(c(rep(1,10),rep(0,10),rep(1,20)), c(rep(0,10),rep(1,10),rep(1,20)))

# (3) model: within item dimensionality and 2PL estimation  
mod5 <- tam.mml.2pl(resp, Q=Q, irtmodel="2PL" , 
		control=list(maxiter=20) )
summary(mod5)

# item difficulties
mod5$item
# item loadings
mod5$B

##############################################
# Example 6: ordered data - Generalized partial credit model
data( data.gpcm )

# Ex6.1: 2PL model
mod6_1 <- tam.mml.2pl( resp=data.gpcm , irtmodel="2PL" , control=list( maxiter=200) )
mod6_1$item # item intercepts
mod6_1$B    # for every category a separate slope parameter is estimated

# Ex6.2: Generalized partial credit model
mod6_2 <- tam.mml.2pl( resp= data.gpcm , irtmodel="GPCM" , control=list( maxiter=200) )
mod6_2$B    # joint slope parameter for all categories

# Ex6.3: some fixed entries of slope matrix B
# B: nitems x maxK x ndim  
#   ( number of items x maximum number of categories x number of dimensions)
# set two constraints
B.fixed <- matrix( 0 , 2 , 4 )
# set second item, second category, first dimension to 2.3
B.fixed[1,] <- c(2,3,1,2.3)
# set third item, first category, first dimension to 1.4
B.fixed[2,] <- c(3,2,1,1.4)

# estimate item parameter with variance fixed (by default)
mod6_3 <- tam.mml.2pl( resp=data.gpcm , irtmodel="2PL" , B.fixed = B.fixed , 
            control=list( maxiter=200) )
mod6_3$B

# Ex 6.4: estimate the same model, but estimate variance
mod6_4 <- tam.mml.2pl( resp=data.gpcm , irtmodel="2PL" , B.fixed = B.fixed , 
            est.variance = TRUE , control=list( maxiter=350) )
mod6_4$B

# Ex 6.5: partial credit model
mod6_5 <- tam.mml( resp=data.gpcm ,control=list( maxiter=200) )
mod6_5$B

# Ex 6.6: partial credit model with Conquest parametrization
A1 <- .A.PCM2( resp=data.gpcm)	# create design matrix
mod6_6 <- tam.mml( resp=data.gpcm , A=A1 )
summary(mod6_6)

##############################
# Example 7: design matrix for slopes for the
#            generalized partial credit model

# (1) simulate data from a model with a design matrix
set.seed(789)
I <- 42
b <- seq( -2 , 2 , len=I)
# create design matrix for loadings
E <- matrix( 0 , I , 5 )
E[ seq(1,I,3) , 1 ] <- 1
E[ seq(2,I,3) , 2 ] <- 1
E[ seq(3,I,3) , 3 ] <- 1
ind <- seq( 1, I , 2 ) ; E[ ind , 4 ] <- rep( c( .3 , -.2 ) , I )[ 1:length(ind) ]
ind <- seq( 2, I , 4 ) ; E[ ind , 5 ] <- rep( .15 , I )[ 1:length(ind) ]
E

# true basis slope parameters
lambda <- c( 1 , 1.2 , 0.8 , 1 , 1.1 )
# calculate item slopes
a <- apply( E , 1 , FUN=function(ll){ sum(ll*lambda) } )

N <- 4000
theta <- rnorm( N )
aM <- outer( rep(1,N) , a )
bM <- outer( rep(1,N) , b )
pM <- plogis( aM * ( matrix( theta , nrow=N , ncol=I  ) - bM ) )
dat <- 1 * ( pM > runif( N*I ) )
colnames(dat) <- paste("I" , 1:I , sep="")

mod7 <- tam.mml.2pl( resp = dat , irtmodel="GPCM.design" , E=E )
mod7$B

# recalculate estimated basis parameters
lm( mod7$B[,2,1] ~ 0+ as.matrix(E ) )
##   Call:
##   lm(formula = mod7$B[, 2, 1] ~ 0 + as.matrix(E))
##   Coefficients:
##   as.matrix(E)1  as.matrix(E)2  as.matrix(E)3  as.matrix(E)4  as.matrix(E)5  
##          0.9904         1.1896         0.7817         0.9601         1.2132  

###################################################
# Example 8: Differential item functioning
#  A first example of a Multifaceted Rasch Model
data(data.ex08)
# Facet is only female; 10 items are studied

formulaA <- ~ item+female+item*female
resp <- data.ex08[["resp"]]
facets <- as.data.frame( data.ex08[["facets"]] )

mod8 <- tam.mml.mfr( resp= resp , facets= facets , 
    formulaA = formulaA )    
summary(mod8)          

####################################################
# Example 9: Multifaceted Rasch Model
data(sim.mfr) ; data(sim.facets)

# two way interaction item and rater
formulaA <- ~item+item:step + item*rater
mod9a <- tam.mml.mfr( resp=sim.mfr , facets=sim.facets , formulaA = formulaA )
mod9a$xsi.facets
summary(mod9a)

# three way interaction item, female and rater
formulaA <- ~item+item:step + female*rater + female*item*step
mod9b <- tam.mml.mfr( resp=sim.mfr , facets=sim.facets , formulaA = formulaA )
summary(mod9b)

####################################################
# Example 10: Model with raters.
#   Persons are arranged in multiple rows which is indicated
#   by multiple person identifiers.
data(data.ex10)
dat <- data.ex10
head(dat)
##     pid rater I0001 I0002 I0003 I0004 I0005
## 1     1     1     0     1     1     0     0
## 451   1     2     1     1     1     1     0
## 901   1     3     1     1     1     0     1
## 452   2     2     1     1     1     0     1
## 902   2     3     1     1     0     1     1

facets <- dat[ , "rater" , drop=FALSE ] # define facet (rater)
pid <- dat$pid      # define person identifier (a person occurs multiple times)
resp <- dat[ , -c(1:2) ]        # item response data
formulaA <- ~ item * rater      # formula

mod10 <- tam.mml.mfr( resp=resp , facets=facets , formulaA = formulaA , pid=dat$pid )
summary(mod10)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marginal maximum likelihood estimation (MML)}
%% \keyword{Joint maximum likelihood estimation (JML)}% __ONLY ONE__ keyword per line
\keyword{Multidimensional item response model}