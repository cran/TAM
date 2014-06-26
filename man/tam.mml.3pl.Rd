\name{tam.mml.3pl}
\alias{tam.mml.3pl}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
3PL Structured Item Response Model in \pkg{TAM}
}

\description{
This estimates a 3PL model with design matrices for item slopes and
item intercepts. Discrete distributions of the latent variables are
allowed which can be log-linearly smoothed.		
}


\usage{
tam.mml.3pl(resp, Y = NULL, group = NULL, formulaY = NULL, dataY = NULL, ndim = 1, 
  pid = NULL, xsi.fixed = NULL, xsi.inits = NULL, xsi.prior = NULL, 
  beta.fixed = NULL, beta.inits = NULL, variance.fixed = NULL, variance.inits = NULL, 
  est.variance = TRUE, A = NULL, notA=FALSE , Q = NULL, E = NULL,gammaslope.des = "2PL", 
  gammaslope = NULL, gammaslope.fixed = NULL, 
  est.some.slopes = TRUE, gammaslope.constr.V = NULL, gammaslope.constr.c = NULL, 
  gammaslope.prior = NULL, est.guess = NULL, guess = rep(0, ncol(resp)), 
  guess.prior = NULL, skillspace = "normal", theta.k = NULL, delta.designmatrix = NULL, 
  delta.fixed = NULL, pweights = NULL, control = list() )
}
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
     A matrix with two columns (in the same way defined as in
     \code{xsi.fixed} with initial value for \eqn{\xi}.
}
  \item{xsi.prior}{
     Item-specific prior distribution for \eqn{\xi} parameters. It is
     assumed that \eqn{\xi_k \sim N( \mu_k , \sigma_k^2 )}. The first column
     in \code{xsi.prior} is \eqn{\mu_k}, the second is \eqn{\sigma_k}.     
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
  \item{A}{
     An optional array of dimension \eqn{ I \times (K+1) \times N_\xi}.
     Only \eqn{\xi} parameters are estimated, entries in \eqn{A}
     only correspond to the design.
}
\item{notA}{An optional logical indicating whether all entries in
	the \eqn{A} matrix are set to zero and no item intercept
	\eqn{\xi} should be estimated.}	
  \item{Q}{
  An optional \eqn{I \times D} matrix (the Q-matrix) which specifies the 
  loading structure of items on dimensions.  
}
  \item{E}{
Optional design array for item slopes \eqn{\gamma}. It is a four 
dimensional array of size \eqn{I \times (K+1) \times D \times N_\gamma}
containing items, categories, dimensions, \eqn{\gamma} parameter.
}
\item{gammaslope.des}{
Optional string indicating type of item slope parameter to be estimated.
\code{gammaslope.des="2PL"} estimates a slope parameter for an item,
\code{gammaslope.des="2PLcat"} for an item and a 
}
  \item{gammaslope}{
Initial or fixed vector of \eqn{\gamma} parameters
}
  \item{gammaslope.fixed}{
An optional matrix containing fixed values in the \eqn{\gamma} vector.
First column: parameter index; second colunmn: fixed value.
}
  \item{est.some.slopes}{
An optional logical indicating whether some item slopes should be estimated.
}
  \item{gammaslope.constr.V}{
An optional constraint matrix \eqn{V} for item slope parameters \eqn{\gamma}
}
  \item{gammaslope.constr.c}{
An optional constraint vector \eqn{c} for item slope parameters \eqn{\gamma}
}
  \item{gammaslope.prior}{
     Item-specific prior distribution for \eqn{\gamma} parameters. It is
     assumed that \eqn{\gamma_k \sim N( \mu_k , \sigma_k^2 )}. The first column
     in \code{gammaslope.prior} is \eqn{\mu_k}, the second is \eqn{\sigma_k}. 
}
  \item{est.guess}{
An optional vector of integers indicating for which items a guessing
parameter should be estimated. Same integers correspond to same estimated 
guessing parameters. A value of 0 denotes an item for which no guessing
parameter is estimated.
}
  \item{guess}{
Fixed or initial guessing parameters
}
  \item{guess.prior}{
     Item-specific prior distribution for guessing parameters \eqn{c_i}. It is
     assumed that \eqn{c_i \sim N( a_i , b_i)}. The first column
     in \code{gammaslope.prior} is \eqn{a_i}, the second is \eqn{b_i}.
}
  \item{skillspace}{
Skill space: normal distribution (\code{"normal"}) or discrete
distribution (\code{"discrete"}).
}
  \item{theta.k}{
A matrix of the \eqn{\bold{\theta}} skill space in case of a discrete
distribution (\code{skillspace="discrete"}).
}
  \item{delta.designmatrix}{
A design matrix of the log-linear model for the skill space in case of a discrete
distribution (\code{skillspace="discrete"}).

}
  \item{delta.fixed}{
Fixed \eqn{\delta} values of the log-linear skill space.
}
  \item{pweights}{
Optional vector of person weights.
}
  \item{control}{
See \code{\link{tam.mml}} for more details.
}
}

\details{
The item response model for item \eqn{i} and category \eqn{h} and no guessing
parameters can be written as
\deqn{ P( X_{i} = h | \bold{\theta} ) \propto exp( \sum_d b_{ihd} \theta_d +  
 \sum_k a_{ih} \xi_k ) }
For item slopes, a linear decomposition is allowed 
\deqn{ b_{ihd} =  \sum_k   e_{ihdk} \gamma_k }
In case of a guessing parameter, the item response function is
\deqn{ P( X_{i} = h | \bold{\theta} ) = c_i + ( 1 - c_i ) \cdot
	( 1 + exp( - \sum_d b_{ihd} \theta_d - \sum_k a_{ih} \xi_k ) )^{-1} 
		}
For possibilities of definitions of the design matrix \eqn{E=(e_{ihdk})}
see Formann (2007), for the strongly related linear logistic latent
class model see Formann (1992). Linear constraints on \eqn{\gamma} 
can be specified by equations \eqn{V \gamma = c} and using the arguments
\code{gammaslope.constr.V} and \code{gammaslope.constr.c}.

Like in \code{\link{tam.mml}}, a multivariate linear regression 
\deqn{ \bold{\theta} = Y \beta + \bold{\epsilon}}
assuming a multivariate normal distribution on the residuals \eqn{\bold{\epsilon}}
can be specified (\code{skillspace="normal"}).

Alternatively, a log-linear distribution of the skill classes \eqn{P(\theta)}
(\code{skillspace="discrete"}) 
is performed \deqn{log P(\theta ) = D_{ \delta } \delta } See Xu and
von Davier (2008). The design matrix \eqn{D_{\delta}} can be specified in
\code{delta.designmatrix}. The theta grid \eqn{\bold{\theta}} of the skill space 
can be defined in \code{theta.k}.
}


\value{
The same output as in \code{\link{tam.mml}} and additional entries		
 \item{delta}{Parameter vector \eqn{\delta}}
 \item{gammaslope}{Estimated \eqn{\gamma} item slope parameters}
 \item{E}{Used design matrix \eqn{E}}
}
% \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}


\references{
Formann, A. K. (1992). Linear logistic latent class analysis for polytomous data. 
\emph{Journal of the American Statistical Association}, \bold{87}, 476-486.	

Formann, A. K. (2007). (Almost) Equivalence between conditional and mixture maximum
likelihood estimates for some models of the Rasch type. In M. von Davier & C. H. Carstensen
(Eds.), \emph{Multivariate and mixture distribution Rasch models} (pp. 177-189). 
New York: Springer.
	
Xu, X., & von Davier, M. (2008). \emph{Fitting the structured general diagnostic
model to NAEP data}. ETS Research Report ETS RR-08-27. Princeton, ETS.
}


%\author{
%%  ~~who you are~~
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{tam.mml}}.

See the \code{\link[CDM:slca]{slca}} function in the \pkg{CDM} package
for a similar method.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Dichotomous data | sim.rasch
#############################################################################

data(sim.rasch)
dat <- sim.rasch
# some control arguments
ctl.list <- list(maxiter= 20 ) # increase the number of iterations in applications!

#*** Model 1: Rasch model, normal trait distribution
mod1 <- tam.mml.3pl(resp=dat , skillspace= "normal" , est.some.slopes=FALSE , 
        control=ctl.list)
summary(mod1)

#*** Model 2: Rasch model, discrete trait distribution
#  choose theta grid
theta.k <- seq( -3 , 3 , len=7 )	# discrete theta grid distribution
# define symmetric trait distribution
delta.designmatrix <- matrix( 0 , nrow=7 , ncol= 4 )
delta.designmatrix[4,1] <- 1
delta.designmatrix[c(3,5),2] <- 1
delta.designmatrix[c(2,6),3] <- 1
delta.designmatrix[c(1,7),4] <- 1
mod2 <- tam.mml.3pl(resp=dat , skillspace= "discrete" , est.some.slopes=FALSE , 
           theta.k=theta.k, delta.designmatrix=delta.designmatrix , control=ctl.list)
summary(mod2)

#*** Model 3: 2PL model
mod3 <- tam.mml.3pl(resp=dat , skillspace= "normal" , gammaslope.des="2PL" ,control=ctl.list)
summary(mod3)

#*** Model 4: 3PL model
# estimate guessing parameters for items 3,7,9 and 12
I <- ncol(dat)
est.guess <- rep(0, I )
# set parameters 9 and 12 equal -> same integers
est.guess[ c(3,7,9,12) ] <- c( 1 , 3 , 2 , 2 )
# starting values guessing parameter
guess <- .2*(est.guess > 0)
# estimate model
mod4 <- tam.mml.3pl(resp=dat , skillspace= "normal" , gammaslope.des="2PL" ,
        control=ctl.list , est.guess = est.guess , guess=guess)
summary(mod4)

#*** Model 5: 3PL model, add some prior Beta distribution
guess.prior <- matrix( 0 , nrow=I , ncol=2 )
guess.prior[ est.guess  > 0 , 1] <- 5
guess.prior[ est.guess  > 0 , 2] <- 17
mod5 <- tam.mml.3pl(resp=dat , skillspace= "normal" , gammaslope.des="2PL" ,
        control=ctl.list , est.guess = est.guess , guess=guess , guess.prior=guess.prior)
summary(mod5)

#*** Model 6: 2PL model with design matrix for item slopes
I <- 40	     # number of items
D <- 1       # dimension
maxK <- 2    # maximum number of categories
Ngam <- 13   # number of different slope parameters
E <- array( 0 , dim=c(I,maxK,D,Ngam) )
# joint slope parameters for items 1 to 10, 11 to 20, 21 to 30
E[ 1:10 , 2 , 1 , 2 ] <- 1
E[ 11:20 , 2 , 1 , 1 ] <- 1
E[ 21:30 , 2 , 1 , 3 ] <- 1
for (ii in 31:40){   E[ii,2,1,ii - 27 ] <- 1 }
# estimate model
mod6 <- tam.mml.3pl(resp= dat , control= ctl.list ,   E=E , est.variance=FALSE  )
summary(mod6)

#*** Model 7: 2PL model with design matrix of slopes and slope constraints
Ngam <- dim(E)[4]   # number of gamma parameters
# define two constraint equations
gammaslope.constr.V <- matrix( 0 , nrow=Ngam , ncol=2 )
gammaslope.constr.c <- rep(0,2)
# set sum of first two xlambda entries to 1.8
gammaslope.constr.V[1:2,1] <- 1
gammaslope.constr.c[1] <- 1.8
# set sum of entries 4, 5 and 6 to 3
gammaslope.constr.V[4:6,2] <- 1
gammaslope.constr.c[2] <- 3
mod7 <- tam.mml.3pl(resp= dat, control= ctl.list ,  E=E ,  est.variance=FALSE , 
   gammaslope.constr.V=gammaslope.constr.V , gammaslope.constr.c=gammaslope.constr.c)
summary(mod7)

#**** Model 8: Located latent class Rasch model with estimated three skill points

# three classes of theta's are estimated
TP <- 3
theta.k <- diag(TP)
# because item difficulties are unrestricted, we define the sum of the estimated
# theta points equal to zero
Ngam <- TP  # estimate three gamma loading parameters which are discrete theta points
E <- array( 0 , dim=c(I,2,TP,Ngam) )
E[,2,1,1] <- E[,2,2,2] <- E[,2,3,3] <- 1
gammaslope.constr.V <- matrix( 1 , nrow=3 , ncol=1 )
gammaslope.constr.c <- c(0)
# initial gamma values
gammaslope <- c( -2 , 0 , 2 )
# estimate model
mod8 <- tam.mml.3pl(resp= dat, control= ctl.list ,  E=E ,  skillspace="discrete" ,
     theta.k=theta.k , gammaslope=gammaslope , gammaslope.constr.V=gammaslope.constr.V , 
     gammaslope.constr.c=gammaslope.constr.c )            
summary(mod8)

#############################################################################
# EXAMPLE 2: Polytomous data
#############################################################################

data( data.mg , package="CDM")
dat <- data.mg[1:1000, paste0("I",1:11)]
               
#*** Model 1: 1-dimensional 1PL estimation, normal skill distribution
mod1 <- tam.mml.3pl(resp=dat , skillspace="normal" , 
         control= list( maxiter= 30 )  , gammaslope.des = "2PL" ,  
         est.some.slopes=FALSE ,  est.variance=TRUE  )
summary(mod1)

#*** Model 2: 1-dimensional 2PL estimation, discrete skill distribution
# define skill space
theta.k <- matrix( seq(-5,5,len=21) )
# allow skew skill distribution
delta.designmatrix <- cbind( 1 , theta.k , theta.k^2 , theta.k^3 )
# fix 13th xsi item parameter to zero
xsi.fixed <- cbind( 13 , 0 )
# fix 10th slope paremeter to one
gammaslope.fixed <- cbind( 10 , 1 )
# estimate model
mod2 <- tam.mml.3pl(resp=dat , skillspace="discrete" , theta.k=theta.k ,
          delta.designmatrix=delta.designmatrix,  control= list(maxiter=30),
          gammaslope.des = "2PL" , xsi.fixed=xsi.fixed , 
          gammaslope.fixed=gammaslope.fixed)
summary(mod2)

#*** Model 3: 2-dimensional 2PL estimation, normal skill distribution

# define loading matrix
Q <- matrix(0,11,2)
Q[1:6,1] <- 1   # items 1 to 6 load on dimension 1
Q[7:11,2] <- 1  # items 7 to 11 load on dimension 2
# estimate model
mod3 <- tam.mml.3pl(resp=dat , control= list(maxiter=30),
          gammaslope.des = "2PL" , Q=Q )                                                           
summary(mod3)

#############################################################################
# SIMULATED EXAMPLE 3: dichotomous data with guessing
#############################################################################

#*** simulate data
set.seed(9765)
N <- 4000   # number of persons
I <- 20     # number of items
b <- seq( -1.5 , 1.5 , len=I )
guess <- rep(0 , I )
guess.items <- c(6,11,16)
guess[ guess.items ] <- .33
library(sirt)
dat <- sirt::sim.raschtype( rnorm(N) , b=b , fixed.c=guess )

#*** Model 1: Difficulty + guessing model, i.e. fix slopes to 1
est.guess <- rep(0,I)
est.guess[ guess.items ] <- seq(1, length(guess.items) )
# define prior distribution
guess.prior <- matrix( cbind( 5 , 17 ) , I , 2 , byrow=TRUE )
guess.prior[ ! est.guess , ] <- 0
# estimate model
mod1 <- tam.mml.3pl(resp= dat , guess = guess , est.guess = est.guess , 
        guess.prior=guess.prior , control= ctl.list ,est.variance=TRUE  ,
        est.some.slopes=FALSE )          
summary(mod1)

#*** Model 2: estimate a joint guessing parameter
est.guess <- rep(0,I)
est.guess[ guess.items ] <- 1
# estimate model
mod2 <- tam.mml.3pl(resp= dat , guess = guess , est.guess = est.guess , 
        guess.prior=guess.prior , control= ctl.list ,est.variance=TRUE  ,
        est.some.slopes=FALSE )          
summary(mod2)

#############################################################################
# SIMULATED EXAMPLE 4: Latent class model with two classes
#      See slca Simulated Example 2 in the CDM package
#############################################################################

#*** simulate data
set.seed(9876)
I <- 7  # number of items
# simulate response probabilities
a1 <- round(runif(I,0 , .4 ),4)
a2 <- round(runif(I , .6 , 1 ),4)
N <- 1000   # sample size
# simulate data in two classes of proportions .3 and .7
N1 <- round(.3*N)
dat1 <- 1 * ( matrix(a1,N1,I,byrow=TRUE) > matrix( runif( N1 * I) , N1 , I ) )
N2 <- round(.7*N)
dat2 <- 1 * ( matrix(a2,N2,I,byrow=TRUE) > matrix( runif( N2 * I) , N2 , I ) )
dat <- rbind( dat1 , dat2 )
colnames(dat) <- paste0("I" , 1:I)

# define design matrices
TP <- 2     # two classes
theta.k <- diag(TP)     # there are theta dimensions -> two classes
# The idea is that latent classes refer to two different "dimensions".
# Items load on latent class indicators 1 and 2, see below.
E <- array(0 , dim=c(I,2,2,2*I) )
items <- colnames(dat)
dimnames(E)[[4]] <- c(paste0( colnames(dat) , "Class" , 1),
          paste0( colnames(dat) , "Class" , 2) )
# items, categories , classes , parameters
# probabilities for correct solution
for (ii in 1:I){
    E[ ii , 2 , 1 , ii ] <- 1    # probabilities class 1
    E[ ii , 2 , 2 , ii+I ] <- 1  # probabilities class 2
                    }
# define A matrix for xsi parameters which does not contain any parameter contributions
A <- array( 0 , dim=c(I,2,2) )
xsi.fixed <- cbind( cbind(1,2) , 0 )                  
# estimate model
mod1 <- tam.mml.3pl(resp= dat , E=E , A=A , control= list(maxit=20) , skillspace="discrete" , 
          theta.k=theta.k )          
summary(mod1)         
# compare simulated and estimated data
cbind( mod1$rprobs[,2,1] , a2 )  # Simulated class 2
cbind( mod1$rprobs[,2,2] , a1 )  # Simulated class 1

#############################################################################
# SIMULATED EXAMPLE 5: Located latent class model, Rasch model
#      See slca Simulated Example 4 in the CDM package
#############################################################################

#*** simulate data
set.seed(487)
I <- 15  # I items
b1 <- seq( -2 , 2 , len=I)      # item difficulties
N <- 2000    # number of persons
# simulate 4 theta classes
theta0 <- c( -2.5 , -1 , 0.3 , 1.3 )  # skill classes
probs0 <- c( .1 , .4 , .2 , .3 )      # skill class probabilities
TP <- length(theta0)
theta <- theta0[ rep(1:TP, round(probs0*N)  ) ]
library(sirt)
dat <- sirt::sim.raschtype( theta , b1 )

#*** Model 1: Located latent class model with 4 classes
maxK <- 2
Ngam <- TP
E <- array( 0 , dim=c(I , maxK , TP ,  Ngam ) )
dimnames(E)[[1]] <- colnames(dat)
dimnames(E)[[2]] <- paste0("Cat", 1:(maxK) )
dimnames(E)[[3]] <- paste0("Class", 1:TP)
dimnames(E)[[4]] <- paste0("theta", 1:TP) 
# theta design
for (tt in 1:TP){   E[1:I , 2 , tt ,  tt] <- 1       }
theta.k <- diag(TP)
# set eighth item difficulty to zero
xsi.fixed <- cbind( 8 , 0 )
# initial gamma parameter
gammaslope <- seq( -1.5 , 1.5 , len=TP)
# estimate model
mod1 <- tam.mml.3pl(resp= dat , E=E , xsi.fixed=xsi.fixed , 
           control= list(maxiter=100) , skillspace="discrete" , 
           theta.k=theta.k , gammaslope=gammaslope)          
summary(mod1)
# compare estimated and simulated theta class locations
cbind( mod1$gammaslope , theta0 )
# compare estimated and simulated latent class proportions
cbind( mod1$pi.k , probs0 )

#############################################################################
# SIMULATED EXAMPLE 6: DINA model with two skills
#      See slca Simulated Example 5 in the CDM package
#############################################################################

#*** simulate data
set.seed(487)
N <- 3000   # number of persons
# define Q-matrix
I <- 9  # 9 items
NS <- 2 # 2 skills 
TP <- 4 # number of skill classes
Q <- scan(nlines=3, text=
  "1 0   1 0   1 0
   0 1   0 1   0 1
   1 1   1 1   1 1",
   )
Q <- matrix(Q , I, ncol=NS, byrow=TRUE)
# define skill distribution
alpha0 <- matrix( c(0,0,1,0,0,1,1,1) , nrow=4,ncol=2,byrow=TRUE)
prob0 <- c( .2 , .4 , .1 , .3 )
alpha <- alpha0[ rep( 1:TP , prob0*N) ,]
# define guessing and slipping parameters
guess <- round( runif(I, 0 , .4 ) , 2 )
slip <- round( runif(I , 0 , .3 ) , 2 )
# simulate data according to the DINA model
dat <- CDM::sim.din( q.matrix=Q , alpha=alpha , slip=slip , guess=guess )$dat

#*** Model 1: Estimate DINA model
# define E matrix which contains the anti-slipping parameters
maxK <- 2
Ngam <- I
E <- array( 0 , dim=c(I , maxK , TP ,  Ngam ) )
dimnames(E)[[1]] <- colnames(dat)
dimnames(E)[[2]] <- paste0("Cat" , 1:(maxK) )
dimnames(E)[[3]] <- c("S00","S10","S01","S11")
dimnames(E)[[4]] <- paste0( "antislip" , 1:I ) 
# define anti-slipping parameters in E
for (ii in 1:I){  
        # define latent responses
        latresp <- 1*( alpha0 \%*\% Q[ii,]  == sum(Q[ii,]) )[,1]
        # model slipping parameters
        E[ii , 2 , latresp == 1, ii ] <- 1
                 }
# skill space definition
theta.k <- diag(TP)
gammaslope <- rep( qlogis( .8 ) , I ) 

# estimate model
mod1 <- tam.mml.3pl(resp= dat , E=E , control= list(maxiter=100) , skillspace="discrete" , 
          theta.k=theta.k , gammaslope=gammaslope)          
summary(mod1)
# compare estimated and simulated latent class proportions
cbind( mod1$pi.k , probs0 )
# compare estimated and simulated guessing parameters
cbind( mod1$rprobs[,2,1] , guess )
# compare estimated and simulated slipping parameters
cbind( 1 - mod1$rprobs[,2,4] , slip )

#############################################################################
# SIMULATED EXAMPLE 7: Mixed Rasch model with two classes
#      See slca Simulated Example 3 in the CDM package
#############################################################################

#*** simulate data
set.seed(987)
library(sirt)
# simulate two latent classes of Rasch populations
I <- 15  # 6 items
b1 <- seq( -1.5 , 1.5 , len=I)      # difficulties latent class 1
b2 <- b1    # difficulties latent class 2
b2[ c(4,7, 9 , 11 , 12, 13) ] <- c(1 , -.5 , -.5 , .33 , .33 , -.66 )
b2 <- b2 - mean(b2)
N <- 3000    # number of persons
wgt <- .25       # class probability for class 1
# class 1
dat1 <- sirt::sim.raschtype( rnorm( wgt*N ) , - b1 )
# class 2
dat2 <- sirt::sim.raschtype( rnorm( (1-wgt)*N , mean= 1 , sd=1.7) , - b2 )
dat <- rbind( dat1 , dat2 )

# The idea is that each grid point class x theta is defined as new
# dimension. If we approximate the trait distribution by 7 theta points
# and are interested in estimating 3 latent classes, then we need
# 7*3=21 dimensions.

#*** Model 1: Rasch model
# theta grid
theta.k1 <- seq( -3 , 3 , len=7 )
TT <- length(theta.k1)
#-- define theta design matrix
theta.k <- diag(TT)
#-- delta designmatrix
delta.designmatrix <- matrix( 0 , TT , ncol=3 )
delta.designmatrix[ , 1] <- 1
delta.designmatrix[, 2:3] <- cbind( theta.k1 , theta.k1^2 )

#-- define loading matrix E
E <- array( 0 , dim=c(I,2,TT,I + 1) )  # last parameter is constant 1
for (ii in 1:I){
    E[ ii , 2 , 1:TT , ii ] <- -1   # '-b' in '1*theta - b'
    E[ ii , 2 , 1:TT , I+1] <- theta.k1  # '1*theta' in '1*theta - b'    
                }                               
# initial gammaslope parameters
par1 <- qlogis( colMeans( dat ) )
gammaslope <- c( par1 , 1 )
# sum constraint of zero on item difficulties
gammaslope.constr.V <- matrix( 0 , I+1 , 1 )
gammaslope.constr.V[ 1:I , 1] <- 1  # Class 1
gammaslope.constr.c <- c(0)
# fixed gammaslope parameter
gammaslope.fixed <- cbind( I+1 , 1 )
# estimate model
mod1 <- tam.mml.3pl(resp= dat1 , E=E , control= list(maxiter=10) , skillspace="discrete" , 
      theta.k=theta.k , delta.designmatrix=delta.designmatrix ,
      gammaslope=gammaslope , gammaslope.constr.V=gammaslope.constr.V ,
      gammaslope.constr.c=gammaslope.constr.c , gammaslope.fixed=gammaslope.fixed ,
      notA=TRUE , est.variance=FALSE)                                            
summary(mod1)

#*** Model 2: Mixed Rasch model with two latent classes
# theta grid
theta.k1 <- seq( -4 , 4 , len=7 )
TT <- length(theta.k1)
#-- define theta design matrix
theta.k <- diag(2*TT)   # 2*7=14 classes
#-- delta designmatrix
delta.designmatrix <- matrix( 0 , 2*TT , ncol=6 )
# Class 1
delta.designmatrix[1:TT , 1] <- 1
delta.designmatrix[1:TT , 2:3] <- cbind( theta.k1 , theta.k1^2 )
# Class 2
delta.designmatrix[TT+1:TT , 4] <- 1
delta.designmatrix[TT+1:TT , 5:6] <- cbind( theta.k1 , theta.k1^2 )

#-- define loading matrix E
E <- array( 0 , dim=c(I,2,2*TT,2*I + 1) )  # last parameter is constant 1
dimnames(E)[[1]] <- colnames(dat)
dimnames(E)[[2]] <- c("Cat0","Cat1")
dimnames(E)[[3]] <- c( paste0("Class1_theta" , 1:TT) , paste0("Class1_theta" , 1:TT) )
dimnames(E)[[4]] <- c( paste0("b_Class1_" , colnames(dat)) , 
       paste0("b_Class2_" , colnames(dat)) , "One")
for (ii in 1:I){
  # Class 1 item parameters
    E[ ii , 2 , 1:TT , ii ] <- -1   # '-b' in '1*theta - b'
    E[ ii , 2 , 1:TT , 2*I+1] <- theta.k1  # '1*theta' in '1*theta - b'
  # Class 2 item parameters
    E[ ii , 2 , TT + 1:TT , I + ii ] <- -1 
    E[ ii , 2 , TT + 1:TT , 2*I+1] <- theta.k1 
                }
# initial gammaslope parameters
par1 <- qlogis( colMeans( dat ) )
gammaslope <- c( par1 , par1 + runif(I, -2  ,2 ) , 1 )
# sum constraint of zero on item difficulties within a class
gammaslope.constr.V <- matrix( 0 , 2*I+1 , 2 )
gammaslope.constr.V[ 1:I , 1] <- 1  # Class 1
gammaslope.constr.V[ I + 1:I , 2] <- 1  # Class 2
gammaslope.constr.c <- c(0,0)
# fixed gammaslope parameter
gammaslope.fixed <- cbind( 2*I+1 , 1 )

# estimate model
mod1 <- tam.mml.3pl(resp= dat , E=E , control= list(maxiter=30) , skillspace="discrete" , 
      theta.k=theta.k , delta.designmatrix=delta.designmatrix ,
      gammaslope=gammaslope , gammaslope.constr.V=gammaslope.constr.V ,
      gammaslope.constr.c=gammaslope.constr.c , gammaslope.fixed=gammaslope.fixed ,
      notA=TRUE)          
summary(mod1)
# latent class proportions
aggregate( mod1$pi.k , list( rep(1:2, each=TT)) , sum )
# compare simulated and estimated item parameters
cbind( b1 , b2 , - mod1$gammaslope[1:I] , - mod1$gammaslope[I + 1:I ] ) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{3PL model}
\keyword{Structured item response model}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
