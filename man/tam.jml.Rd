\name{tam.jml}

\alias{tam.jml}
\alias{tam.jml2}
\alias{summary.tam.jml}

\title{
Joint Maximum Likelihood Estimation
}

\description{
This function estimate unidimensional
item response models with joint maximum likelihood (JML,
 see e.g. Linacre, 1994).
}

\usage{
tam.jml(resp, group = NULL, adj = .3, disattenuate = FALSE, bias = TRUE, 
    xsi.fixed = NULL, xsi.inits = NULL, A = NULL, B = NULL, Q = NULL, 
    ndim = 1, pweights = NULL, control = list())

tam.jml2(resp, group = NULL, adj = .3 , disattenuate = FALSE, bias = TRUE, 
    xsi.fixed = NULL, xsi.inits = NULL, A = NULL, B = NULL, Q = NULL, 
    ndim = 1, pweights = NULL, control = list())

\method{summary}{tam.jml}(object,file=NULL,\dots)   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
     A matrix of item responses. Missing responses must be declared
     as \code{NA}.
}
  \item{group}{
     An optional vector of group identifier
}
  \item{disattenuate}{
    An optional logical indicating whether the person parameters 
    should be disattenuated due to unreliability? 
    The disattenuation is conducted by applying the Kelley formula.
}
\item{adj}{Adjustment constant which is subtracted or added to extreme
scores (score of zero or maximum score. The default is a value of 0.3}
  \item{bias}{
     A logical which indicates if JML bias should be reduced
     by multiplying item parameters by the adjustment factor
     of \eqn{(I-1)/I}
}
  \item{xsi.fixed}{
     An optional matrix with two columns for fixing some of the
     basis parameters \eqn{\xi} of item intercepts.
     1st column: Index of \eqn{\xi} parameter, 2nd column:
     Fixed value of \eqn{\xi} parameter
}
  \item{xsi.inits}{
     An optional vector of initial \eqn{\xi} parameters. Note that
     all parameters must be specified and the vector is not of the
     same format as \code{xsi.fixed}.
}
  \item{A}{
     A design array \eqn{A} for item category intercepts.
     For item \eqn{i} and category \eqn{k}, the threshold is
     specified as \eqn{ \sum _j a_{ikj}   \xi_j}.
}
  \item{B}{
     A design array for scoring item category responses.
     Entries in \eqn{B} represent item loadings on
     abilities \eqn{\theta}.
}
  \item{Q}{
     A Q-matrix which defines loadings of items on dimensions.
}

  \item{ndim}{
   Number of dimensions in the model. The default is 1.
}
  \item{pweights}{
   An optional vector of person weights.
}
  \item{control}{
     A list of control arguments. See \code{\link{tam.mml}}
     for more details.
}
  \item{object}{
Object of class \code{tam.jml} (only for \code{summary.tam}
function)
}
\item{file}{
A file name in which the summary output will be written 
(only for \code{summary.tam.jml} function)
    }
\item{\dots}{
	Further arguments to be passed
	}   

}

\details{
The function \code{tam.jml2} is just a bit more efficient (for most
applications) implementation than \code{tam.jml}. 
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries
\item{item}{Data frame with item parameters}
\item{xsi}{Vector of item parameters \eqn{\xi}}
\item{errorP}{Standard error of item parameters \eqn{\xi}}
\item{theta}{MLE in final step}
\item{errorWLE}{Standard error of WLE}
\item{WLE}{WLE in last iteration}
\item{WLEreliability}{WLE reliability}
\item{PersonScores}{Scores for each person (sufficient statistic)}
\item{ItemScore}{Sufficient statistic for each item parameter}
\item{PersonMax}{Maximum person score}
\item{ItemMax}{Maximum item score}
%\item{outfitPerson}{Person outfit statistic}
%\item{outfitItem}{Item outfit statistic}
%\item{infitPerson}{Person infit statistic}
%\item{infitItem}{Item infit statistic}
%\item{outfitPerson_t}{The \eqn{t} value of person outfit}
%\item{outfitItem_t}{The \eqn{t} value of item outfit}
%\item{infitPerson_t}{The \eqn{t} value of person infit}
%\item{infitItem_t}{The \eqn{t} value of item infit}
\item{deviance}{Deviance}
\item{deviance.history}{Deviance history in iterations}
\item{resp}{Original data frame}
\item{resp.ind}{Response indicator matrix}
\item{group}{Vector of group identifiers (if provided as an argument)}
\item{pweights}{Vector of person weights}
\item{A}{Design matrix \eqn{A} of item intercepts}
\item{B}{Loading (or scoring) matrix \eqn{B}}
\item{nitems}{Number of items}
\item{maxK}{Maximum number of categories}
\item{nstud}{Number of persons in \code{resp}}
\item{resp.ind.list}{Like \code{resp.ind}, only in the format of a list}
\item{xsi.fixed}{Fixed \eqn{\xi} item parameters}
\item{control}{Control list}
\item{\dots}{}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\references{
Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}. 
Chicago: MESA Press.
}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For estimating the same class of models with marginal
maximum likelihood estimation see \code{\link{tam.mml}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Dichotomous data
#############################################################################
data(sim.rasch)
resp <- sim.rasch[1:700 , seq( 1 , 40 , len=10)  ]  # subsample
# estimate the Rasch model with JML (function 'tam.jml')
mod1a <- tam.jml(resp=resp) 
summary(mod1a)
itemfit <- tam.fit(mod1a)$fit.item

# the same model but using function 'tam.jml2'
mod1b <- tam.jml2(resp=resp)
summary(mod1b)
itemfit <- tam.fit(mod1b)$fit.item

# compare results with Rasch model estimated by MML
mod1c <- tam.mml(resp=resp )
# plot estimated parameters
plot( mod1b$xsi , mod1c$xsi$xsi , pch=16 ,
    xlab= expression( paste( xi[i] , " (JML)" )) ,
    ylab= expression( paste( xi[i] , " (MML)" )) ,
    main="Item Parameter Estimate Comparison")
lines( c(-5,5) , c(-5,5) , col="gray" )

# Now, the adjustment pf .05 instead of the default .3 is used.
mod1d <- tam.jml2(resp=resp , adj=.05)
# compare item parameters
round( rbind( mod1b$xsi , mod1d$xsi ) , 3 )
##          [,1]   [,2]   [,3]   [,4]   [,5]  [,6]  [,7]  [,8]  [,9] [,10]
##   [1,] -2.076 -1.743 -1.217 -0.733 -0.338 0.147 0.593 1.158 1.570 2.091
##   [2,] -2.105 -1.766 -1.233 -0.746 -0.349 0.139 0.587 1.156 1.574 2.108

# person parameters for persons with a score 0, 5 and 10
pers1 <- data.frame( "score_adj0.3"= mod1b$PersonScore , "theta_adj0.3"= mod1b$theta , 
           "score_adj0.05"= mod1d$PersonScore , "theta_adj0.05"= mod1d$theta  )
round( pers1[ c(698, 683, 608) , ] ,3  )
##       score_adj0.3 theta_adj0.3 score_adj0.05 theta_adj0.05
##   698          0.3       -4.404          0.05        -6.283
##   683          5.0       -0.070          5.00        -0.081
##   608          9.7        4.315          9.95         6.179

\dontrun{
#############################################################################
# EXAMPLE 2: Partial credit model
#############################################################################

data(data.gpcm)
# the same model but using function 'tam.jml2'
mod2 <- tam.jml2(resp=data.gpcm)
mod2$xsi     # extract item parameters
summary(mod2)

#############################################################################
# EXAMPLE 3: Facet model estimation using joint maximum likelihood
#            data.ex10; see also Example 10 in tam.mml
#############################################################################

data(data.ex10)
dat <- data.ex10
## > head(dat)
##  pid rater I0001 I0002 I0003 I0004 I0005
##    1     1     0     1     1     0     0
##    1     2     1     1     1     1     0
##    1     3     1     1     1     0     1
##    2     2     1     1     1     0     1
##    2     3     1     1     0     1     1

facets <- dat[ , "rater" , drop=FALSE ] # define facet (rater)
pid <- dat$pid      # define person identifier (a person occurs multiple times)
resp <- dat[ , -c(1:2) ]        # item response data
formulaA <- ~ item * rater      # formula

# use MML function only to restructure data and input obtained design matrices
# and processed response data to tam.jml (-> therefore use only 2 iterations) 
mod3a <- tam.mml.mfr( resp=resp , facets=facets , formulaA = formulaA , 
             pid=dat$pid ,  control=list(maxiter=2) )

# use modified response data mod3a$resp and design matrix mod3a$A
resp1 <- mod3a$resp
# JML with 'tam.jml'
mod3b <- tam.jml( resp=resp1 , A=mod3a$A , control=list(maxiter=200) )
# 'tam.jml2'
mod3c <- tam.jml2( resp=resp1 , A=mod3a$A , control=list(maxiter=200) )	
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Joint maximum likelihood estimation (JML)}
\keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
