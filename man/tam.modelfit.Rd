\name{tam.modelfit}
\alias{tam.modelfit}
\alias{summary.tam.modelfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Q3 and Adjusted Q3 Statistic and Further Model Fit Statistics
}
\description{
This function computes the Q3 statistic (Yen, 1984) and an
adjusted variant of it (see Details). In addition,
some effect sizes of model fit (\code{MADaQ3}, \eqn{MADRESIDCOV},
\eqn{SRMR}) are implemneted.
}


\usage{
tam.modelfit(tamobj, progress = TRUE)

\method{summary}{tam.modelfit}(object,...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tamobj}{
Object of class \code{tam}
}
  \item{progress}{
An optional logical indicating whether progress should
be displayed
}

\item{object}{Object of class \code{tam.modelfit}}
\item{\dots}{Further arguments to be passed}
}

\details{
For each item \eqn{i} and each person \eqn{n}, residuals
\eqn{e_{ni}=X_{ni}-E(X_{ni})} are computed. The expected
value \eqn{E(X_{ni})} is obtained by integrating the individual
posterior distribution.

The Q3 statistic of item pairs \eqn{i} and \eqn{j} is defined as the
correlation \eqn{Q3_{ij} = Cor( e_{ni} , e_{nj} )}.
It is known that under local independence,
the expected value is slightly smaller than zero. Therefore,
an adjusted Q3 statistic (\code{aQ3}; \eqn{aQ3_{ij}})
is computed by subtracting the average of all Q3 statistics from 
Q3. To control for multiple testing, 
a p value adjustment by the method of
Holm (\code{p.holm}) is employed (see Chen, de la Torre & Zhang, 2013).

An effect size of model fit (\code{MADaQ3}) is defined as
the average of absolute values of \eqn{aQ3} statistics. 

The SRMSR (standardized root mean square root of squared residuals,
Maydeu-Olivaras, 2013) is based on comparing residual correlations 
of item pairs
\deqn{ SRMSR = \sqrt{ \frac{1}{ J(J-1)/2 } \sum_{i < j} 
	( r_{ij} - \hat{r}_{ij} )^2 } }

The \eqn{MADRESIDCOV} statistic (McDonald & Mok, 1995) is based on comparing 
residual covariances of item pairs
\deqn{ MADRESIDCOV = \frac{1}{ J(J-1)/2 } \sum_{i < j} 
	| c_{ij} - \hat{c}_{ij} |  }	
This statistic is just multiplied by 100 in the output of this function.	
}


\value{
A list with following entries
  \item{stat.MADaQ3}{Global fit statistic \code{MADaQ3} and
  global model test with \code{p} value obtained
  by Holm adjustment}
  \item{chi2.stat}{Data frame with chi square tests of condional independence
     for every item pair (Chen & Thissen, 1997)}
  \item{fitstat}{Model fit statistics \eqn{100 \cdot MADRESIDCOV} and
      \eqn{SRMR}}
  \item{modelfit.test}{Test statistic of global fit based on multiple
  testing correction of \eqn{\chi^2} statistics
  		}
  \item{stat.itempair}{Q3 and adjusted Q3 statistic for all item pairs}
  \item{residuals}{Residuals}
  \item{Q3.matr}{Matrix of \eqn{Q_3} statistics}
  \item{aQ3.matr}{Matrix of adjusted \eqn{Q_3} statistics}  
}

%  \item{chisquare.itemfit}{Item-wise \eqn{\chi ^2} local independence statistics}

%	res <- list( "stat.MADaQ3"=stat.MADaQ3 , "chi2.stat"=chi2.stat ,
%		   "fitstat" = fitstat , "modelfit.test"= modelfit.test , 
%		   "stat.itempair"=dfr , "residuals" = residM )



\references{
Chen, J., de la Torre, J., & Zhang, Z. (2013). 
Relative and absolute fit evaluation in cognitive diagnosis modeling. 
\emph{Journal of Educational Measurement}, \bold{50}, 123-140.	

Chen, W., & Thissen, D. (1997). Local dependence indexes for item pairs 
using item response theory. \emph{Journal of Educational and Behavioral Statistics},
\bold{22}, 265-289.		

Maydeu-Olivares, A. (2013). Goodness-of-fit assessment of item response
theory models (with discussion). 
\emph{Measurement: Interdisciplinary Research and Perspectives},
\bold{11}, 71-137.

McDonald, R. P., & Mok, M. M.-C. (1995). Goodness of fit in item response models.
\emph{Multivariate Behavioral Research}, \bold{30}, 23-40.

Yen, W. M. (1984). Effects of local item dependence on the fit and equating 
performance of the three-parameter logistic model. 
\emph{Applied Psychological Measurement}, \bold{8}, 125-145.
}

%\author{
%%  ~~who you are~~
%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#############################################################################
# EXAMPLE 1: data.cqc01
#############################################################################

data(data.cqc01)
dat <- data.cqc01

#*** estimation
mod1 <- tam.mml( dat )
#*** model fit
res1 <- tam.modelfit( tamobj = mod1 )
summary(res1)
# display item pairs with five largest adjusted Q3 statistics
res1$stat.itempair[1:5,c("item1","item2","aQ3","p","p.holm")]

\dontrun{
#############################################################################
# SIMULATED EXAMPLE 2: Rasch model
#############################################################################

set.seed(8766)
N <- 1000    # number of persons
I <- 20      # number of items
# simulate responses
library(sirt)
dat <- sirt::sim.raschtype( rnorm(N) , b=seq(-1.5,1.5,len=I) )
#*** estimation
mod1 <- tam.mml( dat )
summary(dat)
#*** model fit
res1 <- tam.modelfit( tamobj = mod1)
summary(res1)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Q3 statistic}
\keyword{Model fit}
\keyword{summary}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
