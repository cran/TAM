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
     Should individual posterior distributions be approximated by a
     normal distribution? The default is \code{FALSE}. In this case,
     the number of ability nodes \code{ntheta} can be substantially
     smaller than 2000, say 200 or 500.
}
  \item{samp.regr}{
     Should regression coefficients be fixed in the imputation or
     also sampled from their posterior distribution? The default is \code{FALSE}.
}
  \item{np.adj}{
     This parameter defines the spread of the random theta values
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
%% ...
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
data(sim.rasch)
# estimate Rasch model
mod <- tam.mml(sim.rasch)
# draw 5 plausible values (without a normality
# assumption of the posterior) and 2000 ability nodes
pv1 <- tam.pv( mod , nplausible=5 , ntheta=2000 )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plausible value imputation}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
