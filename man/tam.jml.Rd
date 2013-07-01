\name{tam.jml}

\alias{tam.jml}
\alias{tam.jml2}

\title{
Joint Maximum Likelihood Estimation
}

\description{
This function estimate unidimensional
item response models with joint maximum likelihood (JML,
 see e.g. Linacre, 1994).
}

\usage{
tam.jml(resp, group = NULL, disattenuate = FALSE, bias = TRUE, 
    xsi.fixed = NULL, xsi.inits = NULL, A = NULL, B = NULL, Q = NULL, 
    R = NULL, ndim = 1, pweights = NULL, control = list())

tam.jml2(resp, group = NULL, disattenuate = FALSE, bias = TRUE, 
    xsi.fixed = NULL, xsi.inits = NULL, A = NULL, B = NULL, Q = NULL, 
    R = NULL, ndim = 1, pweights = NULL, control = list())
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
    Should the person parameters disattenuated due to
    unreliability? The disattenuation is conducted by applying
    the Kelley formula.
}
  \item{bias}{
     A logical which indicates if JML bias shold be reduced
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
     An optional vector of initial \eqn{\xi} parameters
}
  \item{A}{
     A design array \eqn{A} for item category intercepts.
     For item \eqn{i} and category \eqn{k}, the threshold is
     specified is \eqn{ \sum _j a_{ikj}   \xi_j}.
}
  \item{B}{
     A design array for scoring item category responses.
     Entries in \eqn{B} represent item loadings on
     abilities \eqn{\theta}.
}
  \item{Q}{
     A \eqn{Q} matrix which defines loadings of items on dimensions.
}
  \item{R}{
     This argument is not yet used.
}
  \item{ndim}{
   Number of dimensions in the model. The default is 1.
}
  \item{pweights}{
   An optional vector of person weights.
}
  \item{control}{
     A list of control arguments. See \code{\link{tam}}
     for more details.
}
}

\details{
The function \code{tam.jml2} is just a faster
implementation of \code{tam.jml}. It is in general
recommended to use \code{tam.jml2}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries
\item{xsi}{Vector of item parameters \eqn{\xi}}
\item{errorP}{Standard error of item parameters \eqn{\xi}}
\item{theta}{WLE in final step}
\item{errorWLE}{Standard error of WLE}
\item{WLE}{WLE in last iteration}
\item{WLEreliability}{WLE reliability}
\item{PersonScores}{Scores for each person (sufficient statistic)}
\item{ItemScore}{Sufficient statistic for each item parameter}
\item{PersonMax}{Maximum person score}
\item{ItemMax}{Maximum item score}
\item{outfitPerson}{Person outfit statistic}
\item{outfitItem}{Item outfit statistic}
\item{infitPerson}{Person infit statistic}
\item{infitItem}{Item infit statistic}
\item{outfitPerson_t}{The \eqn{t} value of person outfit}
\item{outfitItem_t}{The \eqn{t} value of item outfit}
\item{infitPerson_t}{The \eqn{t} value of person infit}
\item{infitItem_t}{The \eqn{t} value of item infit}
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
maximum likelihood estimation see \code{\link{tam}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
data(sim.rasch)
# estimate the Rasch model with JML
mod <- tam.jml2(resp=sim.rasch) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Joint maximum likelihood estimation (JML)}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
