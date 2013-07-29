\name{anova.tam}
\alias{anova.tam}
\alias{anova.tam.mml}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Test for Model Comparisons
}
\description{
This function compares two models estimated \code{tam}
or \code{tam.mml} using a likelihood ratio test.
}
\usage{
\method{anova}{tam}(object,object1,\dots)

\method{anova}{tam.mml}(object,object1,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{tam} or \code{tam.mml} 
}
  \item{object1}{
Object of class \code{tam} or \code{tam.mml} 
}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%% Alexander Robitzsch
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
# 1PL estimation
mod1 <- tam.mml(resp=sim.rasch) 
# 2PL estimation
mod2 <- tam.mml.2pl(resp=sim.rasch , irtmodel="2PL") 
# Model comparison
anova( mod1 , mod2 )
## > anova( mod1 , mod2 )
##     Model   loglike Deviance Npars      AIC      BIC    Chisq df       p
## 1 Model 1 -42077.88 84155.77    41 84237.77 84467.40 54.05078 39 0.05508
## 2 Model 2 -42050.86 84101.72    80 84261.72 84709.79       NA NA      NA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Likelihood ratio test}
\keyword{anova}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
