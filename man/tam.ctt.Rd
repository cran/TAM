\name{tam.ctt}
\alias{tam.ctt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classical Test Theory Based Statistics
}
\description{
This function computes some item statistics based on classical test
theory.
}
\usage{
tam.ctt(resp, wlescore=NULL, pvscores=NULL, group=NULL , progress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
A data frame with unscored or scored item responses
}
  \item{wlescore}{
A vector with person parameter estimates, e.g. weighted likelihood
estimates obtained from \code{tam.wle}
}
  \item{pvscores}{
A matrix with plausible values, e.g. obtained from \code{tam.pv}
}
\item{group}{Vector of group identifiers if descriptive statistics shall
be groupwise calculated}
  \item{progress}{
Display progress?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with following columns:	
\item{index}{Index variable in this data frame}
\item{group}{Group identifier}
\item{itemno}{Item number}
\item{item}{Item}
\item{N}{Number of students responding to this item}
\item{Categ}{Category label}
\item{AbsFreq}{Absolute frequency of category}
\item{RelFreq}{Relative frequency of category}
\item{rpb.WLE}{Point biserial correlation of an item category and the WLE}
\item{M.WLE}{Mean of the WLE of students in this item category}
\item{SD.WLE}{Standard deviation of the WLE of students in this item category}
\item{rpb.PV}{Point biserial correlation of an item category and the PV}
\item{M.PV}{Mean of the PV of students in this item category}
\item{SD.PV}{Standard deviation of the PV of students in this item category}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
For dichotomously scored data, \code{rpb.WLE} is the ordinary point biserial
correlation of an item and a test score (here the WLE).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#################################################
# Example data.mc
data(data.mc)
# estimate Rasch model for scored data.mc data
mod <- tam.mml( resp=data.mc$scored )
# estimate WLE
w1 <- tam.wle( mod )
# estimate plausible values
set.seed(789)
p1 <- tam.pv( mod , ntheta=500 , normal.approx=TRUE )$pv

# CTT results for raw data
stat1 <- tam.ctt( resp=data.mc$raw , wlescore=w1$theta , pvscores=p1[,-1] )

# CTT results for scored data
stat2 <- tam.ctt( resp=data.mc$scored , wlescore=w1$theta , pvscores=p1[,-1] )

# descriptive statistics for different groups
# define group identifier
group <- c( rep(1,70) , rep(2,73) )
stat3 <- tam.ctt( resp=data.mc$raw , wlescore=w1$theta , pvscores=p1[,-1] , group=group)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Classical test theory statistics}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
