\name{tam.threshold}
\alias{tam.threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Thurstonian Thresholds
}
\description{
This function estimates Thurstonian thresholds for item category
parameters of (generalized) partial credit models (see Details).
}
\usage{
tam.threshold(tamobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tamobj}{
Object of class \code{tam}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with Thurstonian thresholds. Rows correspond to items and
columns to item steps.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%% \pkg{TAM} authors
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##############################################
# Example 1: ordered data - Partial credit model
data( data.gpcm )

# Model 1: partial credit model
mod1 <- tam.mml( resp=data.gpcm ,control=list( maxiter=200) )
summary(mod1)
##   Item Parameters -A*Xsi
##        item   N     M AXsi_.Cat1 AXsi_.Cat2 AXsi_.Cat3 B.Cat1.Dim1 B.Cat2.Dim1 B.Cat3.Dim1
##   1 Comfort 392 0.880     -1.302      1.154      3.881           1           2           3
##   2    Work 392 1.278     -1.706     -0.847      0.833           1           2           3
##   3 Benefit 392 1.163     -1.233     -0.404      1.806           1           2           3

# Calculation of Thurstonian thresholds
tam.threshold(mod1)
##                Cat1      Cat2     Cat3
##   Comfort -1.325226 2.0717468 3.139801
##   Work    -1.777679 0.6459045 1.971222
##   Benefit -1.343536 0.7491760 2.403168
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Thurstonian thresholds}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
