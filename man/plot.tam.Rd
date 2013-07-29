\name{plot.tam}
\alias{plot.tam}
\alias{plot.tam.mml}
\alias{plot.tam.jml}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Expected Response Curves for Unidimensional Item Response Models
}
\description{
S3 plot method for objects of class \code{tam}, \code{tam.mml}
or \code{tam.mml}.
}

\usage{
\method{plot}{tam}(x, low = -3, high = 3, ngroups = 6, wle = NULL, 
     export = TRUE, observed = TRUE, overlay = FALSE, ask = FALSE, ...)
	
\method{plot}{tam.mml}(x, low = -3, high = 3, ngroups = 6, wle = NULL, 
     export = TRUE, observed = TRUE, overlay = FALSE, ask = FALSE, ...)
     
\method{plot}{tam.jml}(x, low = -3, high = 3, ngroups = 6, wle = NULL,      
     export = TRUE, observed = TRUE, overlay = FALSE, ask = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{
Object of class \code{tam}, \code{tam.mml}
or \code{tam.mml}.
}
  \item{low}{
Lowest \eqn{\theta} value to be displayed
}
  \item{high}{
Highest \eqn{\theta} value to be displayed
}
  \item{ngroups}{
Number of score groups to be displayed. The default are
six groups.
}
  \item{wle}{
Use WLE estimate for displaying observed scores.
}
  \item{export}{
Export all graphics separately in png format in
the subfolder \code{'Plots'} of the working directory?
}
  \item{observed}{
A logical which indicates if observed response curve should
be displayed?
}
  \item{overlay}{
Overlay expected score functions
}
  \item{ask}{
A logical which asks for changing the graphic from item to item.
The default is \code{FALSE}.
}
  \item{\dots}{
Further arguments to be passed
}
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
mod <- tam.mml(sim.rasch)
plot(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line