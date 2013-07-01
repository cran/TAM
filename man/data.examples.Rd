\name{data.examples}


\alias{data.ex08}
\alias{data.ex10}

\docType{data}

\title{
Data Examples in \pkg{TAM} Package
}
\description{
Data examples included in the \pkg{TAM} package
}
\usage{
data(data.ex08)
data(data.ex10)
}

\format{
Data \code{data.ex08} for Example 8 in \code{\link{tam}} 
has the following format: \cr
\code{List of 2} \cr
\code{ $ facets:'data.frame': 1000 obs. of  1 variable:} \cr
\code{  ..$ female: int [1:1000] 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ resp  : num [1:1000, 1:10] 1 1 1 0 1 0 1 1 0 1 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:10] "I0001" "I0002" "I0003" "I0004" ...} \cr

Data \code{data.ex10} for Example 10 in \code{\link{tam}} 
has the following format: \cr
\code{'data.frame':   675 obs. of  7 variables:} \cr
\code{ $ pid  : int  1 1 1 2 2 3 3 4 4 5 ...} \cr
\code{ $ rater: int  1 2 3 2 3 1 2 1 3 1 ...} \cr
\code{ $ I0001: num  0 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ I0002: num  1 1 1 1 1 0 1 1 1 1 ...} \cr
\code{ $ I0003: num  1 1 1 1 0 0 0 1 0 1 ...} \cr
\code{ $ I0004: num  0 1 0 0 1 0 1 0 1 0 ...} \cr
\code{ $ I0005: num  0 0 1 1 1 0 0 1 0 1 ...} \cr


}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%The dataset is copied from the \pkg{ltm} package.
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{
These examples are used in the \code{\link{tam}} Examples.
	}
%\examples{
%data(data.gpcm)
%summary(data.gpcm)
%}
\keyword{datasets}
