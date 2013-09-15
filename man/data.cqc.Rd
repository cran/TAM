\name{data.cqc}


\alias{data.cqc01}
\alias{data.cqc02}

\docType{data}
\title{
More Datasets and Examples (Similar to ConQuest Examples)
}


\description{
Datasets and examples similar to the ones in the ConQuest manual
(Wu, Adams, Wilson, & Haldane, 2007).
}


\usage{
data(data.cqc01)
data(data.cqc02)
}


\format{

\itemize{
\item \code{data.cqc01} contains 512 persons on 
12 dichotomous items of following format

\code{'data.frame':   512 obs. of  12 variables:} \cr
\code{ $ BSMMA01: int  1 1 0 1 1 1 1 1 0 0 ...} \cr
\code{ $ BSMMA02: int  1 0 1 1 0 1 1 1 0 0 ...} \cr
\code{ $ BSMMA03: int  1 1 0 1 1 1 1 1 1 0 ...} \cr
\code{ [...]} \cr
\code{ $ BSMSA12: int  0 0 0 0 1 0 1 1 0 0 ...} \cr

\item \code{data.cqc02} contains 431 persons on 8 variables of following
format

\code{'data.frame':   431 obs. of  8 variables:} \cr
\code{ $ It1: int  1 1 2 0 2 1 2 2 2 1 ...} \cr
\code{ $ It2: int  3 0 1 2 2 3 2 2 1 1 ...} \cr
\code{ $ It3: int  1 1 1 0 1 1 0 0 1 0 ...} \cr
\code{ [...]} \cr
\code{ $ It8: int  3 1 0 0 3 1 3 0 3 0 ...} \cr

}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}



%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Wu, M. L., Adams, R. J., Wilson, M. R. & Haldane, S. (2007).
\emph{ACER ConQuest Version 2.0}. Mulgrave.
\url{https://shop.acer.edu.au/acer-shop/group/CON3}
}


\examples{
\dontrun{	
#############################################################################
# EXAMPLE 01: Rasch model data.cqc01
#############################################################################
data(data.cqc01)

#*** Model 01: Estimate Rasch model
mod01 <- tam.mml( dat )
summary( mod01 )

#############################################################################
# EXAMPLE 02: Partial credit model and rating scale model data.cqc02
#############################################################################
data(data.cqc02)

# Model 02a: Partial credit model in ConQuest parametrization 'item+item*step'
mod02a <- tam.mml( dat , irtmodel="PCM2" )
summary( mod02a  )	# summary
tam.fit( mod02a )	# item fit

# Model 02b: Rating scale model
mod02b <- tam.mml( dat , irtmodel="RSM" )
summary( mod02b )
}
}

\keyword{datasets}
