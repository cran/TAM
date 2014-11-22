\name{lavaanify.IRT}
\alias{lavaanify.IRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Slight Extension of the \code{lavaan} Syntax, with Focus on Item Response Models 
}

\description{
This functions slightly extends the \code{lavaan}
syntax implemented in the \pkg{lavaan} package
(see \code{\link[lavaan:lavaanify]{lavaanify}}).
Guessing and slipping parameters can be specified 
by using the operators \code{?=g1} and \code{?=s1},
respectively. The operator \code{__} can be used for a convenient
specification for groups of items. For example, \code{I1__I5} refers
to items \code{I1,...,I5}. The operator \code{__} can also be used for
item labels (see Example 2).
}

\usage{
lavaanify.IRT(lavmodel, items=NULL, data=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lavmodel}{
A model in \code{lavaan} syntax plus the additional
operators \code{?=g1} and \code{?=s1}.
}
\item{items}{Optional vector of item names}
\item{data}{Optional data frame with item responses}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with following entries	
\item{lavpartable}{A \code{lavaan} parameter table}
\item{lavaan.syntax}{Processed syntax for \code{lavaan} package}
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[lavaan:lavaanify]{lavaanify}} (\pkg{lavaan})


See \code{\link[sirt:tam2mirt]{tam2mirt}} (\pkg{sirt}) for
converting objects of class \code{tam} into \code{mirt}
objects.

See \code{\link[sirt:lavaan2mirt]{lavaan2mirt}} (\pkg{sirt}) 
for estimating models in the \pkg{mirt} package using \code{lavaan} syntax.
}

\examples{
#############################################################################
# EXAMPLE 1: lavaan syntax with guessing and slipping parameters
#############################################################################

# define model in lavaan
lavmodel <- "
    F =~ A1+c*A2+A3+A4
    # define slipping parameters for A1 and A2
    A1 + A2 ?= s1
    # joint guessing parameter for A1 and A2
    A1+A2 ?= c1*g1
    A3 | 0.75*t1
    # fix guessing parameter to .25 and 
    # slipping parameter to .01 for item A3
    A3 ?= .25*g1+.01*s1
    A4 ?= c2*g1
    A1 | a*t1
    A2 | b*t1
      "

# process lavaan syntax
lavpartable <- lavaanify.IRT(lavmodel)$lavpartable
  ##     id lhs op rhs user group free ustart exo label eq.id unco
  ##  1   1   F =~  A1    1     1    1     NA   0           0    1
  ##  2   2   F =~  A2    1     1    2     NA   0     c     0    2
  ##  3   3   F =~  A3    1     1    3     NA   0           0    3
  ##  4   4   F =~  A4    1     1    4     NA   0           0    4
  ##  5   5  A3  |  t1    1     1    0   0.75   0           0    0
  ##  6   6  A1  |  t1    1     1    5     NA   0     a     0    5
  ##  7   7  A2  |  t1    1     1    6     NA   0     b     0    6
  ##  8   8  A1 ?=  s1    1     1    7     NA   0           0    7
  ##  9   9  A2 ?=  s1    1     1    8     NA   0           0    8
  ##  10 10  A1 ?=  g1    1     1    9     NA   0    c1     1    9
  ##  11 11  A2 ?=  g1    1     1    9     NA   0    c1     1   10
  ##  12 12  A3 ?=  g1    1     1    0   0.25   0           0    0
  ##  13 13  A3 ?=  s1    1     1    0   0.01   0           0    0
  ##  14 14  A4 ?=  g1    1     1   10     NA   0    c2     0   11

\dontrun{
#############################################################################
# EXAMPLE 2: Usage of "__" and "?=" operators
#############################################################################  

library(sirt)
data(data.read, package="sirt")
dat <- data.read
items <- colnames(dat)

lavmodel <- "
   F1 =~ A1+A2+ A3+lam4*A4
   # equal item loadings for items B1 to B4
   F2 =~ lam5*B1__B4
   # different labelled item loadings of items C1 to C4
   F3 =~ lam9__lam12*C1__C4
   # item intercepts
   B1__B2 | -0.5*t1
   B3__C1 | int6*t1
   # guessing parameters
   C1__C3 ?= g1
   C4 + B1__B3 ?= 0.2*g1
   # slipping parameters
   A1__B1 + B3__C2 ?= slip1*s1
   # residual variances
   B1__B3 ~~ errB*B1__B3 
   A2__A4 ~~ erra1__erra3*A2__A4
    " 
lav2 <- lavaanify.IRT( lavmodel , data=dat)
lav2$lavpartable
cat( lav2$lavaan.syntax )

#** simplified example
lavmodel <- "
   F1 =~ A1+lam4*A2+A3+lam4*A4
   F2 =~ lam5__lam8*B1__B4
   F1 ~~ F2
   F1 ~~ 1*F1
   F2 ~~ 1*F2
    "   
lav3 <- lavaanify.IRT( lavmodel , data=dat)
lav3$lavpartable
cat( lav3$lavaan.syntax )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{lavaan}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
