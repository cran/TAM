\name{TAM-internal}

%%----------------------------------------
 \alias{calc_exp_TK}
 \alias{calc_posterior.v2}
 \alias{calc_prob.v5}
 \alias{stud_prior.v2}

 \alias{mstep.regression}
 \alias{Mstep_slope.v2}

 \alias{resp.pattern3}
 \alias{theta.sq}
 \alias{theta.sq2}
 
 \alias{rowcumsums}
 \alias{rowMaxs}

 \alias{tam.jml.xsi}
 \alias{tam.jml.xsi2}

 \alias{add.lead}
%%----------------------------------------

\title{
Computation of intermediate information
}
\description{
Listing of all the intermediate calculation methods currently 
active in TAM. Not to be called by user.
}
\usage{
calc_exp_TK(rprobs, A, itemwt, p, ip, nitems, resp.ind)
calc_posterior.v2(rprobs, gwt, resp, nitems, resp.ind.list, 
    normalization = TRUE, thetasamp.density = NULL, snodes = 0, resp.ind=NULL,
    avoid.zerosum=FALSE , logprobs = FALSE)
calc_prob.v5(iIndex, A, AXsi, B, xsi, theta, nnodes, maxK, 
    recalc = TRUE)
stud_prior.v2(theta, Y, beta, variance, nstud, nnodes, ndim, YSD, unidim_simplify)
tam.jml.xsi(resp, resp.ind, A, B, nstud, nitems, maxK, convM, 
    ItemScore, theta, xsi, Msteps, pweightsM, est.xsi.index)
tam.jml.xsi2(resp, resp.ind, A, A.0, B, nstud, nitems, maxK, convM, 
    ItemScore, theta, xsi, Msteps, pweightsM, est.xsi.index, rp3, 
    rp3.sel, rp3.pweightsM)

mstep.regression(resp, hwt, resp.ind, pweights, pweightsM, Y, 
    theta, theta2, YYinv, ndim, nstud, beta.fixed, variance, Variance.fixed, 
    group, G, snodes = 0 , thetasamp.density = NULL , nomiss=FALSE)
Mstep_slope.v2(B_orig, B, B_obs, B.fixed, max.increment, nitems, 
    A, AXsi, xsi, theta, nnodes, maxK, itemwt, Msteps, ndim, 
    convM, irtmodel, progress, est.slopegroups, E, basispar,se.B,
    equal.categ)

resp.pattern3(x)
theta.sq(theta)     # R version
theta.sq2(theta)    # Rcpp version

rowcumsums(m1)
rowMaxs(mat, na.rm = FALSE)

add.lead(x, width=max(nchar(x)))
}

%%\arguments{
%%  \item{rprobs}{
%%     ~~Describe \code{rprobs} here~~
%%}
%%  \item{A}{
%%     ~~Describe \code{A} here~~
%%}
%%  \item{itemwt}{
%%     ~~Describe \code{itemwt} here~~
%%}
%%  \item{p}{
%%     ~~Describe \code{p} here~~
%%}
%%  \item{ip}{
%%     ~~Describe \code{ip} here~~
%%}
%%  \item{nitems}{
%%     ~~Describe \code{nitems} here~~
%%}
%%  \item{resp.ind}{
%%     ~~Describe \code{resp.ind} here~~
%%}
%%}

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

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
