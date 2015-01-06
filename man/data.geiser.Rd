\name{data.geiser}
\alias{data.geiser}
\docType{data}

\title{
Dataset from Geiser et al. (2006)
}

\description{
This is a subsample of the dataset used in Geiser et al. (2006) and
Geiser and Eid (2010).
}

\usage{
data(data.geiser)
}

\format{
A data frame with 519 observations on the following 24 variables
  
\code{'data.frame':   519 obs. of  24 variables:} \cr
\code{ $ mrt1 : num  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ mrt2 : num  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ mrt3 : num  0 0 0 0 0 0 0 0 1 0 ...} \cr
\code{ $ mrt4 : num  0 0 0 0 0 1 0 0 0 0 ...} \cr
\code{[...]} \cr
\code{ $ mrt23: num  0 0 0 0 0 0 0 1 0 0 ...} \cr
\code{ $ mrt24: num  0 0 0 0 0 0 0 0 0 0 ...} \cr
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}


\references{
Geiser, C., & Eid, M. (2010). Item-Response-Theorie.
In C. Wolf & H. Best (Hrsg.). \emph{Handbuch der sozialwissenschaftlichen
Datenanalyse} (S. 311-332). VS Verlag fuer Sozialwissenschaften. 	
	
Geiser, C., Lehmann, W., & Eid, M. (2006). Separating rotators from
nonrotators in the mental rotations test: A multigroup latent class analysis. 
\emph{Multivariate Behavioral Research}, \bold{41(3)}, 261-293.
}



\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Latent trait and latent class models (Geiser et al., 2006, MBR)
#############################################################################

data(data.geiser)
dat <- data.geiser

#**********************************************
# Model 1: Rasch model
tammodel <- "
  LAVAAN MODEL:
    F =~ 1*mrt1__mrt12
    F ~~ F
  ITEM TYPE:
    ALL(Rasch)
    "
mod1 <- tamaan( tammodel , dat)    
summary(mod1)

#**********************************************
# Model 2: 2PL model
tammodel <- "
  LAVAAN MODEL:
    F =~ mrt1__mrt12
    F ~~ 1*F
    "
mod2 <- tamaan( tammodel , dat)    
summary(mod2)

# model comparison Rasch vs. 2PL
anova(mod1,mod2)

#*********************************************************************
#*** Model 3: Latent class analysis with four classes  

tammodel <- "
ANALYSIS:
  TYPE=LCA;
  NCLASSES(4);   # 4 classes
  NSTARTS(10,20); # 10 random starts with 20 iterations
LAVAAN MODEL:
  F =~ mrt1__mrt12
    "    
mod3 <- tamaan( tammodel , resp=dat  )
summary(mod3)

# extract item response functions
imod2 <- IRT.irfprob(mod3)[,2,]
# plot class specific probabilities
matplot( imod2 , type="o" , pch=1:4 , xlab="Item" , ylab="Probability" )
legend( 10 ,1 , paste0("Class",1:4) , lty=1:4 , col=1:4 , pch=1:4 )

#*********************************************************************
#*** Model 4: Latent class analysis with five classes  

tammodel <- "
ANALYSIS:
  TYPE=LCA;
  NCLASSES(5); 
  NSTARTS(10,20);
LAVAAN MODEL:
  F =~ mrt1__mrt12
    "    
mod4 <- tamaan( tammodel , resp=dat  )
summary(mod4)

# compare different models
AIC(mod1); AIC(mod2); AIC(mod3); AIC(mod4)
BIC(mod1); BIC(mod2); BIC(mod3); BIC(mod4)

#############################################################################
# EXAMPLE 2: Rasch model and mixture Rasch model (Geiser & Eid, 2010)
#############################################################################

### TO BE INCLUDED
}
}
\keyword{datasets}
