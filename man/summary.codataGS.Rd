\name{summary.CodataGS}
\alias{summary}
\alias{summary.CodataGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary method for CodataGS objects}
\description{
A summary method for the object class \code{CodataGS}
}
\usage{
\method{summary}{CodataGS}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A CodataGS object }
  \item{\dots}{ arguments not used }
}
\details{
  Provides a concise summary of CodataGS objects.
}
\examples{
#######
#Simulation part
set.seed(1234)
N <- 200 #Number of individuals
k <- 300 #Number of SNPs with all marker positions including a QTL
Z1 <- matrix(0, N, k )
Z2 <- matrix(0, N, k )
Z1[1:N, 1] <- rbinom(N, 1, 0.5) #Simulated phased SNP matrices
Z2[1:N, 1] <- rbinom(N, 1, 0.5) 
LD.par <- 0.2 #A parameter to simulate LD. 0 gives full LD, and 0.5 no LD
for (j in 2:k) {
  Z1[1:N, j] <- abs( Z1[1:N, j-1] - rbinom(N, 1, LD.par) )
  Z2[1:N, j] <- abs( Z2[1:N, j-1] - rbinom(N, 1, LD.par) )
} 
Z <- Z1 + Z2 #Genotypic SNP matrix
x1 <- c(rep(1,k/2), rep(0,k/2)) #An indicator for the SNPs. 
#The first k/2 SNPs and the last k/2 have different variances
#Simulate linear predictor for the random effect variance
lin.pred <- 0 + 2*x1 
X.snp <- model.matrix( ~ x1 ) #Corresponding design matrix
u <- rnorm(k, 0 , sqrt( exp(lin.pred) )) 
#Took the square root here because it is the SD that is specified.
#and exp() because we are modelling a log link.
u.scaled <- u/as.numeric( sqrt( var( crossprod(t(Z), u) )) ) 
#Scaled by the variance of the breeding values
e <- rnorm(N) #A residual variance
mu <- 0
y <- mu + crossprod(t(Z),u.scaled) + e
######
#Estimation part
mod1 <- genomicEBV.w.codata(y = as.numeric(y), 
          X = matrix(1, N, 1), Z = Z, X.SNPcodata = X.snp)
summary(mod1)
}
