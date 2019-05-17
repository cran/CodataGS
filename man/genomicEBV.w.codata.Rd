\name{genomicEBV.w.codata}
\alias{genomicEBV.w.codata}
\title{
Performs genomic prediction based on SNP codata. 
}
\description{
The main function of the package. The input includes response values, a design matrix for the fixed effects, a matrix with SNP genotype data and a design matrix for the SNP codata.
}
\usage{
genomicEBV.w.codata(y, X, Z, X.SNPcodata, Z.test = NULL, max.iter = 100, conv.crit = 1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Response values
}
  \item{X}{
Design matrix for the fixed effects
}
  \item{Z}{
Genotype matrix with element values of 0, 1 or 2
}
  \item{X.SNPcodata}{
Design matrix for the linear predictor of the SNP variances.
}
  \item{Z.test}{
An optional genotype matrix for a test data set.
}
  \item{max.iter}{
The maximum number of iterations
}
  \item{conv.crit}{
The value of the convergence criterion.
}
}
\details{
By specifying the matrix \code{Z.test} in the input, the function computes predicted genomic breeding values for an out-of-sample data set.
}
\value{
\item{gEBV}{Genomic breeding values}
\item{predicted.gEBV}{Genomic breeding values based on the genotypes in \code{Z.test} }
\item{w}{Computed SNP weights}
\item{u}{Fitted SNP effects}
\item{beta}{Fitted fixed effects}
\item{disp.beta}{Fitted coefficients in the linear predictor for the SNP variance model}
\item{Converge}{Shows whether the algorithm has converged or not}
\item{iter}{The number of iterations used}
}
\author{
Lars Ronnegard
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
#To fit gBLUP just specify X.SNPcodata = matrix(1, k, 1)
cat("Correlation between true and estimated BV for the codata model:")
cat(cor(crossprod(t(Z),u.scaled), mod1$gEBV), "\n")
}
