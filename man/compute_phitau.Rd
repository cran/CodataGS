\name{compute_phitau}
\alias{compute_phitau}

\title{
Computes models for the variance components
}
\description{
This function computes the residual variance, the SNP variances and the linear predictor for the SNP variance model.
}
\usage{
compute_phitau(dev, hv, devu, hvu, X.rand.disp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dev}{
Deviance values
}
  \item{hv}{
Hat values for the observed response values
}
  \item{devu}{
Deviance values computed for the random effects
}
  \item{hvu}{
Hat values for the random effects
}
  \item{X.rand.disp}{
Design matrix used in the linear predictor for the SNP variance model.
}
}

\value{
\item{var.e }{Residual variance}
\item{phi }{Vector of SNP variances}
\item{coef }{Fitted coefficients for the linear predictor in the SNP variance model}
}

\author{
Lars Ronnegard
}

\examples{
set.seed(1234)
N <- 20 #Number of individuals
k <- 30 #Number of SNPs with all marker positions including a QTL
#Simulated deviances and hat values
dev <- rnorm(N)^2
hv <- runif(N, 0.1, 0.5)
devu <- rnorm(k)^2
hvu <- runif(k, 0.1, 0.85)
X.rand.disp <- matrix(1, k, 1)
sim.res <- compute_phitau(dev, hv, devu, hvu, X.rand.disp)
}
