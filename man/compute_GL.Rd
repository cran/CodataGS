\name{compute_GL}
\alias{compute_GL}
\title{
Computes genomic relationship matrix
}
\description{
This function computes the genomic relationship matrix, G, together with its matrix square root, L.
}
\usage{
compute_GL(Z, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Scaled matrix with genotype information
}
  \item{w}{
weights
}
}

\value{

  \item{L}{Square root matrix of G }
  \item{svdVec}{Vectors in the Single Value Decomposition of G}
  \item{svdD}{Diagonal elements in the Single Value Decomposition of G}
  \item{wZt}{weights times the transpose of Z}
}

\author{
Lars Ronnegard
}
\examples{
set.seed(1234)
N <- 20 #Number of individuals
k <- 30 #Number of SNPs with all marker positions including a QTL
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
sim.res <- compute_GL(Z, w = rep(1,k))
}