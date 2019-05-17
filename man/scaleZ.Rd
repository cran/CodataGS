\name{scaleZ}
\alias{scaleZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scales the genotype matrix.
}
\description{
Scales the genotype matrix so that ZZ' gives the genomic relationship matrix.
}
\usage{
scaleZ(Z, freq1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Genotype matrix with element values 0, 1 and 2
}
  \item{freq1}{
Optional input parameter with allele frequencies. A vector of length equal to the number of columns in Z.
}
}
\value{
\item{Z}{Scaled genotype matrix}
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
sim.res <- scaleZ(Z)
}
