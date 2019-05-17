\name{Transform}
\alias{Transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforms hat values
}
\description{
The function calls the hat.transf function.
}
\usage{
Transform(X, L, var.e, var.u, v, svdVec, svdD, wZt, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Design matrix for the fixed effects
}
  \item{L}{
Square root matrix of the genomic relationship matrix, G
}
  \item{var.e}{
Residual variance
}
  \item{var.u}{
Genetic variance
}
  \item{v}{
Random effects
}
  \item{svdVec}{
Vector from the Single Value Decomposition of G
}
  \item{svdD}{
Diagonal elements of the Single Value Decomposition of G
}
  \item{wZt}{
Weights times the transpose of the scaled genotype matrix
}
  \item{w}{
Fitted SNP weights
}
}
\value{
\item{u}{SNP effects}
\item{qu}{Hat values for the SNP effects}
}
\author{
Lars Ronnegard
}

