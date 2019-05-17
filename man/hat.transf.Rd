\name{hat.transf}
\alias{hat.transf}
\title{
Transforms hat values
}
\description{
Transforms hat values between the SNP-BLUP model and the gBLUP model.
}
\usage{
hat.transf(C22, transf, vc, k, N, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C22}{
Submatrix of the inverse of the LHS in the MME
}
  \item{transf}{
A transformation matrix.
}
  \item{vc}{
Genetic variance
}
  \item{k}{
Number of SNPs
}
  \item{N}{
Number of individuals
}
  \item{w}{
SNP weights
}
}
\value{
Transformed hat values
}
\author{
Lars Ronnegard
}
