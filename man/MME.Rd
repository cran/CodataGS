\name{MME}
\alias{MME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixed model equations
}
\description{
A fast version of the Henderson's mixed model equations (MME)
}
\usage{
MME(y, X, Z, var.e, var.u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Response
}
  \item{X}{
Design matrix for fixed effects
}
  \item{Z}{
Design matrix for the random effects
}
  \item{var.e}{
Residual variance
}
  \item{var.u}{
Genetic variance
}
}
\value{
 \item{beta}{Estimates of fixed effects}
 \item{v}{Fitted random effects}
 \item{hv}{Hat values}
 \item{dev}{Deviances}
}
\author{
Lars Ronnegard
}
