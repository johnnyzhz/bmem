\name{bmem.sem}
\alias{bmem.sem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a mediaiton model using SEM technique
}
\description{
Estimate a mediaiton model using SEM technique
}
\usage{
bmem.sem(x, ram, N, indirect, moment=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A covariance matrix}
  \item{ram}{A path diagram from \code{specify.model}}
  \item{N}{Sample size}
  \item{indirect}{A vector of indirect effects}
  \item{moment}{Whether mean strucuture is used. The default is FALSE}
  \item{\dots}{Options that can be supplied to function \code{sem}.}
}
\seealso{
\code{\link{bmem.list.cov}}, \code{\link{bmem.pair.cov}}, \code{\link{bmem.mi.cov}}, \code{\link{bmem.em.cov}}
}