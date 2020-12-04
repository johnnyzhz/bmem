\name{bmem.pair.boot}
\alias{bmem.pair.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap for pairwise deletion}
\description{Bootstrap for pairwise deletion}
\usage{
bmem.pair.boot(x, ram, boot = 1000, parallel=FALSE,ncore=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{boot}{Number of bootstraps. Default is 1000.}
  \item{parallel}{Whether to use parallel method to calculate.}
  \item{ncore}{Numbers of core for parallel method.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
