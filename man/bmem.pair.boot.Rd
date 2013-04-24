\name{bmem.pair.boot}
\alias{bmem.pair.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap for pairwise deletion}
\description{Bootstrap for pairwise deletion}
\usage{
bmem.pair.boot(x, ram, indirect, boot = 1000, moment = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{boot}{Number of bootstraps. Default is 1000.}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
