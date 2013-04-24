\name{bmem.list.boot}
\alias{bmem.list.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap for listwise deletion method
}
\description{
Bootstrap for listwise deletion method
}
\usage{
bmem.list.boot(x, ram, indirect, boot = 1000, moment = FALSE, ...)
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
