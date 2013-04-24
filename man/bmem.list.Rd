\name{bmem.list}
\alias{bmem.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a mediaiton model based on listwise deletion
}
\description{
Estimate a mediaiton model based on listwise deletion
}
\usage{
bmem.list(x, ram, indirect, moment = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
