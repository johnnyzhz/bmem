\name{bmem.pair.jack}
\alias{bmem.pair.jack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jackknife for pairwise deletion}
\description{Jackknife for pairwise deletion}
\usage{
bmem.pair.jack(x, ram, indirect, moment = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}  
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
