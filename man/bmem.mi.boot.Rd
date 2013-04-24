\name{bmem.mi.boot}
\alias{bmem.mi.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap for multiple imputation
}
\description{
Bootstrap for multiple imputation
}
\usage{
bmem.mi.boot(x, ram, indirect, v, m = 10, boot = 1000, 
             moment = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{m}{Number of imputations}
  \item{boot}{Number of bootstraps. Default is 1000.}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
