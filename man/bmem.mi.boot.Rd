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
bmem.mi.boot(x, ram, v, m = 10, boot = 1000, parallel=FALSE, ncore=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{m}{Number of imputations}
  \item{boot}{Number of bootstraps. Default is 1000.}
  \item{parallel}{Whether to use parallel method to calculate.}
  \item{ncore}{Numbers of core for parallel method.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
