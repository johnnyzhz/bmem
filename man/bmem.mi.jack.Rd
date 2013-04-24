\name{bmem.mi.jack}
\alias{bmem.mi.jack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jackknife for multiple imputation}
\description{Jackknife for multiple imputation}
\usage{
bmem.mi.jack(x, ram, indirect, v, m = 10, moment = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{m}{Number of imputations.}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
