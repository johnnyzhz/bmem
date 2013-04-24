\name{bmem.em}
\alias{bmem.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a mediation model based on EM covariance matrix
}
\description{
Estimate a mediation model based on EM covariance matrix
}
\usage{
bmem.em(x, ram, indirect, v, robust = FALSE, 
        varphi = 0.1, st= "i", moment = FALSE, 
        max_it = 500, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{robust}{Roubst method}
  \item{varphi}{Percent of data to be downweighted}
  \item{st}{Starting values}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{max_it}{Maximum number of iterations in EM}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
