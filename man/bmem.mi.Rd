\name{bmem.mi}
\alias{bmem.mi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a mediation model based on multiple imputation
}
\description{
Estimate a mediation model based on multiple imputation
}
\usage{
bmem.mi(x, ram, v, m = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{m}{Number of imputations.}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
