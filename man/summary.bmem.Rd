\name{summary.bmem}
\alias{summary.bmem}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate bootstrap confidence intervals}
\description{Calculate bootstrap confidence intervals}
\usage{
\method{summary}{bmem}(object, ci='bc', cl=.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An output object from the function \code{\link{bmem}}}
  \item{ci}{\code{norm}: normal approximation CI, \code{perc}: percentile CI, \code{bc}: bias-corrected CI, \code{bca}: BCa}
  \item{cl}{Confidence level. Can be a vector.}
  \item{...}{other options can be used for the generic summary function.}
}
\details{
The other type of confidence intervals can be constructed from the output of the function \code{\link{bmem}}. Note if the BCa is required, the \code{ci='BCa'} should have been specified in the function \code{\link{bmem}}.
}
\value{
The on-screen output includes the parameter estimates, bootstrap standard errors, and CIs.
}
