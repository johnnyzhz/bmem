\name{bmem.bs}
\alias{bmem.bs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap but using the Bollen-Stine method
}
\description{
The same as \code{\link{bmem}} but using the Bollen-Stine method
}
\usage{
bmem.bs(x, ram, indirect, v, ci='bc', cl=.95, 
        boot=1000, max_it=500, ...)
}
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{ci}{\code{norm}: normal approximation CI, \code{perc}: percentile CI, \code{bc}: bias-corrected CI, \code{bca}: BCa}
  \item{cl}{Confidence level. Can be a vector.}
  \item{boot}{Number of bootstraps}
  \item{max_it}{Maximum number of iterations in EM}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}

\value{
The on-screen output includes the parameter estimates, bootstrap standard errors, and CIs.
}

\author{
Zhiyong Zhang and Lijuan Wang
}
\references{
Zhang, Z. & Wang, L. (2011) Four methods for mediation analysis with missing data.

Zhang, Z. (2011) Robust mediation analysis with missing data and auxiliary variables.}
\seealso{
\code{\link{bmem}}, \code{\link{bmem.sobel}}, \code{\link{bmem.plot}}
}
