\name{bmem.ci.bc1}
\alias{bmem.ci.bc1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias-corrected confidence intervals (for a single variable)
}
\description{
Bias-corrected confidence intervals (for a single variable)
}
\usage{
bmem.ci.bc1(x, b, cl = 0.95)
}
\arguments{
  \item{x}{A vector from a bootstrap output.}
  \item{b}{Parameter estimate from the original sample}
  \item{cl}{Confidence level. Default 0.95.}
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\seealso{
\code{\link{bmem.ci.norm}}, \code{\link{bmem.ci.p}}, \code{\link{bmem.ci.bca}}
}
