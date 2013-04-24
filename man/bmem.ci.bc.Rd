\name{bmem.ci.bc}
\alias{bmem.ci.bc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias-corrected confidence intervals
}
\description{
Bias-corrected confidence intervals
}
\usage{
bmem.ci.bc(par.boot, par0, cl=.95)
}
\arguments{
  \item{par.boot}{A bootstrap object.}
  \item{par0}{Original estimate}
  \item{cl}{Confidence level. Default 0.95.}
}

\value{
BC confidence intervals. The output includes - estimates, bootstrap standard errors, and confidence intervals.
}

\author{
Zhiyong Zhang and Lijuan Wang
}
\seealso{
\code{\link{bmem.ci.norm}}, \code{\link{bmem.ci.p}}, \code{\link{bmem.ci.bca}}
}
