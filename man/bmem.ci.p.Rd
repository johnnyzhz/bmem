\name{bmem.ci.p}
\alias{bmem.ci.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Percentile confidence interval
}
\description{
Percentile confidence interval
}
\usage{
bmem.ci.p(par.boot, par0, cl = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par.boot}{A bootstrap object.}
  \item{par0}{Original estimate}
  \item{cl}{Confidence level. Default 0.95.}
}
\value{
Percentile confidence intervals. The output includes - estimates, bootstrap standard errors, and confidence intervals.
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\seealso{
\code{\link{bmem.ci.bca}}, \code{\link{bmem.ci.norm}}, \code{\link{bmem.ci.bc}}
}