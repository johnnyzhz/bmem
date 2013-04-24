\name{bmem.ci.bca}
\alias{bmem.ci.bca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias-corrected and accelerated confidence intervals
}
\description{
Bias-corrected and accelerated confidence intervals
}
\usage{
bmem.ci.bca(par.boot, par0, jack, cl = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par.boot}{A bootstrap object.}
  \item{par0}{Original estimate}
  \item{jack}{A Jackknife object.}
  \item{cl}{Confidence level. Default 0.95.}  
}
\value{
BCa confidence intervals. The output includes - estimates, bootstrap standard errors, and confidence intervals.
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\seealso{
\code{\link{bmem.ci.norm}}, \code{\link{bmem.ci.p}}, \code{\link{bmem.ci.bc}}, \code{\link{bmem.list.jack}}, \code{\link{bmem.pair.jack}}, \code{\link{bmem.mi.jack}}, \code{\link{bmem.em.jack}},
}
