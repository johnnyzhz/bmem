\name{bmem.ci.norm}
\alias{bmem.ci.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence interval based on normal approximation
}
\description{
Confidence interval based on normal approximation
}
\usage{
bmem.ci.norm(par.boot, par0, cl = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par.boot}{A bootstrap object.}
  \item{par0}{Original estimate}
  \item{cl}{Confidence level. Default 0.95.}
}
\value{
Normal confidence intervals. The output includes - estimates, bootstrap standard errors, and confidence intervals.
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\seealso{
\code{\link{bmem.ci.bca}}, \code{\link{bmem.ci.p}}, \code{\link{bmem.ci.bc}}
}