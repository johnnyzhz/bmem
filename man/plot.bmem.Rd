\name{plot.bmem}
\alias{plot}
\alias{plot.bmem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the bootstrap distribution
}
\description{
Plot of the bootstrap distribution
}
\usage{
\method{plot}{bmem}(x, par, ...)
}
\arguments{
  \item{x}{A bmem object}
  \item{par}{Name of parameter to be plotted.}
  \item{...}{Options used for the generic plot function.}
}

\value{
Generate the bootstrap histogram for a chosen parameter.}
\author{
Zhiyong Zhang and Lijuan Wang
}
\references{
Zhang, Z. & Wang, L. (2011) Four methods for mediation analysis with missing data.

Zhang, Z. (2011) Robust mediation analysis with missing data and auxiliary variables.}
\seealso{
\code{\link{bmem}}, \code{\link{bmem.sobel}}, \code{\link{bmem.plot}}
}
