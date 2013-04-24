\name{bmem.sobel}
\alias{bmem.sobel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mediation analysis using sobel test (for complete data only)
}
\description{
Mediation analysis using sobel test (for complete data only)
}
\usage{
bmem.sobel(x, ram, indirect, moment=FALSE, ...)
}
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{moment}{Covariance or moment analysis}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}

\value{The on-screen output includes the parameter estimates and sobel standard errors.
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\references{
Zhang, Z. & Wang, L. (2011) Four methods for mediation analysis with missing data.

Zhang, Z. (2011) Robust mediation analysis with missing data and auxiliary variables.
}
\seealso{
\code{\link{bmem}}, \code{\link{bmem.sobel}}, \code{\link{bmem.plot}}
}
