\name{bmem.em.boot}
\alias{bmem.em.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap for EM
}
\description{
Bootstrap for EM
}
\usage{
bmem.em.boot(x, ram, indirect, v, robust = FALSE, 
             varphi = 0.1, st= "i", boot = 1000, 
             moment = FALSE, max_it = 500, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{robust}{Roubst method}
  \item{varphi}{Percent of data to be downweighted}
  \item{st}{Starting values}
  \item{boot}{Number of bootstraps. Default is 1000.}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{max_it}{Maximum number of iterations in EM}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
\details{
The indirect effect can be specified using equations such as \code{a*b}, \code{a*b+c}, and \code{a*b*c+d*e+f}. A vector of indirect effects can be used \code{indirect=c('a*b', 'a*b+c')}.
}
\value{
  \item{par.boot}{Parameter estimates from bootstrap samples}
  \item{par0}{Parameter estimates from the orignal samples}
}
\author{Zhiyong Zhang and Lijuan Wang}
