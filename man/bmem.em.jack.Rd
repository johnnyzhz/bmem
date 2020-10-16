\name{bmem.em.jack}
\alias{bmem.em.jack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jackknife estimate using EM
}
\description{
Jackknife estimate using EM
}
\usage{
bmem.em.jack(x, ram, v, robust = FALSE, 
             varphi = 0.1, st= "i", 
             max_it = 500, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{robust}{Roubst method}
  \item{varphi}{Percent of data to be downweighted}
  \item{st}{Starting values}
  \item{max_it}{Maximum number of iterations in EM}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
