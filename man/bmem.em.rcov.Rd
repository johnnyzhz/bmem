\name{bmem.em.rcov}
\alias{bmem.em.rcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of robust covariance matrix
}
\description{
Estimation of robust covariance matrix
}
\usage{
bmem.em.rcov(xmis, varphi=.1, moment=FALSE, max_it=1000, st='i')
}
\arguments{
  \item{xmis}{Missing data pattern}
  \item{varphi}{Percent of data to be downweighted}
  \item{moment}{Moment analysis if TRUE}
  \item{max_it}{Maximum number of iteration}
  \item{st}{Starting values}
}

\value{
An interval function to calculate the robust covaraince matrix}

\author{
Zhiyong Zhang and Lijuan Wang
}

