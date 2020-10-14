\name{bmem}
\alias{bmem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mediation analysis based on bootstrap}
\description{Mediation analysis based on bootstrap}
\usage{
bmem(x, ram, indirect, v, method='tsml', ci='bc', cl=.95, 
     boot=1000, m=10, varphi=.1, st='i', robust=FALSE, 
     max_it=500, moment=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{ram}{RAM path for the mediaiton model}
  \item{indirect}{A vector of indirect effec}
  \item{v}{Indices of variables used in the mediation model. If omitted, all variables are used.}
  \item{method}{\code{list}: listwise deletion, \code{pair}: pairwise deletion, \code{mi}: multiple imputation, \code{em}: EM algorithm.}
  \item{ci}{\code{norm}: normal approximation CI, \code{perc}: percentile CI, \code{bc}: bias-corrected CI, \code{bca}: BCa}
  \item{cl}{Confidence level. Can be a vector.}
  \item{boot}{Number of bootstraps}
  \item{m}{Number of imputations}
  \item{varphi}{Percent of data to be downweighted}
  \item{st}{Starting values}
  \item{robust}{Robust method}
  \item{moment}{Select mean structure or covariance analysis. moment=FALSE, covariance analysis. moment=TRUE, mean and covariance analysis.}
  \item{max_it}{Maximum number of iterations in EM}
  \item{\dots}{Other options for \code{\link{sem}} function can be used.}
}
\details{
The indirect effect can be specified using equations such as \code{a*b}, \code{a*b+c}, and \code{a*b*c+d*e+f}. A vector of indirect effects can be used \code{indirect=c('a*b', 'a*b+c')}.
}
\value{
The on-screen output includes the parameter estimates, bootstrap standard errors, and CIs.
}
\author{
Zhiyong Zhang and Lijuan Wang
}
\references{
Zhang, Z., & Wang, L. (2013). Methods for mediation analysis with missing data. Psychometrika, 78(1), 154-184.
}
