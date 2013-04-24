\name{bmem.moments}
\alias{bmem.moments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the moments of a data set
}
\description{
Calculate the moments of a data set using either listwise deletion or pairwise deletion
}
\usage{
bmem.moments(x, type=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data set}
  \item{type}{How to deal with missing data. 0: listwise deletion; 1: pairwise deletion}
}