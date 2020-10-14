\name{power.basic}
\alias{power.basic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conducting power analysis based on Sobel test
}
\description{
Different from \code{\link{power.boot}}, this function conduct power analysis based on the Sobel test.
}
\usage{
power.basic(model, indirect = NULL, nobs, nrep = 1000, alpha = 0.95, 
skewness = NULL, kurtosis = NULL, ovnames = NULL, se = "default", 
estimator = "default", parallel = "no", ncore = 1,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A model specified using lavaan notation and above. See \code{\link{model.syntax}} for basic model specification. 

For the power analysis, the population parameter values should be provided in the following way. For example, the coefficient between math and HE is .39. Then it is specified as start(.39). If the parameter will be referred in the mediation effect, a label should be given as a modifier as b*HE+start(.39)*HE.

model<-'
math ~ c*ME+start(0)*ME + b*HE+start(.39)*HE
HE ~ a*ME+start(.39)*ME
'

}
  \item{indirect}{
The indirect or other composite effects are specified in the following way

indirect<-'
ab: = a*b
abc := a*b + c
'
}
  \item{nobs}{
Number of observations for power analysis. If it is a vector, multiple group analysis will be conducted.
}
  \item{nrep}{
Number of replications for Monte Carlo simulation. At least 1,000 is recommended.
}
  \item{alpha}{
The alpha level is used to obtain the confidence interval for model parameters.
}
  \item{skewness}{
A vector to give the skewness for the observed variables.
}
  \item{kurtosis}{
A vector to give the kurtosis for the observed variables.
}
  \item{ovnames}{
A vector to give the variable names for the observed variables. This is only needed when the skewness and kurtosis are provided. The skewness, kurtosis and variable names should be in the same order.
}
  \item{se}{
How to calculate the standard error, for example, robust standard error can be specified using se="robust".
}
  \item{estimator}{
Estimation methods to be used here.
}
  \item{parallel}{
Parallel methods, snow or multicore, can be used here. 
}
  \item{ncore}{
Number of cores to be used in parallel. By defautl, the maximum number of cores are used.
}
  \item{\dots}{
Other named arguments for lavaan can be passed here.
}
}

\value{
\item{power}{power for all parameters and required ones in the model}
\item{coverage}{coverage probability}
\item{pop.value}{Population parameter values}
\item{results}{A list to give all intermediate results}
\item{data}{The last data set generated for checking purpose}
}

\examples{
\dontrun{
ex1model<-'
math ~ c*ME+start(0)*ME + b*HE+start(0.39)*HE
HE ~ a*ME+start(0.39)*ME
'

indirect<-'ab:=a*b'

N<-50

system.time(non.normal<-power.basic(ex1model, indirect, N, 
 nrep=2000,parallel='multicore', skewness=c(-.3, -.7, 1.3), 
 kurtosis=c(1.5, 0, 5), ovnames=c('ME', 'HE', 'math'), ncore=8))
 
summary(non.normal)

}
}
