\name{jackGroupSE}
\alias{jackGroupSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delete-a-Petri-dish jackknife
}
\description{
This function estimates standard errors for regression model parameters by using the fully-iterated delete-a-group jackknife (Yu and Peng, 2008). This is asymptotically equivalent to the cluster-robust 'sandwich' variance estimator (Lipsiz, 1994).
}
\usage{
jackGroupSE(obj, data, units)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a fitted model object}
  \item{data}{dataset containing the variables used in obj}
  \item{units}{a vector containing the code for Petri dishes}
}
\details{
Care should be taken to ensure that the assumptions for using a fully-iterated jackknife are valid.
}
\value{
This function returns a data.frame containing the estimated parameters, the SEs and the robust SEs
}
\references{
Yu, B., Peng, Y., 2008. Mixture cure models for multivariate survival data. Computational Statistics and Data Analysis 52, 1524–1532.
Lipsitz, S.R., Dear, K.B.G., Zhao, L., 1994. Jackknife Estimators of Variance for Parameter Estimates from Estimating Equations with Applications to Clustered Survival Data. Biometrics 50, 842–846.
}
\author{
Andrea Onofri
}
\examples{
data(rape)
\donttest{
modHTE <- drmte( nSeeds ~ timeBef + timeAf + Psi,
               data=rape, fct=HTE1())
robust <- jackGroupSE(modHTE, rape, rape$Dish)
robust
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed germination }% use one of  RShowDoc("KEYWORDS")
