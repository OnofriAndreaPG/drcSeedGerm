\name{HaloTnorm}
\alias{HaloTnorm}
\alias{HaloTnorm.fun}

\title{
Halotime model with normal distribution of base salt concentration (derived from Bakhshandeh et al., 2020)
}

\description{
This model relates the time-course of the proportion of germinated seeds to salt concentration (C) in the substrate. It is based on a normal distribution of base salt concentration within the seed lot. The equation is:

\deqn{ P(t) = \Phi \left\{ \frac{C + (\theta_{halo} / t) - C_b }{\sigma_{C_b}} \right\} }

where \eqn{\Phi} is a gaussian cumulative distribution function for base water potential. This models describes the distribution of base salt concentration within the population, while the distribution of germination times is indirectly modelled, but it is not, in itself, gaussian (you see that 't' is at the denominator).
The 'HaloTNorm.fun()' is a generic function, which can be used for plotting or other applications, while the 'HaloTNorm()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
This model is derived from the halo-thermal time model described in Bakhshandeh et al. (2020).
}

\usage{
HaloTnorm()
HaloTnorm.fun(time, SConc, ThetaHalo, SConcb50, sigmaSConcb)
}
\arguments{
The 'HaloTNorm()' function has no arguments. The 'HaloTNorm.fun()' function has the following arguments:
  \item{time}{time}
  \item{Conc}{salt concentration in the substrate}
  \item{ThetaHalo}{halotime constant}
  \item{SConcb50}{median base water potential}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
}
\value{
The 'HTNorm.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate. The 'HTNorm()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Bradford, K.J., 2002. Applications of hydrothermal time to quantifying and modeling seed germination and dormancy. Weed Science 50, 248–260.
}
\author{
Andrea Onofri
}
\note{
This function is for use with the R function 'drmte()'
}

\examples{
# Read data
datasetOr <- read.csv("https://www.casaonofri.it/_datasets/Halotime4genotypes.csv",
                      check.names = F)
datasetOr[,4:9] <- lapply(datasetOr[,4:9], as.numeric)
dataset <- melt_te(datasetOr, count_cols = 4:9, treat_cols = 1:3,
            monitimes = c(3,4,5,7,10,14), n.subjects = 50)
head(dataset)
dataset <- dataset[dataset$CV == "Saturnin",]

# Model fitting
modg <- drmte(count ~ timeBef + timeAf + SConc, fct=HaloTnorm(),
              data=dataset)
summary(modg)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
