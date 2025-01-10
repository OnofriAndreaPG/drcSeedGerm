\name{HaloTnorm}
\alias{HaloTnorm}
\alias{HaloTnorm.fun}

\title{
Halotime model with normal distribution of base salt concentration (derived from Bakhshandeh et al., 2020)
}

\description{
This model relates the time-course of the proportion of germinated seeds to salt concentration (C) in the substrate. It is based on a normal distribution of base salt concentration within the seed lot. The equation is:

\deqn{ P(t) = 1 - \Phi \left\{ \frac{C + (\theta_{halo} / t) - C_b }{\sigma_{C_b}} \right\} }

where \eqn{\Phi} is a gaussian cumulative distribution. This models describes the distribution of base salt concentration within the population, while the distribution of germination times is indirectly modelled, but it is not, in itself, gaussian (you see that 't' is at the denominator).
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
  \item{SConc}{salt concentration in the substrate}
  \item{ThetaHalo}{halotime constant}
  \item{SConcb50}{median base salt concentration}
  \item{sigmaSConcb}{standard deviation for the base salt concentration within the seed lot}
}
\value{
The 'HaloTNorm.fun()' function returns the proportion of germinated seeds, for any given values of time and salt concentration in the substrate. The 'HaloTNorm()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Bakhshandeh, E., Bradford, K.J., Pirdashti, H., Vahabinia, F., Abdellaoui, R., 2020. A new halothermal time model describes seed germination responses to salinity across both sub- and supra-optimal temperatures. Acta Physiologia Plantarum 42, 137.
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
                      check.names = FALSE)
datasetOr[,4:9] <- lapply(datasetOr[,4:9], as.numeric)
dataset <- melt_te(datasetOr, count_cols = 4:9, treat_cols = 1:3,
            monitimes = c(3,4,5,7,10,14), n.subjects = 50)
dataset <- dataset[dataset$CV == "Saturnin",]

# Model fitting
modg <- drmte(count ~ timeBef + timeAf + SConc, fct=HaloTnorm(),
              data=dataset)
summary(modg)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
