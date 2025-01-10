\name{HaloTLL}
\alias{HaloTLL}
\alias{HaloTLL.fun}

\title{
Halotime model with log-logistic distribution of base salt concentration (derived from Bakhshandeh et al., 2020)
}

\description{
This model relates the time-course of the proportion of germinated seeds to salt concentration (C) in the substrate. It is based on a log-logistic distribution of base salt concentration within the seed lot. The equation is:

\deqn{ P(t) = 1 - \frac{1} {1 + \exp \left\{ - \frac{1}{\sigma} \log \left[ \frac{\left(C + \theta_{halo} / t \right)}{C_b} \right] \right\} } }

This models describes the distribution of base salt concentration within the population, while the distribution of germination times is indirectly modelled, but it is not, in itself, logistic (you see that 't' is at the denominator).
The 'HaloTLL.fun()' is a generic function, which can be used for plotting or other applications, while the 'HaloTLL()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
This model is derived from the halo-thermal time model described in Bakhshandeh et al. (2020).
}

\usage{
HaloTLL()
HaloTLL.fun(time, SConc, thetaHalo, SConcb50, sigma)
}
\arguments{
The 'HaloTLL()' function has no arguments. The 'HaloTLL.fun()' function has the following arguments:
  \item{time}{time}
  \item{SConc}{salt concentration in the substrate}
  \item{thetaHalo}{halotime constant}
  \item{SConcb50}{median base salt concentration}
  \item{sigma}{shape parameter for the base salt concentration within the seed lot}
}
\value{
The 'HaloTLL.fun()' function returns the proportion of germinated seeds, for any given values of time and salt concentration in the substrate. The 'HaloTLL()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
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
modg <- drmte(count ~ timeBef + timeAf + SConc, fct=HaloTLL(),
              data=dataset)
summary(modg)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
