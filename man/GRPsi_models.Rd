\name{GR_Psi_models}
\alias{GRPsiLin.fun}
\alias{GRPsiLin}
\alias{GRPsiPol.fun}
\alias{GRPsiPol}
\alias{GRPsiPol2.fun}
\alias{GRPsiPol2}
\alias{GRPsi.Lin}
\alias{GRPsi.Pol}
\alias{GRPsi.Pol2}
\title{Hydro-time models based on the germination rate}
\description{
These models describe the response of germination rate to water potential in the substrate.
}
\usage{
GRPsiLin()
GRPsiPol()
GRPsiPol2()
GRPsi.Lin()
GRPsi.Pol()
GRPsi.Pol2()
GRPsiLin.fun(Psi, Psib, thetaH)
GRPsiPol.fun(Psi, Psib, thetaH)
GRPsiPol2.fun(Psi, Psib, thetaH)
}

\arguments{
The functions 'GRPsiLin()', 'GRPsiPol()', 'GRPsiPol2()', GRPsi.Lin(), GRPsi.Pol(), GRPsi.Pol2() have no arguments. The functions 'GRPsiLin.fun()', 'GRPsiPol.fun()' and 'GRPsiPol2.fun()' have the following arguments:
  \item{Psi}{Water potential in the substrate}
  \item{Psib}{Base water potential within the population}
  \item{thetaH}{Hydro-time parameter}
}

\details{
The functions 'GRPsiLin()', 'GRPsiPol()', 'GRPsiPol2()', GRPsi.Lin(), GRPsi.Pol(), GRPsi.Pol2() are meant to be used with the 'drm()' function in the 'drc' package ('GRPsiLin()' and 'GRPsi.Lin()', 'GRPsiPol()' and 'GRPsi.Pol()', 'GRPsiPol2()' and 'GRPsiPol2()' are totally equivalent, apart from the names). The functions 'GRPsiLin.fun()', 'GRPsiPol.fun()' and 'GRPsiPol2.fun()' are used for general purposes (plotting and other uses). Details about these functions and the meaning of parameters are described in Bradford (2002) and in the package ducumentation (see references below).
}

\value{
The 'GRPsiLin.fun()', 'GRPsiPol.fun()' and 'GRPsiPol2.fun()' functions return the germination rate for any given values of water potential in the substrate. The 'GRPsiLin()', 'GRPsiPol()' and 'GRPsiPol2()' (and 'GRPsi.Lin()', 'GRPsi.Pol()' and 'GRPsi.Pol2()') functions return a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drc()' function.
}

\references{
    See package documentation at: \url{https://www.statforbiology.com/_seedtutorial/}
}
\author{
Andrea Onofri
}
\examples{
library(drcte)
# Observed data
Psi <- c(-2, -1.5, -1.2, -1, -0.8, -0.6, -0.4, -0.25,
         -0.12, -0.06, -0.03, 0)
GR <- c(0, 0, 0, 0, 0.0585, 0.094, 0.1231, 0.1351,
        0.1418, 0.1453, 0.1458, 0.1459)
Psi2 <- c(-0.5, -0.6, -0.7, -0.8, -0.9, -1, -1.1, -1.2,
          -1.5)
GR2 <- c(1.4018, 1.0071, 0.5614, 0.3546, 0.2293, 0, 0,
         0, 0)

# Model fitting
modHT1 <- drm(GR ~ Psi, fct = GRPsiLin())
modHT2 <- drm(GR ~ Psi, fct = GRPsiPol())
modHT3 <- drm(GR2 ~ Psi2, fct = GRPsiPol2())

summary(modHT1)
summary(modHT2)
summary(modHT2)
}
