\name{GR_Psi_models}
\alias{GRPsiLin.fun}
\alias{GRPsiLin}
\alias{GRPsiPol.fun}
\alias{GRPsiPol}
\alias{GRPsiPol2.fun}
\alias{GRPsiPol2}
\title{Hydro-time models for seed germination}
\description{
These models describe the response of germination rate to water potential in the substrate.
}
\usage{
GRPsiLin()
GRPsiPol()
GRPsiPol2()
}

\arguments{
  These functions have no arguments
}

\details{
These functions are meant to be used with the 'drm()' function in the 'drc' package. The function 'GRPsiLin()' assumes that the germination rate increases linearly with water potential in the substrate. 'GRPsiPol()' and 'GRPsiPol2()' assume a curvilinear increase with concave and convex shapes, respectively.
}
\value{
The value returned is a list containing the nonlinear function, the self starter function and the parameter names.
}
\references{
\url{https://www.statforbiology.com/2021/stat_seedgermination_htt2step/#one-threshold-polynomial}
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
