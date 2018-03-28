\name{lotusCum}
\alias{lotusCum}
\docType{data}
\title{
Germination assay with Lotus ornithopodioides
}
\description{
Seeds of L. ornithopodioides were harvested in three stages after anthesis and the time course of their germination was evaluated in Petri dishes. This is the dataset as required for nonlinear regression analysis and needs to be transformed into a form that is suitable for time-to-event analysis in R.
}
\usage{data("lotusCum")}
\format{
  A data frame with 180 observations on the following 6 variables.
  \describe{
    \item{\code{Stage}}{a factor with levels \code{A} \code{B} \code{C}}
    \item{\code{Dish}}{a numeric vector: code for Petri dishes}
    \item{\code{Time}}{a numeric vector: inspection time in days}
    \item{\code{nSeeds}}{a numeric vector: number of germinated seeds at each time}
    \item{\code{nCum}}{a numeric vector: cumulative number of germinated seeds at each time}
    \item{\code{Prop}}{a numeric vector: cumulative proportion of germinated seeds at each time}
  }
}
\details{
There were 25 seeds per Petri dish. The columns represent the characteristics of each Petri dish. The columns from 1 to 15 represent the number of germinated seeds counted in each assessment time. Seeds of Lotus ornithopodioides were harvested at the following stages (A) Green fleshy pods with brilliant gree sedds (27 Days After Anthesis; DAA), (B) green pergamenous pods with brown soft seeds (35 DAA) and (C) green pergamenous pods with brown seeds and moderately hard (41 DAA).
}
\source{
Data are kindly provided by Fabio Gresta (University of Reggio Calabria, Italy).
}
\references{
Gresta, F., Avola, G., Onofri, A., Anastasi, U., Cristaudo, A., 2011. When Does Hard Coat Impose Dormancy in Legume Seeds? Lotus and Scorpiurus Case Study. Crop Science 51, 1739–1747.
}
\examples{
data(lotusCum)
head(lotusCum)
moniTime <- lotusCum$Time
count <- lotusCum$nCum
nViable <- rep(25, length(lotusCum[,1]))
Dish <- as.factor(lotusCum$Dish)
treatGroups <- lotusCum[,1]
dataset_sd <- makeDrm2(moniTime, count, nViable, Dish, treatGroups)
head(dataset_sd)
count <- lotusCum$nSeeds
dataset_sd <- makeDrm2(moniTime, count, nViable, Dish, treatGroups, cumulative = FALSE)
head(dataset_sd)
}
\keyword{datasets}
