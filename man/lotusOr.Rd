\name{lotusOr}
\alias{lotusOr}
\docType{data}
\title{
Germination assay with Lotus ornithopodioides
}
\description{
Seeds of L. ornithopodioides were harvested in three stages after anthesis and the time course of their germination was evaluated in Petri dishes. This is the original dataset, as collected by the technicians and needs to be transformed into a form that is suitable for time-to-event analysis in R.
}
\usage{data("lotusOr")}
\format{
  A data frame with 12 observations on the following 17 variables.
  \describe{
    \item{\code{Stage}}{a factor with levels \code{A} \code{B} \code{C}}
    \item{\code{Dish}}{a numeric vector: code of Petri dish}
    \item{\samp{1}}{a numeric vector: counts at day 1}
    \item{\samp{2}}{a numeric vector: counts at day 2}
    \item{\samp{3}}{a numeric vector: counts at day 3}
    \item{\samp{4}}{a numeric vector: counts at day 4}
    \item{\samp{5}}{a numeric vector: counts at day 5}
    \item{\samp{6}}{a numeric vector: counts at day 6}
    \item{\samp{7}}{a numeric vector: counts at day 7}
    \item{\samp{8}}{a numeric vector: counts at day 8}
    \item{\samp{9}}{a numeric vector: counts at day 9}
    \item{\samp{10}}{a numeric vector: counts at day 10}
    \item{\samp{11}}{a numeric vector: counts at day 11}
    \item{\samp{12}}{a numeric vector: counts at day 12}
    \item{\samp{13}}{a numeric vector: counts at day 13}
    \item{\samp{14}}{a numeric vector: counts at day 14}
    \item{\samp{15}}{a numeric vector: counts at day 15}
  }
}
\details{
Every line of data represents a Petri dish. There were 25 seeds per Petri dish. The columns represent the characteristics of each Petri dish. The columns from 1 to 15 represent the number of germinated seeds counted in each assessment time. Seeds of Lotus ornithopodioides were harvested at the following stages (A) Green fleshy pods with brilliant gree sedds (27 Days After Anthesis; DAA), (B) green pergamenous pods with brown soft seeds (35 DAA) and (C) green pergamenous pods with brown seeds and moderately hard (41 DAA).
}
\source{
Data are kindly provided by Fabio Gresta (University of Reggio Calabria, Italy).
}
\references{
Gresta, F., Avola, G., Onofri, A., Anastasi, U., Cristaudo, A., 2011. When Does Hard Coat Impose Dormancy in Legume Seeds? Lotus and Scorpiurus Case Study. Crop Science 51, 1739–1747. }
\examples{
data(lotusOr)
counts <- lotusOr[,3:length(lotusOr[1,])]
treat <- data.frame(tratt=lotusOr[,1])
nViable <- rep(25,12)
moniTimes <- c(1:15)
datasetDrc <- makeDrm(counts=counts, treat=treat, nViable=nViable, moniTimes)
head(datasetDrc, 30)
}
\keyword{time-to-event-analysis, seed germination}
