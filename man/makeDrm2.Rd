\name{makeDrm2}
\alias{makeDrm2}
\title{
Reshape a seed germination datasets for time-to-event model fitting.
}
\description{
This function reshapes a datasets organised as necessary for nonlinear regression into the kind of dataset required by the drmte() function in the 'drcte' package. It works with either the counts of germinated seeds at each monitoring time or the cumulative counts at each monitoring time.
}
\usage{
makeDrm2(counts, treat, nViable, moniTimes, Dish,  cumulative = T)
}

\arguments{
\item{counts}{
Vector listing the counts of germinated seeds in each Petri dish at each assessment time.
}
\item{treat}{
Dataframe listing, for each row of data, the corresponding level of experimental factors (one factor per column)
}
\item{nViable}{
A vector listing the number of viable seeds, at the beginning of the assay. This number is the same for all observations belonging to the same dish.
}
\item{moniTimes}{
Vector of monitoring times.
}
\item{Dish}{
Vector of codes for dishes.
}
\item{cumulative}{
Logical: True if counts are cumulative, False if they are not.
}
}
\value{
Returns a dataframe
}
\author{
Andrea Onofri
}
\examples{
# makeDrm2 (deprecated, example not run)
# data(lotusCum)
# moniTime <- lotusCum$Time
# count <- lotusCum$nCum
# nViable <- rep(25, length(lotusCum[,1]))
# Dish <- as.factor(lotusCum$Dish)
# treatGroups <- lotusCum[,1]
# dataset_sd <- makeDrm2(count, treatGroups, nViable, moniTime, Dish)
# head(dataset_sd)
# count <- lotusCum$nSeeds
# dataset_sd <- makeDrm2(count, treatGroups, nViable, moniTime, Dish,  # cumulative = FALSE)
# head(dataset_sd)
}
\keyword{seed germination}
