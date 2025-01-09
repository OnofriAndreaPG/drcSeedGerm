\name{makeDrm}
\alias{makeDrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshape a seed germination datasets for time-to-event analyses (deprecated)
}
\description{
MakeDrm: this function reshapes a common field book for germination assays into the form that is required for time-to-event analyses with drm() in the drc package and drmte() in the drcte package. The common field book has one row per each Petri dish and the counts of germinated seeds at each assessment time are listed in different columns.
}
\usage{
makeDrm(counts, treat, nViable, moniTimes)
}

\arguments{
\item{counts}{
Dataframe listing the counts of germinated seeds in each Petri dish (rows) at each assessment time (columns).
}
\item{treat}{
Dataframe listing for each dish the levels of each treatment. It has as many columns as there are treatments in the experiment.
}
\item{nViable}{
A vector with the number of viable seeds per dish, at the beginning of the assay
}
\item{moniTimes}{
Vector of monitoring times. It needs to be of the same length as the number of columns in 'counts'.
}
}
\value{
This function returns a dataframe
}
\author{
Andrea Onofri
}
\keyword{seed germination}% use one of  RShowDoc("KEYWORDS")
