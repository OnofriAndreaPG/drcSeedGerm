# Make drm, now superseded

library(drcSeedGerm)
data(lotusOr)
lotusOr
counts <- lotusOr[,3:length(lotusOr[1,])]
treat <- data.frame(tratt=lotusOr[,1])
nViable <- rep(25,12)
moniTimes <- c(1:15)
devtools::load_all()
makeDrm(counts=counts, treat=treat, nViable=nViable, moniTimes)
head(datasetDrc, 30)

datasetDrc <- makeDrm(counts = counts, treat=treat, nViable=nViable, moniTimes)
head(datasetDrc, 30)
datasetDrc

ncols_counts <- length(counts[1,])
ncols_treat <- length(treat[1,])
sprova <- data.frame(treat, counts)

melt_te(sprova, count_cols = (ncols_treat + 1):(ncols_counts+ncols_treat),
        treat_cols = 1:ncols_treat,
            monitimes = moniTimes, n.subjects = nViable)
head(datasetG, 16)
