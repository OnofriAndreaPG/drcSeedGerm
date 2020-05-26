# This function organises a dataset to be submitted to nonlinear regression
# Start: field book
# Arrival: dataset with cumulative proportions of germinated seeds

#dataset is a dataframe listing for each Petri dish (rows) the number of germinated seeds
#at each assessment time (culumns);
#treat is a dataframe listing for each Petri dish (rows) the levels of each treatment (culumns)
#seeds = vector listing the number of viable seeds for each Petri dish (same length as
#number of rows in dataset and treat.
#moniTimes = vector of monitoring times. Seme length as number of columns in dataset

makeNlin <- function(counts, treat, nViable, moniTimes){
  dataset <- makeDrm(counts, treat, nViable, moniTimes)
  dataset <- dataset[is.finite(dataset$timeAf)==T, ]
  output <- data.frame(dataset[,1:length(treat[1,])], Time = dataset$timeAf,
    propCum = dataset$propCum, row.names = 1:length(dataset[,1]))
  output
}
