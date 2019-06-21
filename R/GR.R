GRate <- function(mod, respLev, type="absolute", vcov. = sandwich){
  #mod <- obj; x <- Psi; respLev <- g
  vcMat <- vcov.
  if(is.null(mod$fct[["name"]]) == T) {
    choice <- "noName"
  }else{choice <- mod$fct[["name"]]}
  temp <- ED2.drc(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcMat, display=F)
  if( choice == "LL.4" |  choice == "LL.3" | choice == "LL.2" |
      choice == "LN.4" |  choice == "LN.3" | choice== "LN.2" |
      choice == "W1.4" |  choice == "W1.3" | choice == "W1.2" |
      choice == "W2.4" |  choice == "W2.3" | choice == "W2.2" ) {
     GRval <- 1/temp[,1]
     GRes <- (temp[,2] * 1/temp[,1]^2)
  }else{
     GRval <- temp[,1]
     GRes <- temp[,2]
  }
  GRval[is.nan(GRval)==T] <- 0
  GRes[is.nan(GRes)==T] <- 0
  GR <- data.frame(Estimate=GRval, SE=GRes)
  row.names(GR) <- gsub("e:", "GR:", row.names(GR), fixed=T)
  GR
}

GTime <- function(mod, respLev, type="absolute", vcov. = sandwich){
  vcMat <- vcov.
  if(is.null(mod$fct[["name"]]) == T) {
    choice <- "noName"
  }else{choice <- mod$fct[["name"]]}
  temp <- ED2.drc(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcMat, display=F)
  if( choice == "LL.4" |  choice == "LL.3" | choice == "LL.2" |
      choice == "LN.4" |  choice == "LN.3" | choice== "LN.2" |
      choice == "W1.4" |  choice == "W1.3" | choice == "W1.2" |
      choice == "W2.4" |  choice == "W2.3" | choice == "W2.2" ) {
     GRval <- temp[,1]
     GRes <- temp[,2]
  }else{
     GRval <- 1/temp[,1]
     GRes <- (temp[,2] * 1/temp[,1]^2)
  }
  GRval[is.nan(GRval)==T] <- Inf
  GRes[is.nan(GRes)==T] <- NA
  GR <- data.frame(Estimate=GRval, SE=GRes)
  row.names(GR) <- gsub("e:", "T:", row.names(GR), fixed=T)
  GR
  }
