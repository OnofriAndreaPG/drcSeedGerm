# Halotime model with normal distribution of base salt concentration
# 9/1/2024
HaloTnorm.fun <- function(time, SConc, ThetaHalo, SConcb50, sigmaSConcb){
  pnorm((SConc + (ThetaHalo/time) - SConcb50)/sigmaSConcb, lower.tail = F) }
"HaloTnorm" <- function(){
  fct <- function(x, parm){
    HaloTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
  names <- c("ThetaHalo", "SConcb50", "sigmaSConcb")
  name <- "HaloTnorm"
  text <- "Halotime model with normal distribution of SConcb"
  ss <- function(data){
    x1 <- data[, 1]
    x2 <- data[, 2]
    y <- data[, 3]
    pseudoY <- qnorm((y+10e-6)*0.99, lower.tail = F)
    mod <- lm(pseudoY ~ I(1/x1) + x2)
    sigmaSConcb <- 1/coef(mod)[3]
    SConcb50 <- -coef(mod)[1]*sigmaSConcb
    ThetaHalo <- coef(mod)[2]*sigmaSConcb
    return(c(ThetaHalo, SConcb50, sigmaSConcb))
  }
  GR <- function(parms, respl, reference="control", type="relative", SConc){
    # print(respl)
    HaloTnorm.gra <- function(p1, p2, p3, SConc, g) {
      # parameters: ThetaHalo, SConcb50, sigmaSConcb
      GR <- (p3 * qnorm(g, lower.tail = F) - SConc + p2 ) / p1 #returns rate
      GR <- ifelse(GR > 0, GR, 0)
      1/GR # returns time
    }
    p1 <- as.numeric(parms[1])
    p2 <- as.numeric(parms[2])
    p3 <- as.numeric(parms[3])
    g <- respl
    if(type=="absolute"){

      EDp <- HaloTnorm.gra(p1, p2, p3, SConc, g)

      #Approximation of derivatives(finite differences)
      d1.1 <- HaloTnorm.gra(p1, p2, p3, SConc, g)
      d1.2 <- HaloTnorm.gra(p1 + 10e-6, p2, p3, SConc, g)
      d1 <- (d1.2 - d1.1)/10e-6

      d2.1 <- HaloTnorm.gra(p1, p2, p3, SConc, g)
      d2.2 <- HaloTnorm.gra(p1, p2  + 10e-6, p3, SConc, g)
      d2 <- (d2.2 - d2.1)/10e-6

      d3.1 <- HaloTnorm.gra(p1, p2, p3, SConc, g)
      d3.2 <- HaloTnorm.gra(p1, p2, p3 + 10e-6, SConc, g)
      d3 <- (d3.2 - d3.1)/10e-6
      EDder <- c(d1, d2, d3)
    } else{ if(type=="relative") {
      .Pmax <- pnorm((SConc - p2 )/p3, lower.tail = FALSE)
      grel <- .Pmax*g
      EDp <- HaloTnorm.gra(p1, p2, p3, SConc, grel)

      #Approximation of derivatives(finite differences)
      d1.1 <- HaloTnorm.gra(p1, p2, p3, SConc, grel)
      d1.2 <- HaloTnorm.gra(p1 + 10e-6, p2, p3, SConc, grel)
      d1 <- (d1.2 - d1.1)/10e-6

      d2.1 <- HaloTnorm.gra(p1, p2, p3, SConc, grel)
      d2.2 <- HaloTnorm.gra(p1, p2  + 10e-6, p3, SConc, grel)
      d2 <- (d2.2 - d2.1)/10e-6

      d3.1 <- HaloTnorm.gra(p1, p2, p3, SConc, grel)
      d3.2 <- HaloTnorm.gra(p1, p2, p3 + 10e-6, SConc, grel)
      d3 <- (d3.2 - d3.1)/10e-6
      EDder <- c(d1, d2, d3)
    } }
    return(list(EDp, EDder))
  }

  deriv1 <- function(x, parm){
    #Approximation by using finite differences
    d1.1 <- HaloTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
    d1.2 <- HaloTnorm.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3])
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HaloTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
    d2.2 <- HaloTnorm.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3])
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HaloTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
    d3.2 <- HaloTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6))
    d3 <- (d3.2 - d3.1)/10e-6

    cbind(d1, d2, d3)
  }
  returnList <- list(fct=fct, ssfct=ss, name = name,
                     names=names, text=text, edfct=GR,
                     deriv1=deriv1)
  class(returnList) <- "drcMean"
  invisible(returnList)
}
