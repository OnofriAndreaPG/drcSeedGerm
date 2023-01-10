\name{HTE3}
\alias{HTE3}
\alias{HTE3.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydro-time model with shifted exponential for Pmax and linear model for GR50
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HTE3()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function ()
{
    HTE1.fct <- function(x, parm) {
        S <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2], parm[,
            3], parm[, 4], parm[, 5])
        return(S)
    }
    HTE1.names <- c("G", "Psib", "sigmaPsib", "thetaH", "b")
    HTE1.ss <- function(data) {
        data <- subset(data, is.finite(data[, 1]) == T)
        result <- c()
        PsiF <- factor(data[, 2])
        for (i in 1:length(levels(PsiF))) {
            temp <- subset(data, data[, 2] == levels(PsiF)[i])
            x <- temp[, 1]
            y <- temp[, 3]
            modT <- try(nls(y ~ NLSLL.3(x, a, b, c)), silent = T)
            if (class(modT) == "try-error") {
                res <- as.numeric(levels(PsiF)[i])
                result <- c(result, res)
            }
        }
        result
        dataset_cum <- subset(data, is.finite(data[, 1]) == T)
        if (is.null(result) != T) {
            for (i in 1:length(result)) dataset_cum <- subset(dataset_cum,
                dataset_cum[, 2] != result[i])
        }
        PsiF <- factor(dataset_cum[, 2])
        x1 <- dataset_cum[, 1]
        x2 <- dataset_cum[, 2]
        y <- dataset_cum[, 3]
        modI <- drm(y ~ x1, fct = LL.3(), curveid = PsiF, pmodels = list(~1,
            ~PsiF - 1, ~PsiF - 1), data = dataset_cum)
        psiLevels <- as.numeric(levels(PsiF))
        b <- -coef(modI)[1]
        Pmax <- coef(modI)[2:(length(psiLevels) + 1)]
        modPmax <- drm(Pmax ~ psiLevels, fct = PmaxPsi1())
        G <- coef(modPmax)[3]
        Psib <- coef(modPmax)[1]
        sigmaPsib <- coef(modPmax)[2]
        GR50 <- 1/coef(modI)[(length(psiLevels) + 2):length(coef(modI))]
        modGR <- drm(GR50 ~ psiLevels, fct = GRPsiLin())
        thetaH <- coef(modGR)[2]
        Psib2 <- coef(modGR)[1]
        psib <- mean(Psib, Psib2)
        return(c(G, psib, sigmaPsib, thetaH, b))
    }
    HTE1.GR <- function(parms, respl, reference = "control",
        type = "relative", Psi) {
        G <- as.numeric(parms[1])
        Psib <- as.numeric(parms[2])
        sigmaPsib <- as.numeric(parms[3])
        thetaH <- as.numeric(parms[4])
        b <- as.numeric(parms[5])
        g <- respl/100
        if (type == "absolute") {
            .Pmax <- Pmax2b.fun(Psi, G, Psib, sigmaPsib)
            .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
            .temp2 <- (.Pmax - g)/g
            .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
            .GR50 <- GR501.fun(Psi, Psib, thetaH)
            .GR50 <- ifelse(.GR50 > 0, .GR50, 0)
            res <- as.numeric(exp(-(1/b) * log(.temp2) + log(1/.GR50)))
            EDp <- 1/res
            .GR <- EDp
            EDder <- c(1/.GR * ((1/b) * ((1/g) * (.Pmax/G)/((1/g) *
                (.Pmax - g))))/(1/.GR^2), -(((1/b) * ((1/g) *
                (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) * (1/sigmaPsib)))/((.Pmax -
                g)/g)) + 1/(Psi - Psib))/(1/.GR)), -(((1/b) *
                ((1/g) * (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) *
                  ((Psi - Psib) * (1/sigmaPsib^2))))/((.Pmax -
                  g)/g)))/(1/.GR)), -(1/(Psi - Psib)/(1/.GR50))/(1/.GR),
                -(1/.GR * (1/b^2 * log((1/g) * (.Pmax - g)))/exp(-(1/b) *
                  log((1/g) * (.Pmax - g)) + log(1/.GR50))^2))
        }
        else {
            if (type == "relative") {
                .Pmax <- Pmax2b.fun(Psi, G, Psib, sigmaPsib)
                .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
                .temp2 <- (1 - g)/g
                .GR50 <- GR501.fun(Psi, Psib, thetaH)
                .GR50 <- ifelse(.GR50 > 0, .GR50, 0)
                res <- as.numeric(exp(-(1/b) * log(.temp2) +
                  log(1/.GR50)))
                EDp <- 1/res
                .GR <- EDp
                EDder <- c(0, -(1/.GR * ((1/.GR50)/(Psi - Psib)/(1/.GR50))/(1/.GR50)^2),
                  0, -(1/.GR * (1/(Psi - Psib)/(1/.GR50))/(1/.GR)^2),
                  -(.GR * log((1 - g)/g))/(b^2))
            }
        }
        return(list(EDp, EDder))
    }
    HTE1.deriv1 <- function(x, parm) {
        d1.1 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], parm[, 5])
        d1.2 <- HTE1.fun(x[, 1], x[, 2], (parm[, 1] + 1e-05),
            parm[, 2], parm[, 3], parm[, 4], parm[, 5])
        d1 <- (d1.2 - d1.1)/1e-05
        d2.1 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], parm[, 5])
        d2.2 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], (parm[, 2] +
            1e-05), parm[, 3], parm[, 4], parm[, 5])
        d2 <- (d2.2 - d2.1)/1e-05
        d3.1 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], parm[, 5])
        d3.2 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            (parm[, 3] + 1e-05), parm[, 4], parm[, 5])
        d3 <- (d3.2 - d3.1)/1e-05
        d4.1 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], parm[, 5])
        d4.2 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], (parm[, 4] + 1e-05), parm[, 5])
        d4 <- (d4.2 - d4.1)/1e-05
        d5.1 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], parm[, 5])
        d5.2 <- HTE1.fun(x[, 1], x[, 2], parm[, 1], parm[, 2],
            parm[, 3], parm[, 4], (parm[, 5] + 1e-05))
        d5 <- (d5.2 - d5.1)/1e-05
        cbind(d1, d2, d3, d4, d5)
    }
    HTE1.text <- "Hydro-time model with shifted exponential for Pmax and linear model for GR50"
    HTE1 <- list(fct = HTE1.fct, ssfct = HTE1.ss, names = HTE1.names,
        text = HTE1.text, edfct = HTE1.GR, deriv1 = HTE1.deriv1)
    class(HTE1) <- "drcMean"
    invisible(HTE1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
