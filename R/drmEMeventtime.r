"drmEMeventtime" <- 
function(dose, resp, multCurves, doseScaling = 1)
{ 
    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
      Fstart <- multCurves(dose[, -2] / doseScaling, c)
      dose2 <- dose[, -1]
      Fend <- multCurves(dose2 / doseScaling, c)

      ifelse(is.matrix(dose2)==T,
             Fend[is.finite(dose2[, 1])==F] <- 1,
             Fend[!is.finite(dose2)] <- 1)
        
      temp <- Fend - Fstart
      temp[temp==0] <- 10e-8
        
      return( -sum(resp * log(temp)) )  
        # minus in front of sum() as maximization is done as minimization
    }    

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )
        
        c(
        -object$"fit"$value,  # oops a constant is missing!
        object$"sumList"$"df.residual"
        )
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFeventtime" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
