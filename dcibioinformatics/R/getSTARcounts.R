#' getSTARcounts
#' The function illustrates how to use the myLS Rcpp function
#' from this package from within an R function to fit the model
#' Y = X beta + epsilon 
#' What is returned is the vector betahat=LSE(beta)
#' @param Y the vector of response
#' @param X the matrix of covariates (does not include intercept)
#' @return \code{betahat}
#' @export
getSTARcounts <- function(Y, X) {
    betahat <- myLSrcpp(Y, X)
    return(betahat)
}
