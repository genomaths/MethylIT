#' @rdname estimateECDF
#'
#' @title A variant of Empirical Cumulative Distribution Function 'ecdf'
#' @description This function is used to reduce the number of points used to
#'     build a ecdf of an experimental variable. When a variable has a very
#'     high amount of experimental values (several millions) could be
#'     computationally time consuming to perform a good-of-fit test and a
#'     non-linear regression estimation for a theoretical CDF based in such a
#'     big number of values.
#' @details The histogram cell midpoints values are used to build estimateECDF.
#'
#' @param x numeric vector
#' @param npoints minimum number of non-missing values
#'
#' @return ecdf of numeric vector
#'
#' @examples
#'     x <- sample(1:500, 50, replace=TRUE)
#'     estimateECDF(x, npoints = 15)
#'
#' @importFrom graphics hist
#' @importFrom stats approxfun
#' @keywords internal
#' @export
estimateECDF <- function(x, npoints = 10) {
    x <- sort(x)
    n <- length(x)
    if (n < 10) 
        stop("'x' must have 10 or more non-missing values")
    
    ## To reduce the number of point used in the
    ## nonlinear fit
    h <- hist(x, breaks = npoints, plot = FALSE)
    vals <- h$mids
    rval <- approxfun(vals, cumsum(h$counts)/n, method = "constant", 
        yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    return(rval)
}
