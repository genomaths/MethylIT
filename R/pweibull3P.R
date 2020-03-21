#' @rdname pweibull3P
#'
#' @title Weibull distribution with three parameters
#' @description Density, distribution function, quantile function and random
#'     generation for the Weibull distribution with three parameters
#'
#' @param q vector of quantiles
#' @param shape shape parameter, or slope, defaulting to 1
#' @param scale scale parameter, or characteristic life,  defaulting to 1
#' @param mu location parameter, or failure free life,  defaulting to 0
#'
#' @return 3 parameters Weibull distribution
#'
#' @importFrom stats pweibull
#'
#' @examples
#' num.samples <- 10000
#' shape <- 0.75
#' scale <- 1
#' x <- rweibull(num.samples, shape = shape, scale = scale)
#' wei.model <- weibull3P(x)
#' y <- pweibull3P(x,
#'                 shape = as.numeric(wei.model$Estimate[1]),
#'                 scale = as.numeric(wei.model$Estimate[2]),
#'                 mu = as.numeric(wei.model$Estimate[3]) )
#'
#' @export
pweibull3P <- function(q, shape = 1, scale = 1, mu = 0) {
    y <- q - mu
    pweibull(y, shape, scale)
}
