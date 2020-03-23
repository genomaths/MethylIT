#' @rdname AICquasiPoisson
#'
#' @title AIC for Quasi-Poisson glm model
#' @description AIC for Quasi-Poisson glm model dpois
#'
#' @param fitObj a fitted model object
#'
#' @return numeric AIC for Quasi-Poisson glm model
#'
#' @importFrom stats predict coef dpois
#' @keywords internal
.AICquasiPoisson <- function(fitObj) {
    LogLike <- sum(dpois(fitObj$data$count, lambda = exp(predict(fitObj)), 
        log = TRUE))
    return(2 * (length(coef(fitObj)) - LogLike))
}
