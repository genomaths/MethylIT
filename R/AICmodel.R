#' @rdname AICmodel
#'
#' @title Akaike's Information Criterion (AIC)
#' @description this function permits the estimation of the AIC for models for
#'     which the function 'AIC' from the 'stats' package does not work.
#' @details if for a given model 'm' AIC(m) works, then AICmodel(m) = AIC(m).
#'
#' @param model if provided, it is an R object from where the residuals and
#'     model parameters can be retrieved using resid(model) and coef(model),
#'     respectively.
#' @param residuals if provided, it is numerical vector with the residuals:
#'     residuals = observed.values - predicted.values, where predicted values
#'     are estimated from the model. If the parameter 'model' is not provided,
#'     then this parameter must be provided.
#' @param np number of model parameters. If the parameter 'model' is not
#'     provided, then 'np' and 'residuals' must be provided.
#'
#' @return AIC numerical value
#'
#' @examples
#' ## Build a set of point holding an exponential decay
#' ## and adding random noise.
#' set.seed(77)
#' x = runif(100, 1, 5)
#' y = 2 * exp(-0.5 * x) + runif(100, 0, 0.1)
#' plot(x, y)
#'
#' ## Non-linear regression
#' nlm <- nls(Y ~ a * exp( b * X), data = data.frame(X=x, Y=y),
#'             start=list(a=1.5, b=-0.7),
#'             control=nls.control(maxiter=10^4, tol=1e-05),
#'             algorithm='port')
#'
#' ## Estimations of Akaike information criteria given by 'AIC' function
#' ## from stats' R package and from 'AICmodel' function are equal.
#' AICmodel(nlm) == AIC(nlm)
#'
#' ## Now, using residuals from the fitted model:
#' res = y - coef(nlm)[1] * exp(coef(nlm)[2] * x)
#'
#' AICmodel(residuals=res, np=2) == AIC(nlm)
#'
#' @importFrom stats resid coef
#'
#' @keywords internal
#' @export
AICmodel <- function(model = NULL, residuals = NULL,
    np = NULL) {
    if (is.null(model) && is.null(residuals)) {
        stop(paste("At least one of the parameters 'model'",
            " or 'residual' must be provided"))
    }
    if (!is.null(model)) {
        if (is.null(np))
            np <- length(coef(model))
        RESID <- resid(model)
    }
    if (is.null(np)) {
        stop("The number of model parameters 'np' must be provided")
    }
    if (!is.null(residuals))
        RESID <- residuals

    sse <- sum(RESID^2, na.rm = TRUE)
    N <- length(RESID)
    return(N * (1 + log(2 * pi) + log(sse/N)) + 2 *
        (1L + np))
}
