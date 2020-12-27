#' @rdname estimateBetaDist
#'
#' @title Select the beta distribution that fit specified quantiles
#' @description This function perform a rough estimation of the shape
#'     parameters of beta distribution
#' @details In order to obtain the estimates for shape parameters of beta
#' distribution, the squared of the difference between the empirical cumulative
#' distribution function (ecdf) & the theoretical cdf is minimized using the
#' Non-Linear Minimization function of 'stats' package.
#'
#' @param q prior probabilities
#' @param init.pars  initial parameter values. Defaults is NULL and an initial
#' guess is estimated using \code{\link[stats]{optim}} function. If the initial
#' guessing fails initial parameter values are to alpha = 1 &
#' beta = 1, which imply the parsimony pseudo-counts greater than zero.
#' @param maxiter,ftol,ptol,gradtol Optional parameters for
#' \code{\link[minpack.lm]{nlsLM}} and \code{\link[stats]{nlm}} functions.
#'
#' @return the estimated values of the shape parameters of the selected beta
#'    distribution.
#'
#' @importFrom stats ecdf optim nlm pbeta
#' @importFrom minpack.lm nls.lm
#' @keywords internal
.estimateBetaDist <- function(
                            q,
                            init.pars = NULL,
                            maxiter = 1024,
                            ftol = 1e-12,
                            ptol = 1e-12,
                            gradtol = 1e-8) {

    ## q: prior probabilities
    ## init.pars: initial parameter values. Defaults to alpha = 1 & beta = 1,
    ## which the parsimony pseudo-counts greater than zero.

    Q <- ecdf(q)  ## Empirical Cumulative Distribution Function
    dat <- data.frame(Q = Q(q), q = q)

    if(is.null(init.pars)) {
        init.pars <- try(suppressWarnings(
                        optim(par = c(1,1),
                              min.RSS, data = dat,
                              method = "BFGS",
                              control = list(maxit = 500,
                                             abstol = (1e-12)))),
                        silent = TRUE)
        init.pars <- c(init.pars$par, min(q, na.rm = TRUE))
        if (inherits(init.pars, "try-error")) {
            init.pars <- c(1,1,min(q, na.rm = TRUE))
        }
    }

    formula <- as.formula("Y ~ pbeta3( X, shape1, shape2, mu)")
    start <- list(shape1 = init.pars[1], shape2 = init.pars[2],
                    mu = init.pars[3])
    fit1 <- suppressWarnings(try(nlsLM(formula,
                                       data = data.frame(X = q, Y = Q(q)),
                                       start = start,
                                       control = list(maxiter = maxiter,
                                                      ftol = ftol,
                                                      ptol = ptol)),
                                 silent = TRUE))

    formula <- as.formula("Y ~ pbeta( X, shape1, shape2)")
    start <- list(shape1 = init.pars[1], shape2 = init.pars[2])
    fit2 <- suppressWarnings(try(nlsLM(formula,
                                       data = data.frame(X = q, Y = Q(q)),
                                       start = start,
                                       control = list(maxiter = maxiter,
                                                      ftol = ftol,
                                                      ptol = ptol)),
                                 silent = TRUE))


    if (!inherits(fit1, "try-error") && !inherits(fit2, "try-error"))
        fit <- list(fit1, fit2)[[ which.min(c(AIC(fit1), AIC(fit2))) ]]

    if (!inherits(fit1, "try-error") && inherits(fit2, "try-error"))
        fit <- fit1

    if (inherits(fit1, "try-error") && !inherits(fit2, "try-error"))
        fit <- fit2

    validfit <- (!inherits(fit1, "try-error") || !inherits(fit2, "try-error"))

    if (validfit)
        fit_summary <- try(summary(fit), silent = TRUE)

    if (inherits(fit_summary, "try-error")) {
        fit <- try(suppressWarnings(
            nlm(min.RSS, p = init.pars,
                data = dat, gradtol = gradtol)),
            silent = TRUE)
        if (!inherits(fit, "try-error")) {
            fit <- try(suppressWarnings(
                optim(par = init.pars,
                    min.RSS, data = dat, method = "BFGS",
                    control = list(maxit = 500,
                                    abstol = (10^-8)))),
                silent = TRUE)

            if (inherits(fit, "try-error")) {
                par <- c(0, 0)
            } else par <- fit$par
        } else par <- fit$estimate
    } else par <- fit_summary$parameters[1:2, 1]
    par
}

### ----------- Auxiliary function --------------------
### # ---- 3P beta CDF --- #
pbeta3 <- function(
                    q,
                    shape1 = 2,
                    shape2 = 3,
                    mu = 0,
                    lower.tail = TRUE,
                    log.p = FALSE ) {

    pbeta( q - mu,
           shape1,
           shape2,
           lower.tail = lower.tail,
           log.p = log.p )
}


optFun <- function(par, y, q, probfun) {
    alpha <- par[1]  ## This corresponds to the initial
    ## starting parameter for alpha
    beta <- par[2]  ## This corresponds to the initial
    ## starting parameter for beta Because optim
    ## minimizes a function, the
    return((y - pbeta(q, alpha, beta))^2)
}

#### Beta fitting In order to obtain the estimates for
#### shape paramaters the squared of the difference
#### between the ecdf & the theoretical cdf is
#### minimized
min.RSS <- function(data, p) {
    alpha <- p[1]  ## This corresponds to the initial
    ## starting parameter for alpha
    beta <- p[2]  ## This corresponds to the initial
    ## starting parameter for beta Because optim
    ## minimizes a function, the
    with(data, sum((Q - pbeta(q, alpha, beta))^2))
}

