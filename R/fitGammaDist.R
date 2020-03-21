#' @rdname fitGammaDist
#'
#' @title Nonlinear fit of Gamma CDF (Gamma)
#' @description This function performs the nonlinear fit of GGamma CDF of a
#'     variable x
#' @details The algorithm tries to fit the two-parameter Gamma CDF
#' ('Gamma2P') or the three-parameter Gamma ('Gamma3P') using a modification of
#' Levenberg-Marquardt algorithm implemented in function 'nls.lm' from
#' 'minpack.lm' package that is used to perform the nonlinear fit.
#' Cross-validations for the nonlinear regressions (R.Cross.val) were performed
#' in each methylome as described in reference (1). In addition, Stein's formula
#' for adjusted R squared (rho) was used as an estimator of the average
#' cross-validation predictive power (1).
#'
#' If the number of values to fit is >10^6, the fitting to a GGamma CDF would be
#' a time consuming task. To reduce the computational time, the data can be
#' 'summarized' into 'npoints' (bins) and used as the new predictors.
#' @param x numerical vector
#' @param probability.x probability vector of x. If not provided, the values
#'     are estimated using the empirical cumulative distribution function
#'     ('ecdf') from 'stats' R package.
#' @param parameter.values initial parameter values for the nonlinear fit. If
#'     the locator parameter is included (mu != 0), this must be given as
#'     parameter.values = list(shape = 'value', scale = 'value', mu = 'value')
#'     or if mu = 0, as: parameter.values = list(shape = 'value',
#'     scale = 'value'). If not provided, then an initial guess is provided.
#' @param location.par whether to consider the fitting to generalized gamma
#'     distribution (Gamma) including the location parameter, i.e., a Gamma
#'     with four parameters (GGamma3P).
#' @param sample.size size of the sample.
#' @param npoints number of points used in the fit.
#' @param maxiter positive integer. Termination occurs when the number of
#'     iterations reaches maxiter. Default value: 1024.
#' @param ftol non-negative numeric. Termination occurs when both the actual
#'     and predicted relative reductions in the sum of squares are at most ftol.
#'     Therefore, ftol measures the relative error desired in the sum of
#'     squares. Default value: 1e-12
#' @param ptol non-negative numeric. Termination occurs when the relative error
#'     between two consecutive iterates is at most ptol. Therefore, ptol
#'     measures the relative error desired in the approximate solution.
#'     Default value: 1e-12.
#' @param maxfev integer; termination occurs when the number of calls to fn has
#'     reached maxfev. Note that nls.lm sets the value of maxfev to
#'     100*(length(par) + 1) if maxfev = integer(), where par is the list or
#'     vector of parameters to be optimized.
#' @param nlms Logical. Whether to return the nonlinear model object
#'     \code{\link[minpack.lm]{nls.lm}}. Default is FALSE.
#' @param verbose if TRUE, prints the function log to stdout
#'
#' @return Model table with coefficients and goodness-of-fit results:
#'     Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
#'     coefficient covariance matrix.
#'
#' @references 1. Stevens JP. Applied Multivariate Statistics for the Social
#'     Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 06/03/2016
#'
#' @examples
#' set.seed(126)
#' x <- rgamma(1000, shape = 1.03, scale = 2.1)
#' fitGammaDist(x)
#'
#' @importFrom stats var stepfun as.formula coef AIC pweibull BIC vcov knots
#'   deviance BIC cor quantile
#' @importFrom nls2 nls2
#' @importFrom utils head
#' @importFrom stats ecdf nls.control
#' @importFrom minpack.lm nlsLM nls.lm nls.lm.control
#' @export

fitGammaDist <- function(x, probability.x, parameter.values,
    location.par = FALSE, sample.size = 20, npoints = NULL,
    maxiter = 1024, ftol = 1e-12, ptol = 1e-12, maxfev = 1e+05,
    nlms = FALSE, verbose = TRUE) {
    ind <- which(x > 0)
    if (length(ind) > sample.size) {
        x <- x[ind]
        if (missing(probability.x))
            Fy <- ecdf(x)
    } else stop("*** Sample size is lower than the set minimun size: ",
        sample.size)

    pgamma3p <- function(q, shape, scale, mu) pgamma(q -
        mu, shape = shape, scale = scale)

    getPreds <- function(par, q) {
        if (length(par) > 2) {
            pred <- pgamma3p(q, shape = par[1], scale = par[2],
                mu = par[3])
        } else pred <- pgamma(q, shape = par[1], scale = par[2])
        return(pred)
    }

    optFun <- function(par, probfun, quantiles, prob,
        eval = FALSE) {
        START <- as.list(par)
        START$q <- quantiles
        EVAL <- try(do.call(probfun, START), silent = TRUE)
        if (inherits(EVAL, "try-error"))
            return(NA)
        EVAL[is.nan(EVAL)] <- 0
        RSS <- (prob - EVAL)^2
        if (eval) {
            return(EVAL)
        } else return(RSS)
    }

    N <- length(x)

    ## To reduce the number of points to be used in the
    ## fit
    if (!is.null(npoints) && npoints < N) {
        x <- pretty(x, n = npoints)
        x <- x[x > 0]
    }

    n <- length(x)  ## size of the sample used in the computation
    pX <- Fy(x)

    if (verbose && !location.par && !is.null(npoints)) {
        message(paste0("*** Trying nonlinear fit to a generalized 2P Gamma ",
            "distribution model (summarized data: ",
            npoints, " values)..."))
    } else {
        if (verbose && location.par && !is.null(npoints)) {
            message(paste0("*** Trying nonlinear fit to a 3P Gamma ",
                "distribution model ", npoints, " values..."))
        }
    }

    ## =============== starting parameter values
    ## =========== #
    if (missing(parameter.values)) {
        MEAN <- mean(x, na.rm = TRUE)
        VAR <- var(x, na.rm = TRUE, use = "everything")
        MIN <- min(x, na.rm = TRUE)

        shape = MEAN^2/VAR
        mu = MIN
        scale = VAR/MEAN

        if (location.par)
            starts <- list(shape = shape, scale = scale,
                mu = mu[1]) else starts <- list(shape = shape, scale = scale)
    } else starts = parameter.values

    ## ============ END starting parameter values
    ## ========== #

    ## ==================== Fitting models =================
    ## Try approach to 'Gamma3P'
    if (location.par) {
        if (verbose)
            message(paste0("*** Trying nonlinear fit to a 3P Gamma ",
                "distribution model ..."))
        FIT <- try(nls.lm(par = starts, fn = optFun,
            probfun = pgamma3p, quantiles = x, prob = pX,
            control = nls.lm.control(maxiter = maxiter,
                ftol = ftol, maxfev = maxfev, ptol = ptol)),
            silent = TRUE)
        if (inherits(FIT, "try-error")) {
            messg = paste0("* The 'Gamma' model did not fit the data.\n",
                "Trying to fit based on approach to 'Gamma2P' model ...")
            message(messg)
            starts <- unname(starts)

            starts <- list(shape = starts[1], scale = starts[2])
            FIT <- try(nls.lm(par = starts, fn = optFun,
                probfun = pgamma, quantiles = x, prob = pX,
                control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                        maxfev = maxfev, ptol = ptol)),
                silent = TRUE)
        }
    }

    # Try 'Gamms2P'
    if (!location.par) {
        if (verbose)
            message(paste0("*** Trying nonlinear fit to a 2P Gamma ",
                "distribution model ..."))
        FIT <- try(nls.lm(par = starts, fn = optFun,
            probfun = pgamma, quantiles = x, prob = pX,
            control = nls.lm.control(maxiter = maxiter,
                ftol = ftol, maxfev = maxfev, ptol = ptol)),
            silent = TRUE)
    }

    if (inherits(FIT, "try-error") && !location.par) {
        starts <- list(shape = 1, scale = scale)
        FIT <- try(nls.lm(par = starts, fn = optFun,
            probfun = pgamma, quantiles = x, prob = pX,
            control = nls.lm.control(maxiter = maxiter,
                ftol = ftol, maxfev = maxfev, ptol = ptol)),
            silent = TRUE)
    }


    if (!inherits(FIT, "try-error")) {
        m <- length(coef(FIT))
        ## **** R squares ****
        Adj.R.Square <- (1 - (deviance(FIT)/((n - m) *
            var(pX, use = "everything"))))
        Adj.R.Square <- ifelse(is.na(Adj.R.Square) ||
                                    Adj.R.Square < 0, 0, Adj.R.Square)

        ## Stein adjusted R square
        if (m > 2)
            rho <- ((n - 1)/(n - 4)) * ((n - 2)/(n - 5)) * ((n + 1)/n)
        else rho <- ((n - 1)/(n - 3)) * ((n - 2)/(n - 4)) * ((n + 1)/n)
        rho <- 1 - rho * (1 - Adj.R.Square)
        rho = ifelse(is.na(rho) | rho < 0, 0, rho)

        ##--Crossvalidation standard model for Nonlinear regression: x versus r
        if (verbose) {
            cat(paste("*** Performing nonlinear regression model ",
                    "crossvalidation...\n"))
        }
        set.seed(123)

        cros.ind.1 <- sample.int(N, size = round(N/2))
        cros.ind.2 <- setdiff(seq_len(N), cros.ind.1)
        starts1 <- as.list(coef(FIT))

        if (length(starts1) > 2) {
            FIT1 <- try(nls.lm(par = starts1, fn = optFun,
                probfun = pgamma3p, quantiles = x[cros.ind.1],
                prob = pX[cros.ind.1],
                control = nls.lm.control(maxiter = maxiter,
                                    ftol = ftol, maxfev = maxfev, ptol = ptol)),
                silent = TRUE)
            if (inherits(FIT1, "try-error")) {
                FIT1 <- try(nls.lm(par = starts, fn = optFun,
                                probfun = pgamma3p, quantiles = x[cros.ind.1],
                                prob = pX[cros.ind.1],
                                control = nls.lm.control(maxiter = maxiter,
                                                        ftol = ftol,
                                                        maxfev = maxfev,
                                                        ptol = ptol)),
                            silent = TRUE)
            }
        } else {
            FIT1 <- try(nls.lm(par = starts1, fn = optFun,
                probfun = pgamma, quantiles = x[cros.ind.1],
                prob = pX[cros.ind.1],
                control = nls.lm.control(maxiter = maxiter,
                                    ftol = ftol, maxfev = maxfev, ptol = ptol)),
                silent = TRUE)
            if (inherits(FIT1, "try-error")) {
                FIT1 <- try(nls.lm(par = starts, fn = optFun,
                                    probfun = pgamma, quantiles = x[cros.ind.1],
                                    prob = pX[cros.ind.1],
                                    control = nls.lm.control(maxiter = maxiter,
                                                            ftol = ftol,
                                                            maxfev = maxfev,
                                                            ptol = ptol)),
                            silent = TRUE)
            }
        }

        if (length(starts1) > 2) {
            FIT2 <- try(nls.lm(par = starts1, fn = optFun, probfun = pgamma3p,
                                quantiles = x[cros.ind.2],
                                prob = pX[cros.ind.2],
                                control = nls.lm.control(maxiter = maxiter,
                                                        ftol = ftol,
                                                        maxfev = maxfev,
                                                        ptol = ptol)),
                        silent = TRUE)
            if (inherits(FIT2, "try-error")) {
                starts <- c(shape = shape, scale = scale, mu = mu[1])
                FIT2 <- try(nls.lm(par = starts, fn = optFun,
                                probfun = pgamma3p, quantiles = x[cros.ind.2],
                                prob = pX[cros.ind.2],
                                control = nls.lm.control(maxiter = maxiter,
                                                        ftol = ftol,
                                                        maxfev = maxfev,
                                                        ptol = ptol)),
                            silent = TRUE)
            }
        } else {
            FIT2 <- try(nls.lm(par = starts1, fn = optFun,
                                probfun = pgamma, quantiles = x[cros.ind.2],
                                prob = pX[cros.ind.2],
                                control = nls.lm.control(maxiter = maxiter,
                                                        ftol = ftol,
                                                        maxfev = maxfev,
                                                        ptol = ptol)),
                        silent = TRUE)
            if (inherits(FIT2, "try-error")) {
                starts <- c(shape = shape, scale = scale)
                FIT2 <- try(nls.lm(par = starts, fn = optFun,
                                    probfun = pgamma, quantiles = x[cros.ind.2],
                                    prob = pX[cros.ind.2],
                                    control = nls.lm.control(maxiter = maxiter,
                                                        ftol = ftol,
                                                        maxfev = maxfev,
                                                        ptol = ptol)),
                            silent = TRUE)
            }
        }

        if (inherits(FIT1, "try-error") && inherits(FIT2, "try-error"))
            R.cross.FIT <- 0 else {
            ## prediction using model 1
            p.FIT1 <- getPreds(coef(FIT1), x[cros.ind.2])
            R.FIT1 <- cor(p.FIT1, pX[cros.ind.2], use = "complete.obs")
            ## prediction using model 2
            p.FIT2 <- getPreds(coef(FIT2), x[cros.ind.1])
            R.FIT2 <- cor(p.FIT2, pX[cros.ind.1], use = "complete.obs")

            R.cross.FIT <- (length(p.FIT1) * R.FIT1 +
                length(p.FIT2) * R.FIT2)/(length(p.FIT1) +
                length(p.FIT2))
        }
        res <- pX - getPreds(coef(FIT), x)

        if (m > 2) {
            COV = try(vcov(FIT), silent = TRUE)
            if (inherits(COV, "try-error"))
                COV = matrix(NA, nrow = 3, ncol = 3)

            stats <- data.frame(summary(FIT)$coefficients,
                Adj.R.Square = c(Adj.R.Square, "", ""), rho = c(rho, "", ""),
                R.Cross.val = c(R.cross.FIT, "", ""),
                DEV = c(deviance(FIT), "", ""),
                AIC = c(AICmodel(FIT, residuals = res, np = 4), "", ""),
                BIC = c(BICmodel(FIT, residuals = res, np = 4), "", ""),
                COV = COV, n = c(N - 3, n - 3, n - 3),
                model = c("Gamma3P", "", ""),
                stringsAsFactors = FALSE)
        } else {
            COV = try(vcov(FIT), silent = TRUE)
            if (inherits(COV, "try-error"))
                COV = matrix(NA, nrow = 2, ncol = 2)
            stats <- data.frame(summary(FIT)$coefficients,
                Adj.R.Square = c(Adj.R.Square, ""),
                rho = c(rho, ""), R.Cross.val = c(R.cross.FIT, ""),
                DEV = c(deviance(FIT), ""),
                AIC = c(AICmodel(FIT, residuals = res, np = 3), ""),
                BIC = c(BICmodel(FIT, residuals = res, np = 3), ""), COV = COV,
                COV.mu = c(NA, NA), n = c(N - 2, n - 2),
                model = c("Gamma2P", ""), stringsAsFactors = FALSE)
        }
    } else {
        warning(paste("Data did not fit to the model.",
            "Returning empty coefficient table."))
        stats <- data.frame(NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    }

    colnames(stats) <- c("Estimate", "Std. Error",
        "t value", "Pr(>|t|))", "Adj.R.Square", "rho",
        "R.Cross.val", "DEV", "AIC", "BIC", "COV.shape",
        "COV.scale", "COV.mu", "df", "model")
    if (nlms)
        stats <- list(stats, nlms = FIT)
    return(stats)
}

