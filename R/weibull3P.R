#' @rdname weibull3P
#'
#' @title Nonlinear fit of Weibull CDF
#' @description This function performs the nonlinear fit of Weibull CDF of a
#'  variable x
#' @details The script algorithm first try to fit the two-parameter Weibull CDF
#'  (Weibull2P). If Weibull2P did not fit, then the algorithm will try to fit
#'  Weibull3P. The Levenberg-Marquardt algorithm implemented in 'minpack.lm' R
#'  package is used to perform the nonlinear fit. Cross-validations for the
#'  nonlinear regressions (R.Cross.val) were performed in each methylome as
#'  described in reference (1). In addition, Stein's formula for adjusted R
#'  squared (rho) was used as an estimator of the average cross-validation
#'  predictive power (1).
#'
#' @param X numerical vector
#' @param sample.size size of the sample
#' @param model Distribution model to fit, two-parameters and three-parameters
#'  Weibull model ('Weibull2P' or simply '2P' and 'Weibull3P' or '3P). Default
#'  is 'all' and the model with the best AIC criterion is reported.
#' @param npoints number of points used in the fit
#' @param maxiter positive integer. Termination occurs when the number of
#'  iterations reaches maxiter. Default value: 1024
#' @param tol A positive numeric value specifying the tolerance level for the
#'  relative offset convergence criterion. Default value: 1e-12,
#' @param ftol non-negative numeric. Termination occurs when both the actual and
#'  predicted relative reductions in the sum of squares are at most ftol.
#'  Therefore, ftol measures the relative error desired in the sum of squares.
#'  Default value: 1e-12,
#' @param ptol non-negative numeric. Termination occurs when the relative error
#'  between two consecutive iterates is at most ptol. Therefore, ptol measures
#'  the relative error desired in the approximate solution. Default value:
#'  1e-12,
#' @param minFactor  A positive numeric value specifying the minimum step-size
#'  factor allowed on any step in the iteration. The increment is calculated
#'  with a Gauss-Newton algorithm and successively halved until the residual sum
#'  of squares has been decreased or until the step-size factor has been reduced
#'  below this limit. Default value: 10^-6
#' @param nlms Logical. Whether to return the nonlinear model object
#'  \code{\link[minpack.lm]{nls.lm}}. Default is FALSE.
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... other parameters
#'
#' @return Model table with coefficients and goodness-of-fit results:
#'  Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
#'  coefficient covariance matrix.
#'
#' @references 1. Stevens JP. Applied Multivariate Statistics for the Social
#'  Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 06/03/2016 <https://github.com/genomaths>
#'
#' @examples
#' x <- rweibull(1000, shape=0.75, scale=1)
#' weibull3P(x, sample.size = 100)
#'
#' @importFrom stats var stepfun as.formula coef AIC pweibull BIC vcov knots
#'  deviance BIC cor quantile
#' @importFrom nls2 nls2
#' @importFrom utils head
#' @importFrom stats ecdf
#' @importFrom minpack.lm nlsLM
#'
#' @export
weibull3P <- function(X, sample.size = 20, model = c("all",
    "2P", "3P", "Weibull2P", "Weibull3P"), npoints = NULL,
    maxiter = 1024, tol = 1e-12, ftol = 1e-12, ptol = 1e-12,
    minFactor = 10^-6, nlms = FALSE, verbose = TRUE, ...) {
    model <- match.arg(model)
    if (model == "Weibull2P") model <- "2P"
    if (model == "Weibull3P") model <- "3P"

    ind <- which(X > 0)
    if (length(ind) > sample.size) {
        X <- X[ind]
        Fy <- ecdf(X)
        N <- length(X)
        MIN <- min(X, na.rm = TRUE)
        MEAN <- mean(X, na.rm = TRUE)
        VAR <- var(X, na.rm = TRUE, use = "na.or.complete")

        ## To reduce the number of points to be used in the
        ## fit
        if (!is.null(npoints) && npoints < N) {
            X <- pretty(X, n = npoints)
            X <- X[X > 0]
        }
        n <- length(X)  ## size of the sample used in the computation
        pX <- Fy(X)

        if (verbose)
            cat("*** Non-linear regression \n")

        ## ------------------ Weibull 2P ---------------------
        if (model == "all" || model == "2P") {
            if (verbose)
                message("*** Trying nonlinear fit a Weibull 2P model...")

            shape <- log(-log(1 - 0.31))/log(quantile(X,0.31)/quantile(X, 0.63))
            starts <- list(shape = unname(shape),
                            scale = unname(quantile(X, 0.63)))
            formula <- as.formula("Y ~ pweibull( X, shape, scale )")

            FIT1 <- suppressWarnings(try(nlsLM(formula,
                            data = data.frame(X = X, Y = pX), start = starts,
                            control = list(maxiter = maxiter, ftol = ftol,
                                            ptol = ptol)),
                silent = TRUE))
        }

        ## ------------------ Weibull 3P ---------------------
        if (model == "all" || model == "3P") {
            starts <- list(shape = log(2), scale = (2/log(2)),  mu = 0)
            if (model == "all") {
                if (!inherits(FIT1, "try-error"))
                    starts <- list(shape = unname(coef(FIT1)[1]),
                    scale = unname(coef(FIT1)[2]),
                    mu = MIN)
            }
            if (verbose)
                message("*** Trying nonlinear fit a Weibull 3P model ...")
            formula <- as.formula("Y ~ pweibull( X - mu, shape, scale )")
            FIT2 <- suppressWarnings(try(nlsLM(formula,
                data = data.frame(X = X, Y = pX), start = starts,
                control = list(maxiter = maxiter, ftol = ftol, ptol = ptol)),
                silent = TRUE))
            if (inherits(FIT2, "try-error")) {
                starts <- list(shape = 1, scale = 1, mu = 0)
                FIT2 <- suppressWarnings(try(nlsLM(formula,
                                            data = data.frame(X = X, Y = pX),
                                            start = starts,
                                            control = list(maxiter = maxiter,
                                                        ftol = ftol,
                                                        ptol = ptol)),
                                            silent = TRUE))
            }

            if (inherits(FIT2, "try-error")) {
                starts <- list(shape = (MEAN^2/VAR), scale = (VAR/MEAN),
                                mu = MIN)
                FIT2 <- suppressWarnings(try(nlsLM(formula,
                                            data = data.frame(X = X, Y = pX),
                                            start = starts,
                                            control = list(maxiter = maxiter,
                                                            ftol = ftol,
                                                            ptol = ptol)),
                                            silent = TRUE))
            }
        }

        if (model == "all" || model == "2P") {
            if (!inherits(FIT1, "try-error")) {
                if (sum(summary(FIT1)$parameters[c(1, 2), 4] > 0.05) == 2 ||
                    coef(FIT1)[2] < 1e-07) PASS1 <- FALSE
                else PASS1 <- TRUE
            }
        } else PASS1 <- FALSE
        if (model == "all" || model == "3P") {
            if (!inherits(FIT2, "try-error")) {
                suma <- sum(summary(FIT2)$parameters[c(1, 2, 3), 4] > 0.05)
                if (suma == 3 || coef(FIT2)[2] < 1e-07 ||
                    coef(FIT2)[3] < 0) {
                    PASS2 <- FALSE
                } else PASS2 <- TRUE
            }
        } else PASS2 <- FALSE

        if (PASS1 && PASS2) {
            if (AIC(FIT1) <= AIC(FIT2)) {
                FIT <- FIT1
                if (verbose) message("*** Weibull-2P is the best fitted model")
            } else {
                FIT <- FIT2
                if (verbose) message("*** Weibull-3P is the best fitted model")
            }
        }
        if (PASS1 && !PASS2) {
            FIT <- FIT1
            if (verbose) message("*** Weibull-2P is the best fitted model")
        }
        if (!PASS1 & PASS2) {
            FIT <- FIT2
            if (verbose) message("*** Weibull-3P is the best fitted model")
        }

        if (model == "all") rm(FIT1, FIT2)
        if (model == "2P") rm(FIT1)
        if (model == "3P") rm(FIT2)

        if (PASS1 || PASS2) {
            m <- length(coef(FIT))
            ## **** R squares ****
            Adj.R.Square <- 1 - (deviance(FIT)/((n - m) *
                                                var(pX, use = "everything")))
            Adj.R.Square <- ifelse(is.na(Adj.R.Square) ||
                Adj.R.Square < 0, 0, Adj.R.Square)

            ## Stain adjusted R square
            if (m > 2)
                rho <- ((n - 1)/(n - 4)) * ((n - 2)/(n - 5)) *
                        ((n + 1)/n) else rho <- ((n - 1)/(n - 3)) *
                        ((n - 2)/(n - 4)) * ((n + 1)/n)
            rho <- 1 - rho * (1 - Adj.R.Square)
            rho <- ifelse(is.na(rho) | rho < 0, 0, rho)

            ##---- Crossvalidation standard model for Nonlinear regression:
            ##---- x versus r
            if (verbose)
                cat(paste("*** Performing nonlinear regression model",
                        "crossvalidation...\n"))
            if (m > 2) {
                getPreds <- function(par, x)
                    pweibull(x - par[3], par[1], par[2])
                formula <- as.formula("Y ~ pweibull( X - mu, shape, scale )")
            } else {
                getPreds <- function(par, x) stats::pweibull(x, par[1], par[2])
                formula <- as.formula("Y ~ pweibull( X, shape, scale )")
            }
            cros.ind.1 <- sample.int(n, size = round(n/2))
            cros.ind.2 <- setdiff(seq_len(n), cros.ind.1)
            starts <- as.list(coef(FIT))

            FIT1 <- try(nlsLM(formula, data = data.frame(X = X[cros.ind.1],
                                                        Y = pX[cros.ind.1]),
                                start = starts,
                                control = list(maxiter = maxiter, ptol = ptol)),
                        silent = TRUE)

            if (inherits(FIT1, "try-error")) {
                FIT1 <- try(nls2(formula, data = data.frame(X = X[cros.ind.1],
                                                            Y = pX[cros.ind.1]),
                                start = starts, algorithm = "plinear-brute",
                                control = list(maxiter = maxiter, tol = tol,
                                                minFactor = minFactor)),
                            silent = TRUE)
            }

            FIT2 <- try(nlsLM(formula, data = data.frame(X = X[cros.ind.2],
                                                        Y = pX[cros.ind.2]),
                            start = starts, control = list(maxiter = maxiter,
                                                            ptol = ptol)),
                        silent = TRUE)

            if (inherits(FIT2, "try-error")) {
                FIT2 <- try(nls2(formula, data = data.frame(X = X[cros.ind.2],
                                                            Y = pX[cros.ind.2]),
                                start = starts, algorithm = "plinear-brute",
                                control = list(maxiter = maxiter, tol = tol,
                                                minFactor = minFactor)),
                    silent = TRUE)
            }

            if (inherits(FIT1, "try-error") && inherits(FIT2, "try-error")) {
                R.cross.FIT <- 0
            } else {
                ## prediction using model 1
                p.FIT1 <- getPreds(coef(FIT1), X[cros.ind.2])
                R.FIT1 <- cor(p.FIT1, pX[cros.ind.2], use = "complete.obs")
                ## prediction using model 2
                p.FIT2 <- getPreds(coef(FIT2), X[cros.ind.1])
                R.FIT2 <- cor(p.FIT2, pX[cros.ind.1], use = "complete.obs")

                R.cross.FIT <- (length(p.FIT1) * R.FIT1 +
                                    length(p.FIT2) * R.FIT2)/(length(p.FIT1) +
                                                                length(p.FIT2))
            }

            fit_summary <- try(summary(FIT), silent = TRUE)

            if (!inherits(fit_summary, "try-error")) {
                if (m > 2) {
                    stats <- data.frame(fit_summary$parameters[c(1, 2, 3), ],
                                        Adj.R.Square = c(Adj.R.Square,"", ""),
                                        rho = c(rho, "", ""),
                                        R.Cross.val = c(R.cross.FIT, "", ""),
                                        DEV = c(deviance(FIT), "", ""),
                                        AIC = c(AIC(FIT), "", ""),
                                        BIC = c(BIC(FIT), "", ""),
                                        COV = vcov(FIT), n = c(N, n, n),
                                        model = c("Weibull3P", "", ""),
                                        stringsAsFactors = FALSE)
                } else {
                    stats = data.frame(fit_summary$parameters[c(1, 2),],
                                       Adj.R.Squar = c(Adj.R.Square, ""),
                                       rho = c(rho, ""),
                                       R.Cross.val = c(R.cross.FIT, ""),
                                       DEV = c(deviance(FIT), ""),
                                       AIC = c(AIC(FIT), ""),
                                       BIC = c(BIC(FIT), ""),
                                       COV = vcov(FIT), COV.mu = c(NA,NA),
                                       n = c(N, n), model = c("Weibull2P", ""),
                                       stringsAsFactors = FALSE)
                }
            } else {
                warning(paste("\nThe estimation of the regression statistics",
                    "was not possible.\n Probably the Choleski decomposition",
                    "of the covariance matrix was not possible.\n",
                    "Returning empty coefficient table."))
                stats <- data.frame(NA, NA, NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
            }
        } else {
            ## TODO log error in screen / file
            warning(paste("Data did not fit to the model.",
                "Returning empty coefficient table."))
            stats <- data.frame(NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        }
    } else {
        warning(paste("Insufficient data to fit the model.",
            "Returning empty coefficient table."))
        stats <- data.frame(NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA, NA, NA, NA)
    }
    colnames(stats) <- c("Estimate", "Std. Error",
        "t value", "Pr(>|t|))", "Adj.R.Square", "rho",
        "R.Cross.val", "DEV", "AIC", "BIC", "COV.shape",
        "COV.scale", "COV.mu", "n", "model")
    if (nlms) stats <- list(stats, nlms = FIT)
    return(stats)
}
