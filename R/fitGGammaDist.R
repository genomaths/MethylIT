#' @rdname fitGGammaDist
#'
#' @title Nonlinear fit of Generalized Gamma CDF (GGamma)
#' @description This function performs the nonlinear fit of GGamma CDF of a
#'     variable x
#' @details The script algorithm tries to fit the three-parameter GGamma CDF
#'     ("GGamma3P") or the four-parameter GGamma ("GGamma4P") using a
#'     modification of Levenberg-Marquardt algorithm implemented in function
#'     'nls.lm' from 'minpack.lm' package that is used to perform the nonlinear
#'     fit. Cross-validations for the nonlinear regressions (R.Cross.val) were
#'     performed in each methylome as described in reference [1]. In addition,
#'     Stein's formula for adjusted R squared (rho) was used as an estimator of
#'     the average cross-validation predictive power [1].
#'
#'     If the number of values to fit is >10^6, the fitting to a GGamma CDF
#'     would be a time consuming task. To reduce the computational time, the
#'     option summarized.data' can be set 'TRUE'. If npoint != NULL, the
#'     original variable values are summarized into 'npoint' bins and their
#'     midpoints are used as the new predictors. In this case, only the
#'     goodness-of-fit indicators AIC and R.Cross.val are estimated based on all
#'     the original variable x values.
#'
#' @param x numerical vector
#' @param parameter.values initial parameter values for the nonlinear fit. If
#'     the locator paramter is included (mu != 0), this must be given as
#'     parameter.values = list(alpha = 'value', scale = 'value', mu = 'value',
#'     psi = 'value') or if mu = 0, as: parameter.values =list(alpha = 'value',
#'     scale = 'value', psi = 'value'). If not provided, then an initial guess
#'     is provided.
#' @param location.par whether to consider the fitting to generalized gamma
#'     distribution (GGamma) including the location parameter, i.e., a GGamma
#'     with four parameters (GGamam4P).
#' @param sample.size size of the sample.
#' @param npoints number of points used in the fit. If the number of points if
#'     greater than 10^6, then the fit is automatically set to npoints = 999999.
#'     However, the reported values for R.Cross.val, AIC, and BIC are
#'     computed taking into account the whole set of points.
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
#' @param ... arguments passed to or from other methods.
#'
#' @return Model table with coefficients and goodness-of-fit results:
#'     Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
#'     coefficient covariance matrix. If \strong{nlms = TRUE}, then a list with
#'     nonlinear model object \code{\link[minpack.lm]{nls.lm}} is returned.
#'
#' @references 1. Stevens JP. Applied Multivariate Statistics for the Social
#'     Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 06/03/2016
#'
#' @examples
#' set.seed(123)
#' x <- rggamma(2000, alpha = 1.03, psi = 0.75, scale = 1.1)
#' fitGGammaDist(x)
#'
#' @importFrom stats var stepfun as.formula coef AIC pweibull BIC vcov knots
#'   deviance BIC cor quantile
#' @importFrom nls2 nls2
#' @importFrom utils head
#' @importFrom stats ecdf nls.control coef
#' @importFrom minpack.lm nlsLM nls.lm nls.lm.control
#' @export

fitGGammaDist <- function(x, parameter.values, location.par = FALSE,
                          sample.size=20, npoints=NULL, maxiter=1024,
                          ftol=1e-12, ptol=1e-12, maxfev = 1e+5,
                          nlms = FALSE, verbose=TRUE, ...) {
  ind <- which(x > 0)
  if (length(ind) > sample.size) {
       x <- x[ind]
       Fy <- ecdf(x)
  } else stop("*** Sample size is lower than the set minimun size: ",
              sample.size)

   if (location.par) {
       formula <- as.formula("Y ~ pggamma(X, alpha, scale, mu, psi)")
       getPreds <- function(par, q) pggamma(q, alpha = par[1], scale = par[2],
                                         mu = par[3], psi = par[4])
   } else {
      formula <- as.formula("Y ~ pggamma(X, alpha, scale, mu = 0, psi)")
      getPreds <- function(par, q) pggamma(q, alpha = par[1], scale = par[2],
                                         psi = par[3])
   }

   optFun <- function(par, probfun, quantiles, prob, eval = FALSE) {
       START <- as.list(par)
       START$q <- quantiles
       EVAL <- try(do.call(probfun, START), silent = TRUE)
       if (inherits(EVAL, "try-error")) return(NA)
       EVAL[is.nan(EVAL)] <- 0
       RSS <- (prob - EVAL) ^ 2
       if (eval) {
           return(EVAL)
       } else return(RSS)
   }

   N <- length(x)
   if (N > 10 ^ 6) npoints = 999999

   if (!is.null(npoints) && npoints < N) {
       DENS <- hist(x, breaks = npoints, plot = FALSE)
       X <- DENS$mids
       Y <- Fy(X)
       probFun <- pggamma
       N <- length(X)

       if (verbose && location.par) {
           message(paste0("*** Trying nonlinear fit to a generalized 4P Gamma ",
                       "distribution model (summarized data: ",
                       npoints," values)..."))
       }
       if (verbose && !location.par) {
           message(paste0("*** Trying nonlinear fit to a generalized 3P Gamma ",
                       "distribution model (summarized data: ", npoints,
                       " values)..."))
       }
   } else {
       if (verbose && location.par) {
           message(paste0("*** Trying nonlinear fit to a generalized 4P Gamma ",
                       "distribution model ", npoints," values..."))
       }
       if (verbose && !location.par) {
           message(paste0("*** Trying nonlinear fit to a generalized 3P Gamma ",
                     "distribution model ", npoints," values..."))
       }
   }

   if (is.null(npoints)) {
       X = x
       Y = Fy(X)
       probFun <- pggamma
   }

   ## =============== starting parameter values =========== #
   if (missing(parameter.values)) {
       MEAN <- mean(X, na.rm = TRUE)
       VAR <- var(X, na.rm = TRUE)
       MIN <- min( X, na.rm = TRUE)

       alpha = MEAN^2/VAR
       mu = MIN
       scale = VAR/MEAN
       psi = 0.75

       if (location.par) {
           starts <- c(alpha = alpha, scale = scale, mu = mu[1], psi = 1)
       } else starts <- c(alpha = alpha, scale = scale, psi = 1)
   } else starts = parameter.values

   ## ============ END starting parameter values ========== #

   ## ==================== Fitting models ================= #
   FIT <- try(nls.lm(par = starts, fn = optFun, probfun = probFun,
                   quantiles = X, prob = Y,
                   control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                           maxfev = maxfev, ptol = 1e-12)),
               silent = TRUE)
   if (inherits( FIT, "try-error") && location.par) {
       starts <- c(alpha = alpha, scale = scale, psi = 1)
       FIT <- try(nls.lm(par = starts, fn = optFun, probfun = probFun,
                       quantiles = X, prob = Y,
                       control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = 1e-12)),
                   silent = TRUE)
   }
   if (inherits( FIT, "try-error") && location.par && probFun != dggamma) {
       DENS <- hist(x, breaks = npoints, plot = FALSE, ...)
       X <- DENS$mids
       Y <- DENS$density
       probFun <- dggamma
       FIT <- try(nls.lm(par = starts, fn = optFun, probfun = probFun,
                       quantiles = X, prob = Y,
                       control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = 1e-12)),
                   silent = TRUE)
   }
   if (inherits( FIT, "try-error") && !location.par && probFun != dggamma) {
       starts <- c(alpha = alpha, scale = scale, psi = 1)
       DENS <- hist(x, breaks = npoints, plot = FALSE, ...)
       X <- DENS$mids
       Y <- DENS$density
       probFun <- dggamma
       FIT <- try(nls.lm(par = starts, fn = optFun, probfun = probFun,
                       quantiles = X, prob = Y,
                       control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = 1e-12)),
                   silent = TRUE)
   }

   if (!inherits( FIT, "try-error" )) {
       ## **** R squares ****
       Adj.R.Square <- (1 - (deviance(FIT) / ((N - length(coef(FIT))) *
                                             var(Y, use="everything"))))
       Adj.R.Square <- ifelse(is.na(Adj.R.Square) || Adj.R.Square < 0,
                           0, Adj.R.Square)

       ## Stein adjusted R square
       rho = (1 - ((N - 2) / (N - 3)) * ((N + 1) / (N)) * (1 - Adj.R.Square))
       rho = ifelse( is.na( rho ) | rho < 0, 0, rho )

       ##--- Crossvalidation standard model for Nonlinear regression: x versus r
       if (verbose) {
           cat(paste("*** Performing nonlinear regression model ",
                   "crossvalidation...\n" ))
       }
       set.seed(123)

       cros.ind.1 <- sample.int(N, size=round(N / 2))
       cros.ind.2 <- setdiff(1:N, cros.ind.1)
       starts1 <- as.list(coef(FIT))

       FIT1 <- try(nls.lm(par=starts1, fn=optFun, probfun=probFun,
                           quantiles=X[ cros.ind.1 ], prob=Y[cros.ind.1],
                           control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)

       if (inherits( FIT1, "try-error")) {
           FIT1 <- try(nls.lm(par=starts, fn=optFun, probfun=probFun,
                           quantiles=X[ cros.ind.1 ], prob=Y[cros.ind.1],
                           control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       }

       FIT2 <- try(nls.lm(par=starts1, fn=optFun, probfun=probFun,
                       quantiles=X[cros.ind.2], prob=Y[cros.ind.2],
                       control=nls.lm.control(maxiter=maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent=TRUE)

       if (inherits( FIT2, "try-error")) {
           FIT2 <- try(nls.lm(par=starts, fn=optFun, probfun=probFun,
                           quantiles=X[ cros.ind.2 ], prob = Y[cros.ind.2],
                           control=nls.lm.control(maxiter = maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       }

       if (inherits(FIT1, "try-error") || inherits(FIT2, "try-error"))
                       R.cross.FIT <- 0
       else {
           ## ---- prediction using model 1 & 2 ----
           n <- length(x)
           Y <- Fy(x)
           cros.ind.1 <- sample.int(n, size = round(n / 2))
           cros.ind.2 <- setdiff(1:n, cros.ind.1)

           p.FIT1 <- getPreds(coef(FIT1), x[cros.ind.2])
           p.FIT2 <- getPreds(coef(FIT2), x[cros.ind.1])

           if (!all(is.na(p.FIT1)) && !all(is.na(p.FIT2))) {
               R.FIT1 <- cor(p.FIT1, Y[cros.ind.2], use="complete.obs")
               R.FIT2 <- cor(p.FIT2, Y[cros.ind.1], use="complete.obs")
               R.cross.FIT <- (length(p.FIT1) * R.FIT1 + length(p.FIT2)* R.FIT2)
               R.cross.FIT <- R.cross.FIT/(length(p.FIT1) + length(p.FIT2))
           } else R.cross.FIT <- 0
       }
       res <- Y - getPreds(coef(FIT), x)
       options(stringsAsFactors = FALSE)

       if (length( coef(FIT)) > 3) {
           COV = try(vcov(FIT), silent = TRUE)
           if (inherits(COV, "try-error")) COV = matrix(NA, nrow = 4, ncol = 4)

           stats <- data.frame(summary(FIT)$coefficients,
                           Adj.R.Square=c( Adj.R.Square, "", "", ""),
                           rho=c(rho, "", "", ""),
                           R.Cross.val=c(R.cross.FIT, "", "", ""),
                           DEV=c(deviance(FIT), "", "", "" ),
                           AIC=c(AICmodel(FIT, residuals=res, np=4), "", "",""),
                           BIC=c(BICmodel(FIT, residuals=res, np=4), "","", ""),
                           COV=COV, n=c(N, n, n, n),
                           stringsAsFactors = FALSE)
       }
       else {
           COV = try(vcov(FIT), silent = TRUE)
           if (inherits(COV, "try-error")) COV = matrix(NA, nrow = 3, ncol = 3)
           stats <- data.frame(summary(FIT)$coefficients,
                           Adj.R.Square=c(Adj.R.Square, "", ""),
                           rho=c(rho, "", ""),
                           R.Cross.val=c(R.cross.FIT, "", ""),
                           DEV=c(deviance(FIT), "", ""),
                           AIC=c(AICmodel(FIT, residuals=res, np=3), "", ""),
                           BIC=c(BICmodel(FIT, residuals=res, np=3), "", ""),
                           COV=COV,
                           COV.mu=c(NA, NA, NA),
                           n=c(N, n, n),
                           stringsAsFactors = FALSE)
       }
   } else {
       warning(paste("Data did not fit to the model.",
                   "Returning empty coefficient table."))
       stats <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                        NA, NA, NA)
   }

   colnames(stats) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|))",
                        "Adj.R.Square", "rho", "R.Cross.val", "DEV", "AIC",
                        "BIC", "COV.alpha", "COV.scale", "COV.psi",
                        "COV.mu", "N")
   if (nlms) stats <- list(stats, nlms = FIT)
   return(stats)
}

