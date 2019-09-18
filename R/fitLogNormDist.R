#' @rdname fitLogNormDist
#'
#' @title Nonlinear fit of Log-Normal CDF (LogNorm)
#' @description This function performs the nonlinear fit of GGamma CDF of a
#'     variable x
#' @details The algorithm tries to fit the two-parameter LogNorm CDF
#'     using a modification of Levenberg-Marquardt algorithm implemented in
#'     function 'nls.lm' from 'minpack.lm' package that is used to perform the
#'     nonlinear fit. Cross-validations for the nonlinear regressions
#'     (R.Cross.val) were performed in each methylome as described in reference
#'     [1]. In addition, Stein's formula for adjusted R squared (rho) was used
#'     as an estimator of the average cross-validation predictive power [1].
#'
#'     If the number of values to fit is >10^6, the fitting to a LogNorm CDF
#'     would be a time consuming task. To reduce the computational time, the
#'     option 'summarized.data' can be set 'TRUE'. If summarized.data = TRUE,
#'     the original variable values are summarized into 'npoint' bins and their
#'     midpoints are used as the new predictors. In this case, only the
#'     goodness-of-fit indicators AIC and R.Cross.val are estimated based on all
#'     the original variable x values.
#'
#' @param x numerical vector
#' @param probability.x probability vector of x. If not provided, the values
#'     are estimated using the empirical cumulative distribution function
#'     ('ecdf') from 'stats' R package.
#' @param parameter.values initial parameter values for the nonlinear fit. If
#'     the locator paramter is included (mu != 0), this must be given as
#'     parameter.values = list(alpha = 'value', scale = 'value', mu = 'value')
#'     or if mu = 0, as: parameter.values = list(alpha = 'value',
#'     scale = 'value'). If not provided, then an initial guess is provided.
#' @param summarized.data Logic value. If TRUE (default: FALSE), summarized
#'     data based on 'npoints' are used to perform the nonlinear fit.
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
#'   reached maxfev. Note that nls.lm sets the value of maxfev to
#'   100*(length(par) + 1) if maxfev = integer(), where par is the list or
#'   vector of parameters to be optimized.
#' @param verbose if TRUE, prints the function log to stdout
#'
#' @return Model table with coefficients and goodness-of-fit results:
#'     Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
#'     coefficient covariance matrix.
#'
#' @references 1. Stevens JP. Applied Multivariate Statistics for the Social
#'     Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 04/09/2019
#'
#' @examples
#' set.seed(126)
#' x <- rlnorm(1000, meanlog = 1.03, sdlog = 2.1)
#' fitLogNormDist(x)
#'
#' @importFrom stats var stepfun as.formula coef AIC pweibull BIC vcov knots
#'   deviance BIC cor quantile
#' @importFrom nls2 nls2
#' @importFrom utils head
#' @importFrom stats ecdf nls.control plnorm
#' @importFrom minpack.lm nlsLM nls.lm nls.lm.control
#' @export

fitLogNormDist <- function(x, probability.x, parameter.values,
                           summarized.data=FALSE, sample.size=20, npoints=NULL,
                           maxiter=1024, ftol=1e-12, ptol=1e-12, maxfev = 1e+5,
                           verbose=TRUE) {
   ind <- which(x > 0)
   if (length(ind) > sample.size) {
       x <- x[ind]
       if (missing(probability.x)) Fy <- ecdf(x)
   } else stop("*** Sample size is lower than the set minimun size: ",
               sample.size)

   getPreds <- function(par, q) plnorm(q, meanlog=par[1], sdlog=par[2])

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
   if (summarized.data && is.null(npoints)) {
       npoints = min(10 ^ 6, N)
   }

   if (!is.null(npoints) && npoints < N) {
       F0 <- estimateECDF(x, npoints = npoints)
       X0 <- knots(F0)
       pX0 <- F0(X0)
       N <- length(X0)

       if (verbose && !is.null(npoints)) {
           message(paste0("*** Trying nonlinear fit to a generalized 2P Gamma ",
                       "distribution model (summarized data: ", npoints,
                       " values)..."))
       }
   }

   if (summarized.data) {X = X0; pX = pX0} else {X = x; pX = Fy(x)}

   ## =============== starting parameter values =========== #
   if (missing(parameter.values)) {
       starts <- c(meanlog = mean( log1p( X ), na.rm = TRUE ),
                   sdlog = sd(log1p( X ), na.rm = TRUE ) )
   } else starts <- parameter.values
   ## ============ END starting parameter values ========== #

   ## ==================== Fitting models ================= #
   if (verbose)
       message(paste0("*** Trying nonlinear fit to a Log-Normal Dist.",
                   "distribution model ..."))
   FIT <- try(nls.lm(par = starts, fn = optFun, probfun = plnorm,
                   quantiles = X, prob = pX,
                   control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                           maxfev = maxfev, ptol = 1e-12)),
               silent = TRUE)

   if (!inherits( FIT, "try-error" )) {
       ## **** R squares ****
       Adj.R.Square <- (1 - (deviance(FIT) / ((N - length(coef(FIT))) *
                                               var(pX, use="everything"))))
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

       FIT1 <- try(nls.lm(par = starts1, fn = optFun, probfun = plnorm,
                         quantiles=X[ cros.ind.1 ], prob=pX[cros.ind.1],
                         control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                                 maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       if (inherits(FIT1, "try-error")) {
           FIT1 <- try(nls.lm(par=starts, fn=optFun, probfun=plnorm,
                           quantiles=X[ cros.ind.1 ], prob=pX[cros.ind.1],
                           control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                       silent = TRUE)
       }

       FIT2 <- try(nls.lm(par=starts1, fn=optFun, probfun=plnorm,
                         quantiles=X[ cros.ind.2 ], prob=pX[cros.ind.2],
                         control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       if (inherits( FIT2, "try-error")) {
           starts <- c(meanlog = mean( log1p( X ), na.rm = TRUE ),
                       sdlog = sd(log1p( X ), na.rm = TRUE ) )
           FIT2 <- try(nls.lm(par=starts, fn=optFun, probfun=plnorm,
                           quantiles=X[ cros.ind.2 ], prob=pX[cros.ind.2],
                           control=nls.lm.control(maxiter=maxiter, ftol=ftol,
                                               maxfev = maxfev, ptol = ptol)),
                       silent = TRUE)
       }

       if (inherits(FIT1, "try-error") && inherits(FIT2, "try-error"))
           R.cross.FIT <- 0
       else {
           if (summarized.data) {
               n <- length(x)
               pX <- Fy(x)
               cros.ind.1 <- sample.int(n, size = round(n / 2))
               cros.ind.2 <- setdiff(1:n, cros.ind.1)
           }
           n <- length(x)

           ## prediction using model 1
           p.FIT1 <- getPreds(coef(FIT1), x[cros.ind.2])
           R.FIT1 <- cor(p.FIT1, pX[cros.ind.2], use="complete.obs")
           ## prediction using model 2
           p.FIT2 <- getPreds(coef(FIT2), x[cros.ind.1])
           R.FIT2 <- cor(p.FIT2, pX[cros.ind.1], use="complete.obs")

           R.cross.FIT <- (length(p.FIT1) * R.FIT1 + length(p.FIT2) * R.FIT2) /
                           (length(p.FIT1) + length(p.FIT2))
       }
       res <- pX - getPreds(coef(FIT), x)

       COV = try(vcov(FIT), silent = TRUE)
       if (inherits(COV, "try-error")) COV = matrix(NA, nrow = 2, ncol = 2)
       stats <- data.frame(summary(FIT)$coefficients,
                           Adj.R.Square=c(Adj.R.Square, ""),
                           rho=c(rho, ""),
                           R.Cross.val=c(R.cross.FIT, ""),
                           DEV=c(deviance(FIT), ""),
                           AIC=c(AICmodel(FIT, residuals=res, np=3), ""),
                           BIC=c(BICmodel(FIT, residuals=res, np=3), ""),
                           COV=COV,
                           n=c(N - 2, n - 2))
   } else {
       warning(paste("Data did not fit to the model.",
                   "Returning empty coefficient table."))
       stats <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
   }

   colnames(stats) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|))",
                       "Adj.R.Square", "rho", "R.Cross.val", "DEV", "AIC",
                       "BIC", "COV.mean", "COV.sd", "df")
   return(stats)
}

