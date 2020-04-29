#' @rdname nonlinearFitDist
#'
#' @title Nonlinear fit of Information divergences distribution
#' @description A wrapper to call functions 'Weibull3P' and
#'     'fitGGammaDist' to operate on list of GRanges.
#' @details The algorithm prepares the information divergence variable
#' to try fitting Weibull or generalized gamma distribution model to the data.
#' If Weibull distribution is selected (default: 'Weibull'), function
#' 'Weibull2P' first attempts fitting to the two-parameter Weibull CDF
#' (Weibull2P). If Weibull2P did not fit, then the algorithm will try to fit
#' Weibull3P. The Levenberg-Marquardt algorithm implemented in R package
#' 'minpack.lm' is used to perform the nonlinear fit. Cross-validations for the
#' nonlinear regressions (R.Cross.val) are performed in each methylome as
#' described in reference (1-2). In addition, Stein's formula for adjusted R
#' squared (rho) is used as an estimator of the average cross-validation
#' predictive power (2).
#'
#' If 'GGamma3P' is selected the call to function 'fitGGammaDist' permits the
#' fitting to the three-parameter GGamma CDF ('GGamma3P'). The fit to the
#' four-parameter GGamma ('GGamma4P') is also available. GGamma distribution are
#' fitted using a modification of Levenberg-Marquardt algorithm implemented in
#' function 'nls.lm' from the 'minpack.lm' R package. Notice that the fit to
#' GGamma distribution is computationally time consuming (see ?fitGGammaDist for
#' additional information).
#'
#' @param LR A list of GRanges objects with information divergence values in
#'     their meta-columns.
#'@param column An integer number denoting the index of the GRanges column
#'     where the information divergence is given. Default column = 1
#' @param dist.name Name(s) of the distribution to fit. A single character
#'     string or character vector naming the distribution(s): 'Weibull'
#'     (default), gamma with three-parameter (Gamma3P), gamma with two-parameter
#'     (Gamma2P), generalized gamma with three-parameter ('GGamma3P') or
#'     four-parameter ('GGamma4P'), and Log-Normal (LogNorm).
#' @param sample.size size of the sample
#' @param location.par whether to consider the fitting to generalized gamma
#'     distribution (GGamma) including the location parameter, i.e., a GGamma
#'     with four parameters (GGamma4P).
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MethylIT
#'     analysis. If 'absolute = TRUE', then TV is transformed into |TV|, which
#'     is an information divergence that can be fitted to Weibull or to
#'     Generalized Gamma distribution.
#' @param npoints number of points to be used in the fit. Default is NULL.
#' @param model Optional. Only when dist.name = 'Weibull'. A selection of the
#'     distribution model, two-parameters and three-parameters Weibull model
#'     ('2P' and '3P'). Default is 'all' and the model with the best AIC
#'     criterion is reported. Alternatively, just use dist.name = 'Weibull2P' or
#'     dist.name = 'Weibull3P'.
#' @param maxiter positive integer. Termination occurs when the number of
#'     iterations reaches maxiter. Default value: 1024
#' @param tol A positive numeric value specifying the tolerance level for the
#'     relative offset convergence criterion. Default value: 1e-12,
#' @param ftol non-negative numeric. Termination occurs when both the actual
#'     and predicted relative reductions in the sum of squares are at most ftol.
#'     Therefore, ftol measures the relative error desired in the sum of
#'     squares. Default value: 1e-12
#' @param ptol non-negative numeric. Termination occurs when the relative error
#'     between two consecutive iterates is at most ptol. Therefore, ptol
#'     measures the relative error desired in the approximate solution. Default
#'     value: 1e-12,
#' @param minFactor A positive numeric value specifying the minimum step-size
#'     factor allowed on any step in the iteration. The increment is calculated
#'     with a Gauss-Newton algorithm and successively halved until the residual
#'     sum of squares has been decreased or until the step-size factor has been
#'     reduced below this limit. Default value: 10^-6.
#' @param maxfev integer; termination occurs when the number of calls to fn has
#'     reached maxfev. Note that nls.lm sets the value of maxfev to
#'     100*(length(par) + 1) if maxfev = integer(), where par is the list or
#'     vector of parameters to be optimized.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see
#'     \code{\link[BiocParallel]{bplapply}} function from BiocParallel package).
#' @param tasks integer. The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call to
#'     a function, such as bplapply, bpmapply etc. A task is the division of the
#'     X argument into chunks. When tasks == 0 (default), X is divided as evenly
#'     as possible over the number of workers (see
#'     \code{\link[BiocParallel]{MulticoreParam-class}} from BiocParallel
#'     package).
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... other parameters
#'
#' @return Model table with coeficients and goodness-of-fit results:
#'     Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
#'     coefficient covariance matrix.
#' @references
#' \enumerate{
#'     \item R. Sanchez and S. A. Mackenzie, “Information Thermodynamics of
#'           Cytosine DNA Methylation,” PLoS One, vol. 11, no. 3, p. e0150427,
#'           Mar. 2016.
#'     \item Stevens JP. Applied Multivariate Statistics for the Social
#'             Sciences. Fifth Edit. Routledge Academic; 2009.
#' }
#' @seealso \code{\link{gofReport}}
#' @author Robersy Sanchez 01/31/2018 <https://github.com/genomaths>
#'
#' @examples
#' ## Load a dataset with Hellinger Divergence of methylation levels on it.
#' data(HD)
#'
#' ## The nonlinear fit based on three-parameter GGamma distribution
#' nlms2 <- nonlinearFitDist(HD, npoints = 100, dist.name = 'GGamma3P',
#'                             verbose = FALSE)
#'
#' ## Weilbull distribution is a particular case of GGamma.
#' nlms <- nonlinearFitDist(HD, npoints = 100, verbose = FALSE)
#'
#' ## The goodness-of-fit indicators AIC suggests that the best fitted model
#' ## is obtained with GGamma distribution (in this example).
#' res <- mapply(function(m1,m2) as.numeric(c(Weibull = m1$AIC[1],
#'                                         GGamma = m2$AIC[1])),
#'                 nlms, nlms2)
#' rownames(res) <-c('Weibull', 'GGamma')
#' res
#'
#' ## However, the Cross-validations correlation coefficient is saying that
#' ## the Weibull distribution would be a little better probability
#' ## predictor.
#' res <- mapply(function(m1,m2) as.numeric(c(Weibull = m1$R.Cross.val[1],
#'                                         GGamma = m2$R.Cross.val[1])),
#'                 nlms, nlms2)
#' rownames(res) <-c('Weibull', 'GGamma')
#' res
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom GenomicRanges GRanges GRangesList mcols
#'
#' @export
nonlinearFitDist <- function(LR, column = 9, dist.name = "Weibull",
    sample.size = 20, location.par = FALSE, absolute = FALSE,
    npoints = NULL, model = "all", maxiter = 1024,
    tol = 1e-12, ftol = 1e-12, ptol = 1e-12, minFactor = 10^-6,
    num.cores = NULL, tasks = 0L, maxfev = 1e+05, verbose = TRUE,
    ...) {

    ## ---------------------- valid 'InfDiv' object ------------------------ ##
    validateClass(LR)
    ## --------------------------------------------------------------------- ##
    ##

    sn <- names(LR)
    toFit <- function(k, dist.name, column, sample.size, npoints,
        maxiter, tol, ftol, ptol, minFactor, verbose) {
        if (verbose)
            message("* Processing sample #", k, " ",
                sn[k])
        x <- LR[[k]]
        x <- mcols(x[, column])[, 1]
        if (absolute) x = abs(x)
        x <- x[x > 0]
        x <- switch(dist.name,
                    LogNorm = fitLogNormDist(x, sample.size = sample.size,
                                            npoints = npoints,
                                            maxiter = maxiter,
                                            ftol = ftol, ptol = ptol,
                                            verbose = verbose),
                    Weibull = weibull3P(x, sample.size = sample.size,
                                        model = model, npoints = npoints,
                                        maxiter = maxiter, tol = tol,
                                        ftol = ftol, ptol = ptol,
                                        minFactor = minFactor,
                                        verbose = verbose, ...),
                    Weibull2P = weibull3P(x, sample.size = sample.size,
                                        model = "2P", npoints = npoints,
                                        maxiter = maxiter, tol = tol,
                                        ftol = ftol, ptol = ptol,
                                        minFactor = minFactor,
                                        verbose = verbose, ...),
                    Weibull3P = weibull3P(x, sample.size = sample.size,
                                        model = "3P", npoints = npoints,
                                        maxiter = maxiter, tol = tol,
                                        ftol = ftol, ptol = ptol,
                                        minFactor = minFactor,
                                        verbose = verbose, ...),
                    Gamma2P = fitGammaDist(x, location.par = FALSE,
                                        sample.size = sample.size,
                                        npoints = npoints, maxiter = maxiter,
                                        ftol = ftol, ptol = ptol,
                                        verbose = verbose, ...),
                    Gamma3P = fitGammaDist(x, location.par = TRUE,
                                        sample.size = sample.size,
                                        npoints = npoints, maxiter = maxiter,
                                        ftol = ftol, ptol = ptol,
                                        verbose = verbose, ...),
                    GGamma3P = fitGGammaDist(x, location.par = FALSE,
                                            sample.size = sample.size,
                                            npoints = npoints,
                                            maxiter = maxiter, ftol = ftol,
                                            ptol = ptol,
                                            verbose = verbose, ...),
                    GGamma4P = fitGGammaDist(x, location.par = TRUE,
                                            sample.size = sample.size,
                                            npoints = npoints,
                                            maxiter = maxiter, ftol = ftol,
                                            ptol = ptol,
                                            verbose = verbose, ...))
        if (dist.name == "GGamma4P" && sum(is.na(x)) == 15) {
            x <- fitGGammaDist(x, location.par = FALSE,
                            sample.size = sample.size, npoints = npoints,
                            maxiter = maxiter, ftol = ftol, ptol = ptol,
                            verbose = verbose)
        }
        x <- structure(x, class = c("ProbDistr", "data.frame"))
        return(x)
    }

    if (length(dist.name) == 1 && length(LR) > 1)
        dist.name <- rep(dist.name, length(LR))
    if (is.null(num.cores)) {
        x <- mapply(toFit, seq_along(LR), dist.name,
            MoreArgs = list(sample.size = sample.size,
                        columnn = column, npoints = npoints,
                        maxiter = maxiter, tol = tol,
                        ftol = ftol, ptol = ptol,
                        minFactor = minFactor,
                        verbose = verbose),
            SIMPLIFY = FALSE)
    } else {
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(workers = num.cores,
                tasks = tasks)
        } else {
            bpparam <- SnowParam(workers = num.cores,
                type = "SOCK")
        }
        x <- bpmapply(toFit, seq_along(LR), dist.name,
            MoreArgs = list(sample.size = sample.size,
                        column = column, npoints = npoints,
                        maxiter = maxiter, tol = tol,
                        ftol = ftol, ptol = ptol,
                        minFactor = minFactor,
                        verbose = verbose),
            SIMPLIFY = FALSE, BPPARAM = bpparam)
    }
    names(x) <- sn
    x <- structure(x, class = c("ProbDistrList", "list"))
    return(x)
}


