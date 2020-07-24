#' @rdname meth_levels
#' @title Compute methylation levels
#' @description This function computes the
#' @param GR,x A \code{\link[GenomicRanges]{GRanges-class}} object
#' (\strong{\emph{GR}}) or 'data.frame' (\strong{\emph{x}}) with a matrix of
#' counts in the meta-columns (methylated mC and unmethylated uC cytosines) or a
#' list of \code{\link[GenomicRanges]{GRanges-class}} objects.
#' @param columns Vector of one or two integer numbers denoting the indexes of
#' the columns where the methylated and unmethylated read counts are found.
#' Unless specified in the parameter 'columns', the methylation counts must be
#' given in the first four columns: 'mC1' and 'uC1' methylated and unmethylated
#' counts for control sample, and 'mC2' and 'uC2' methylated and unmethylated
#' counts for treatment sample, respectively.
#' @param Bayesian logical(1). Whether to perform the estimations based on
#' posterior estimations of methylation levels.
#' @param min.coverage An integer or an integer vector of length 2. Cytosine
#' sites where the coverage in both samples, 'x' and 'y', are less than
#' min.coverage' are discarded. The cytosine site is preserved, however, if the
#' coverage is greater than 'min.coverage' in at least one sample. If
#' 'min.coverage' is an integer vector, then the corresponding min coverage is
#' applied to each sample.
#' @param tv logical(1). Whether to compute the total variation distance at each
#' cytosine site. That is, the difference of methylation levels.
#' @param bay.tv logical(1). Whether to compute the total variation distance at
#' each cytosine site based on Bayesian estimation of methylation levels.
#' @param preserve.dt logical(1). Option of whether to preserve all
#' the metadata from the original 'data.frame' or
#' \code{\link[GenomicRanges]{GRanges-class}} object.
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to use,
#' i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS). These parameters will be passed to
#' \code{\link{uniqueGRanges}}.
#' @param verbose if TRUE, prints the function log to stdout

#' @aliases meth_levels
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com})
#' @examples
#' ## The read count data are created
#' num.samples <- 250
#' s <- 1:num.samples
#' gr <- data.frame(chr = 'chr1', start = s, end = s,
#'                 strand = sample(c("+", "-"), num.samples, replace = TRUE),
#'                 mCc = rnbinom(size = num.samples, mu = 4, n = 500),
#'                 uCc = rnbinom(size = num.samples, mu = 4, n = 500),
#'                 mCt = rnbinom(size = num.samples, mu = 4, n = 500),
#'                 uCt = rnbinom(size = num.samples, mu = 4, n = 500))
#'
#' gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE)
#'
#' gr <- meth_levels(GR = gr,
#'                 columns = c(mC1 = 1, uC1 = 2,
#'                             mC2 = 3, uC2 = 4),
#'                 preserve.dt = TRUE,
#'                 Bayesian = TRUE, tv = TRUE, bay.tv = TRUE,
#'                 num.cores = 1)
setGeneric("meth_levels",
           function(GR,
                    x,
                    columns = c(mC1 = 1, uC1 = 2, mC2 = NULL, uC2 = NULL),
                    Bayesian = FALSE,
                    min.coverage = 4,
                    tv = FALSE,
                    bay.tv = FALSE,
                    filter = FALSE,
                    preserve.dt = FALSE,
                    num.cores = 1,
                    tasks = 0L,
                    verbose = TRUE, ...) standardGeneric("meth_levels"))

#' @aliases meth_levels
#' @rdname meth_levels
#' @export
setMethod("meth_levels", signature(x = "data.frame"),
          function(GR = NULL,
                   x,
                   columns = c(mC1 = 1, uC1 = 2, mC2 = NULL, uC2 = NULL),
                   Bayesian = FALSE,
                   min.coverage = 4,
                   tv = FALSE,
                   bay.tv = FALSE,
                   filter = FALSE,
                   preserve.dt = FALSE,
                   verbose = TRUE, ...) {

    colm <- is.element(c("mC1", "uC1", "mC2", "uC2"), names(columns))

    if (preserve.dt) y <- x
    x <- data.matrix(x[, columns])
    if (ncol(x) < 2) {
              stop("If counts are provided, then at least two columns must ",
                    "carry numerical values")
    }

    if (filter) {
        r1 <- rowSums(x[, 1:2])
        ind1 <- which(r1 > min.coverage)

        if (colm[3] && colm[4]) {
            r2 <- rowSums(x[, 3:4])
            ind2 <- which(r2 > min.coverage)
            ind <- union(ind1, ind2)
            rm(ind1, ind2, r1, r2)
            x <- x[ind, ]
            if (preserve.dt) y <- y[ind, ]
        }
        else {
            x <- x[ind1, ]
            if (preserve.dt) y <- y[ind1, ]
        }
    }

    if (colm[3] && colm[4]) {
        if (tv || !Bayesian) {
            p1 <- meth.level(x[, 1:2], Bayesian = FALSE, verbose = verbose)
            p2 <- meth.level(x[, 3:4], Bayesian = FALSE, verbose = verbose)
            if (tv) TV <- p2 - p1
        }
        if (Bayesian) {
            p1 <- meth.level(x[, 1:2], Bayesian = Bayesian, verbose = verbose)
            p2 <- meth.level(x[, 3:4], Bayesian = Bayesian, verbose = verbose)
            if (bay.tv) bay.TV <- p2 - p1
        }
        if (preserve.dt) x <- data.frame(y, p1, p2)
        else x <- data.frame(p1, p2)
    } else {
        p1 <- meth.level(x[, 1:2], Bayesian = Bayesian, verbose = verbose)
        if (preserve.dt) x <- data.frame(y, p1)
        else x <- data.frame(p1)
    }

    if (tv && colm[3] && colm[4]) x$TV <- TV
    if (bay.tv && colm[3] && colm[4]) x$bay.TV <- bay.TV

    return(x)
})

#' @aliases meth_levels
#' @rdname meth_levels
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @export
setMethod("meth_levels", signature(GR = "GRanges"),
          function(
            GR,
            x = NULL,
            columns = c(mC1 = 1, uC1 = 2, mC2 = NULL, uC2 = NULL),
            Bayesian = FALSE,
            min.coverage = 4,
            tv = FALSE,
            bay.tv = FALSE,
            filter = FALSE,
            preserve.dt = FALSE,
            verbose = TRUE, ...) {

            colm <- is.element(c("mC1", "uC1", "mC2", "uC2"), names(columns))
            x <- data.matrix(mcols(GR)[, columns])
            if (filter) {
                r1 <- rowSums(x[, 1:2])
                ind1 <- which(r1 > min.coverage)

                if (colm[3] && colm[4]) {
                    r2 <- rowSums(x[, 3:4])
                    ind2 <- which(r2 > min.coverage)
                    ind <- union(ind1, ind2)
                    rm(ind1, ind2, r1, r2, x); gc()
                    GR <- GR[ ind ]
                } else GR <- GR[ ind1 ]
            }

            if (preserve.dt) {
                mcols(GR) <- data.frame(mcols(GR),
                                    meth_levels(GR = NULL,
                                        x = data.frame(mcols(GR)[, columns]),
                                        Bayesian = Bayesian,
                                        columns = columns,
                                        min.coverage = 4,
                                        tv = tv,
                                        bay.tv = bay.tv,
                                        filter = FALSE,
                                        preserve.dt = FALSE,
                                        verbose = verbose, ...))
            }
            else
                mcols(GR) <- meth_levels(GR = NULL,
                                        x = data.frame(mcols(GR)[, columns]),
                                        Bayesian = Bayesian,
                                        columns = columns,
                                        min.coverage = 4,
                                        tv = tv,
                                        bay.tv = bay.tv,
                                        filter = FALSE,
                                        preserve.dt = FALSE,
                                        verbose = verbose, ...)
            return(GR)
})

#' @aliases meth_levels
#' @rdname meth_levels
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @export
setMethod("meth_levels", signature(GR = "list"),
          function(
            GR,
            x = NULL,
            columns = c(mC1 = 1, uC1 = 2, mC2 = 0, uC2 = 0),
            Bayesian = FALSE,
            min.coverage = 4,
            tv = FALSE,
            bay.tv = FALSE,
            filter = FALSE,
            preserve.dt = FALSE,
            num.cores = 1,
            tasks = 0L,
            verbose = TRUE, ...) {

            progressbar <- FALSE

            if (verbose) progressbar <- TRUE
            if (Sys.info()["sysname"] == "Linux") {
              bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                        progressbar = progressbar)
            } else {
              bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                   progressbar = progressbar)
            }
            GR <- bplapply(GR, meth_levels, Bayesian = Bayesian,
                           columns = columns, min.coverage = min.coverage,
                           tv = tv, bay.tv = bay.tv, filter = filter,
                           preserve.dt = preserve.dt, num.cores = num.cores,
                           tasks = tasks, verbose, BPPARAM = bpparam)
            return(GR)
})


### ======================== Auxiliary function ============================

meth.level <- function(x, Bayesian, verbose) {

    n <- rowSums(x)

    if (Bayesian) {
        if (nrow(x) < 10)
            stop(paste("*** You must provide at least 10 cytosine sites ",
                        "to apply a Bayessian approach \n",
                        "using beta distributed priors"))
        if (verbose)
            cat("*** Estimating betaBinomial-posteriors... \n")

        ## Naive distribution q (methylation levels).  In a
        ## Bayesian framework with uniform priors, the
        ## methylation level can be defined as: meth_levels =
        ## ( mC + 1 )/( mC + uC + 2 ).
        q <- (x[, 1] + 1)/(n + 2)

        ## The shape parameters estimated with 'nlm'
        beta <- .estimateBetaDist(q)
        ## Assuming beta priors
        n[n == 0] <- 2
        p <- .betaBinPosteriors(x[, 1], n, a = beta[1], b = beta[2])
        return(p)
    }

    return(x[, 1]/n)
}
