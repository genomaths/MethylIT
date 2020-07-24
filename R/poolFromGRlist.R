#' @rdname poolFromGRlist
#'
#' @title Methylation pool from a list of GRanges objects with methylation read
#'     counts
#' @description This function will build a GRanges methylation pool from a list
#'     of GRanges objects
#' @details The list of GRanges objects (LR) provided to build a virtual
#'     methylome should be an output of the function 'readCounts2GRangesList' or
#'     at least each GRanges must have the columns named 'mC' and 'uC', for the
#'     read counts of methylated and unmethylated cytosines, respectively.
#'
#' @param LR  List of GRanges objects to build a virtual individual (methylation
#'     pool). It is assumed that the list of GRanges was obtained with
#'     \code{\link{readCounts2GRangesList}}. That is, the metacolumn from each
#'     GRanges object must contain the columns named 'mC' (number of reads
#'     signaling methylated cytosine) and 'uC' (number of reads signaling
#'     non-methylated cytosine). If more than two columns are carried on each
#'     GRanges object, then the parameter "columns" denoting the column numbers
#'     where "uC" and "mC" are located must be passed to
#'     \code{\link{uniqueGRanges}} function.
#'
#' @param stat statistic used to estimate the methylation pool: row 'mean', row
#'     'median', row 'sum', or Jacknife row mean ('jackmean') of methylated and
#'     unmethylated read counts across individuals. Notice that, for only two
#'     samples, 'jackmean' makes not sense. Since the centrality statistics are
#'     sensitive to extreme values, stat = 'sum' is an attractive option.
#'     However, in this last case, a further correction for the minimum coverage
#'     for the reference sample must be taken into account in a furhter
#'     estimation of the Hellinger divergence of methylation levels, which is
#'     explained in the detail section from the help of function
#'     \code{\link{estimateDivergence}}. A conservative option is 'mean', which
#'     will return the group centroid.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bplapply function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. Value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param prob Logic. Whether the variable for pooling is between 0 and 1 (a
#'     probability), e.g., methylation levels. If TRUE, then Fisher's
#'     transformation is applied, the row mean is computed for each cytosine
#'     site and returned in the original measurement scale between 0 and 1 by
#'     using the inverse of Fisher's transformation.
#' @param column If prob == TRUE, then the 'column' from the LR metacolumns
#'     where the prob values are found must be provided. Otherwise, column = 1L.
#' @param jstat If stat = 'jackmean', then any of the 'stat' possible values:
#'     'sum', 'mean', or 'median' can be used to compute, for each cytosine
#'     site, the Jacknife vector of the selected statistics and then to
#'     compute the corresponding mean. Default is jstat = 'sum'.
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Additional parameters for \code{\link{uniqueGRanges}} function.
#'
#' @return A GRanges object
#'
#' @examples
#' gr1 <- makeGRangesFromDataFrame(
#' data.frame(chr = 'chr1', start = 11:15, end = 11:15, strand = '*',
#' mC = 1, uC = 1:5), keep.extra.columns = TRUE)
#' gr2 <- makeGRangesFromDataFrame(
#' data.frame(chr = 'chr1', start = 11:15, end = 11:15,
#' strand = '*', mC = 1, uC = 1:5), keep.extra.columns = TRUE)
#'
#' answer <- poolFromGRlist(list(gr1, gr2), stat = 'mean', verbose = FALSE)
#'
#' @importFrom matrixStats rowMedians
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom S4Vectors mcols<-
#' @importFrom methods as
#'
#' @export
poolFromGRlist <- function(LR, stat = c("mean", "median", "jackmean", "sum"),
                        num.cores = 1, tasks = 0L, prob = FALSE, column = 1L,
                        jstat = c("sum", "mean", "median"),
                        verbose = TRUE, ...) {
    stat <- match.arg(stat)
    jstat <- match.arg(jstat)
    jstat <- eval(parse(text = jstat))

    if (verbose)
        message("*** Building a unique GRanges object from the list...")
    if (inherits(LR, "list")) {
        LR <- try(as(LR, "GRangesList"))
    }
    if (inherits(LR, "GRangesList")) {
        if (prob) {
            ## Apply Fisher transformation
            message("* prob == TRUE. Applying Fisher transformation at column ",
                column)
            LR = lapply(LR, function(GR) {
                mcols(GR) <- atanh(as.vector(mcols(GR)[, column]))
                return(GR)
            })
        }
        x0 <- uniqueGRanges(LR, num.cores = num.cores,
                            tasks = tasks, verbose = verbose, ...)
    } else {
        if (inherits(LR, "GRanges")) {
            x0 <- LR
            if (prob) {
                ## Apply Fisher transformation
                mcols(x0) <- atanh(as.vector(mcols(x0)[, column]))
            }
        } else {
            stop(paste0("Object LR is neither a list of GRanges objects,",
                " a 'GRangesList' object or a GRanges object"))
        }
    }

    if (verbose)
        message("*** Building a virtual methylome...")
    x1 <- as.matrix(mcols(x0))
    cn <- colnames(mcols(x0))
    statist <- function(x, stat) {
        x <- switch(stat,
                    sum = rowSums(x),
                    mean = round(rowMeans(x)),
                    median = round(rowMedians(x)), jackmean = round(rowJMean(x,
                stat = jstat)))
    }

    if (prob) {
        ## Apply inverse of Fisher transformation
        prob <- statist(x1, stat = "mean")
        mcols(x0) <- data.frame(prob = tanh(prob))
        x0$prob[is.na(x0$prob)] <- 0
    } else {
        idx <- grep("mC", cn)
        mC <- statist(x1[, idx], stat = stat)
        idx <- grep("uC", cn)
        uC <- statist(x1[, idx], stat = stat)
        mcols(x0) <- data.frame(mC, uC)
    }
    return(x0)
}

# ============== Auxiliary function for Jacknife mean estimation ============= #

jackstat <- function(x, stat = mean) {
    ## x is a vector
    n <- length(x)
    vls <- numeric(length = 1)
    return(vapply(seq_along(x), function(k) stat(x[-k],
        na.rm = TRUE), vls))
}

# --- x is a matrix Function 'apply' returns column
# vectors !
rowJMean <- function(x, stat = mean) colMeans(apply(x,
    1, jackstat, stat = stat), na.rm = TRUE)

## colJMean <- function(x, stat = mean)
## colMeans(apply(x, 2, jackstat, stat = stat),
## na.rm = TRUE) rowJMean1 <- function(x) apply(x,
## 1, jackstat) colJMean1 <- function(x) apply(x, 2,
## jackstat) x = matrix(1:9, nrow = 3) rowJMean1(x)
## jackstat(x[2,]) mean(jackstat(x[2,]))
## colMeans(rowJMean1(x)) rowJMean(x)



