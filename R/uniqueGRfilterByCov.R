#' @rdname uniqueGRfilterByCov
#' @title  Unique GRanges of methylation read counts filtered by coverages
#' @description Given two GRanges objects, samples '1' and '2', this function
#'     will filter by coverage each cytosine site from each GRanges object.
#' @details Cytosine sites with 'coverage' > 'min.coverage' in at least one of
#'     the samples are preserved. Positions with 'coverage' < 'min.coverage' in
#'     both samples, 'x' and 'y', are removed. Positions with 'coverage' >
#'     'percentile' (e.g., 99.9 percentile) are removed as well. It is expected
#'     that the columns of methylated and unmethylated counts are given.
#' @param x An object from the classes 'GRanges', 'InfDiv', or 'pDMP' with
#'     methylated and unmethylated counts in its meta-column. If the argument
#'     'y' is not given, then it is assumed that the first four columns of the
#'     GRanges metadata 'x' are counts: methylated and unmethylated counts for
#'     samples '1' and '2'.
#' @param y A GRanges object with methylated and unmethylated counts in its
#'     meta-column. Default is NULL. If x is a 'InfDiv', or 'pDMP', then 'y' is
#'     not needed, since samples '1' and '2' are the first four columns of
#'     these objects.
#' @param min.coverage An integer or an integer vector of length 2. Cytosine
#'  sites where the coverage in both samples, 'x' and 'y', are less than
#'  'min.coverage' are discarded. The cytosine site is preserved, however, if
#'  the coverage is greater than 'min.coverage' in at least one sample. If
#'  'min.coverage' is an integer vector, then the corresponding min coverage is
#'  applied to each sample.
#' @param min.meth An integer or an integer vector of length 2. Cytosine sites
#'     where the numbers of read counts of methylated cytosine in both samples,
#'     '1' and '2', are less than 'min.meth' are discarded. If 'min.meth' is an
#'     integer vector, then the corresponding min number of reads is applied to
#'     each sample.
#' @param min.umeth An integer or an integer vector of length 2. Min number of
#'     reads to consider cytosine position. Specifically cytosine positions
#'     where (uC <= min.umeth) & (mC > 0) & (mC < min.meth) hold will be
#'     removed, where mC and uC stand for the numbers of methylated and
#'     unmethylated reads. Default is min.umeth = 0.
#' @param min.sitecov An integer. The minimum total coverage. Only sites where
#'  the total coverage (cov1 + cov2) is greater than 'min.sitecov' are
#'  considered for downstream analysis, where cov1 and cov2 are the coverages
#'  for samples 1 and 2, respectively.
#' @param percentile Threshold to remove the outliers (PCR bias) from each file
#'     and all files stacked. If 'high.coverage = NULL', then the threshold
#'     \eqn{q} will be computed as:
#'     \deqn{q1 = quantile(cov1, probs=percentile)}
#'     \deqn{q2 = quantile(cov2, probs=percentile)}
#'     \deqn{q = min(q1, q2)}
#'
#'     where \eqn{cov1} and \eqn{cov2} are the coverage vectors from samples 1
#'     and 2, respectively.
#' @param high.coverage An integer for read counts. Cytosine sites having
#' higher coverage than this are discarded. Default is NULL. If
#' \strong{high.coverage} is not NULL, then the \strong{percentile} argument is
#' disregarded and \strong{high.coverage} is used as threshold to remove the
#' PCR bias.
#' @param columns Vector of integer numbers of the columns (from each GRanges
#'     meta-column) where the methylated and unmethylated counts are provided.
#'     If not provided, then the methylated and unmethylated counts are assumed
#'     to be at columns 1 and 2, respectively.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations. Default value: TRUE
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bplapply function from
#'     BiocParallel package).
#' @param tasks Integer(1). The number of tasks per job. value must be a scalar
#'  integer >= 0L. In this documentation a job is defined as a single call to a
#'  function, such as bplapply, bpmapply etc. A task is the division of the X
#'  argument into chunks. When tasks == 0 (default), X is divided as evenly as
#'  possible over the number of workers (see MulticoreParam from BiocParallel
#'  package).
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Additional parameters for 'uniqueGRanges' function.
#'
#' @return A GRanges object with the columns of methylated and unmethylated
#'   counts filtered for each cytosine position.
#'
#' @examples
#' ### Create new data
#' df1 <- data.frame(chr = 'chr1', start = 11:16, end = 11:16,
#'                 mC = c(2,10,7,9,1,10), uC = c(30,20,4,8,0,10))
#'
#' df2 <- data.frame(chr = 'chr1', start = 12:18, end = 12:18,
#'                 mC2 = 1:7, uC2 = 0:6)
#'
#' gr1 <- makeGRangesFromDataFrame(df1, keep.extra.columns = TRUE)
#' gr2 <- makeGRangesFromDataFrame(df2, keep.extra.columns = TRUE)
#'
#' ## Filtering
#' r1 <- uniqueGRfilterByCov(gr1, gr2, ignore.strand = TRUE)
#' r1
#'
#' ## Cytosine position with coordinate 15 (rows #2) can pass the
#' ## filtering conditions of min.coverage = 4 and lead to meaningless
#' ## situations with methylation levels p = 1/(1 + 0) = 1
#' r1[2]
#'
#' ## The last situation can be prevent, in this case, by setting
#' ## min.meth = 1:
#' r1 <- uniqueGRfilterByCov(gr1, gr2, min.meth = 1, ignore.strand = TRUE)
#' r1
#' @export
uniqueGRfilterByCov <- function(x, y = NULL, min.coverage = 4,
    min.meth = 0, min.umeth = 0, min.sitecov = 4, percentile = 0.9999,
    high.coverage = NULL, columns = c(mC = 1, uC = 2),
    num.cores = 1L, ignore.strand = FALSE, tasks = 0L,
    verbose = TRUE, ...) {

    if (!is.null(y)) {
        x <- x[, columns]
        y <- y[, columns]
        x <- uniqueGRanges(list(x, y), num.cores = num.cores,
            tasks = tasks, ignore.strand = ignore.strand,
            verbose = verbose, ...)
    } else {
        if (!is(x, "GRanges")) {
            # -----------------valid 'pDMP' or 'InfDiv'  object ---------------
            validateClass(x)
            # --------------------------------------------------------------- #
        }
    }

    if (length(min.coverage) == 1)
        min.coverage <- c(min.coverage, min.coverage)
    if (length(min.meth) == 1)
        min.meth <- c(min.meth, min.meth)
    if (length(min.umeth) == 1)
        min.umeth <- c(min.umeth, min.umeth)

    cov1 <- rowSums(as.matrix(mcols(x[, c(1, 2)])))
    cov2 <- rowSums(as.matrix(mcols(x[, c(3, 4)])))
    total_cov <- cov1 + cov2

    if (is.null(high.coverage)) {
        q1 <- quantile(cov1, probs = percentile)
        q2 <- quantile(cov2, probs = percentile)
        q <- min(q1, q2)
    } else q <- high.coverage

    idx1 <- which((cov1 >= min.coverage[1]) | (cov2 >= min.coverage[2]))
    if (!(length(idx1) > 0))
        stop("*** Some filtering condition from min.coverage = c(",
            paste(min.coverage, collapse = ","), ") is not hold by the sample")
    idx2 <- which((cov1 <= q) & (cov2 <= q))
    idx <- intersect(idx1, idx2)
    idx <- intersect(idx, which(total_cov >= min.sitecov))
    if (!(length(idx) > 0))
        stop("*** Some filtering condition is not hold by the sample")

    x <- x[idx]
    if (max(min.meth) > 0) {
        c1 <- mcols(x[, 1])[, 1]
        c2 <- mcols(x[, 3])[, 1]

        idx <- which((c1 >= min.meth[1]) | (c2 >= min.meth[2]))
        if (!(length(idx) > 0))
            stop("*** Some filtering condition from min.meth = c(",
                paste(min.meth, collapse = ","), ") is not hold by the sample")

        x <- x[idx]

        # To remove positions similar to, e.g., c1 = 20,
        # 40, c2 = 1 & t2 = 0, not captured on the above
        # filtering conditions (see example).
        if (length(x) == 1) {
            cov1 <- sum(as.matrix(mcols(x[, c(1, 2)])))
            cov2 <- sum(as.matrix(mcols(x[, c(3, 4)])))
        } else {
            cov1 <- rowSums(as.matrix(mcols(x[, c(1,
                2)])))
            cov2 <- rowSums(as.matrix(mcols(x[, c(3,
                4)])))
        }

        c1 <- mcols(x[, 1])[, 1]
        c2 <- mcols(x[, 3])[, 1]

        t1 <- mcols(x[, 2])[, 1]
        t2 <- mcols(x[, 4])[, 1]

        idx <- which((t1 <= min.umeth[1]) & (c1 > 0) &
            (c1 <= min.meth[1]) & cov1 < min.coverage[1])
        idx1 <- which((t2 <= min.umeth[2]) & (c2 >
            0) & (c2 <= min.meth[2]) & cov2 < min.coverage[2])
        idx <- union(idx, idx1)
        if (length(idx) > 0)
            x <- x[-idx]
    }
    return(x)
}
