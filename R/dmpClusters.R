#' @rdname dmpClusters
#' @name dmpClusters
#' @title DMP clustering
#' @description Given a 'pDMP' (or 'InfDiv') object carrying DMPs (methylated
#' cytosines) detected in Methyl-IT downstream analysis, function
#' \strong{\emph{'dmpClusters'}} build clusters of DMPs, which can be further
#' tested to identify differentially methylated regions (DMRs) with
#' \code{\link{countTest2}} function.
#'
#' @details Two algorithmic approaches are implemented, named: "relaxed" and
#'   "fixed.int" (see the description of parameter 'method'). The "fixed.int" is
#'   mostly addressed to find specific methylation patterns, but the price is
#'   the number of DMRs found is lower.
#'
#'  The number of DMPs reported in each cluster corresponds to the numbers of
#'  sites inside the cluster where DMPs were found in at least one of the
#'  samples (from control or from treatment). That is, \strong{dmpClusters} is
#'  only a tool to locate regions with high density of DMPs from all the
#'  samples. It does not detect DMRs. It is assumed that only DMP coordinates
#'  are given in the \emph{'GR'} object. That is, all the sites provided are
#'  considered in the computation.
#'
#' @param GR An object from \strong{\emph{'pDMP'}} class, which is returned by
#'     \code{\link{selectDIMP}} function.
#' @param maxDist maximum distance at which two reported bases sites from the
#'     same cluster can be separated. Default: \eqn{maxDist = 3}.
#' @param minNumDMPs Minimum number of DMPs inside of each cluster.
#'     Default: \eqn{minNumDMPs = 1}.
#' @param maxClustDist Clusters separated by a distance lesser than
#'   'maxClustDist' positions are merged. Default: \eqn{maxClustDist = NULL}. If
#'   \eqn{0 < maxClustDist < maxDist}, then maxClustDist will be recalculated as
#'   \eqn{maxClustDist = maxDist + 1}.
#' @param method Two different approaches are implemented to clustering DMPs:
#' \describe{
#'   \item{\strong{"relaxed":}}{DMP ranges which are separated by a distance
#'          less than \emph{'maxClustDist'} are merged and ranges with less
#'          than \emph{'minNumDMPs'} are removed.}
#'   \item{\strong{"fixed.int":}}{It will generate a partition where the
#'          distance between consecutive DMPs is not greater than
#'          \emph{'maxDist'}. Ranges with less than \emph{'minNumDMPs'} are
#'          removed. If, additionally, a value \emph{maxClustDist > 0 } is
#'          provided, then the "relaxed" approach is applied to the ranges from
#'          the first step.}
#' }
#'
#' @param chromosomes vector of characters labeling the chromosomes included in
#'     the analysis. Default: chromosomes = NULL (all chromosomes are included).
#' @param ignore.strand Same as in
#'  \code{\link[GenomicRanges]{findOverlaps-methods}.}
#' @param num.cores,tasks integer(1). The number of cores to use, i.e. at most
#'     how many child processes will be run simultaneously (see
#'     \code{\link[BiocParallel]{bplapply}} function from BiocParallel
#'     package).The number of tasks per job. value must be a scalar integer >=
#'     0L (see MulticoreParam from BiocParallel package).
#' @param verbose if TRUE, prints the function log to stdout.
#' @param ... Further parameters for uniqueGRanges function.
#'
#' @return A GRanges object with the numbers of positions inside each cluster,
#'     where DMPs were reported in at least one of the samples.
#' @references
#' \enumerate{
#'     \item Sanchez R, Mackenzie SA (2016) Information Thermodynamics of
#'         Cytosine DNA Methylation. PLOS ONE 11(3): e0150427.
#'         https://doi.org/10.1371/journal.pone.0150427
#' }
#'
#' @importFrom GenomicRanges findOverlaps makeGRangesFromDataFrame GRanges
#' @importFrom IRanges gaps IRanges disjoin
#' @importFrom BiocGenerics unique start end strand strand<-
#' @importFrom S4Vectors subjectHits  mcols mcols<-
#' @importFrom BiocParallel MulticoreParam bpmapply SnowParam
#' @importFrom data.table data.table
#' @importFrom IRanges reduce
#' @importFrom methods show
#' @importFrom stats na.omit
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @export
#' @examples
#' ## Get a dataset of dmps from the package
#' data(dmps)
#'
#' ## Build clusters of DMPs taking into account the DNA strand
#' x1 = dmpClusters2(GR = dmps, maxDist = 7, minNumDMPs = 6,
#'                  method = "fixed.int", ignore.strand = FALSE,
#'                  verbose = FALSE)
#' data.frame(x1)
#'
#'\donttest{
#' ## Build clusters of DMPs ignoring DNA strand and maxClustDist = 7
#' x2 = dmpClusters2(GR = dmps, maxDist = 7, minNumDMPs = 6,
#'                   maxClustDist = 7, method = "fixed.int",
#'                   num.cores=2L, ignore.strand = TRUE,
#'                   verbose = FALSE)
#' DataFrame(data.frame(x2))
#'
#' ## The relaxed approach with method = "relaxed"
#' x3 = dmpClusters2(GR = dmps, minNumDMPs = 6, method = "relaxed",
#'                   maxClustDist = 10, ignore.strand = TRUE,
#'                   verbose = FALSE)
#' DataFrame(data.frame(x3))
#'
#' ## ==== Setting up the experiment design to test for DMRs ===
#' nams <- names(dmps)
#' dmps_at_clusters <- getDMPatRegions(GR = dmps, regions = x3,
#'                                     ignore.strand = TRUE)
#' dmps_at_clusters <- uniqueGRanges(dmps_at_clusters, columns = 2L,
#'                                 ignore.strand = TRUE, type = 'equal',
#'                                 verbose = FALSE)
#' colnames(mcols(dmps_at_clusters)) <- nams
#'
#' colData <- data.frame(condition = factor(c('CT', 'CT', 'CT',
#'                                         'TT', 'TT', 'TT'),
#'                                         levels = c('CT', 'TT')),
#'                     nams, row.names = 2)
#'
#' ## Build a RangedGlmDataSet object
#' ds <- glmDataSet(GR = dmps_at_clusters, colData = colData)
#'
#' ## ================ Testing to detect DMRs ===========
#' dmrs <- countTest2(DS = ds, num.cores = 4L, minCountPerIndv = 4,
#'                 maxGrpCV = c(1, 1), Minlog2FC = 0.5,
#'                 CountPerBp = 0.001, test = 'LRT',
#'                 verbose = TRUE)
#' dmrs
#'}
#'

dmpClusters <- function(GR, maxDist = 3, minNumDMPs = 1,
                        maxClustDist = NULL,
                        method = c("relaxed", "fixed.int"),
                        chromosomes = NULL, ignore.strand = TRUE,
                        num.cores = 1L, tasks = 0L,
                        verbose = TRUE, ...) {
    validateClass(GR)
    if (!inherits(GR, "pDMP") || !inherits(GR, "InfDiv"))
        stop("*** GR object must inherits from 'pDMP' class",
            " which is returned by calling 'selectDMP' function. \n",
            "Or it must inherits from 'InfDiv' class")
    GR <- structure(GR, class = "list")
    GR <- GRangesList(GR)
    if (length(GR) == 1) GR <- unlist(GR, use.names = FALSE)

    method <- match.arg(method)

    if (method == "fixed.int") {

        GR <- meth_status(gr = GR, chromosomes = chromosomes,
                          ignore.strand = ignore.strand,
                          num.cores = num.cores, tasks = tasks,
                          verbose = verbose, ...)
        GR <- unlist(GR, use.names = FALSE)
        GR <- GR[GR$signal > 0, ]
    }
    if (method == "relaxed") {
        if (is.null(maxClustDist))
            stop('\nIf method = "relaxed", then a value for "maxClustDist"',
                 ' must be provided')
        # if gr is a GRangesList object
        if (inherits(GR, "GRangesList")) {
            if (verbose)
                message("* Building a unique GRange objects ... \n")
            GR <- uniqueGRanges(GR, chromosomes = chromosomes, columns = 9L,
                                missing = 0, type = "equal",
                                ignore.strand = ignore.strand,
                                num.cores = num.cores, tasks = tasks,
                                verbose = FALSE, ...)
        }
        mcols(GR) <- rowSums(as.matrix(mcols(GR)), na.rm = TRUE)
        colnames(mcols(GR)) <- "signal"
    }

    if (verbose) message("\n *** Building clusters ...")

    GR = sortBySeqnameAndStart(GR)
    signals <- GR[, "signal"]
    signals$cluster <- NULL

    # To find the word frameworks 'cl'
    if (method == "fixed.int") {
        GR <- reduce(GR, min.gapwidth = maxDist + 1,
                    ignore.strand = ignore.strand)
        GR <- getDMPatRegions(GR = signals, regions = GR,
                            ignore.strand = ignore.strand)
        if (minNumDMPs > 0) GR <- GR[GR$dmps >= minNumDMPs]
    }

    if (verbose) message("*** Counting DMPs in clusters ...")

    if (!is.null(maxClustDist) && is.numeric(maxClustDist)) {

        if (maxClustDist < maxDist) maxClustDist <- maxDist + 1
        if (verbose)
            message("\n *** Joining clusters separated by a distance < ",
                    maxClustDist, " bp")
        GR <- reduce(GR, min.gapwidth = maxClustDist + 1,
                     ignore.strand = ignore.strand)

        GR <- getDMPatRegions(GR = signals, regions = GR,
                            ignore.strand = ignore.strand)
        if (minNumDMPs > 0) GR <- GR[GR$dmps >= minNumDMPs]
    }
    return(GR[, "dmps"])
}

# ========================= Auxiliary functions ========================= #
#
## Methylation status at each cytosine base
meth_status <- function(gr, chromosomes = NULL, ignore.strand = TRUE,
                    num.cores = 1L, tasks = 0L, verbose = TRUE, ...) {
    ## Set parallel computation
    progressbar = FALSE
    if (verbose)
        progressbar = TRUE

    if (num.cores > 1) {
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(workers = num.cores,
                                    tasks = tasks, progressbar = progressbar)
        } else bpparam <- SnowParam(workers = num.cores,
                                    type = "SOCK", progressbar = progressbar)
    }

    status <- function(y, status, CHR, verbose) {
        if (verbose)
            message("\n *** Processing chromosome:  ",
                CHR, "\n")
        seqlevels(y, pruning.mode = "coarse") <- CHR
        min.start <- min(start(y))
        max.end <- max(end(y))
        if (status == 0) {
            y <- disjoin(gaps(y))
            starts <- start(y)
            ends <- end(y)
            idx <- which(starts > min.start)
            starts <- starts[idx]
            ends <- ends[idx]
            idx <- which(ends < max.end)
            starts <- starts[idx]
            ends <- ends[idx]

            if (num.cores > 1) {
                post <- bpmapply(function(s, e) s:e, starts,
                                ends, USE.NAMES = FALSE, BPPARAM = bpparam)
            } else {
                post <- mapply(function(s, e) s:e, starts,
                                ends, USE.NAMES = FALSE)
            }

            post <- sort(unique(unlist(post)))
            y <- GRanges(seqnames = CHR,
                        ranges = IRanges(start = post, end = post))
        }

        if (status == 1) {
            starts <- start(y)
            ends <- end(y)
            if (!all(starts == ends)) {
                post <- mapply(function(s, e) s:e, starts, ends,
                                USE.NAMES = FALSE)
                post <- sort(unique(unlist(post)))
                y <- GRanges(seqnames = CHR, ranges = IRanges(start = post,
                                                            end = post))
            }
        }
        mcols(y) <- data.frame(signal = status)
        return(y)
    }

    # if gr is a GRangesList object
    if (inherits(gr, "GRangesList")) {
        if (verbose)
            message("* Building a unique GRange objects ... \n")
        gr <- uniqueGRanges(gr, chromosomes = chromosomes, missing = 0,
                            type = "equal", ignore.strand = ignore.strand,
                            num.cores = num.cores, tasks = tasks,
                            verbose = FALSE, ...)
        mcols(gr) <- NULL
    }

    if (is.null(chromosomes))
        CHRs <- unique(as.character(seqnames(gr))) else {
        CHRs <- chromosomes
    }

    kCHR <- length(CHRs)
    if (verbose)
        cat("* Computing the genomic signal for status '1' ...\n")
    r <- lapply(CHRs, function(chr) status(y = gr,
                                status = 1, CHR = chr, verbose = verbose))
    if (verbose)
        cat("* Computing the genomic signal for status '0' ...\n")
    y <- lapply(CHRs, function(chr) status(y = gr,
                                status = 0, CHR = chr, verbose = verbose))
    if (verbose)
        cat("* Building the methylation signal ...\n")

    if (num.cores > 1) {
        gr <- bpmapply(function(x, y) sortBySeqnameAndStart(c(x, y)),
                        r, y, BPPARAM = bpparam)
    }
    else gr <- mapply(function(x, y) sortBySeqnameAndStart(c(x, y)), r, y)

    rm(r, y)
    gc()
    gr <- lapply(gr, unique)
    names(gr) <- CHRs
    return(GRangesList(gr))
}

