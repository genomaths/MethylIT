#' @rdname dmpClusters
#' @name dmpClusters
#' @title DMP clustering
#' @description Given a 'pDMP' object carrying DMPs obtained in Methyl-IT
#'     downstream analysis, function \strong{\emph{"dmpClusters"}} build
#'     clusters of DMPs, which can be further tested to identify differentially
#'     methylated regions (DMRs) with \code{\link{countTest2}} function.
#' @details DNA base coordinates provided in the \strong{\emph{GR}} object are
#' used to build a binary string of zeros and ones . The binary string of zeros
#' and ones is used in the detection of cluster of DMPs. Postulating that a
#' reported DNA base is found at the beginning and the end of a cluster,
#' genome-wide screening can be performed where two consecutive targeting base
#' positions are separated by less than a given threshold
#' \strong{\emph{maxDist}} bases (1). DMPs from all the samples included in
#' \strong{\emph{GR}} are considered to build the clusters.
#'
#' The number of DMPs reported in each cluster corresponds to the numbers of
#' positions inside the cluster where DMPs were reported in at least one
#' of the samples.
#'
#' @param GR An object from \strong{\emph{'pDMP'}} class, which is returned by
#'     \code{\link{selectDIMP}} function.
#' @param maxDist maximum distance at which two reported bases sites from the
#'     same cluster can be separated. Default: \eqn{maxDist = 3}.
#' @param minNumDMPs minimum number of marked bases inside of each cluster.
#'     Default: \eqn{minNumDMPs = 1}.
#' @param chromosomes vector of characters labeling the chromosomes included in
#'     the analysis. Default: chromosomes = NULL (all chromosomes are included).
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
#' ## Get a dataset of potential signals and the estimated cutpoint from the
#' ## package
#' data(PS, cutpoint)
#'
#' ## The estimated cutpoints are used to discriminate signals from the noise.
#' ## That is, DMPs are selected using the cupoints
#' dmps <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)
#'
#' ## Build clusters of DMPs
#' x1 = dmpClusters(GR = dmps, maxDist = 7, minNumDMPs = 6, num.cores=2L,
#'                 verbose = FALSE)
#' x1
#'
#' ## ==== Setting up the experiment design to test for DMRs ===
#' nams <- names(dmps)
#' dmps_at_clusters <- getDIMPatGenes(GR = dmps, GENES = x1,
#'                                     ignore.strand = TRUE)
#' dmps_at_clusters <- uniqueGRanges(dmps_at_clusters, columns = 2L,
#'                                     ignore.strand = TRUE, type = "equal",
#'                                     verbose = FALSE)
#' colnames(mcols(dmps_at_clusters)) <- nams
#'
#' colData <- data.frame(condition = factor(c("CT", "CT", "CT",
#'                                            "TT", "TT", "TT"),
#'                                          levels = c("CT", "TT")),
#'                     nams,
#'                     row.names = 2)
#'
#' ## Build a RangedGlmDataSet object
#' ds <- glmDataSet(GR = dmps_at_clusters, colData = colData)
#'
#' ## ================ Testing for DMRs ===========
#' dmrs <- countTest2(DS = ds, num.cores = 4L,
#'                    minCountPerIndv = 4,
#'                    maxGrpCV = c(1, 1),
#'                    Minlog2FC = 0.5,
#'                    CountPerBp = 0.001,
#'                    test = "LRT",
#'                    verbose = TRUE)
#' dmrs
dmpClusters <- function(GR, maxDist = 3,
                        minNumDMPs = 1,
                        chromosomes = NULL,
                        num.cores = 1L, tasks = 0L,
                        verbose = TRUE, ...) {
   validateClass(GR)
   if (!inherits(GR, "pDMP"))
       stop("*** GR object must inherits from 'pDMP' class", " which is
           returned by calling 'selectDMP' function.")
   GR <- structure(GR, class = "list")
   GR <- GRangesList(GR)
   if (length(GR) == 1) GR <- unlist(GR, use.names =FALSE)

   signal <- cluster <- NULL

   GR <- meth_status(gr = GR, chromosomes = chromosomes, num.cores = num.cores,
                       tasks = tasks,  verbose = verbose)
   GR <- unlist(GR, use.names = FALSE)

   if (verbose) message("\n *** Building clusters ...")

   GR = sortBySeqnameAndStart(GR)
   strand(GR) <- rep("*", length(GR))
   signals <- GR[, "signal"]
   GR <- GR[GR$signal > 0, ]
   signals$cluster <- NULL

   # To find the word frameworks "cl"
   GR <- reduce(GR, min.gapwidth = maxDist)
   GR$cluster <- paste(seqnames(GR), start(GR), end(GR), sep = "_")

   # To assign the cl ids to the signal
   Hits <- findOverlaps(signals, GR, ignore.strand = TRUE)
   signals$cluster[queryHits(Hits)] <- GR$cluster[subjectHits(Hits)]
   signals <- unique(signals)

   # To count the number of sites "1" inside each cluster
   if (verbose) message("*** Counting DMPs in clusters ...")
   cl <- data.table(as.data.frame(signals))
   cl <- cl[!is.na(cl$cluster), list(seqnames = unique(seqnames),
                                   start = min(start), end = max(end),
                                   dmps = sum(signal)),
                       by = cluster]
   cl <- makeGRangesFromDataFrame(data.frame(cl), keep.extra.columns = TRUE)
   if(minNumDMPs > 0) cl <- cl[cl$dmps >= minNumDMPs]
   return(cl[, "dmps"])
}

# ========================= Auxiliary signal ================================= #
meth_status <- function(gr, chromosomes = NULL, num.cores = 1L,
                       tasks = 0L, verbose = TRUE) {
   ## Set parallel computation
   progressbar = FALSE
   if (verbose) progressbar = TRUE

   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks = tasks,
                               progressbar = progressbar)
   } else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                               progressbar = progressbar)
   status <- function(y, status, CHR, verbose) {
       if (verbose) message("\n *** Processing chromosome:  ", CHR, "\n")
       seqlevels(y, pruning.mode="coarse") <- CHR
       min.start <- min(start(y))
       max.end <- max(end(y))
       if(status == 0) {
           y <- disjoin(gaps(y))
           starts <- start(y)
           ends <- end(y)
           idx <- which(starts > min.start)
           starts <- starts[idx]
           ends <- ends[idx]
           idx <- which(ends < max.end)
           starts <- starts[idx]
           ends <- ends[idx]
           post<- bpmapply(function(s,e) s:e, starts, ends, USE.NAMES=FALSE,
                           BPPARAM = bpparam)

           post <- sort(unique(unlist(post)))
           y <- GRanges(seqnames=CHR, ranges=IRanges(start=post, end=post))
       }

       if(status == 1) {
           starts <- start(y)
           ends <- end(y)
           if (!all(starts == ends)) {
               post <- mapply(function(s,e) s:e, starts, ends, USE.NAMES=FALSE)
               post <- sort(unique(unlist(post)))
               y <- GRanges(seqnames=CHR, ranges=IRanges(start=post, end=post))
           }
       }
       mcols(y) <- data.frame(signal=status)
       return(y)
   }

   # if gr is a GRangesList object
   if (inherits(gr, "GRangesList")) {
       if (verbose)
           message("* Building a unique GRange objects ... \n")
       gr <- uniqueGRanges(gr, chromosomes = chromosomes, missing = 0,
                           type = "equal", num.cores = num.cores, tasks = tasks,
                           verbose = FALSE)
       mcols(gr) <- NULL
   }

   if(is.null(chromosomes)) CHRs <- unique(as.character(seqnames(gr))) else {
       CHRs <- chromosomes
   }

   kCHR <- length(CHRs)
   if (verbose) cat("* Computing the genomic signal for status '1' ...\n")
   r <- lapply(CHRs, function(chr) status(y = gr, status = 1, CHR = chr,
                                       verbose = verbose))
   if (verbose) cat("* Computing the genomic signal for status '0' ...\n")
   y <- lapply(CHRs, function(chr) status(y = gr, status = 0, CHR = chr,
                                       verbose = verbose))
   if (verbose) cat("* Building the methylation signal ...\n")
   gr <- bpmapply(function(x,y) sortBySeqnameAndStart(c(x,y)), r, y,
                   BPPARAM = bpparam)
   rm(r,y); gc()
   gr <- lapply(gr, unique)
   names(gr) <- CHRs
   return(GRangesList(gr))
}


