#' @rdname meth_level
#' @title Compute methylation levels
#' @description This function computes the
#' @param GR A matrix of counts or \code{\link{GenomicRanges}{"GRanges"}}
#' object with the table of counts in the meta-columns (methylated mC and
#' unmethylated uC cytosines) or list of matrices of GRanges objects . Unless
#' specified in the parameter 'columns', the methylation counts must be given in
#' the first four columns: 'mC1' and 'uC1' methylated and unmethylated counts
#' for control sample, and 'mC2' and 'uC2' methylated and unmethylated counts
#' for treatment sample, respectively.
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
#' @param preserve.gr logical(1). Option of whether to preserve all
#' the metadata from the original GRanges object.
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to use,
#' i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS). These parameters will be passed to
#' \code{\link[GenomicWordFramework]{uniqueGRanges}}.
#' @param verbose if TRUE, prints the function log to stdout

#' @aliases meth_level
setGeneric("meth_level",
           function(GR,
                    Bayesian = FALSE,
                    columns = c(mC1 = 1, uC1 = 2, mC2 = 3, uC2 = 4),
                    min.coverage = 4,
                    tv = FALSE,
                    bay.tv = FALSE,
                    preserve.gr = FALSE,
                    num.cores = 1,
                    tasks = 0L,
                    verbose = TRUE, ...) standardGeneric("meth_level"))

## Set class or "matrix
# setClassUnion("GRanges_OR_matrix", c("GRanges", "matrix"))

#' @aliases meth_level
#' @rdname meth_level
setMethod("meth_level", signature(GR = "GRanges"),
          function(GR,
                   Bayesian = FALSE,
                   columns = c(mC1 = 1, uC1 = 2, mC2 = 3, uC2 = 4),
                   min.coverage = 4,
                   tv = FALSE,
                   bay.tv = FALSE,
                   preserve.gr = FALSE,
                   verbose = TRUE, ...) {

            if (inherits(GR, "GRanges")) x <- as.matrix(mcols(GR))
            x <- x[, columns]
            if (ncol(x) < 4) {
              stop("If counts are provided, then 'length(columns) = 4'",
                   " is expected\n If you are providing methylation level, ",
                   "please set 'meth.level = TRUE'")
            }

            r0 <- rowSums(x)
            ind <- which(r0 > min.coverage); rm(r0)
            x <- x[ind, ]
            GR <- GR[ ind ]
            n1 <- x[, 1] + x[, 2]
            n2 <- x[, 3] + x[, 4]
            n <- cbind(n1, n2)
            x0 <- x
            p1 <- x[, 1]/n1
            p2 <- x[, 3]/n2

            ## if the coverage is zero in control or the
            ## reference individual, then p is NaN (NA). By
            ## definition these cytosine sites has methylation
            ## levels p = 0.
            p1[is.na(p1)] <- 0
            p2[is.na(p2)] <- 0

            ## Ordinary TV
            if (tv) TV <- p2 - p1

            if (Bayesian) {
              if (nrow(x) < 10)
                stop(paste("*** You must provide at least 10 cytosine sites ",
                           "to apply a Bayessian approach \n",
                           "using beta distributed priors"))
              if (verbose)
                cat("*** Estimating betaBinomial-posteriors... \n")
              ## Naive distribution q (methylation levels).  In a
              ## Bayesian framework with uniform priors, the
              ## methylation level can be defined as: meth_level =
              ## ( mC + 1 )/( mC + uC + 2 ).
              q1 <- (x[, 1] + 1)/(n1 + 2)
              q2 <- (x[, 3] + 1)/(n2 + 2)

              ## The shape parameters estimated with 'nlm'
              beta1 <- .estimateBetaDist(q1)
              beta2 <- .estimateBetaDist(q2)
              ## Assuming beta priors
              n1[n1 == 0] <- 2
              n2[n2 == 0] <- 2

              p1 <- .betaBinPosteriors(x[, 1], n1, a = beta1[1], b = beta1[2])
              p2 <- .betaBinPosteriors(x[, 3], n2, a = beta2[1], b = beta2[2])
              if (bay.tv) bay.TV <- p2 - p1
            }

            if (preserve.gr) mcols(GR) <- data.frame(mcols(GR), p1, p2)
            else mcols(GR) <- data.frame(p1, p2)
            if (tv) GR$TV <- TV
            if (bay.tv) GR$bay.TV <- bay.TV

            return(GR)
          })

#' @aliases meth_level
#' @rdname meth_level
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod("meth_level", signature(GR = "list"),
          function(
            GR,
            Bayesian = FALSE,
            columns = c(mC1 = 1, uC1 = 2, mC2 = 3, uC2 = 4),
            min.coverage = 4,
            tv = FALSE,
            bay.tv = FALSE,
            preserve.gr = FALSE,
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
            GR <- bplapply(GR, meth_level, Bayesian, columns, min.coverage,
                           tv, bay.tv, preserve.gr, num.cores, tasks, verbose,
                           BPPARAM = bpparam)
            return(GR)
})
