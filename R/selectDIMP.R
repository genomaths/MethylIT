#' @rdname selectDIMP
#'
#' @title Selection of DIMPs
#' @description For a given cutpoint (previously estimated with the function
#'     estimateCutPoint), 'selectDIMP' will return the differentially
#'     informative methyated positions (DIMPs). DIMPs are cytosine positions for
#'     which the divergence is greater than the cutpoint.
#' @details Theoretically a DIMP denotes a cytosine position with high
#'     probability to be differentially methylated. That is, in the statistical
#'     molecular-biophysics context, a DIMP must be considered only in a
#'     probabilistic term and not as an absolute deterministic experimental
#'     output.
#'
#'     The uncertainty and dynamics of the DNA methylation process, the
#'     continuous action of the omnipresent thermal fluctuations, as well as,
#'     the inherent stochasticity of the biochemical reactions make it
#'     impossible to ensure whether a specific cytosine position is methylated
#'     in an absolutely deterministic sense. Notice that the concept of DIMP is
#'     not applicable to a single cell (if we use an instrumentation/protocol
#'     to directly measure methylation at the molecular level, and not via PCR),
#'     since a concrete, single DNA cytosine position in a single cell is
#'     methylated or not methylated.
#'
#'     However, when pooling DNA extracted from a tissue, the previous
#'     reasonings about uncertainty hold plus an additional uncertainty factor:
#'     cells from the same tissue are not synchronized but are found in
#'     different stages of their ontogenetic developments. Hence, the DIMP
#'     concept holds in the mentioned circumstances where the uncertainty of
#'     methylation is present.
#'
#' @param LR An object from "pDMP" class.
#' @param div.col Number of the column where the divergence variable (i.e.,
#'     Hellinger divergence or total variation) is located in the GRanges
#'     meta-columns.
#' @param pval.col If the cutpoints is a p-value, then the column number for
#'     p-values should be provided. Default: NULL. Notice that one of the
#'     parameter values div.col or pval.col must be given.
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MethylIT
#'     analysis. If 'absolute = TRUE', then TV is tranformed into |TV|, which is
#'     an information divergence that can be fitted to Weibull or to Generalized
#'     Gamma distribution. So, if the nonlinear fit was performed for |TV|, then
#'     here absolute must be set to TRUE.
#' @param cutpoint Cutpoint to select DIMPs. Cytosine positions with divergence
#'     greater than 'cutpoint' will selected as DIMPs. Cutpoints are estimated
#'     with the function 'estimateCutPoint'.
#' @param tv.col Column number for the total variation to be used for filtering
#'     cytosine positions (if provided).
#' @param tv.cut If tv.cut and tv.col are provided, then cytosine sites k with
#'     abs(TV_k) < tv.cut are removed.
#'
#' @return An object from "pDMP" class containing only differentially
#'     informative position (DIMPs).
#'
#' @examples
#' ## Get a dataset of potential signals and the estimated cutpoint from the
#' ## package
#' data(PS, cutpoint)
#'
#' ## The estimated cutpoints are used to discriminate signals from the noise.
#' ## That is, DMPs are selected using the cupoints
#' DMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)
#'
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom S4Vectors mcols
#' @export
selectDIMP <- function(LR, div.col=NULL, pval.col=NULL, absolute=FALSE,
                       cutpoint, tv.col=NULL, tv.cut=NULL) {
   if (!inherits(LR, what=c("pDMP", "InfDiv")))
       stop("*** LR object must be from 'pDMP' or 'InfDiv' class")

  # -------------------------- valid "pDMP" object ----------------------------
  validateClass(LR)
  # ------------------------------------------------------------------------- #

   if (is.null(div.col) && is.null(pval.col))
       stop("*** One of the parameters 'div.col'
           or 'pval.col' must be not NULL")
   if (is.null(div.col)) target.col = pval.col else target.col = div.col

   for (k in 1:length(LR)) {
       x <- LR[[k]]
       if (!is.null(tv.cut) && !is.null(tv.col))
           x <- x[abs(mcols(x[, tv.col])[, 1]) > tv.cut]
       if (is.null(div.col)) {
           LR[[k]] <- x[mcols(x[, target.col])[, 1] < cutpoint]
       } else {
           if (absolute)
               LR[[k]] <- x[abs(mcols(x[, target.col])[, 1]) > cutpoint]
           else LR[[k]] <- x[mcols(x[, target.col])[, 1] > cutpoint]
       }
   }
   return(LR)
}
