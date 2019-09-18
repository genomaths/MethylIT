#' @rdname getPotentialDIMP
#'
#' @title Potential methylation signal
#' @description This function perform a selection of the cytosine sites
#'     carrying the potential methylation signal. The potential signals from
#'     controls and treatments are used as prior classification in further step
#'     of signal detection.
#' @details The potential signals are cytosine sites k with information
#'     divergence (DIV_k) values greater than the DIV(alpha = 0.05). The value
#'     of alpha can be specified. For example, potential signals with
#'     DIV_k > DIV(alpha = 0.01) can be selected. For each sample, cytosine
#'     sites are selected based on the corresponding nonlinear fitted
#'     distribution model that has been supplied.
#' @param LR An object from 'InfDiv' or "testDMP" class. These objects are
#'     previously obtained with function \code{\link{estimateDivergence}} or
#'     \code{\link{FisherTest}}.
#' @param nlms A list of distribution fitted models (output of
#'     'fitNonlinearWeibullDist' function) or NULL. If NULL, then empirical
#'     cumulative distribution function is used to get the potential DIMPs.
#' @param div.col Column number for divergence variable is located in the
#'     meta-column.
#' @param dist.name name of the distribution to fit: Weibull2P (default:
#'     "Weibull2P"), Weibull three-parameters (Weibull3P), gamma with
#'     three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
#'     generalized gamma with three-parameter ("GGamma3P") or four-parameter
#'     ("GGamma4P"), the empirical cumulative distribution function (ECDF) or
#'     "None". If \strong{dist.name != "None"}, and \strong{nlms != NULL}, then
#'     a column named "wprob" with a probability vector derived from the
#'     application of model "nlms" will be returned.
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MethylIT
#'     analysis. If 'absolute = TRUE', then TV is transformed into |TV|, which
#'     is an information divergence that can be fitted to Weibull or to
#'     Generalized Gamma distribution. So, if the nonlinear fit was performed
#'     for |TV|, then absolute must be set to TRUE.
#' @param alpha A numerical value (usually alpha < 0.05) used to select
#'     cytosine sites k with information divergence (DIV_k) for which Weibull
#'     probability P[DIV_k > DIV(alpha)].
#' @param pval.col An integer denoting a column from each GRanges object from
#'     LR where p-values are provided when \strong{dist.name == "None"} and
#'     \strong{nlms == NULL}. Default is NULL. If NUll and
#'     \strong{dist.name == "None"} and \strong{nlms == NULL}, then a column
#'     named \strong{adj.pval} will used to select the potential DMPs.
#' @param tv.col Column number for the total variation to be used for filtering
#'     cytosine positions (if provided).
#' @param tv.cut If tv.cut and tv.col are provided, then cytosine sites k with
#'     abs(TV_k) < tv.cut are removed before to perform the ROC analysis.
#' @param hdiv.col Optional. A column number for the Hellinger distance to be
#'     used for filtering cytosine positions. Default is NULL.
#' @param hdiv.cut If hdiv.cut and hdiv.col are provided, then cytosine sites k
#'     with hdiv < hdiv.cut are removed.
#' @param min.coverage Cytosine sites with coverage less than min.coverage are
#'     discarded. Default: 0
#' @param pAdjustMethod method used to adjust the p-values from other
#'     approaches like Fisher's exact test, which involve multiple comparisons
#'     Default is NULL. Do not apply it when a probability distribution model
#'     is used (\strong{when nlms is given}), since it makes not sense.
#' @return A list of GRanges objects, each GRanges object carrying the selected
#'     cytosine sites and and the Weibull probability P[DIV_k > DIV(alpha)].
#' @importFrom stats p.adjust.methods p.adjust plnorm pweibull pgamma
#' @importFrom S4Vectors mcols mcols<-
#'
#' @export
#'
#' @examples
#' ## Get a dataset of Hellinger divergency of methylation levels and their
#' ## corresponding best nonlinear fit distribution models from the package
#' data(HD, nlms)
#' PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 9L, alpha = 0.05)
#'
getPotentialDIMP <- function(LR, nlms=NULL, div.col, dist.name = "Weibull2P",
                           absolute=FALSE, alpha=0.05, pval.col = NULL,
                           tv.col=NULL, tv.cut=NULL, min.coverage=NULL,
                           hdiv.col = NULL, hdiv.cut = NULL,
                           pAdjustMethod = NULL) {

   # -------------------------- valid "InfDiv" object ------------------------ #
   validateClass(LR)
   # ------------------------------------------------------------------------- #

   cl <- inherits(LR, "testDMP")
   model <- (!is.null(nlms) && dist.name != "None")

   if (!is.null(hdiv.cut) && is.null(hdiv.col)) {
       cat("\n")
       stop("You set hdiv.cut = ", hdiv.cut, ".",
            " You must provide 'hdiv.col' as well")
   }

   if (is.null(hdiv.cut) && !is.null(hdiv.col)) {
       cat("\n")
       stop("You set hdiv.col = ", hdiv.col, ".",
           " You must provide 'hdiv.cut' as well")
   }

   if (!is.null(tv.cut) && is.null(tv.col)) {
       cat("\n")
       stop("You set tv.cut = ", tv.cut, ".",
           " You must provide 'tv.col' as well")
   }

   if (is.null(hdiv.cut) && !is.null(hdiv.col)) {
       cat("\n")
       stop("You set tv.col = ", tv.col, ".",
           " You must provide 'tv.cut' as well")
   }

   P <- function(k) {
       d <- LR[[k]]
       if (length(d) > 0) {
           if (!is.null(min.coverage)) {
               cov1 <- d$c1 + d$t1
               cov2 <- d$c2 + d$t2
               idx <- which((cov1 >= min.coverage) | (cov2 >= min.coverage))
               d <- d[ idx ]
           }
           q <- mcols(d[, div.col])[, 1]

           if (dist.name == "ECDF") ECDF <- ecdf(q)

           if (!is.null(tv.col) && !is.null(tv.cut)) {
               idx <- which( abs(mcols(d[, tv.col])[, 1]) > tv.cut)
               d <- d[ idx ]
               q <- q[ idx ]
           }

           if (absolute) q = abs(q)
           if (!is.null(nlms)) {
               m <- nlms[[k]]
               m <- m[, 1]
           } else  {
               if (!cl || !model) {
                   dist.name <- "ECDF"
                   ECDF <- ecdf(q)
               }
           }

           if (dist.name != "ECDF" && model) {
               p <- switch(dist.name,
                           LogNorm=plnorm(q, meanlog=m[1], sdlog=m[2],
                                           lower.tail=FALSE),
                           Weibull2P=pweibull(q, shape=m[1], scale=m[2],
                                               lower.tail=FALSE),
                           Weibull3P=pweibull(q - m[3], shape=m[1], scale=m[2],
                                               lower.tail = FALSE),
                           Gamma2P=pgamma(q, shape=m[1], scale=m[2],
                                           lower.tail = FALSE),
                           Gamma3P=pgamma(q - m[3], shape=m[1], scale=m[2],
                                           lower.tail = FALSE),
                           GGamma3P=pggamma(q, alpha=m[1], scale=m[2], psi=m[3],
                                           lower.tail = FALSE),
                           GGamma4P=pggamma(q, alpha=m[1], scale=m[2], mu=m[3],
                                           psi=m[4], lower.tail = FALSE)
               )
           }
           if (dist.name == "ECDF") p <- (1 - ECDF(q))
           if (!model && cl && is.null(pval.col)) p <- d$adj.pval
           else if (!model && is.numeric(pval.col)) p <- mcols(d)[, pval.col]
           if (!model && !cl) p <- (1 - ECDF(q))

           if (!is.null(pAdjustMethod)) {
               pAdjustMethod <- match.arg(pAdjustMethod, p.adjust.methods)
               p <- p.adjust(p, method = pAdjustMethod)
           }

           idx <- which(p < alpha)
           p <- p[ idx ]
           d <- d[ idx ]
           if (!is.null(hdiv.cut) && !is.null(hdiv.col)) {
               idx <- which(mcols(d[, hdiv.col])[, 1] > hdiv.cut)
               d <- d[ idx ]
               p <- p[ idx ]
           }
           mcols(d) <- data.frame(mcols(d), wprob = p)
       } else mcols(d) <- data.frame(mcols(d), wprob = numeric(0))
       return(d)
   }
   sn <- names(LR)
   LR <- lapply(1:length(LR), P, keep.attr = TRUE)
   names(LR) <- sn
   if (cl) {
       LR <- structure(LR, class = c("pDMP", "InfDiv", "testDMP", "list"))
   } else LR <- structure(LR, class = c("pDMP", "InfDiv", "list"))
   return(LR)
}
