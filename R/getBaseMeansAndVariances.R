#' @rdname getBaseMeansAndVariances
#'
#' @title Calculate the means and variances by row
#' @description Calculate the means and variances in each row
#'
#' @param counts matrix of methylation counts
#' @param sizeFactors sizeFactors
#' @return data.frame with means and variances by row
#'
#' @importFrom genefilter rowVars
#' @keywords internal
.getBaseMeansAndVariances <- function(counts, sizeFactors) {
    data.frame(baseMean = rowMeans(t(t(counts)/sizeFactors)), 
        baseVar = rowVars(t(t(counts)/sizeFactors)))
}

