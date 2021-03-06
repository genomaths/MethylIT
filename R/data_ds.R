#' Simulated dataset of RangedGlmDataSet class object (DMPs counts)
#'
#' \emph{RangedGlmDataSet} and \emph{glmDataSet} are objects carrying the
#' information on the experimental design for two group comparison of DMP
#' counts by applying generalized linear regression.
#'
#' @format \strong{\emph{ds}} is an \emph{RangedGlmDataSet} with 125 regions
#' and 4 columns (individuals) with factor levels 'A' and 'B'. The accessible
#' objects in the dataset are:
#' \describe{
#'     \item{\strong{GR}}{A \code{\link[GenomicRanges]{GRanges-class}} object of
#'     length(GR) = 125 with the count matrix of DMPs in the metacolumns.}
#'     \item{\strong{counts}}{Count matrix of DMPs with minimal
#'     dim(ds$counts) = c(125, 4).}
#'     \item{\strong{colData}}{A data frame with one columnn named 'condition'.}
#'     \item{\strong{sampleNames}}{Samples names: "A1" "A2" "B1" "B2"}
#'     \item{\strong{levels}}{Design factor levels: "A" "B"}
#' }
#'
#'\strong{\emph{ds}} is an object from class \emph{RangedGlmDataSet} carrying
#'the information to perform DMG analysis with \code{\link{countTest2}}
#'function.
"ds"
