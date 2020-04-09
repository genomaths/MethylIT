#' Classification PCA+QDA model for simulated dataset of DMPs used in examples
#'
#' This data/object carries the information about the classification performance
#' of the combined models of Principal Components \emph{PCA} and Quadratic
#' Discriminant (\emph{"QDA"})) analyses on the set of \code{\link{dmps}}.
#'
#' @format \emph{pcaQda_perf} is an object from class \emph{"pcaQDA"},
#' consisting of a list with the following elements:
#' \describe{
#'     \item{\strong{Performance}}{Classification performance of the
#'     \emph{"PCA+QDA"} model on the set of \code{\link{dmps}}}.
#'     \item{\strong{FDR}}{False discovery rate estimated for the
#'     \emph{"PCA+QDA"} model on the set of \code{\link{dmps}}}.
#'     \item{\strong{model}}{The \emph{"pcaQDA"} model fitted on the set
#'     of \code{\link{dmps}}, carrying the \code{\link[MASS]{qda}} and the
#'     \emph{"PCA"} models. The \emph{"PCA"} is fitted with
#'     \code{\link[stats]{prcomp}}} function.
#' }
#'
#'\strong{\emph{pcaQda_perf}} is an object from "pcaQDA" class, which was
#'obtained applying \code{\link{evaluateDIMPclass}} function on the set of
#'\code{\link{dmps}}.
"pcaQda_perf"
