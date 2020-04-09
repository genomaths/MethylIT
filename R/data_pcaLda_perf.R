#' Classification PCA+LDA model for simulated dataset of DMPs used in examples
#'
#' This data/object carries the information about the classification performance
#' of the combined models of Principal Components \emph{PCA} and Linear
#' Discriminant (\emph{"LDA"})) analyses on the set of \code{\link{dmps}}.
#'
#' @format \emph{pcaLda_perf} is an object from class \emph{"pcaLDA"},
#' consisting of a list with the following elements:
#' \describe{
#'     \item{\strong{Performance}}{Classification performance of the
#'     \emph{"PCA+LDA"} model on the set of \code{\link{dmps}}}.
#'     \item{\strong{FDR}}{False discovery rate estimated for the
#'     \emph{"PCA+LDA"} model on the set of \code{\link{dmps}}}.
#'     \item{\strong{model}}{The \emph{"PCA+LDA"} model fitted on the set
#'     of \code{\link{dmps}}, carrying the \code{\link[MASS]{lda}} and the
#'     \emph{"PCA"} models. The \emph{"PCA"} is fitted with
#'     \code{\link[stats]{prcomp}}} function.
#' }
#'
#'\strong{\emph{pcaLda_perf}} is an object from "pcaLDA" class, which was
#'obtained applying \code{\link{evaluateDIMPclass}} function on the set of
#'\code{\link{dmps}}.
"pcaLda_perf"
