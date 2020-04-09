#' Classification LDA model for simulated dataset of DMPs used in examples
#'
#' This data/object carries the information about the classification performance
#' of a Quadratic Discriminant (\emph{"QDA"})) model on the set of
#' \code{\link{dmps}}.
#'
#' @format \emph{qda_perf} is list object consisting of the following elements:
#' \describe{
#'     \item{\strong{Performance}}{Classification performance of the
#'     \emph{"LDA"} model on the set of \code{\link{dmps}}}.
#'     \item{\strong{FDR}}{False discovery rate estimated for the \emph{"QDA"}
#'     model on the set of \code{\link{dmps}}}
#'     \item{\strong{model}}{The \code{\link[MASS]{qda}} model fitted on the set
#'     of \code{\link{dmps}}}.
#' }
#'
#'\strong{\emph{qda_perf}} is a list object obtained applying
#'\code{\link{evaluateDIMPclass}} function on the set of
#'\code{\link{dmps}}.
"qda_perf"
