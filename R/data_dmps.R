#' Simulated dataset of DMPs used in examples
#'
#' Each individuals sample includes 10000 cytosine positions
#'
#' @format dmps is an object from class 'pDMP' carrying in the meta-columns the
#' following variables:
#' \describe{
#'     \item{\strong{p1}}{methylation level from the reference sample}
#'     \item{\strong{p2}}{methylation level from the treatment sample}
#'     \item{\strong{TV}}{the total variation distance (difference of
#'               methylation levels)}
#'     \item{\strong{hdiv}}{Hellinger divergence}
#'     \item{\strong{wprob}}{the probabilities:
#'     \eqn{wprob = 1 - CDF probability}}
#' }
#'
#'\strong{\emph{dmps}} is an object from class 'pDMP' carrying the same
#'meta-columns as 'HD' (dataset) plus the probabilities: \eqn{wprob = 1 - CDF
#'probability}. \strong{\emph{DMPs}} ('dmps') are obtained from potential DMPs
#'(see \code{\link{PS}}) with \code{\link{selectDIMP}} function.
"dmps"
