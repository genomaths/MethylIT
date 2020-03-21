#' Simulated dataset of Hellinger divergences used in the examples
#'
#' Each individuals sample includes 10000 cytosine positions
#'
#' @format HD is an object from class 'InfDiv' carrying in the meta-columns the
#'     following variables:
#'     \describe{
#'         \item{p1}{methylation level from the reference sample}
#'         \item{p2}{methylation level from the treatment sample}
#'         \item{TV}{the total variation distance (difference of
#'                 methylation levels)}
#'         \item{hdiv}{Hellinger divergence}
#'     }
#'
#'     'HD' was obtained with function \code{\link{estimateDivergence}}.
"HD"
