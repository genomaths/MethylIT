#' @rdname estimateHellingerDiv
#'
#' @title Hellinger divergence of methylation levels
#' @description Given a the methylation levels of two individual, the function
#'     computes the information divergence between methylation levels.
#' @details Each methylation level j for cytosine site i corresponds to a
#'     probability vector p_j = c(p_ij, 1 - p_ij). Then, the information
#'     divergence between methylation levels p1 and p2 is the divergence
#'     between the vectors p1 = c(p_i1, 1 - p_i1) and p2 = c(p_i2, 1 - p_i2).
#'     If the vector of coverage is supplied, then the information divergence
#'     is estimated according to the formula:
#'
#'     hdiv = 2*(n[1] + 1)*(n[2] + 1)*((sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 -
#'     p[1]) - sqrt(1 - p[2]))^2)/(n[1] + n[2] + 2) This formula corresponds to
#'     Hellinger divergence as given in the first formula from Theorem 1 from
#'     reference 1. Otherwise: hdiv = (sqrt(p[1]) - sqrt(p[2]))^2 + (sqrt(1 -
#'     p[1]) - sqrt(1 - p[2]))^2
#'
#' @param p A numerical vector of the methylation levels p = c(p1, p2) of
#'     individuals 1 and 2.
#' @param n if supplied, it is a vector of integers denoting the coverages used
#'     in the estimation of the methylation levels.
#' @return The Hellinger divergence value for the given methylation levels is
#'     returned
#'
#' @examples
#'     p <- c(0.5, 0.5)
#'     estimateHellingerDiv(p)
#'
#' @references ' 1. Basu  A., Mandal  A., Pardo L (2010) Hypothesis testing for
#'     two discrete populations based on the Hellinger distance. Stat Probab
#'     Lett 80: 206-214.
#'
#' @export
estimateHellingerDiv <- function(p, n=NULL) {
   hdiv <- 0
   if (!is.na(sum(p)) && (sum(p) > 0)) {
       if (is.null(n)) {
           hdiv <- (2 * ((sqrt(p[1]) - sqrt(p[2])) ^ 2 + (sqrt(1 - p[1]) -
                   sqrt(1 - p[2])) ^ 2))
       } else {
           ## Hellinger Chi-squared (Hellinger divergence)
           hdiv <- (2 * (n[1] + 1) * (n[2] + 1) * ((sqrt(p[1]) -
                   sqrt(p[2])) ^ 2 + (sqrt(1 - p[1]) - sqrt(1 - p[2])) ^ 2) /
                   (n[1] + n[2] + 2))
       }
   }
   return(hdiv)
}
