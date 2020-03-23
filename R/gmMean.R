#' @rdname gmMean
#'
#' @title Calculate the geometric mean
#' @description Calculate the geometric mean
#'
#' @param x numeric vector
#' @param na.rm Logical value. If TRUE, the NA values will be removed
#'
#' @return geometric mean
#' @keywords internal
.gmMean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}
