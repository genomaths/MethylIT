#' @rdname lapply
#' @title Apply a function over a list-like object preserving its attributes
#' @param x A list-like or vector-like object
#' @param FUN,... See ?base::\code{\link[base]{lapply}} for a description of
#'     these arguments.
#' @param keep.attr Logic. If TRUE, then the original attributes from 'x' are
#'     preserved in the returned list. Default is FALSE.
#' @description lapply returns a list of the same length as 'x', each element of
#'     which is the result of applying FUN to the corresponding element of 'x'.
#' @return Same as in ?base::\code{\link[base]{lapply}} if keep.attr = FALSE.
#'     Otherwise same values preserving original attributes from 'x'.
#' @seealso base::\code{\link[base]{lapply}}
#' @examples
#' # Create a list
#' x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
#' class(x) <- 'nice'
#'
#' # compute the list mean for each list element using 'base::lapply'
#' class(lapply(x, mean))
#'
#' # To preserve attributes
#' class(lapply(x, mean, keep.attr = TRUE))
#' @keywords internal
#' @export
lapply <- function(x, FUN, ...) UseMethod("lapply", x)

#' @name lapply.default
#' @rdname lapply
#' @export
lapply.default <- function(x, FUN, keep.attr = FALSE, ...) {
    if (keep.attr) {
        cl <- class(x)
        nm <- names(x)
        x <- base::lapply(x, FUN, ...)
        x <- structure(x, class = cl, names = nm)
    } else x <- base::lapply(x, FUN, ...)
    return(x)
}
