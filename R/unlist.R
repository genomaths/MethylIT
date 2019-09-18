#' @rdname unlist
#' @title Flatten lists extended to any class
#' @description Given a list 'x' of R objects from the same class and same
#'     format, unlist simplifies it to produce a new R object which contains all
#'     the initial components which in 'x' object.
#' @param x Any list R object.
#' @param recursive Logical. Should unlisting be applied to list components of
#'     x?
#' @param use.names Logical. Should names be preserved?
#' @details This is a method to extend unlist generic function to handle
#'     any list of objects from the same class.
#' @export
#' @examples
#' gr1 <-GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
#'               strand = "+", score = 5L, GC = 0.45)
#' gr2 <-
#'   GRanges(seqnames = c("chr1", "chr1"),
#'           ranges = IRanges(c(7,13), width = 3),
#'           strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
#' gr3 <-
#'   GRanges(seqnames = c("chr1", "chr2"),
#'           ranges = IRanges(c(1, 4), c(3, 9)),
#'           strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
#' grl <- list("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
#' base::unlist(grl) # The default unlist does not work
#' unlist(grl)

setGeneric("unlist", signature = "x")
unlist <- function(x, recursive = TRUE, use.names = TRUE) UseMethod("unlist")

#' @aliases unlist.default
#' @export
#' @keywords internal
unlist.default <- function(x, recursive = TRUE, use.names = TRUE) {
   n <- length(x)
   x <- base::unlist(x, recursive = recursive, use.names = use.names)

   if (length(x) == n) {
     x0 <- try(suppressWarnings(do.call("c", unname(x))), silent = TRUE)
     if (!inherits(x0, "try-error")) x <- x0
   }
   return(x)
}
