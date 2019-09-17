#' @rdname print.qdaDMP
#' @aliases print.qdaDMP
#' @title Printing object from \emph{print.qdaDMP} class by simple
#'     print methods.
#' @param x Object from class \emph{glmDataSet} or from class
#'     \emph{"RangedGlmDataSet"}.
#' @param digits Number of significant digits to be used.
#' @keywords internal
#' @exportMethod print.qdaDMP
print.qdaDMP <- function(x, digits = getOption("digits")) {
   if (!is.null(cl <- x$call)) {
       names(cl)[2L] <- ""
       cat("Call:\n")
       dput(cl, control = NULL)
   }
   cat("\nPrior probabilities of groups:\n")
   print(x$prior, digits)
   cat("\nGroup means:\n")
   print(x$means, digits)
   invisible(x)
}

