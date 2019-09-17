#' @rdname print.ldaDMP
#' @aliases print.ldaDMP
#' @title Printing object from \emph{print.ldaDMP} class by simple
#'     print methods.
#' @param x Object from class \emph{glmDataSet} or from class
#'     \emph{"RangedGlmDataSet"}.
#' @param digits Number of significant digits to be used.
#' @keywords internal
#' @exportMethod print.ldaDMP
print.ldaDMP <- function(x, digits = getOption("digits")) {
   if (!is.null(cl <- x$call)) {
      names(cl)[2L] <- ""
      cat("Call:\n")
      dput(cl, control = NULL)
   }
   cat("\nPrior probabilities of groups:\n")
   print(x$prior, digits)
   cat("\nGroup means:\n")
   print(x$means, digits)
   cat("\nCoefficients of linear discriminants:\n")
   print(x$scaling, digits)
   svd <- x$svd
   names(svd) <- dimnames(x$scaling)[[2L]]
   if (length(svd) > 1L) {
       cat("\nProportion of trace:\n")
       print(round(svd^2/sum(svd^2), 4L), digits)
   }
   invisible(x)
}

