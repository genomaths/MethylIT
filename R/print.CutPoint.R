#' @rdname print.CutPoint
#' @title Printing object from 'CutPoint' class by simple print methods
#' @param x Object of class \code{"CutPoint"}
#' @param digits Number of significant digits to be used.
#' @keywords internal
#' @exportMethod print.CutPoint
print.CutPoint <- function(x, digits = getOption("digits")) {

   postProbCut <- format(signif(x$postProbCut, max(1L, digits - 2L)))
   cutpoint <- format(signif(x$cutpoint, max(1L, digits - 2L)))

   if (x$initModel != "Youden Index" ) {
     cat("Cutpoint estimation with '", x$initModel, "' classifier \n",
         sep = "")
   } else {
     cat("Cutpoint estimation with '", x$initModel, "' \n", sep = "")
   }
   if (!is.na(x$statistic)) {

       sta.val <- format(signif(x$optStatVal, max(1L, digits - 2L)))

       cat("Cutpoint search performed using model posterior probabilities \n")
       cat("\n")
       cat("Posterior probability used to get the cutpoint =",
           postProbCut, "\n")
       cat("Cytosine sites with treatment PostProbCut >=",
           postProbCut, "have a \n")
       cat("divergence value >=", x$postCut, "\n")
       cat("\n")

       cat("Optimized statistic:", x$statistic, "=", sta.val, "\n")
       cat("Cutpoint =", cutpoint, "\n")
       cat("\n")

       cat("Model classifier '", x$classifier, "' \n", sep = "")
       cat("\n")
   } else {
       if (x$initModel != "Youden Index") {
           cat("Posterior probability used to get the cutpoint =",
               postProbCut, "\n")
           cat("Cutpoint =", cutpoint, "\n")
           cat("\n")
           cat("Cytosine sites with treatment PostProbCut >=",
               postProbCut, "have a \n")
           cat("divergence value >=", cutpoint, "\n")
           cat("\n")
           cat("Model classifier", x$classifier, "\n")
           cat("\n")
       }
       if (x$initModel == "Youden Index") {
           cat("Simple cutpoint estimation \n")
           cat("Cutpoint =", cutpoint, "\n")
           cat("\n")
           cat("Cytosine sites from treatment have divergence values >=",
               cutpoint, "\n")
           cat("\n")
       }
   }
  cat("The accessible objects in the output list are: \n")
  print(summary(x))
}
