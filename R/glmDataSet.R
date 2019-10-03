#' @rdname glmDataSet
#' @title Data set constructor for class glmDataSet
#' @description This function is used to build a object suitable to be used
#'     with Methyl-IT \code{link{countTest2}} function.
#' @details Data set constructor for class glmDataSet also validate the object
#' @param GR A GRanges object with the count matrix of DMPs in the metacolumns
#'     (see \emph{'counts'}). If provided, then leave paramater
#'     \emph{'counts = NULL'}.
#' @param counts Count matrix of DMPs with minimal dimmensions 1 (row) x 4
#'     (columns). Column names must corresponds to the rownames from parameter
#'      'colData'.
#' @param colData A data frame with one columnn named 'condition', which must be
#'     a factor with exactly two levels. The rownames of \emph{colData}
#'     individual samples. The row names of \emph{colData} must correspond to
#'     th column names of the count matrix.
#'
#' @return A RangedGlmDataSet object, good for downstream use with
#'     Methyl-IT \code{link{countTest2}} function.
#'
#' @author Robersy Sanchez
#' @examples
#' DMPs <- data.frame(chr = "chr1", start = 1:10,
#'                    end = 1:10,strand = '*',
#'                    treat1 = rnbinom(size = 10, mu = 4, n = 500),
#'                    treat2 = rnbinom(size = 10, mu = 4, n = 500),
#'                    cntrl1 = rnbinom(size = 10, mu = 4, n = 500),
#'                    cntrl2 = rnbinom(size = 10, mu = 4, n = 500))
#' DMPs <- makeGRangesFromDataFrame(DMPs, keep.extra.columns = TRUE)
#' condition <- data.frame(condition = factor(c("TT","TT","CT","CT"),
#'                         levels = c("CT", "TT")))
#' rownames(condition) <- names(mcols(DMPs))
#' DIMR <- glmDataSet(GR = DMPs, colData = condition)
#'
#' @export
glmDataSet <- function(GR = NULL, counts = NULL, colData = NULL) {
   if (is.null(GR) && is.null(counts)) {
       cat("\n")
       stop("'GR' or 'counts' must be provided")
   }
   if (is.null(colData$condition)) {
      cat("\n")
      stop("In 'colData', 'condition' must be provided")
   }
    if (class(colData$condition) != "factor") {
       cat("\n")
       stop("In 'colData','condition' must be a 'factor'")
   }
   if (!is.null(colData$condition)) {
      if (length(levels(colData$condition)) != 2) {
         cat("\n")
         stop("In 'colData', factor 'condition' must have only two levels")
      }
   }
   if (!is.null(GR) && class(GR) != "GRanges") {
       cat("\n")
       stop("'GR' must be a GRanges")
   }
   if (!is.null(counts) && class(counts) != "matrix") {
       counts <- try(as.matrix(counts), silent = TRUE)
       if (inherits(counts, "try-error")) {
           cat("\n")
           stop("'counts' must be a matrix or an object coercible as a matrix")
       }
   }
   if (!is.null(counts)) {
       if (is.null(colnames(counts))) {
           cat("\n")
           stop("The 'counts' matrix must have column names, which \n",
                "must correspond to 'colData' rownames")
       }
       if (any(colnames(counts) != rownames(colData))) {
           cat("\n")
           stop("Individual names in 'counts' must correspond to ",
              "'colData' rownames")
       }
   }
   if (is.null(counts)) {
       if (ncol(mcols(GR)) > 0) {
           if (any(rownames(colData) != colnames(mcols(GR)))) {
               cat("\n")
               stop("Individual names in 'GR' must correspond to \n",
                   "'colData' rownames")
           }
       } else stop("Since 'counts = NULL', GR must have a count matrix in the",
                   " metacolumns or provide 'counts'")
   }
   if (is.null(GR)) {
       if (any(rownames(colData) != colnames(counts))) {
           cat("\n")
           stop("Individual names in 'counts' must correspond to \n",
               "condition rownames")
       }
   }

   if (is.null(GR)) {
       x <- list(counts = counts, colData = colData,
                 sampleNames = rownames(colData),
                 levels = levels(colData$condition),
                 optionData = NULL)
       x <- structure(x, class = c("glmDataSet"))
   }
   else {
       if (is.null(counts)) {
           counts <- as.matrix(mcols(GR))
           mcols(GR) <- NULL
       }
       # To make sure information is not redundant
       if (!is.null(counts) && !is.null(GR)) mcols(GR) <- NULL

       # The final output
       x <- list(GR = GR, counts = counts, colData = colData,
                 sampleNames = rownames(colData),
                 levels = levels(colData$condition),
                 optionData = NULL)
       x <- structure(x, class = c("RangedGlmDataSet"))
   }
   return(x)
}


