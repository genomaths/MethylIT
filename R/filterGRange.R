#' @rdname filterGRange
#'
#' @title Filter methylation counts by coverage in a GRange object
#' @description The function is used to discard the cytosine positions with
#'     coverage values less than 'min.coverage' read counts or values greater
#'     than the specified 'percentile'.
#' @details The input must be a GRanges object with a coverage column in the
#'     metacolumn table or the columns with methylated (mC) and unmethylated
#'     counts (uC).
#'
#' @param x GRanges object
#' @param min.coverage Cytosine sites with coverage less than min.coverage are
#'     discarded. Default: 0
#' @param max.coverage Cytosine sites with coverage greater than max.coverage
#'     are discarded. Default: Inf
#' @param percentile Threshold to remove the outliers from each file and all
#'     files stacked. If percentile is 1, all the outliers stay
#' @param col.names The number of the 'coverage' column. Since no specific table
#'     format for the count data is specified, at least the number of the
#'     'coverage' column must be given, or the number of the columns with
#'     methylated (mC) and unmethylated counts (uC). Then coverage = mC + uC.
#' @param sample.name Name of the sample
#' @param verbose If TRUE, prints the function log to stdout
#'
#' @return  The input GRanges object or list of GRanges objects after filtering
#'     it.
#'
#' @examples
#'     gr1 <- makeGRangesFromDataFrame(
#'         data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                 strand = c("+","-","+","*","."), mC = 1, uC = 1:5),
#'         keep.extra.columns = TRUE)
#'     filterGRange(gr1, min.coverage = 1, max.coverage = 4,
#'                 col.names = c(mC = 1, uC = 2), verbose = FALSE)
#'
#' @importFrom S4Vectors mcols
#'
#' @export
filterGRange <- function(x, min.coverage=4, max.coverage=Inf,
                         percentile=.999,
                         col.names=c(coverage=NULL, mC=NULL, uC=NULL ),
                         sample.name='', verbose=TRUE) {

   if (!(class(x) == 'GRanges')) {
       stop("* Unable to process a non-GRanges object!")
   }
   if (verbose)
       message( "*** Processing the GRanges sample ", sample.name, " ..." )
   if (is.element("coverage", names(col.names))) {
       cov <- mcols(x[ ,col.names['coverage']])[ ,1]
   } else {
       cov <- rowSums(as.matrix(mcols(x[ ,col.names[c('mC','uC')]])))
   }
   idx <- which((cov >= min.coverage) & (cov <= max.coverage))
   cov <- cov[idx]
   x <- x[idx]
   x <- x[cov < quantile(cov, probs=percentile) ]
   return(x)
}
