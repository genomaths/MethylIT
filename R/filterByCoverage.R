#' @rdname filterByCoverage
#'
#' @title Filter methylation counts by coverage
#' @description The function is used to discard the cytosine positions with
#'     coverage values less than 'min.coverage' read counts or values greater
#'     than the specified 'percentile'.
#' @details The input must be a GRanges object or list of GRanges objects with a
#'     coverage column in the meta-column table or the columns with methylated
#'     (mC) and unmethylated counts (uC).
#'
#' @param x GRanges object or list of GRanges
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
#' @param verbose If TRUE, prints the function log to stdout
#'
#' @return The input GRanges object or list of GRanges objects after filtering
#'     them.
#'
#' @examples
#' gr1 <- makeGRangesFromDataFrame(data.frame(chr = 'chr1', start = 11:15,
#' end = 11:15, strand = c('+','-','+','*','.'), mC = 1, uC = 1:5),
#' keep.extra.columns = TRUE)
#'
#' filterByCoverage(gr1, min.coverage = 1, max.coverage = 4,
#' col.names = c(mC = 1, uC = 2), verbose = FALSE)
#'
#' @importFrom S4Vectors mcols
#'
#' @export
filterByCoverage <- function(x, min.coverage = 4, max.coverage = Inf,
    percentile = 0.999, col.names = c(coverage = NULL,
        mC = NULL, uC = NULL), verbose = TRUE) {

    cn <- names(col.names)
    if ((!is.element("mC", cn) && !is.element("uC",
        cn)) && !is.element("coverage", cn))
        stop(paste("* Provide the number of the 'coverage' column or the ",
            "numbers for columns 'mC' & 'uC'"))
    if ((inherits(x, "list") || inherits(x, "GRangesList")) &&
        inherits(x[[1]], "GRanges")) {
        sn <- names(x)
        x <- lapply(seq_len(length(x)), function(k) {
            filterGRange(x = x[[k]], min.coverage = min.coverage,
                max.coverage = max.coverage, percentile = percentile,
                col.names = col.names, sample.name = sn[[k]],
                verbose = verbose)
        })
        names(x) <- sn
    } else {
        if (inherits(x, "GRanges")) {
            x <- filterGRange(x, min.coverage = min.coverage,
                max.coverage = max.coverage, percentile = percentile,
                col.names = col.names, sample.name = "",
                verbose = verbose)
        } else {
            mgs = paste0(" The object provided is not a list of GRanges ",
                "objects or a GRangesList")
            stop(mgs)
        }
    }
    return(x)
}
