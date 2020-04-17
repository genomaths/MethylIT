#' @rdname sortBySeqnameAndStart
#'
#' @title Sorting 'GRanges' objects
#' @description Sorts a GRanges object by seqname and start position
#'
#' @param gr GRanges object
#'
#' @return GRanges object
#'
#' @examples
#' GR <- as(c('chr2:1-1', 'chr1:1-1'), 'GRanges')
#' GR <- sortBySeqnameAndStart(GR)
#'
#' @importFrom BiocGenerics sort start
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqnames
#'
#' @aliases sortBySeqnameAndStart
#' @export
sortBySeqnameAndStart <- function(gr) {
    seqlevels(gr) <- sort(seqlevels(gr))
    return(gr[order(as.factor(seqnames(gr)), start(gr)),])
}

#' @rdname sortBySeqnameAndStart
#' @aliases sortBySeqnameAndEnd
#' @export
sortBySeqnameAndEnd <- function(gr) {
  seqlevels(gr) <- sort(seqlevels(gr))
  return(gr[order(as.factor(seqnames(gr)), end(gr)),])
}
