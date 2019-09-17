#' @rdname getDIMPatGenes
#' @aliases getDIMPatGenes.default
#' @title Count DIMPs at gene-body
#' @description The function counts DIMPs overlapping with gene-body. In fact,
#'     this function also can be used to count DIMPs overlapping with any set of
#'     regions given in a GRanges object.
#'
#' @param GR An obects object from the any of the classes: 'pDMP', 'InfDiv',
#'     GRangesList, GRanges or a list of GRanges.
#' @param GENES A GRanges object with gene coordinates and gene IDs. A
#'     meta-column named 'gene_id' carying the gene ids should be included. If
#'     the meta-column named 'gene_id' is not provided, then gene (region) ids
#'     will be created using the gene (region) coordinates.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the calculations. Default value: TRUE
#' @param ... Not in use.
#' @return A GRanges object
#'
#' @examples
#' num.points <- 10000 # Number of cytosine position with methylation call
#' ## Gene annotation
#' genes <- GRanges(seqnames = "1",
#'                 ranges = IRanges(start = c(3631, 6788, 11649),
#'                                  end = c(5899, 9130, 13714)),
#'                 strand = c("+", "-", "-"))
#' mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
#'                                        "AT1G01030"))
#'
#' set.seed(123) ## To set a seed for random number generation
#' ## GRanges object of the reference with methylation levels in
#' ## its meta-column
#' Ref <- makeGRangesFromDataFrame(
#'     data.frame(chr = '1',
#'                 start = 1:num.points,
#'                 end = 1:num.points,
#'                 strand = '*',
#'                 p2 = rbeta(num.points, shape1 = 1, shape2 = 1.5)),
#'     keep.extra.columns = TRUE)
#'
#' ## List of Granges objects of individuals methylation levels
#' Indiv <- GRangesList(
#'     sample11 = makeGRangesFromDataFrame(
#'         data.frame(chr = '1',
#'                 start = 1:num.points,
#'                 end = 1:num.points,
#'                 strand = '*',
#'                 p2 = rbeta(num.points, shape1 = 1.5, shape2 = 2)),
#'         keep.extra.columns = TRUE),
#'     sample12 = makeGRangesFromDataFrame(
#'         data.frame(chr = '1',
#'                 start = 1:num.points,
#'                 end = 1:num.points,
#'                 strand = '*',
#'                 p2 = rbeta(num.points, shape1 = 1.6, shape2 = 2.1)),
#'         keep.extra.columns = TRUE),
#'     sample21 = makeGRangesFromDataFrame(
#'         data.frame(chr = '1',
#'                 start = 1:num.points,
#'                 end = 1:num.points,
#'                 strand = '*',
#'                 p2 = rbeta(num.points, shape1 = 10, shape2 = 4)),
#'         keep.extra.columns = TRUE),
#'     sample22 = makeGRangesFromDataFrame(
#'         data.frame(chr = '1',
#'                 start = 1:num.points,
#'                 end = 1:num.points,
#'                 strand = '*',
#'                 p2 = rbeta(num.points, shape1 = 11, shape2 = 4)),
#'         keep.extra.columns = TRUE))
#' ## To estimate Hellinger divergence using only the methylation levels.
#' HD <- estimateDivergence(ref = Ref, indiv = Indiv, meth.level = TRUE,
#'                         columns = 1)
#' ## To perform the nonlinear regression analysis
#' nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)
#'
#' ## Next, the potential signal can be estimated
#' PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 4, alpha = 0.05)
#'
#' ## The cutpoint estimation used to discriminate the signal from the noise
#' cutpoints <- estimateCutPoint(PS, control.names = c("sample11", "sample12"),
#'                             treatment.names = c("sample21", "sample22"),
#'                             div.col = 4, verbose = TRUE)
#' ## DIMPs are selected using the cupoints
#' DIMPs <- selectDIMP(PS, div.col = 4, cutpoint = min(cutpoints$cutpoint))
#'
#' ## Finally DIMPs found on genes
#' DIMR <- getDIMPatGenes(GR = DIMPs$sample21, GENES = genes)
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors subjectHits mcols<-
#' @importFrom data.table data.table
#' @importFrom rtracklayer import
#'
#' @export
getDIMPatGenes <- function(GR, ...) UseMethod("getDIMPatGenes")

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.default
getDIMPatGenes.default <- function(GR, GENES, ignore.strand=TRUE) {
   gene_id <- GENES$gene_id
   if(any(is.na(gene_id))) {
      warnings("At least one gene ID is NA. Using gene coordinates as IDs")
      gene_id <- NULL
   }
   if (is.null(gene_id)) {
       chr = seqnames(GENES)
       starts = start(GENES)
       ends = end(GENES)
       strands = strand(GENES)
       GENES$gene_id <- paste(chr, starts, ends, strands, sep = "_")
   }
   Hits <- findOverlaps(GR, GENES, type="within", ignore.strand=ignore.strand)
   DIMP <- GENES[subjectHits(Hits)]
   DIMP <- data.table( as.data.frame( DIMP ) )
   DIMP <- DIMP[ ,list(seqnames=unique(seqnames),
                       start=min(start), end=max(end),
                       DIMPs=length(start)), by=gene_id]
   DIMP <- data.frame(DIMP)
   DIMP <- makeGRangesFromDataFrame(DIMP, keep.extra.columns=TRUE )
   Hits <- findOverlaps(DIMP, GENES, type="within", ignore.strand=ignore.strand)
   GENES <- GENES[subjectHits(Hits), ]
   DIMP <- as.data.frame(DIMP[queryHits(Hits), ])
   mcols(GENES) <- data.frame(GeneID=GENES$gene_id, DIMPs=DIMP$DIMPs)
   return(unique(GENES))
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.pDMP
getDIMPatGenes.pDMP <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES, keep.attr = TRUE))
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.InfDiv
getDIMPatGenes.InfDiv <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES, keep.attr = TRUE))
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.list
getDIMPatGenes.list <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES ))
}


