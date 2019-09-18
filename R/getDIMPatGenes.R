#' @rdname getDIMPatGenes
#' @name getDIMPatGene
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
#' ## Gene annotation
#' genes <- GRanges(seqnames = "1",
#'                 ranges = IRanges(start = c(3631, 6788, 11649),
#'                                  end = c(5899, 9130, 13714)),
#'                 strand = c("+", "-", "-"))
#' mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
#'                                        "AT1G01030"))
#'
#' ## Get a dataset of potential signals and the estimated cutpoint from the
#' ## package
#' data(PS, cutpoint)
#'
#' ## The estimated cutpoints are used to discriminate signals from the noise.
#' ## That is, DMPs are selected using the cupoints
#' DIMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)
#'
#' ## Finally DMPs found on genes
#' DIMR <- getDIMPatGenes(GR = DIMPs$T1, GENES = genes)
#'
#' @importFrom GenomicRanges GRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors subjectHits mcols mcols<-
#' @importFrom data.table data.table
#' @importFrom rtracklayer import
#' @export
getDIMPatGenes <- function(GR, GENES, ignore.strand=TRUE)
   UseMethod("getDIMPatGenes")

#' @rdname getDIMPatGenes
#' @importFrom S4Vectors mcols
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
   if (length(Hits) > 0) {
       DIMP <- GENES[subjectHits(Hits)]
       DIMP <- data.table( as.data.frame( DIMP ) )
       DIMP <- DIMP[ ,list(seqnames=unique(seqnames),
                           start=min(start), end=max(end),
                           DIMPs=length(start)), by=gene_id]
       DIMP <- data.frame(DIMP)
       DIMP <- makeGRangesFromDataFrame(DIMP, keep.extra.columns=TRUE )
       Hits <- findOverlaps(DIMP, GENES, type="within",
                           ignore.strand=ignore.strand)
       GENES <- GENES[subjectHits(Hits), ]
       DIMP <- as.data.frame(DIMP[queryHits(Hits), ])
       mcols(GENES) <- data.frame(GeneID=GENES$gene_id, DIMPs=DIMP$DIMPs)
   } else {
       GENES <- GRanges()
       mcols(GENES) <- data.frame(GeneID = factor(), DIMPs = integer())
   }
   return(unique(GENES))
}

#' @rdname getDIMPatGenes
#' @importFrom S4Vectors mcols
#' @exportMethod getDIMPatGenes.GRanges
#' @export
getDIMPatGenes.GRanges <- function(GR, GENES, ignore.strand=TRUE) {
   vn <- c("hdiv", "TV", "wprob")
   ns <- colnames(mcols(GR))
   nams <- sum(is.element(vn, ns))
   if (nams != 3) {
      warning("'GRanges' metacolumn has incorrect column names")
      cat("\n")
      stop("A 'GRanges' metacolumn must have 'hdiv', 'TV', and 'wprob'",
           " named columns to be used as an argument for 'getDIMPatGenes'")
   }

   GR <- getDIMPatGenes.default(GR, GENES = GENES,
                               ignore.strand = ignore.strand)
   return(GR)
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.pDMP
getDIMPatGenes.pDMP <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES,
               ignore.strand = ignore.strand, keep.attr = TRUE))
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.InfDiv
getDIMPatGenes.InfDiv <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES,
               ignore.strand = ignore.strand, keep.attr = TRUE))
}

#' @rdname getDIMPatGenes
#' @exportMethod getDIMPatGenes.list
getDIMPatGenes.list <- function(GR, GENES, ignore.strand=TRUE) {
   return(lapply(GR, getDIMPatGenes.default, GENES = GENES,
               ignore.strand = ignore.strand, keep.attr = TRUE))
}


