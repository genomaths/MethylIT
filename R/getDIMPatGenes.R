#' @rdname getDIMPatGenes
#' @name getDIMPatGenes
#' @title Count DMPs at gene-body
#' @description The function counts DMPs overlapping with gene-body. In fact,
#'     this function also can be used to count DMPs overlapping with any set of
#'     regions given in a GRanges object.
#'
#' @param GR An objects object from the any of the classes: 'pDMP', 'InfDiv',
#'     GRangesList, GRanges or a list of GRanges.
#' @param GENES A GRanges object with gene coordinates and gene IDs. A
#'     column named \strong{'gene_id'} carrying the gene ids should be
#'     included in the metacolumns. If the meta-column named 'gene_id'
#' is not provided, then gene (region) ids will be created using the
#' gene (region) coordinates.
#' @param output Class of the object to be returned, a "list", or a "GRanges"
#'  object.
#' @param ignore.strand,type Same as for
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}.
#' @param ... optional arguments for
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}. Users must evaluate
#'  whether specific setting makes sense on each particular context.
#' @seealso \code{\link{getDMPatRegions}}
#' @return A a list GRanges object.
#'
#' @examples
#' ## Gene annotation
#' genes <- GRanges(seqnames = '1',
#' ranges = IRanges(start = c(3631, 6788, 11649), end = c(5899, 9130, 13714)),
#' strand = c('+', '-', '-'))
#' mcols(genes) <- data.frame(gene_id = c('AT1G01010', 'AT1G01020',
#' 'AT1G01030'))
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
getDIMPatGenes <- function(GR, GENES, type = "within",
                            ignore.strand = TRUE,
                            output = c("list", "GRanges"), ...)
                        UseMethod("getDIMPatGenes")

#' @rdname getDIMPatGenes
#' @importFrom S4Vectors mcols
#' @export
getDIMPatGenes.default <- function(GR, GENES, type = "within",
                                    ignore.strand = TRUE, ...) {
    gene_id <- GENES$gene_id
    if (any(is.na(gene_id))) {
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
    Hits <- findOverlaps(GR, GENES, type = type,
                        ignore.strand = ignore.strand, ...)
    if (length(Hits) > 0) {
        DIMP <- GENES[subjectHits(Hits)]
        DIMP <- data.table(as.data.frame(DIMP))
        DIMP <- DIMP[, list(seqnames = unique(seqnames), start = min(start),
                            end = max(end), DIMPs = length(start)),
                    by = gene_id]
        DIMP <- data.frame(DIMP)
        DIMP <- makeGRangesFromDataFrame(DIMP, keep.extra.columns = TRUE)
        Hits <- findOverlaps(DIMP, GENES, type = type,
                            ignore.strand = ignore.strand, ...)
        GENES <- GENES[subjectHits(Hits), ]
        DIMP <- as.data.frame(DIMP[queryHits(Hits), ])
        mcols(GENES) <- data.frame(GeneID = GENES$gene_id, DIMPs = DIMP$DIMPs)
    } else {
        GENES <- GRanges()
        mcols(GENES) <- data.frame(GeneID = factor(),
            DIMPs = integer())
    }

    return(unique(GENES))
}

#' @rdname getDIMPatGenes
#' @importFrom S4Vectors mcols
#' @export
getDIMPatGenes.GRanges <- function(GR, GENES, type = "within",
                                    ignore.strand = TRUE, ...) {
    vn <- c("hdiv", "TV", "wprob")
    ns <- colnames(mcols(GR))
    nams <- sum(is.element(vn, ns))
    if (nams != 3) {
        warning("'GRanges' metacolumn has incorrect column names")
        cat("\n")
        stop("A 'GRanges' metacolumn must have 'hdiv', 'TV', and 'wprob'",
            " named columns to be used as an argument for 'getDIMPatGenes'")
    }

    GR <- getDIMPatGenes.default(GR, GENES = GENES, type = type,
                                ignore.strand = ignore.strand, ...)
    return(GR)
}

#' @rdname getDIMPatGenes
#' @export
getDIMPatGenes.pDMP <- function(GR, GENES, type = "within",
                                ignore.strand = TRUE,
                                output = c("list", "GRanges"), ...) {
    output <- match.arg(output)
    gene_id <- GENES$gene_id
    if (any(is.na(gene_id))) gene_id <- NULL

    dmps <- lapply(GR, getDIMPatGenes.default, GENES = GENES, type = type,
                    ignore.strand = ignore.strand, ...)

    if (output == "GRanges") {
        nams <- names(dmps)
        dmps <- uniqueGRanges(dmps, columns = 2L,
                              ignore.strand = ignore.strand,
                              type = 'equal', verbose = FALSE)
        colnames(mcols(dmps)) <- nams

        if (!is.null(gene_id)) {
            hits <- findOverlaps(dmps, GENES, ignore.strand = ignore.strand,
                                 type = 'equal')
            dmps <- dmps[queryHits(hits)]
            names(dmps) <- GENES$gene_id[subjectHits(hits)]
        }
    }

    return(dmps)
}

#' @rdname getDIMPatGenes
#' @export
getDIMPatGenes.InfDiv <- function(GR, GENES, type = "within",
                                ignore.strand = TRUE,
                                output = c("list", "GRanges"), ...) {
    output <- match.arg(output)
    gene_id <- GENES$gene_id
    if (any(is.na(gene_id))) gene_id <- NULL

    dmps <- lapply(GR, getDIMPatGenes.default, GENES = GENES, type = type,
                    ignore.strand = ignore.strand, ...)

    if (output == "GRanges") {
        nams <- names(dmps)
        dmps <- uniqueGRanges(dmps, columns = 2L,
                              ignore.strand = ignore.strand,
                              type = 'equal', verbose = FALSE)
        colnames(mcols(dmps)) <- nams

        if (!is.null(gene_id)) {
            hits <- findOverlaps(dmps, GENES, ignore.strand = ignore.strand,
                                 type = 'equal')
            dmps <- dmps[queryHits(hits)]
            names(dmps) <- GENES$gene_id[subjectHits(hits)]
        }
    }
    return(dmps)
}

#' @rdname getDIMPatGenes
#' @export
getDIMPatGenes.list <- function(GR, GENES, type = "within",
                                ignore.strand = TRUE,
                                output = c("list", "GRanges"), ...) {
    output <- match.arg(output)
    gene_id <- GENES$gene_id
    if (any(is.na(gene_id))) gene_id <- NULL

    dmps <- lapply(GR, getDIMPatGenes.default, GENES = GENES, type = type,
                    ignore.strand = ignore.strand, ...)

    if (output == "GRanges") {
        nams <- names(dmps)
        dmps <- uniqueGRanges(dmps, columns = 2L,
                            ignore.strand = ignore.strand,
                            type = 'equal', verbose = FALSE)
        colnames(mcols(dmps)) <- nams

        if (!is.null(gene_id)) {
            hits <- findOverlaps(dmps, GENES, ignore.strand = ignore.strand,
                                type = 'equal')
            dmps <- dmps[queryHits(hits)]
            names(dmps) <- GENES$gene_id[subjectHits(hits)]
        }
    }

    return(unique(dmps))
}


