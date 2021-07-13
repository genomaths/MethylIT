#' @rdname getDMPatRegions
#' @name getDMPatRegions
#' @title Count DMPs at Genomic Regions
#' @description The function counts DMPs overlapping with genomic regions.
#' In fact, this function operates as \code{\link{getDIMPatGenes}} function,
#' but without the restrictions set for GRanges objects derived from MethylIT
#' pipeline.
#'
#' @param GR An objects object from the any of the classes: 'pDMP', 'InfDiv',
#'     GRangesList, GRanges or a list of GRanges.
#' @param regions A GRanges object with gene coordinates and gene IDs. A
#'     meta-column named 'gene_id' carying the gene ids should be included. If
#'     the meta-column named 'gene_id' is not provided, then gene (region) ids
#'     will be created using the gene (region) coordinates.
#' @param only.hypo,only.hyper logical(1). Whether to select only
#' hypo-methylated or hyper-methylated cytosine sites.
#' @param ignore.strand,type Same as for
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}
#' @param ... optional arguments for
#'  \code{\link[GenomicRanges]{findOverlaps-methods}}. Users must evaluate
#'  whether specific setting makes sense on each particular context.
#' @seealso \code{\link{getDIMPatGenes}}
#' @return A GRanges object
#'
#' @examples
#' ## Gene annotation
#' genes <- GRanges(seqnames = '1',
#'                 ranges = IRanges(start = c(3631, 6788, 11649),
#'                 end = c(5899, 9130, 13714)),
#'                 strand = c('+', '-', '-'))
#' mcols(genes) <- data.frame(gene_id = c('AT1G01010', 'AT1G01020',
#'                                         'AT1G01030'))
#'
#' ## Get a dataset of dmps from the package
#' data(dmps)
#'
#' ## Finally DMPs found on genes
#' dmrs <- getDMPatRegions(GR = dmps, regions = genes)
#'
#' @importFrom GenomicRanges GRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors subjectHits mcols mcols<-
#' @importFrom data.table data.table
#' @importFrom rtracklayer import
#' @export
getDMPatRegions <- function(
                            GR,
                            regions,
                            only.hypo = FALSE,
                            only.hyper = FALSE,
                            type = "within",
                            ignore.strand = TRUE,
                            ...)
                        UseMethod("getDMPatRegions")

#' @rdname getDMPatRegions
#' @importFrom S4Vectors mcols
#' @export
getDMPatRegions.default <- function(
                                    GR,
                                    regions,
                                    only.hypo = FALSE,
                                    only.hyper = FALSE,
                                    type = "within",
                                    ignore.strand = TRUE,
                                    ...) {

    if (!inherits(GR, "GRanges"))
        stop("*** The 'GR' argument must inherit from 'GRanges-class'.")

    if (!is.null(GR$TV) && is.numeric(GR$TV)) {
        if (only.hypo)
            GR <- GR[ GR$TV < 0]
        if (only.hyper)
            GR <- GR[ GR$TV > 0]
    }

    gene_id <- regions$gene_id
    if (any(is.na(gene_id))) {
        warnings("At least one gene ID is NA. Using gene coordinates as IDs")
        gene_id <- NULL
    }
    if (is.null(gene_id)) {
        chr = seqnames(regions)
        starts = start(regions)
        ends = end(regions)
        strands = strand(regions)
        regions$gene_id <- paste(chr, starts, ends, strands, sep = "_")
    }
    Hits <- findOverlaps(GR, regions, type = type,
                        ignore.strand = ignore.strand, ...)
    if (length(Hits) > 0) {
        DMP <- regions[subjectHits(Hits)]
        DMP <- data.table(as.data.frame(DMP))
        DMP <- DMP[, list(seqnames = unique(seqnames),
                        start = min(start), end = max(end),
                        strand = unique(strand),
                        DMPs = length(start)),
            by = gene_id]
        DMP <- data.frame(DMP)
        DMP <- makeGRangesFromDataFrame(DMP, keep.extra.columns = TRUE)
        Hits <- findOverlaps(DMP, regions, type = type,
                            ignore.strand = ignore.strand, ...)
        regions <- regions[subjectHits(Hits), ]
        DMP <- as.data.frame(DMP[queryHits(Hits), ])
        mcols(regions) <- data.frame(id = regions$gene_id, dmps = DMP$DMPs)
    } else {
        regions <- GRanges()
        mcols(regions) <- data.frame(id = factor(), dmps = integer())
    }
    return(unique(regions))
}

#' @rdname getDMPatRegions
#' @importFrom S4Vectors mcols
#' @export
getDMPatRegions.GRanges <- function(
                                    GR,
                                    regions,
                                    only.hypo = FALSE,
                                    only.hyper = FALSE,
                                    type = "within",
                                    ignore.strand = TRUE,
                                    ...) {
    GR <- getDMPatRegions.default(
                                    GR,
                                    regions = regions,
                                    only.hypo = only.hypo,
                                    only.hyper = only.hyper,
                                    type = type,
                                    ignore.strand = ignore.strand,
                                    ...)
    return(GR)
}

#' @rdname getDMPatRegions
#' @export
getDMPatRegions.pDMP <- function(
                                GR,
                                regions,
                                only.hypo = FALSE,
                                only.hyper = FALSE,
                                type = "within",
                                ignore.strand = TRUE,
                                ...) {
    return(lapply( GR,
                   getDMPatRegions.default,
                   regions = regions,
                   only.hypo = only.hypo,
                   only.hyper = only.hyper,
                   type = type,
                   ignore.strand = ignore.strand,
                   keep.attr = TRUE,
                   ...))
}

#' @rdname getDMPatRegions
#' @export
getDMPatRegions.InfDiv <- function(
                                    GR,
                                    regions,
                                    only.hypo = FALSE,
                                    only.hyper = FALSE,
                                    type = "within",
                                    ignore.strand = TRUE,
                                    ...) {
    return(lapply(  GR,
                    getDMPatRegions.default,
                    regions = regions,
                    only.hypo = only.hypo,
                    only.hyper = only.hyper,
                    type = type,
                    ignore.strand = ignore.strand,
                    keep.attr = TRUE,
                    ...))
}

#' @rdname getDMPatRegions
#' @export
getDMPatRegions.list <- function(
                                GR,
                                regions,
                                only.hypo = FALSE,
                                only.hyper = FALSE,
                                type = "within",
                                ignore.strand = TRUE,
                                ...) {
    return(lapply(  GR,
                    getDMPatRegions.default,
                    regions = regions,
                    only.hypo = only.hypo,
                    only.hyper = only.hyper,
                    type = type,
                    ignore.strand = ignore.strand,
                    keep.attr = TRUE,
                    ...))
}


