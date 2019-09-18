#' @rdname uniqueGRanges
#'
#' @title Unique genomic ranges from a list of GRanges objects
#' @description Build an unique GRanges object from a list of Granges objects.
#' @details The metadata of each one of these GRanges must have one or more
#'     columns to yield a unique GRanges object with metadata columns from the
#'     original GRanges objects. Otherwise, a unique GRanges object will be
#'     created  without metadata columns. Additionally, all metadata must be the
#'     same class, e.g. all numeric or all characters, or all factor
#' @param ListOfGranges Objects to combine. A list of GRanges object or a
#'     GRangesList object.
#' @param ncols integer. Number of columns to use from the meta-column of each
#'     GRanges object. Default value: NULL. If NULL, all the columns (from
#'     column 1 to ncols) from each GRanges will be present in the uniqueGRanges
#'     output.
#' @param columns interger number(s) corresponding to the specific column(s) to
#'     use from the meta-column of each GRanges. Default value: NULL. if
#'     provided, the metacolumn from the uniqueGRanges output will contain the
#'     specified columns.
#' @param chromosomes Chromosomes used Default value: NULL
#' @param maxgap See GenomicRanges::findOverlaps in the IRanges package for a
#'     description of these arguments Default value: -1L
#' @param minoverlap See GenomicRanges::findOverlaps in the IRanges package for
#'     a description of these arguments Default value: 1L
#' @param missing A numerical value (default 0) or NA to write in ranges with
#'     missing values. For example, suppose that we want to build a
#'     uniqueGRanges object from the GRanges objects X and Y. If a given range
#'     k from the GRanges object X with metacolum value x is missing in the
#'     GRanges object Y, then the metacolum of range k from
#'     uniqueGRanges(list(X,Y)) object will be the row vector (x,0) or (x,NA)
#'     if missing = NA.
#' @param type By default, any overlap is accepted. By specifying the type
#'     parameter, one can select for specific types of overlap. The types
#'     correspond to operations in Allen's Interval Algebra (see references).
#'     If type is start or end, the intervals are required to have matching
#'     starts or ends, respectively. While this operation seems trivial, the
#'     naive implementation using outer would be much less efficient. Specifying
#'     equal as the type returns the intersection of the start and end matches.
#'     If  type is within, the query interval must be wholly contained within
#'     the subject interval. Note that all matches must additionally satisfy the
#'     minoverlap constraint described above. The maxgap parameter has special
#'     meaning with the special overlap types. For start, end, and equal, it
#'     specifies the maximum difference in the starts, ends or both,
#'     respectively. For within, it is the maximum amount by which the subject
#'     may be wider than the query.
#' @param select When select is "all" (the default), the results are returned as
#'     a Hits object. Otherwise the returned value is an integer vector parallel
#'     to query (i.e. same length) containing the first, last, or arbitrary
#'     overlapping interval in subject, with NA indicating intervals that did
#'     not overlap any intervals in subject.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#'     the overlap calculations. Default value: TRUE
#' @param keep.strand When set to TRUE, the strand information is preserved on
#'     the objects even if ignore.strand is set to TRUE.  This makes it possible
#'     to ignore the strand during overlap calculations but to preserve the
#'     strand information and not overwrite with *. Default value is
#'     keep.strand = !ignore.strand.
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see bplapply function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#' @param verbose if TRUE, prints the function log to stdout
#'
#' @return a GRanges object
#'
#' @examples
#' dfChr1 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                     strand = c("+","-","+","*","."), score = 1:5)
#' dfChr2 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                     strand = c("+","-","+","*","."), score = 1:5)
#' dfChr3 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
#'                     strand = c("+","-","+","*","."), score = 1:5)
#'
#' gr1 <- makeGRangesFromDataFrame(dfChr1, keep.extra.columns = TRUE)
#' gr2 <- makeGRangesFromDataFrame(dfChr2, keep.extra.columns = TRUE)
#' gr3 <- makeGRangesFromDataFrame(dfChr3, keep.extra.columns = TRUE)
#'
#' grList <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
#'
#' uniqueGRanges(grList)
#'
#' @importFrom S4Vectors subjectHits mcols queryHits
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames seqlevels
#' @importFrom BiocGenerics strand end strand<-
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps
#' @importFrom S4Vectors mcols queryHits subjectHits mcols<-
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#'
#' @export
uniqueGRanges <- function(ListOfGranges, ncols=NULL, columns=NULL,
                       chromosomes=NULL, maxgap=-1L, minoverlap=1L, missing=0,
                       type=c("any", "start", "end", "within", "equal"),
                       select=c("all", "first", "last", "arbitrary"),
                       ignore.strand=FALSE, keep.strand=!ignore.strand,
                       num.cores=1, tasks=0L, verbose=TRUE) {

   if (class(ListOfGranges) == "list" && class(ListOfGranges) != "GRangesList")
   {
       GR <- try(as(ListOfGranges, "GRangesList"), silent=TRUE)

       if (inherits(ListOfGranges, "try-error")) {
         numgr <- sum(unlist(lapply(ListOfGranges, function(l)
                                               class(l) == "GRanges")))
           if (numgr != length(ListOfGranges)) {stop(
               "Not all the elements from the list are valid GRanges objects")
           } else {ListOfGranges <- GR; rm(GR); gc()}
       }
   }

   # Set parallel computation
   if (Sys.info()['sysname'] == "Linux") {
       bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
   } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")

   ##samples = names(ListOfGranges)
   unlistfn <- function(x) {
       matrix(unlist(x), ncol=3, byrow=TRUE)
   }
   mzeros <- function(l, nc, seq, hits) {
       if (nc > 1) {
           m <- matrix(missing, l, nc)
           m[subjectHits(hits), 1:nc] <- as.matrix(mcols(seq)[queryHits(hits),
                                                         1:nc])
       } else{
           m <- rep(missing, l)
           seq <- as.vector(mcols(seq)[queryHits(hits), 1])
           if (class(seq) == "factor" || class(seq) == "character") {
               m[subjectHits(hits)] <- as.character(seq)
           } else m[subjectHits(hits)] <- as.numeric(seq)
       }
       m
   }

   if (is.null(chromosomes)) {
       chromosomes <- try(lapply(ListOfGranges, seqlevels), silent=TRUE)
       if (inherits(chromosomes, "try-error")) {
           warning("* Chromosome labels could not be determined. \n",
                   "Arbitrary Arabidopsis thaliana labels are used instead \n")
           chromosomes <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
       }
       chromosomes <- unique(as.character(unlist(chromosomes)))
   }

   cond1 <- all(unlist(lapply(ListOfGranges, function(x) !is.null(mcols(x)))))
   cond2 <- all(unlist(lapply(ListOfGranges, function(x) ncol(mcols(x)) > 0)))

   ## It assigns the data to the GRanges metacolumns
   if (verbose)
       message(" *** Building coordinates for the new GRanges object ..." )

   if (ignore.strand && !keep.strand) {
       if (verbose) {
           message(" *** Setting strand information to * ..." )
       }

       ListOfGranges <- lapply(ListOfGranges, function(GR) {
           strand(GR) <- "*"
           return(GR)
       })
   } else {
       if (verbose) {
           message(" *** Strand information is preserved ..." )
       }
   }

   granges <- unique(unname(unlist(ListOfGranges)))

   if (cond1 && cond2) {
       ## This assigns the data values to the GRanges metacolumns
       l <- length(granges)
       n <- length(ListOfGranges)

       if (verbose) message(" *** Processing GRanges for sample: #", 1 )
       seq <- ListOfGranges[[1]]
       if (is.null(ncols) && is.null(columns)) {
           snames <- colnames(mcols(seq))
           ncol <- length(snames)
       } else {
           if (!is.null(columns)) {
               snames <- colnames(mcols(seq))[columns]
               ncol = length(snames)
               seq = seq[, snames]
           }
           if (!is.null(ncols)) {
               ncol <- ncols
               snames <- colnames(mcols(seq))[1:ncol]
           }
       }

       Hits <- findOverlaps(seq, granges, maxgap=maxgap, minoverlap=minoverlap,
                       type=type, select=select, ignore.strand=ignore.strand)
       x <- mzeros(l, nc=ncol, seq, Hits)
       methl <- data.frame(x, stringsAsFactors=FALSE)
       rm(x, seq, Hits); gc()

       for (k in 2:n) {
           seq <- ListOfGranges[[k]]
           if (is.null(ncols) && is.null(columns)) {
               nam <- colnames(mcols(seq))
               ncol <- length(nam)
           } else {
               if (!is.null(columns)) {
                   nam <- colnames(mcols(seq))[columns]
                   ncol = length(nam)
                   seq = seq[, nam]
               }
               if (!is.null(ncols)) {
                   ncol <- ncols
                   nam <- colnames(mcols(seq))[1:ncol]
               }
           }

           snames <- c(snames, nam)
           if (verbose)
               message(" *** Processing GRanges for sample: #", k, "...")
           Hits <- findOverlaps(seq, granges, maxgap=maxgap,
                           minoverlap=minoverlap, type=type, select=select,
                           ignore.strand=ignore.strand)
           x <- mzeros(l, nc=ncol, seq, Hits)
           methl <- cbind(methl, x)
           rm(x, seq, Hits); gc()
       }
       colnames(methl) <- snames
       mcols(granges) <- methl
       rm(methl); gc()
       colnames(mcols(granges)) <- make.unique(colnames(mcols(granges)))
   } else warnings("At least one GRanges has zero metacolumns")

   if (verbose) message(" *** Sorting by chromosomes and start positions...")
   granges <- sortBySeqnameAndStart(granges)
   return(granges)
}

