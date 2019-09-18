#' @rdname readCounts2GRangesList
#'
#' @title Read files of methylation count tables
#' @description This function is addressed to read files with methylation count
#'     table data commonly generated after the alignment of BS-seq data or found
#'     in GEO database
#' @details Read tables from files with a table methylation count data using
#'     the function fread from the package 'data.table' and and yields a list of
#'     GRanges objects with the information provided.
#'
#' @param filenames Character vector with the file names
#' @param sample.id Character vector with the names of the samples
#'     corresponding to each file
#' @param pattern Chromosome name pattern. Users working on Linux OS can
#'     specify the reading of specific lines from each file by using regular
#'     expressions.
#' @param remove Logic (TRUE). Usually the supplementary files from GEO
#'     datasets are 'gz' compressed. File datasets must be decompressed to be
#'     read. The decompressed files are removed after read if this is set
#'     'TRUE'.
#' @param columns Vector of integer numbers denoting the table columns that
#'     must be read. The numbers for 'seqnames' (chromosomes), 'start', and
#'     'end' (if different from 'start') columns must be given. The possible
#'     fields are: 'seqnames' (chromosomes),'start', 'end', 'strand',
#'     'fraction', percent' (metylation percentage), 'mC' (methylates cytosine),
#'     'uC' (non methylated cytosine), 'coverage', and 'context'
#'     (methylation context). These column headers are not required to be in the
#'     files.
#' @param chromosome.names If provided, for each GRanges object, chromosome
#'     names will be changed to those provided in 'chromosome.names' applying
#'     seqlevels(x) <- chromosome.names'. This option permits to use all the
#'     functionality of the function "seqlevels" defined from package
#'     "GenomeInfoDb", which rename, add, and reorder the seqlevels all at once
#'     (see ?seqlevels).
#' @param chromosomes If provided, it must be a character vector with the names
#'     of the chromosomes that you want to include in the final GRanges objects.
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Additional parameters for 'fread' function from 'data.table'
#'     package
#'
#' @return A list of GRanges objects
#'
#' @examples
#' ## Create a cov file with it's file name including "gz" (tarball extension)
#' filename <- "./file.cov"
#' gr1 <- data.frame(chr = c("chr1", "chr1"), post = c(1,2),
#'                 strand = c("+", "-"), ratio = c(0.9, 0.5),
#'                 context = c("CG", "CG"), CT = c(20, 30))
#' filename <- "./file.cov"
#' write.table(as.data.frame(gr1), file = filename,
#'             col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' ## Read the file. It does not work. Typing mistake: "fractions"
#' LR <- try(readCounts2GRangesList(filenames = filename, remove = FALSE,
#'                             sample.id = "test",
#'                             columns = c(seqnames = 1, start = 2,
#'                                     strand = 3, fractions = 4,
#'                                     context = 5, coverage = 6)),
#'                                     silent = TRUE)
#' file.remove(filename) # Remove the file
#'
#' @importFrom data.table fread
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
readCounts2GRangesList <- function(filenames=NULL, sample.id=NULL, pattern=NULL,
                   remove=FALSE, columns=c(seqnames=NULL, start=NULL, end=NULL,
                   strand=NULL,fraction=NULL, percent=NULL, mC=NULL, uC=NULL,
                   coverage=NULL, context=NULL), chromosome.names=NULL,
                   chromosomes=NULL, verbose=TRUE, ...) {

   colns <- c("seqnames", "start", "end", "strand", "fraction",
            "percent", "mC", "uC", "coverage", "context")
   Check <- ArgumentCheck::newArgCheck()
   if (is.null(filenames)) {
       ArgumentCheck::addError(msg="No file names provided",
                           argcheck=Check)
   }
   for (file in filenames) {
       if (!file.exists(file)) {
           ArgumentCheck::addError(msg=paste0("Unable to find: ",
                                           file), argcheck=Check)
       }
   }
   cn <- names(columns)
   if (!is.element("seqnames", cn) || !is.element("start", cn)) {
       ArgumentCheck::addError(msg=paste0("You must provide the numbers for ",
                   "'seqnames' (chromosomes names), ", "  'start' columns"),
                   argcheck=Check)
   }

   gz.ext = "[.]gz$"
   if (grepl(gz.ext, "", filenames) &&
       sum(is.element(filenames, sub(gz.ext, "", filenames))) > 0) {
       ArgumentCheck::addError(msg=paste0("File duplication. File: \n",
                   filenames[is.element(filenames, sub(gz.ext, "", filenames))],
                   "\n is also provided in compressed format '.gz' \n"),
                   argcheck=Check)
   }

   idx <- is.element(names(columns), colns)
   if (sum(idx) != length(columns)) {
       ArgumentCheck::addError(msg=paste0("Probably you have a typing mistake ",
           "in the column names: ", "'", names(columns)[!idx],"'"),
           argcheck=Check)
   }

   ArgumentCheck::finishArgCheck(Check)

   meth <- lapply(filenames, function(file) {
       if (verbose) message("*** Processing file: ", file)
       gunzippedfile <- sub(gz.ext, "", file)
       if (grepl(pattern=gz.ext, file)) {
           if (file.exists(gunzippedfile)) {
               warning(paste("A file with the name: \n", gunzippedfile,
                   "\n already exists. So, the corresponding file '.gz' was",
                   "not decompressed"), call. = FALSE)
           } else .gunzip(file, remove=FALSE)
       }
       if (.Platform$OS.type == "unix") {
           if (!is.null(pattern)) {
               x <- suppressMessages(
                           fread(paste0("grep ", pattern, " ", gunzippedfile),
                               select=columns, ...))
           } else {
               if (!is.null(chromosomes)) {
                   if (length(chromosomes) > 1) {
                       chrom <- paste0("'", paste0(chromosomes, collapse = "|"),
                                       "'")
                   } else chrom <- chromosomes
                   x <- suppressMessages(
                           fread(paste0("egrep ", chrom, " ", gunzippedfile),
                               select=columns, ...))
               } else x <- suppressMessages(fread(gunzippedfile,
                                               select=columns, ...))
           }
       } else x <- suppressMessages(fread(gunzippedfile, select=columns, ...))

       colns <- colnames(x)
       idx <- is.element(colns, cn)
       if (sum(idx) != length(colns)) {
           colns[!idx] <- cn[!idx]
           colnames(x) <- colns
       }

       if (!is.element("end", cn)) x$end <- x$start
       if (is.element("coverage", cn) && is.element("mC", cn)) {
           x$uC <- x$coverage - x$mC
       }
       if (is.element("fraction", cn) && is.element("coverage", cn)) {
           x$mC <- round(x$fraction * x$coverage)
           x$uC <- round(x$coverage) - x$mC
       }
       if (is.element("percent", cn) && is.element("coverage", cn)) {
           x$mC <- round(x$coverage * x$percent / 100)
           x$uC <- round(x$coverage) - x$mC
       }
       if (!is.null(chromosomes)) {
           idx <- is.element(x$seqnames, chromosomes)
           x <- x[idx]
       }
       if (remove) file.remove(gunzippedfile)
       x <- makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
       ## if requested change chromosomes names
       if (!is.null(chromosome.names)) {
           seqlevels(x) <- chromosome.names
       }

       x <- sortBySeqnameAndStart(x)
       return(x)
       })
   if (!is.null(sample.id)) {
       names(meth) <- sample.id
   }
   return(meth)
   }
