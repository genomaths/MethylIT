#' @rdname gunzip
#'
#' @title Decompress 'gzip' files
#' @description Decompress 'gzip' files.
#' @details This function is originally from the R packages 'R.utils'. Herein,
#'     the version as given by 'GEOquery' is incorporated.
#'
#' @param filename Character vector with the files names to be decompressed
#' @param destname The destination file.
#' @param overwrite Logic (FALSE) indicating whether or not to overwrite a
#'     destfile of the same name.
#' @param remove Logic (TRUE) indicating whether or not to remove the original
#'     file after completion.
#' @param BFR.SIZE The size of the read buffer.
#' @return Invisibly, the number of bytes read.
#'
#' @author Original author: Henrik Bengtsson.
#' @keywords internal
.gunzip <- function(filename, destname = gsub("[.]gz$", 
    "", filename), overwrite = FALSE, remove = TRUE, 
    BFR.SIZE = 1e+07) {
    
    if (filename == destname) 
        stop(sprintf("Argument 'filename' and 'destname' are identical: %s", 
            filename))
    if (!overwrite && file.exists(destname)) 
        stop(sprintf("File already exists: %s", destname))
    inn <- gzfile(filename, "rb")
    on.exit(if (!is.null(inn)) close(inn))
    out <- file(destname, "wb")
    on.exit(close(out), add = TRUE)
    nbytes <- 0
    repeat {
        bfr <- readBin(inn, what = raw(0), size = 1, 
            n = BFR.SIZE)
        n <- length(bfr)
        if (n == 0) 
            break
        nbytes <- nbytes + n
        writeBin(bfr, con = out, size = 1)
    }
    if (remove) {
        close(inn)
        inn <- NULL
        file.remove(filename)
    }
    invisible(nbytes)
}
