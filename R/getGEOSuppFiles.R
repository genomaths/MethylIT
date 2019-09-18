#' @rdname getGEOSuppFiles
#'
#' @title Get Supplemental Files from GEO
#' @description Decompress 'gzip' files.
#' @details Download supplemental files from a specified GEO dataset. This
#'     function is originally provided in the Bioconductor package 'GEOquery'.
#'     The original function download all the supplemental files for a given GEO
#'     accession number. Herein small detail is added to permit only the
#'     download of the specified files and from several GEO accession numbers
#'     with only one call to the function.
#'
#' @param GEO A character vector with GEO accession numbers.
#' @param makeDirectory Logic (FALSE). If GEO accession number is provided,
#'     whether to create a subdirectory for the downloaded files.
#' @param baseDir Directory where files are downloads if GEO accession number
#'     is provided. Default is the current working directory.
#' @param pattern A pattern for the name of the supplementary files from the
#'     GEO dataset. If provided, then only the files with the given pattern are
#'     downloaded. Otherwise, all the supplementary files are downloaded.
#' @param verbose If TRUE, prints the function log to stdout
#' @return A data frame is returned invisibly with rownames representing the
#'     full path of the resulting downloaded files and the records in the
#'     data.frame the output of file.info for each downloaded file.
#' @examples
#' ## Download supplementary files from GEO data set and store "fullpath/name"
#' ## in variable filename. The parameter 'pattern' permits us to download only
#' ## the specified filesCG, in this case, CG and CHG methylation contexts.
#'
#' filenames <- getGEOSuppFiles(GEO = "GSM881757",
#'                 pattern = "G_cytosine.txt.gz")
#'
#' file.remove(filenames) ## Remove the downloaded file
#'
#' @importFrom RCurl getURL
#' @importFrom utils read.table download.file
#' @importFrom XML htmlParse xpathSApply free
#'
#' @export
#'
#' @author Original author: Sean Davis <sdavis2@mail.nih.gov>
getGEOSuppFiles <- function(GEO, makeDirectory=FALSE, baseDir=getwd(),
                            pattern=NULL, verbose=TRUE) {

   getGEO <- function(geo, makeDir, baseDir, pattern, verbose) {
       getDirListing <- function(url, verbose) {
           if (verbose) message(url)
           a <- RCurl::getURL(url)
           if (grepl("<HTML", a, ignore.case=TRUE)) {
               doc <- XML::htmlParse(a)
               links <- XML::xpathSApply(doc, "//a/@href")
               XML::free(doc)
               b <- as.character(links)[-1]
               if (verbose) message("OK")
           } else {
               tmpcon <- textConnection(a, "r")
               b <- read.table(tmpcon)
               close(tmpcon)
           }
           return(b)
       }

       geotype <- toupper(substr(geo, 1, 3))
       storedir <- baseDir
       fileinfo <- c()
       stub <- gsub("\\d{1,3}$", "nnn", geo, perl=TRUE)
       if (geotype == "GSM") {
           url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/",
                               "%s/suppl/"), stub, geo)
       }
       if (geotype == "GSE") {
           url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/",
                               "%s/suppl/"), stub, geo)
       }
       if (geotype == "GPL") {
           url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s",
                                 "/%s/suppl/"), stub, geo)
       }
       fnames <- try(getDirListing(url=url, verbose=verbose), silent=TRUE)
       if (inherits(fnames, "try-error")) {
           message("No supplemental files found.")
           message("Check URL manually if in doubt")
           message(url)
           return(NULL)
       }

       if (!is.null(pattern)) {
           idx <- if (length(pattern) > 1 ) {
                   unlist(lapply(pattern, function(p){
                   grep(p, fnames)}))} else idx = grep(pattern, fnames)
           fnames <- fnames[idx]
       }

       if (makeDirectory) {
           suppressWarnings(dir.create(storedir <- file.path(baseDir, geo)))
       }
       for (i in fnames) {
           download.file(file.path(url, i), destfile = file.path(storedir,i),
                   mode = "wb", quiet = ifelse(verbose, FALSE, TRUE),
                   method = getOption("download.file.method.GEOquery"))
           ## fileinfo[[file.path(storedir, i)]] <-
           ## file.info(file.path(storedir, i))
           fileinfo <- c(fileinfo, file.path(storedir, i))
       }
       return(fileinfo)
   }
   if (verbose) cat("*** Downloading files ...\n")
   x <- lapply(GEO, function(x) {getGEO(geo=x, makeDir=makeDirectory,
                                       baseDir=baseDir, pattern=pattern,
                                       verbose=verbose)})
   return(unlist(x))
}


