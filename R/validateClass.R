#' @rdname validateClass
#' @title Function to validate S3 classes in MethylIT
#' @param LR An object from class 'pDMP' or 'InfDiv'
#' @keywords internal
#' @export
validateClass <- function(LR) UseMethod("validateClass")

#' @rdname validateClass
#' @keywords internal
#' @exportMethod validateClass.default
validateClass.default <- function(LR){
   stop(paste("'validateClass' does not know how to handle object of class",
               class(LR),
               "and can only be used on classes 'pDMP', 'InfDiv'"))

}

#' @rdname validateClass
#' @keywords internal
#' @exportMethod validateClass.pDMP
validateClass.pDMP <- function(LR) {
   vn <- c("hdiv", "TV", "wprob")
   if (any(!unlist(lapply(LR, function(GR) class(GR) == "GRanges")))) {
       warning("At least one element from 'LR' is not a 'GRanges' object")
       cat("\n")
       stop("LR is not a valid 'pDMP' object")
   }
   nams <- unlist(lapply(LR, function(GR) {
       ns <- colnames(mcols(GR))
       sum(is.element(vn, ns))
   }))
   if (any(nams != 3)) {
       warning("At least one element from 'LR' has incorrect column names")
       cat("\n")
       stop("LR is not a valid 'pDMP' object")
   } else invisible(TRUE)
}

#' @rdname validateClass
#' @keywords internal
#' @exportMethod validateClass.InfDiv
validateClass.InfDiv <- function(LR) {
   vn <- c("hdiv", "TV")
   if (any(!unlist(lapply(LR, function(GR) class(GR) == "GRanges")))) {
       warning("At least one element from 'LR' is not a 'GRanges' object")
       cat("\n")
       stop("LR is not a valid 'InfDiv' object")
   }
   nams <- unlist(lapply(LR, function(GR) {
       ns <- colnames(mcols(GR))
       sum(is.element(vn, ns))
   }))
   if (any(nams != 2)) {
       warning("At least one element from 'LR' has incorrect column names")
       cat("\n")
       stop("LR is not a valid 'InfDiv' object")
   } else invisible(TRUE)
}


