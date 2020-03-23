#' @rdname validateClass
#' @title Function to validate S3 classes in MethylIT
#' @description S3 generic used internally by MethylIT to ensure correct object
#'     classes are being used.  validateClass() is called by these functions:
#'     estimateCutPoint evaluateDIMPclass FisherTest getPotentialDIMP
#'     nonlinearFitDist predictDIMPclass predictLogisticR selectDIMP
#'     uniqueGRfilterByCov
#'
#' @param LR An object from class 'pDMP' or 'InfDiv'
#' @return TRUE if object is valid, otherwise warning() and stop() are called
#' @keywords internal
#' @examples
#' data(HD, gof)
#' PS <- getPotentialDIMP(LR = HD, nlms = gof$nlms, div.col = 9L, alpha = 0.05)
#' validateClass(HD)
#' validateClass(PS)
#'
#' @export
validateClass <- function(LR) UseMethod("validateClass")

#' @rdname validateClass
#' @keywords internal
#' @export
validateClass.default <- function(LR) {
    stop(paste("'validateClass' does not know how to handle object of class", 
        class(LR), "and can only be used on classes 'pDMP', 'InfDiv'"))
    
}

#' @rdname validateClass
#' @keywords internal
#' @importFrom S4Vectors mcols
#' @export
validateClass.pDMP <- function(LR) {
    vn <- c("hdiv", "TV", "wprob")
    if (any(!unlist(lapply(LR, function(GR) is(GR, 
        "GRanges"))))) {
        warning("At least one element from 'LR' is not a 'GRanges' object.")
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
        cat("Columns named 'hdiv', 'TV', and 'wprob' must be present in ", 
            "each GRanges metacolumn.")
        cat("\n")
        stop("LR is not a valid 'pDMP' object")
    } else invisible(TRUE)
}

#' @rdname validateClass
#' @keywords internal
#' @importFrom S4Vectors mcols
#' @export
validateClass.InfDiv <- function(LR) {
    vn <- c("hdiv", "TV")
    if (any(!unlist(lapply(LR, function(GR) is(GR, 
        "GRanges"))))) {
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
        cat("Columns named 'hdiv' and 'TV'must be present in ", 
            "each GRanges metacolumn.")
        cat("\n")
        stop("LR is not a valid 'InfDiv' object")
    } else invisible(TRUE)
}


