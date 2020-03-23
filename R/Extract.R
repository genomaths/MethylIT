#' @rdname Extract
#' @name Extract
#' @title Subscript operator for objects from several Methyl-IT classes
#' @param x A list-like or vector-like object
#' @param i Index specifying elements to Extract or replace. Index can be a
#'     numeric or a character vector. See ?base::\code{\link[base]{Extract}} for
#'     a description of these arguments.
#' @description \eqn{x[i]} returns an object of the length as \eqn{i} preserving
#'     the class from the original object \eqn{x}. The following Methyl-IT
#'     classes are considered: 'InfDiv','pDMP', 'glmDataSet', and
#'     'RangedGlmDataSet'.
#' @return Same as in ?base::\code{\link[base]{Extract}}, but preserving
#'     the class from the original object \eqn{x}.
#' @seealso base::\code{\link[base]{Extract}}
#' @importFrom S4Vectors extractROWS DataFrame
#' @keywords internal
#' @export
"[.InfDiv" <- function(x, i) {
    cl <- class(x)
    x <- .subset(x, i)
    attr(x, "class") <- cl
    return(x)
}

#' @rdname Extract
#' @export
#' @keywords internal
"[.pDMP" <- function(x, i) {
    cl <- class(x)
    x <- .subset(x, i)
    attr(x, "class") <- cl
    return(x)
}

#' @rdname Extract
#' @export
#' @keywords internal
"[.glmDataSet" <- function(x, i, j, ...) {
    if (is.null(nrow(x$counts))) 
        X <- t(x$counts) else X <- x$counts
    if (missing(j)) 
        j <- seq_len(ncol(X))
    x$counts <- .subset(X, i, j)
    rn <- rownames(x$colData)
    if (isS4(x$colData)) {
        x$colData <- DataFrame(x$colData[j, ])
        rownames(x$colData) <- rn[j]
    } else {
        x$colData <- data.frame(x$colData[j, ])
        rownames(x$colData) <- rn[j]
    }
    if (!is.null(x$optionData)) {
        x$optionData <- DataFrame(x$optionData[i, ])
    }
    attr(x, "class") <- "glmDataSet"
    return(x)
}

#' @rdname Extract
#' @export
#' @keywords internal
"[.RangedGlmDataSet" <- function(x, i, j, ...) {
    if (is.null(nrow(x$counts))) 
        X <- t(x$counts) else X <- x$counts
    if (missing(j)) 
        j <- seq_len(ncol(X))
    x$GR <- x$GR[i]
    x$counts <- .subset(X, i, j)
    rn <- rownames(x$colData)
    if (isS4(x$colData)) {
        x$colData <- DataFrame(x$colData[j, ])
        rownames(x$colData) <- rn[j]
    } else {
        x$colData <- data.frame(x$colData[j, ])
        rownames(x$colData) <- rn[j]
    }
    if (!is.null(x$optionData)) {
        x$optionData <- DataFrame(x$optionData[i, ])
    }
    attr(x, "class") <- "RangedGlmDataSet"
    return(x)
}

