#' @rdname predict.qdaDMP
#' @name predict.qdaDMP
#' @aliases predict.pcaqda
#' @title Predict function for 'qdaDMP' method
#' @description Predict using a qda model built with function
#'     \code{\link[MethylIT]{evaluateDIMPclass}}
#' @param object To use with function 'predict'. An object from class 'qdaDMP'
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction
#' @param type To use with function 'predict'. The type of prediction
#'     required. The default is 'all' basic predictions: classes and posterior
#'     classification probabilities. Option 'qda.pred' returns the
#'     object given by function 'predict.qda' from MASS package: 'class',
#'     'posterior', 'scores' (cases scores on discriminant variables,
#'     see \code{\link[MASS]{qda}}.
#' @param ... arguments passed to or from other methods.
#' @seealso \code{\link[MethylIT]{estimateCutPoint}}, \code{\link[MASS]{qda}}
#' @keywords internal
#' @export
predict.qdaDMP <- function(object, ...) UseMethod("predict")

predict.qdaDMP <- function(
                        object,
                        newdata,
                        type = c("all", "class", "posterior",
                                 "scores", "qda.pred"),
                        ...) {
    if (!inherits(object, "qdaDMP")) {
        stop("* 'object' must be a model from class 'qdaDMP'")
    }

    type <- match.arg(type)
    vn <- colnames(object$means)

    if (!is.null(newdata) && inherits(newdata, c("pDMP",
        "InfDiv", "GRanges"))) {
        if (inherits(newdata, c("pDMP", "InfDiv")))
            newdata <- unlist(newdata)
        if (is.element("pos", vn)) {
            position <- function(gr) {
                chrs <- split(gr, seqnames(gr))
                gr <- lapply(chrs, function(grc) {
                                        x <- start(grc)
                                        x.min <- min(x)
                                        x.max <- max(x)
                                        delta <- max(c(x.max - x, 1))
                                        return((x - x.min)/(delta))
                })
                return(unlist(gr))
            }
            newdata$pos <- position(newdata)
        }
        if (is.element("logP", vn))
            newdata$logP <- log10(newdata$wprob + 2.2e-308)
        newdata <- mcols(newdata)
        newdata <- newdata[vn]
    } else newdata <- newdata[vn]
    class(object) <- "qda"

    pred <- predict(object, newdata = newdata, prior = object$prior)
    pred <- switch(type,
                qda.pred = pred,
                class = pred$class,
                posterior = pred$posterior,
                scores = pred$x, ## cases scores on discriminant variables
                all = list(class = pred$class, posterior = pred$posterior)
    )
    return(pred)
}
