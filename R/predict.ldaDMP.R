#' @rdname predict.ldaDMP
#' @title Predict function for 'ldaDMP' method
#' @description Predict using a LDA model built with function
#'     \code{\link[MethylIT]{evaluateDIMPclass}}
#' @details NULL
#' @param object To use with function 'predict'. An object from class 'ldaDMP'
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction
#' @param type To use with function 'predict'. The type of prediction
#'     required. The default is "all" given by function 'predict.lda' from MASS
#'     package: 'class', 'posterior', and 'scores' (cases scores on discriminant
#'     variables, see \code{link[MASS]{lda}}).
#' @param ... arguments passed to or from other methods.
#' @seealso \code{link[MethylIT]{estimateCutPoint}}, \code{link[MASS]{lda}}
#' @keywords internal
#' @exportMethod predict.ldaDMP
predict.ldaDMP <- function(object, newdata,
                           type = c("class", "posterior", "scores"), ...) {
   if (!inherits(object, "ldaDMP")) {
       stop("* 'object' must be a model from class 'ldaDMP'")
   }

   vn <- colnames(object$means)

   if (!is.null(newdata) && inherits(newdata, c("pDMP", "InfDiv", "GRanges"))) {
     if (inherits(newdata, c("pDMP", "InfDiv"))) newdata <- unlist(newdata)
     if (is.element("pos", vn)) {
       position <- function(gr) {
         chrs <- split(gr, seqnames(gr))
         gr <- lapply(chrs, function(grc) {
           x <- start(grc)
           x.min <- min(x)
           x.max <- max(x)
           delta <-  max(c(x.max - x, 1))
           return((x - x.min) / (delta))})
         return(unlist(gr))
       }
       newdata$pos <- position(newdata)
     }
     newdata$logP <- log10(newdata$wprob + 2.2e-308)
     newdata <- mcols(newdata)
     newdata <- newdata[vn]
   } else newdata <- newdata[vn]
   class(object) <- "lda"

   pred <- predict(object, newdata = newdata, prior= object$prior)
   pred <- switch(type[1], lda.pred=pred, class=pred$class,
               posterior=pred$posterior,
               scores=pred$x ## cases scores on discriminant variables
               )
  return(pred)
}
