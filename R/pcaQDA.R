#' @rdname pcaQDA
#' @name pcaQDA
#' @title Quadratic Discriminant Analysis (QDA) using Principal Component
#'     Analysis (PCA)
#' @description The principal components (PCs) for predictor variables provided
#'     as input data are estimated and then the individual coordinates in the
#'     selected PCs are used as predictors in the qda
#' @details The principal components (PCs) are obtained using the function
#'     'prcomp' from R pacakage 'stats', while the qda is performed using the
#'     'qda' function from R package 'MASS'. The current application only uses
#'     basic functionalities of mentioned functions. As shown in the example,
#'     'pcaQDA' function can be used in general classification problems.
#'
#' @param formula Same as in 'qda'from pakage 'MASS'.
#' @param data Same as in 'qda'from pakage 'MASS'.
#' @param grouping Same as in 'qda'from pakage 'MASS'.
#' @param n.pc Number of principal components to use in the qda.
#' @param scale Same as in 'prcomp' from pakage 'prcomp'.
#' @param center Same as in 'prcomp' from pakage 'prcomp'.
#' @param tol Same as in 'prcomp' from pakage 'prcomp'v.
#' @param method Same as in 'qda'from pakage 'MASS'.
#' @param max.pc Same as in paramter 'rank.' from prcomp' {prcomp}.
#' @return
#'     Function 'pcaQDA' returns an object ("pcaQDA") consisting of a list with
#'     two objects:
#'         1) 'qda': an object of class 'qda' from package 'MASS'.
#'         2) 'pca': an object of class 'prcomp' from package 'stats'.
#'     For information on how to use these objects see ?qda and ?prcomp.
#'
#' @examples
#' data(iris)
#' qd1 <- pcaQDA(formula = Species ~ Petal.Length + Sepal.Length + Sepal.Width,
#'             data = iris, n.pc = 1, max.pc = 2, scale = TRUE, center = TRUE)
#' ## === Prediction === ##
#' qd2 <- pcaQDA(formula = Species ~., data = iris, n.pc = 1, max.pc = 2,
#'              scale = TRUE, center = TRUE)
#' set.seed(123)
#' idx <- sample.int(150, 40)
#' newdata <- iris[idx, 1:4]
#' newdata.prediction <- predict(qd2, newdata = newdata, type = "all")
#'
#' ## The confusion matrix
#' x <- data.frame(TRUE.class = iris$Species[idx],
#'                PRED.class = newdata.prediction$class)
#' table(x)
#' @importFrom MASS qda
#' @importFrom stats prcomp terms
#' @export
pcaQDA <- function(formula=NULL, data=NULL, grouping=NULL, n.pc=1, scale=FALSE,
                   center=FALSE, tol=1.0e-4, method="moment", max.pc=NULL) {

   Check <- ArgumentCheck::newArgCheck()
   if (!is.null(formula) && class(formula) != "formula") {
       ans <- paste("A formula of the form groups ~ x1 + x2 + ...",
                       "(see ?pcaQDA or ?qda).")
       ArgumentCheck::addError(msg=ans, argcheck=Check)
   }
   if (!is.null(formula) && class(formula) == "formula") {
       vn <- try(attr(terms(formula), "term.labels"), silent=TRUE)
       if (inherits(vn, "try-error")) {
          vn <-  try(setdiff(colnames(data), as.character(formula)[2]),
               silent=TRUE)
       }
       if (inherits(vn, "try-error")) stop("* Error in the formula")
       if (length(vn) < n.pc) {
           ans <- "The number of number predictor variables greater than "
           ans1 <- "the number of principal components: "
           ArgumentCheck::addError(msg=paste0(ans, ans1, n.pc), argcheck=Check)
       }
   }
   if (is.null(formula) && is.null(grouping)) {
       ans <- "A formula or grouping varible must be provided."
       ArgumentCheck::addError(msg=ans, argcheck=Check)
   }
   if (is.null(formula) && !is.null(grouping)) {
       vn <- setdiff(colnames(data), as.character(grouping))
       if (length(vn) < n.pc) {
           ans <- "The number of number predictor variables must be greater than "
           ans1 <- "or equal the number of principal components: "
           ArgumentCheck::addError(msg=paste0(ans, ans1, n.pc, "\n"),
                                   argcheck=Check)
       }
   }
   ArgumentCheck::finishArgCheck(Check)
   if (is.null(formula) && !is.null(grouping)) {
       m <- dim(data)
       if (floor(m[1] / 3) < n.pc) {
           ans <- "The number principal components: "
           ans1 <- " must be lower than the number of individuals N/3 \n"
           warning(paste0(ans, n.pc, ans1))
       }
   }
   if (!is.null(formula)) {
       resp <- as.character(formula)[2]
       pc <- prcomp(x=data[vn], retx=TRUE, center=center, scale.=scale, tol=tol,
                   rank.=max.pc)
       cn <- colnames(pc$x)
       if (ncol(pc$x) > n.pc) {
           ind.coord <- as.data.frame(pc$x[, 1:n.pc])
           colnames(ind.coord) <- cn[1:n.pc]
       } else {
           ind.coord <- as.data.frame(pc$x)
           colnames(ind.coord) <- cn[1:ncol(pc$x)]
       }

       qda.model <- qda(x=ind.coord, grouping=data[resp][,1], tol=tol,
                   method=method)
   } else {
       pc <- prcomp(x=data[vn], retx=TRUE, center=center, scale.=scale,
                tol=tol, rank.=max.pc)
       cn <- colnames(pc$x)
       ind.coord <- as.data.frame(pc$x[, 1:n.pc])
       colnames(ind.coord) <- cn[1:n.pc]
       qda.model <- qda(x=ind.coord, grouping=grouping, tol=tol, method=method)
   }
   model <- structure(list(qda=qda.model, pca=pc), class="pcaQDA")
   return(model)
}

#' @name predict.pcaQDA
#' @rdname pcaQDA
#' @aliases predict.pcaQDA
#' @title Linear Discriminant Analysis (qda) using Principal Component Analysis
#' @description Predict using a PCA-LDA model built with function 'pcaLDA'
#' @details NULL
#' @param object To use with function 'predict'. A 'pcaQDA' object containing a
#'     list of two objects: 1) an object of class inheriting from "qda" and 2)
#'     an object of class inheriting from "prcomp".
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction.
#' @param type To use with function 'predict'. The type of prediction required.
#'     The default is "all" given by function 'predict.qda' from MASS package:
#'     'class', 'posterior', and 'scores' (see ?predict.QDA).
#' @param ... Not in use.
#' @export
predict.pcaQDA <- function(object, newdata,
                           type=c("qda.pred", "class", "posterior",
                               "pca.ind.coord", "all"), ...) {
   if (!inherits(object, "pcaQDA")) {
       stop("* 'object' must be a model from class 'pcaQDA'")
   }
   vn <- rownames(object$pca$rotation)

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
       newnam <- colnames(mcols(newdata))
       newdata$logP <- log10(newdata$wprob + 2.2e-308)
       newdata <- mcols(newdata)
       newdata <- newdata[vn]
       newdata <- as.matrix(newdata)
   } else newdata <- newdata[vn]

   ## Centering and scaling new individuals
   dt.scaled <- scale(newdata, center=object$pca$center, scale=object$pca$scale)
   ## Coordinates of the individuals
   coord_func <- function(ind, loadings){
       x <- loadings * ind
       return(apply(x, 2, sum))
   }
   nc <- ncol(object$qda$means)
   loadings <- object$pca$rotation[, 1:nc]
   if (nc == 1 ) loadings <- as.matrix(loadings)

   ind.coord <- t(apply(dt.scaled, 1, coord_func, loadings))
   if (nc == 1 ) ind.coord <- t(ind.coord)
   pred <- predict(object$qda, newdata=ind.coord, prior=object$qda$prior)
   pred <- switch(type[1], qda.pred=pred, class=pred$class,
               posterior=pred$posterior, pca.ind.coord=ind.coord,
               all=list(class=pred$class, posterior=pred$posterior,
                           pca.ind.coord=ind.coord))
  return(pred)
}

