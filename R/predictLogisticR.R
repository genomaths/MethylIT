#' @rdname predict.LogisticR
#' @name predict.LogisticR
#' @aliases predict.LogisticR
#' @title Predict function for logistic regression model from 'LogisticR' class
#' @description Predict using a logistic model obtained from the output of
#'     function \code{\link{evaluateDIMPclass}}.
#' @details This function is specific for predictions based on a logistic model
#'     given by function \code{\link{evaluateDIMPclass}}. A logistic model is
#'     obtained with 'glm' regression can be used directly with function
#'     'predict' from 'stats' package.
#' @param object To use with function 'predict'. An object from 'LogisticR'
#'     class. A logistic model given by function
#'     \code{\link{evaluateDIMPclass}}.
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction. Optionally, an object from class "GRanges", a list of
#'     GRanges, "pDMP" or "InfDIv", in which to look for variables with which to
#'     predict. If omitted, the fitted linear predictors are used.
#' @param type The type of output required. Possible outputs are: 'class',
#'     "posterior" and "all". The default is "all".
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param ... Not in use.
#'
#' @return If type is set to 'all', then the original 'newdata' with two columns
#'     added, predicted classes and 'posterior' probabilities, in the
#'     meta-columns of each GRange object are given. If 'newdata' is null, then
#'     the predictions given for the model by function
#'     \code{\link[stats]{predict.glm}} are returned. if type is set to 'class'
#'     or to "posterior", then the unlisted predited classification or
#'     posterior classification probabilities are returned.
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @exportMethod predict.LogisticR
predict.LogisticR <- function(object, newdata = NULL,
                              type=c("all", "class", "posterior"),
                              num.cores = 1L, tasks = 0L, ...) {
   type <- type[1]
   if (!is.element(type, c("class", "all", "posterior")))
       stop("The type setting '", type, "' does not exist")
   if (!inherits(object, "LogisticR")) {
       stop("* 'object' must be logistic a model from class 'LogisticR'")
   }

   # ---This builds data frames from newdata from class 'pDMP' or 'InfDiv' ----#
   if (!is.null(newdata) &&
           (inherits(newdata, 'pDMP') || inherits(newdata, 'InfDiv'))) {
       if (!validateClass(newdata)) {
           stop("newdata is not a valid object from class '",
               class(newdata),"'" )
       }
       position <- function(gr) {
           chrs <- split(gr, seqnames(gr))
           gr <- lapply(chrs, function(grc) {
               x <- start(grc)
               x.min <- min(x)
               x.max <- max(x)
               if (x.min == Inf) x.min = 0
               if (x.max == -Inf) x.max = 1
               delta <-  max(c(x.max - x, 1))
               return((x - x.min) / (delta))})
           return(unlist(gr))
       }
       v <- c("hdiv", "TV", "bay.TV", "logP", "pos")
       vn <- setdiff(names(coef(object$modeling)),"(Intercept)")
       v <- v[na.omit(match(vn, v))]
       inter <- unlist(lapply(grep("[:]", vn),
                           function(k) strsplit(vn[k], split = ":")[[1]]))
       vn <- union(v, inter)

       if (is.list(newdata)) {
           dt <- lapply(newdata, function(x) {
                           df <- data.frame(hdiv=x$hdiv, TV=x$TV,
                                            bay.TV=x$bay.TV,
                                           logP=log10(x$wprob + 2.2e-308))
                           if (is.element("pos", vn)) df$pos = position(x)
                           df <- as.matrix(df[, vn])
                           df <- scale(df, center = object$center,
                                       scale = object$scale)
                           return(data.frame(df))
           })
       }
       else {
           dt <- cbind(hdiv=newdata$hdiv, TV=newdata$TV,
                       logP=log10(newdata$wprob + 2.2e-308))
           if (is.element("pos", vn)) dt$pos = position(newdata)
           dt <- as.matrix(dt[, vn])
           dt <- data.frame(scale(df, center = object$center,
                                  scale = object$scale))
       }
   }
   # ---------------------------------------------------------------------#
   object$modeling <- structure(object$modeling, class=c("glm", "lm"))
   if (is.null(newdata)) dt <- NULL

   if (!is.list(dt) && !is.null(newdata)) {
       pred <- predict.glm(object$modeling, newdata = dt, type="response")
       newdata$class <- rep( "CT", length(pred))
       newdata$class[pred > 0.5] <- "TT"
       newdata$posterior <- pred
   }

   if (num.cores > 1 && is.list(dt)) {
       nms <- names(newdata)
       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")
       newdata <- bplapply(1:length(dt),
                       function(k) {
                           p <- predict.glm(object$modeling, newdata = dt[[k]],
                                           type="response")
                           PredClass <- rep( "CT", length(p))
                           PredClass[p > 0.5] <- "TT"
                           newdata[[k]]$class <- PredClass
                           newdata[[k]]$posterior <- p
                           return(newdata[[k]])
                       },
               BPPARAM=bpparam)
       names(newdata) <- nms
   }
   else {
       # === To not depend on BiocParallel package
       if (num.cores == 1 && is.list(dt)) {
           nms <- names(newdata)
           newdata <- lapply(1:length(dt),
                               function(k) {
                                   p <- predict.glm(object$modeling,
                                                   newdata = dt[[k]],
                                                   type="response")
                                   PredClass <- rep( "CT", length(p))
                                   PredClass[p > 0.5] <- "TT"
                                   newdata[[k]]$class <- PredClass
                                   newdata[[k]]$posterior <- p
                                   return(newdata[[k]])
                               }
                       )
           names(newdata) <- nms
       }
   }

   if (type == "class" && !is.list(newdata) && !is.null(newdata)) {
       newdata <- newdata$class
   }

   if (type == "posterior" && !is.list(newdata) && !is.null(newdata)) {
       newdata <- newdata$posterior
   }

   if (type == "class" && inherits(newdata, "list")) {
       newdata <- unlist(lapply(newdata, function(x) x$class))
   }

   if (type == "posterior" && inherits(newdata, "list")) {
       newdata <- unlist(lapply(newdata, function(x) x$posterior))
   }

   if (is.null(newdata)) {
       pred <- predict.glm(object, newdata = NULL,
                             type="response")
       if (type == "class") {
           newdata <- rep( "CT", length(pred))
           newdata[pred > 0.5] <- "TT"
       } else newdata <- pred
   }

   return(newdata)
}
