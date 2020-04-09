#' @rdname simpleCutPoint
#' @title Simple Cutpoint estimation
#' @description Internal function to estimate simple cutpoint according to
#'  Youden Index.
#' @details This function is called by function
#'  \code{\link{estimateCutPoint}}.
#' @param LR,res,control.names,treatment.names,column,div.col Same as in
#'  \code{\link{estimateCutPoint}}
#' @param column,div.col,tv.col,tv.cut,clas.perf Same as in
#'  \code{\link{estimateCutPoint}}
#' @param classifier,n.pc,prop,cutp_data,num.cores,tasks Same as in
#'  \code{\link{estimateCutPoint}}
#'
#' @return Specified in function \code{\link{estimateCutPoint}} for parameter
#' setting \emph{simple = TRUE}
#' @importFrom caret confusionMatrix
#' @importFrom stats addmargins
#' @keywords internal
#' @export
#' @examples
#' ## Get a set of potential DMPS (PS)
#' data(PS, package = 'MethylIT')
#'
#' cutp <- simpleCutPoint(LR = PS, column = c(hdiv = TRUE, TV = TRUE,
#'                                            wprob = TRUE, pos = TRUE),
#'                        classifier = 'qda', n.pc = 4,
#'                        control.names = c('C1', 'C2', 'C3'),
#'                        treatment.names = c('T1', 'T2', 'T3'),
#'                        tv.cut = 0.5, clas.perf = TRUE, prop = 0.6,
#'                        div.col = 9L)
#'
#' cutp
#'
simpleCutPoint <- function(LR, control.names, treatment.names,
                        column, div.col, tv.col = NULL, tv.cut, clas.perf,
                        classifier, n.pc, prop, cutp_data = FALSE, num.cores,
                        tasks, ...) {

    LR = list(unlist(LR[control.names]), unlist(LR[treatment.names]))
    names(LR) <- c("ctrl", "treat")
    LR <- structure(LR, class = "pDMP")
    if (is.null(tv.col)) tv.col = match("TV", colnames(mcols(LR[[1]])))
    classes <- c(rep("CT", length(LR$ctrl)), rep("TT", length(LR$treat)))
    classes <- factor(classes, levels = c("CT", "TT"))

    dt <- .infDiv(LR = LR, div.col = div.col)
    cutp <- .roc(dt = dt)
    cutpoint <- cutp$cutp
    predClasses <- mcols(unlist(LR)[, div.col])[, 1] > cutpoint
    predClasses[predClasses == TRUE] <- "TT"
    predClasses[predClasses == FALSE] <- "CT"
    predClasses <- factor(predClasses, levels = c("CT", "TT"))
    cf.mat <- try(confusionMatrix(data = predClasses,
                                reference = classes, positive = "TT"),
                silent = TRUE)

    if (inherits(cf.mat, "try-error"))
        stop("*** Cutpoint estimation based on Youden index failed. \n",
            "Only ", 100 * sum(dt$treat$idiv > cutpoint)/
        length(dt$treat$idiv),
            " percent of the cytosine positions from treatment are \n",
            "greater than the cutpoint value: ",
            round(cutpoint, 2), "\n\n Other indicators are: \n",
            "sensitivity = ", round(cutp$sens, 3), "\n", "specificity = ",
            round(cutp$spec, 3), "\n", "auc = ", round(cutp$auc, 4),
            "\n\n", " Please try the optimal cutpoint estimation",
            " with machine-learning \n classifiers (simple = FALSE, ",
            "see examples given in \n the help for function:",
            " ?estimateCutPoint)")
    res <- cutpoint_object()
    res$cutpoint <- cutpoint
    if (cutp_data) res$cutpData <- dt

    if (clas.perf) {
        dmps <- selectDIMP(LR, div.col = div.col,
                        cutpoint = cutpoint, tv.col = tv.col,
                        tv.cut = tv.cut)

        conf.mat <- evaluateDIMPclass(LR = dmps,
                                column = column, control.names = "ctrl",
                                treatment.names = "treat",
                                classifier = classifier[1],
                                prop = prop, n.pc = n.pc,
                                output = "conf.mat",
                                num.cores = num.cores, tasks = tasks, ...)

        predClasses <- predict(object = conf.mat$model,
                            newdata = dmps, type = "class")
        predClasses <- factor(predClasses, levels = c("CT", "TT"))
        classes <- c(rep("CT", length(dmps$ctrl)),
                    rep("TT", length(dmps$treat)))
        classes <- factor(classes, levels = c("CT", "TT"))
        conf.matrix <- confusionMatrix(data = predClasses,
                                    reference = classes, positive = "TT")

        res$testSetPerformance <- conf.mat$Performance
        res$testSetModel.FDR <- conf.mat$FDR
        res$model <- conf.mat$model
        res$modelConfMatrix <- conf.matrix
        res$initModel <- "Youden Index"
        res$initModelConfMatrix <- cf.mat
        res$classifier <- classifier[1]
    } else {
        res$modelConfMatrix <- cf.mat
        res$initModel <- "Youden Index"
    }
    return(res)
}

### ================== Auxiliary functions ========================== #

.infDiv <- function(LR, div.col = NULL) {
    ## This builds data frames from the list or ranges
    ## LR to be used for ROC analysis LR: list of sample
    ## GRanges
    if (is.null(div.col))
    stop(paste("* Provide a divergence column"))

    dt <- list(ctrl = data.frame(), treat = data.frame())
    dt$ctrl <- data.frame(idiv = abs(mcols(LR$ctrl[, div.col])[, 1]),
                        treat = "ctrl", stringsAsFactors = TRUE)
    dt$treat <- data.frame(idiv = abs(mcols(LR$treat[, div.col])[, 1]),
                            treat = "treat", stringsAsFactors = TRUE)
    return(dt)
}

### Compute ROC data
.roc <- function(dt) {
    ## This function build the ROC and estimate the
    ## Hellinger divergence cutoff point starting from
    ## which TRUE DMPs are found. dt0 & dt1: data
    ## frames built with HD function dt0: control dt1:
    ## treatment folder: place where the ROC will be
    ## printed
    dt <- do.call(rbind, dt)
    l <- levels(dt$treat)

    dt$status <- as.character(dt$treat)
    dt$status[dt$status == l[1]] <- 0
    dt$status[dt$status == l[2]] <- 1
    dt$status <- as.numeric(dt$status)
    rownames(dt) <- NULL

    m <- as.matrix(base::table(dt$idiv, dt$status))
    m <- addmargins(rbind(0, m), 2)
    nr <- nrow(m)
    m <- apply(m, 2, cumsum)
    sens <- (m[nr, 2] - m[, 2])/m[nr, 2]
    spec <- m[, 1]/m[nr, 1]
    auc <- sum((sens[-1] + sens[-nr])/2 * abs(diff(1 - spec)))
    # Youden Index
    idx <- which.max(sens[-1] + spec[-nr])
    return(list(cutp = as.numeric(names(idx)),
                auc = auc, sens = sens[idx], spec = spec[idx]))
}

## It just build an empty object
cutpoint_object <- function() {
    res <- list(cutpoint = NA, testSetPerformance = NA,
                testSetModel.FDR = NA, model = NA, modelConfMatrix = NA,
                initModel = NA, postProbCut = NA, postCut = NA,
                classifier = NA, statistic = NA, optStatVal = NA,
                cutpData = NA)
    invisible(structure(res, class = c("CutPoint", "list")))
}
