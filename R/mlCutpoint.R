#' @rdname mlCutpoint
#' @title Machine-learning Approach for Cutpoint Estimation
#' @description Internal function to estimate cutpoint following a
#'  machine-learning approach
#' @details This function is called by function \code{\link{estimateCutPoint}}.
#' @param LR,res,control.names,treatment.names,column,div.col Same as in
#'  \code{\link{estimateCutPoint}}
#' @param column,div.col,tv.col,tv.cut,cut.values,stat Same as in
#'  \code{\link{estimateCutPoint}}
#' @param classifier1,classifier2,n.pc,prop,post.cut Same as in
#'  \code{\link{estimateCutPoint}}
#' @param cutp_data,num.cores,tasks Same as in \code{\link{estimateCutPoint}}
#'
#' @return Specified in function \code{\link{estimateCutPoint}} for parameter
#'  setting \emph{simple = FALSE}
#' @importFrom caret confusionMatrix
#' @keywords internal
#' @export
#' @examples
#' ## Get a set of potential DMPS (PS)
#' data(PS, package = 'MethylIT')
#'
#' cutp <- mlCutpoint(LR = PS,
#'                  column = c(hdiv = TRUE, TV = TRUE,
#'                             wprob = TRUE, pos = TRUE),
#'                  classifier1 = 'qda', n.pc = 4,
#'                  control.names = c('C1', 'C2', 'C3'),
#'                  treatment.names = c('T1', 'T2', 'T3'),
#'                  tv.cut = 0.68, prop = 0.6,
#'                  div.col = 9L)
#' cutp

mlCutpoint <- function(LR, control.names, treatment.names,
                    column, div.col, tv.col = NULL, tv.cut, post.cut = 0.5,
                    classifier1, classifier2 = NULL, n.pc, prop = 0.6,
                    stat = 0L,  cut.values = NULL,
                    num.cores = 1L, tasks = 0L, ...) {

    LR = list(unlist(LR[control.names]), unlist(LR[treatment.names]))
    names(LR) <- c("ctrl", "treat")
    LR <- structure(LR, class = "pDMP")
    if (is.null(tv.col)) tv.col = match("TV", colnames(mcols(LR[[1]])))
    classes <- c(rep("CT", length(LR$ctrl)), rep("TT", length(LR$treat)))
    classes <- factor(classes, levels = c("CT", "TT"))

    # ------------- To search for a cutpoint --------------- #

    if (is.null(classifier2[1])) classifier2 <- classifier1[1]
    conf.mat <- evaluateDIMPclass(LR, column = column,
                                control.names = "ctrl",
                                treatment.names = "treat",
                                classifier = classifier1,
                                prop = prop,  output = "conf.mat",
                                n.pc = n.pc, interaction = interaction,
                                num.cores = num.cores, tasks = tasks,
                                ...)
    post <- predict(object = conf.mat$model, newdata = LR, type = "posterior")

    if (classifier1[1] == "logistic" || classifier1[1] == "pca.logistic")
        idx <- which(post > post.cut)
    else idx <- which(post[, 2] > post.cut)
    cutp <- min(mcols(unlist(LR)[idx, div.col])[, 1], na.rm = TRUE)

    res <- cutpoint_object()
    res$postCut <- cutp
    cuts <- cutp
    if (stat %in% seq_len(11)) st <- conf.mat$Performance$byClass[stat]
    if (stat == 0) st <- conf.mat$Performance$overall[1]
    if (stat == 12) st <- conf.mat$FDR

    if (!is.null(cut.values)) {
        if (is.numeric(cut.values)) {
            cuts <- sort(c(cut.values, cutp))
            if (max(cuts) > max(mcols(unlist(LR)[, div.col])[, 1])) {
                warning("Cut values supplied are"," > max Inf Div. Ignored")
                cuts <- cutp
            }
        } else warning("Parameter cut.values must be numeric. Ignored")
    }
    if (length(cuts) > 1) {
        cut_search <- cutpoint_search(LR = LR, column = column,
                                div.col = div.col, cuts = cuts,
                                tv.col = tv.col, tv.cut = tv.cut,
                                classifier = classifier2[1], prop = prop,
                                n.pc = n.pc, st = st, num.cores = num.cores,
                                tasks = tasks, stat = stat, ...)

        dmps <- cut_search$dmps
        cutp <- cut_search$cut
        conf.mat <- cut_search$cfm
        st <- cut_search$st
    } else {
        dmps <- selectDIMP(LR, div.col = div.col,
                            cutpoint = cutp, tv.col = tv.col,
                            tv.cut = tv.cut)
        conf.mat <- evaluateDIMPclass(LR = dmps,
                                    column = column, control.names = "ctrl",
                                    treatment.names = "treat",
                                    classifier = classifier2[1],
                                    prop = prop, output = "conf.mat",
                                    n.pc = n.pc, interaction = interaction,
                                    num.cores = num.cores, tasks = tasks, ...)
        if (stat == 0)  st <- conf.mat$Performance$overall[1]
        if (stat %in% seq_len(11)) st <- conf.mat$Performance$byClass[stat]
    }
    predClasses <- predict(object = conf.mat$model, newdata = dmps,
                            type = "class")
    predClasses <- factor(predClasses, levels = c("CT", "TT"))
    classes <- c(rep("CT", length(dmps$ctrl)),
                rep("TT", length(dmps$treat)))
    classes <- factor(classes, levels = c("CT", "TT"))
    conf.matrix <- confusionMatrix(data = predClasses,
                                    reference = classes, positive = "TT")

    STAT <- c("Accuracy", "Sensitivity", "Specificity",
            "Pos Pred Value", "Neg Pred Value",
            "Precision", "Recall", "F1", "Prevalence",
            "Detection Rate", "Detection Prevalence",
            "Balanced Accuracy", "FDR")
    res$cutpoint <- cutp
    res$testSetPerformance <- conf.mat$Performance
    res$testSetModel.FDR <- conf.mat$FDR
    res$model <- conf.mat$model
    res$modelConfMatrix <- conf.matrix
    res$initModel <- classifier1[1]
    res$classifier <- classifier2[1]
    res$statistic <- STAT[stat + 1]
    res$postProbCut <- post.cut
    res$optStatVal <- st
    return(res)
}

# ========================= Auxiliary function ======================== ###

cutpoint_search <- function(LR, column, div.col, cuts, tv.col, tv.cut,
                            classifier, prop, n.pc, st, num.cores, tasks,
                            stat, ...) {
    k = 1
    opt <- FALSE
    overcut <- FALSE
    while (k < length(cuts) && !opt && !overcut) {
        dmps <- selectDIMP(
                            LR,
                            div.col = div.col,
                            cutpoint = cuts[k],
                            tv.col = tv.col,
                            tv.cut = tv.cut
        )

        min.div <- min(mcols(dmps$treat[, div.col])[, 1], na.rm = TRUE)
        if (k == 1) cutp <- min.div
        if (length(dmps$ctrl) > 0) {
            conf.mat <- evaluateDIMPclass(
                                        LR = dmps,
                                        column = column,
                                        control.names = "ctrl",
                                        treatment.names = "treat",
                                        classifier = classifier[1],
                                        prop = prop,
                                        output = "conf.mat",
                                        n.pc = n.pc,
                                        num.cores = num.cores,
                                        tasks = tasks
            )
            if (stat == 0) {
                st0 <- conf.mat$Performance$overall[1]
                if (st < st0) {
                    st <- st0
                    cutp <- max(cuts[k], min.div, na.rm = TRUE)
                }
                if (st == 1) opt <- TRUE
                k <- k + 1
            }
            if (stat %in% seq_len(11)) {
                st0 <- conf.mat$Performance$byClass[stat]
                if (st < st0) {
                    st <- st0
                    cutp <- max(cuts[k], min.div, na.rm = TRUE)
                }
                if (st == 1) opt <- TRUE
                k <- k + 1
            }
            if (stat == 12) {
                st0 <- conf.mat$FDR
                if (st > st0) {
                    st <- st0
                    cutp <- max(cuts[k], min.div, na.rm = TRUE)
                }
                if (st == 0) opt <- TRUE
                k <- k + 1
            }
        } else {
            overcut <- TRUE
            if (k == 1) {
                st <- conf.mat$Performance$byClass[stat]
                warning("Model classifier ", classifier[1], " is enough")
            }
        }
    }
    return(list(dmps = dmps, cut = cutp, cfm = conf.mat, st = st))
}

