#' @name predictDIMPclass
#' @rdname predictDIMPclass
#' @title Predict DIMP class
#' @description This function classify each DMP as a control or a treatment
#'     DMP
#' @details Predictions only makes sense if the query DMPs belong to same
#'     methylation context and derive from an experiment accomplished under the
#'     same condition set for the DMPs used to build the model.
#' @param LR A list of GRanges objects obtained through the through MethylIT
#'     downstream analysis. Basically, this object is a list of GRanges
#'     containing only differentially methylated position (DMPs). The metacolumn
#'     of each GRanges must contain the column: Hellinger divergence 'hdiv',
#'     total variation 'TV', the probability of potential DMP 'wprob', which
#'     naturally are added in the downstream analysis of MethylIT.
#' @param model A classifier model obtained with the function
#'     'evaluateDIMPclass'.
#' @param conf.matrix Optional. Logic, whether a confusion matrix should be
#'     returned (default, FALSE, see below).
#' @param control.names Optional. Names/IDs of the control samples, which must
#'     be include in the variable LR (default, NULL).
#' @param treatment.names Optional. Names/IDs of the treatment samples, which
#'     must be include in the variable LR (default, NULL).
#' @return The same LR object with tow new columns named 'class' and 'posterior'
#'     added to each GRanges object from LR (default). Based on the model
#'     prediction each DMP is labeled as control 'CT' or as treatment 'TT' in
#'     column 'class'. Column 'posterior' provides, for each DMP, the posterior
#'     probability that the given DMP can be classified as induced by the
#'     'treatment' (a treatment DMP).
#'
#'     Control DMPs classified as 'treatment' are false positives. However, if
#'     the same cytosine position is classified as 'treatment DMP' in both
#'     groups, control and treatment, but with higher posterior probability
#'     in the treatment group, then this would indicate a reinforcement of the
#'     methylation status in such a position induced by the treatment.
#'
#'     If 'conf.matrix' is TRUE and the arguments control.names and
#'     treatment.names are provided, then the overall confusion matrix is
#'     returned.
#' @importFrom caret confusionMatrix
#' @importFrom S4Vectors mcols DataFrame
#' @examples
#'
#' data(cutpoint, PS, package = 'MethylIT')
#'
#' ## DMPs are selected using the cupoints
#' DMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint,
#' tv.cut = 0.92)
#'
#' ## Classification of DMPs into two clases: DMPS from control and DMPs from
#' ## treatment samples and evaluation of the classifier performance (for more
#' ## details see ?evaluateDIMPclass).
#' perf <- evaluateDIMPclass(LR = DMPs, column = c(hdiv = TRUE, TV = TRUE,
#' wprob = TRUE, pos = TRUE), classifier = 'lda', n.pc = 4L,
#' control.names =  c('C1', 'C2', 'C3'), treatment.names = c('T1', 'T2', 'T3'),
#' center = TRUE, scale = TRUE, prop = 0.6)
#'
#' ## Now predictions of DMP for control and treament can be obtained
#' pred = predictDIMPclass(LR = DMPs, model = perf$model, conf.matrix = TRUE,
#' control.names = c('C1', 'C2', 'C3'), treatment.names = c('T1', 'T2', 'T3'))
#' @export
predictDIMPclass <- function(LR, model, conf.matrix = FALSE,
    control.names = NULL, treatment.names = NULL) {
    if (conf.matrix && (is.null(control.names) || is.null(treatment.names))) {
        stop(paste0("* if conf.mat = TRUE, then the character vectors for ",
            "'control.names' and treatment.names must be provided"))
    }

    if (!is.element(class(model), c("LogisticR", "pcaLogisticR",
        "pcaQDA", "pcaLDA", "ldaDMP", "qdaDMP")))
        stop("The ''model' class must be: one the following: \n",
            "'LogisticR', 'pcaLogisticR', 'pcaQDA', 'pcaLDA', 'lda', 'qda'")
    if (inherits(LR, "GRangesList"))
        LR <- as(LR, "list")
    if (class(LR)[1] == "list")
        class(LR) <- "pDMP"
    r1 <- try(validateClass(LR), silent = TRUE)
    if (inherits(r1, "try-error"))
        stop("*** LR is not an object from 'pDMP' class or is not",
            " coercible to 'pDMP' class")

    if (conf.matrix && !is.null(control.names) && !is.null(treatment.names)) {
        sn = names(LR)
        idx.ct = match(control.names, sn)
        if (any(is.na(idx.ct)))
            stop("*** 'control.names' arguments do not match LR names")
        idx.tt = match(treatment.names, sn)
        if (any(is.na(idx.tt)))
            stop("*** 'treatment.names' arguments do not match LR names")
        CT = LR[idx.ct]
        CT = unlist(CT)
        TT = LR[idx.tt]
        TT = unlist(TT)
        classSet = list(CT = CT, TT = TT)
        class(classSet) <- "pDMP"
    } else classSet = LR

    if (inherits(model, "LogisticR")) {
        LR <- predict(object = model, newdata = classSet,
            type = "all")
    } else {
        classifier <- function(GR) {
            pred <- predict(object = model, newdata = GR,
                type = "all")
            GR$class <- pred$class
            if (inherits(model, "pcaLogisticR")) {
                GR$posterior <- pred$posterior
            } else GR$posterior <- pred$posterior[, 2]
            return(GR)
        }
        LR <- lapply(classSet, classifier)
    }

    if (!conf.matrix) {
        return(LR)
    } else {
        TRUE_class <- factor(c(rep("CT", length(CT)),
            rep("TT", length(TT))), levels = c("CT",
            "TT"))
        PRED_class <- factor(c(as.character(LR$CT$class),
            as.character(LR$TT$class)), levels = c("CT",
            "TT"))

        conf.mat <- confusionMatrix(data = PRED_class,
            reference = TRUE_class, positive = "TT")
        m <- conf.mat$table
        conf.mat$FDR <- m[2, 1]/sum(m[2, ])
        return(conf.mat = conf.mat)
    }
}

