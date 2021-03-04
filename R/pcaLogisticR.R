#' @rdname pcaLogisticR
#' @name pcaLogisticR
#' @title Logistic Classification Model using Principal Component
#'     Analysis (PCA)
#' @description Principal components (PCs) are estimated from the predictor
#'     variables provided as input data. Next, the individual coordinates in
#'     the selected PCs are used as predictors in the logistic regression.
#' @details The principal components (PCs) are obtained using the function
#'     \code{\link[stats]{prcomp}}, while the logistic regression is performed
#'     using function \code{\link[stats]{glm}}, both functions from R package
#'     'stats'. The current application only use basic functionalities from the
#'     mentioned functions. As shown in the example, 'pcaLogisticR' function
#'     can be used in general classification problems.
#'
#' @param formula Same as in 'glm' from package 'stats'. One term carrying
#' interaction between two variables can be introduced (with notation as
#' indicated in \code{\link[stats]{formula}} function).
#' @param data Same as in 'glm' from package 'stats'.
#' @param scale Same as in 'prcomp' from package 'prcomp'.
#' @param center Same as in 'prcomp' from package 'prcomp'.
#' @param tol Same as in 'prcomp' from package 'prcomp'.
#' @param n.pc Number of principal components to use in the logistic.
#' @param max.pc Same as in parameter 'rank.' from package 'prcomp'.
#' @return
#' Function 'pcaLogisticR' returns an object ('pcaLogisticR' class)
#' containing a list of two objects:
#' \enumerate{
#'      \item 'logistic': an object of class 'glm' from package 'stats'.
#'      \item 'pca': an object of class 'prcomp' from package 'stats'.
#'      \item reference.level: response level used as reference.
#'      \item positive.level: response level that corresponds to a 'positive'
#'            result. When type = 'response', the probability vector returned
#'            correspond to the probabilities of each individual to be a
#'            result, i.e., the probability to belong to the class of positive
#'            level.
#' }
#'
#' For information on how to use these objects see ?glm and ?prcomp.
#'
#' @examples
#' data(iris)
#' data <- iris[ iris$Species != 'virginica', ]
#' data$Species <- droplevels(data$Species)
#' formula <- Species ~ Petal.Length + Sepal.Length + Petal.Width
#' pca.logistic <- pcaLogisticR(formula = formula,
#'                             data = data, n.pc = 2, scale = TRUE,
#'                             center = TRUE, max.pc = 2)
#' set.seed(123)
#' newdata <- iris[sample.int(150, 40), 1:4]
#' newdata.prediction <- predict(pca.logistic, newdata, type = 'all')
#'
#' @importFrom stats prcomp glm
#'
#' @export
pcaLogisticR <- function(formula = NULL,
                        data = NULL,
                        n.pc = 1,
                        scale = FALSE,
                        center = FALSE,
                        tol = 1e-04,
                        max.pc = NULL) {

    Check <- ArgumentCheck::newArgCheck()
    if (!is.null(formula) && !is(formula, "formula")) {
        ans <- paste("A formula of the form groups ~ x1 + x2 + ... (see",
            "?pcaLogisticR or ?glm).", sep = "")
        ArgumentCheck::addError(msg = ans, argcheck = Check)
    }

    if (!is.null(formula) && is(formula, "formula")) {
        vn <- try(attr(terms(formula), "term.labels"),
            silent = TRUE)
        if (inherits(vn, "try-error")) {
            vn <- try(setdiff(colnames(data), as.character(formula)[2]),
                silent = TRUE)
        }
        if (inherits(vn, "try-error"))
            stop("* Error in the formula")
        if (length(vn) < n.pc) {
            ans <- "The number of number predictor variables greater than "
            ans1 <- "the number of principal components: "
            ArgumentCheck::addError(msg = paste0(ans,
                ans1, n.pc), argcheck = Check)
        }
    }

    if (is.null(formula)) {
        ans <- "A formula or grouping variable must be provided."
        ArgumentCheck::addError(msg = ans, argcheck = Check)
    }

    ArgumentCheck::finishArgCheck(Check)

    m = nrow(data)
    if (floor(m/3) < n.pc) {
        ans <- "The number principal components: "
        ans1 <- " must be lower than the number of individuals N/3"
        warning(paste0(ans, n.pc, ans1))
    }

    idx <- grep(":", vn)[1]
    interactions <- FALSE
    intr <- NA
    if (!is.na(idx)) {
        interactions <- TRUE
        cnam <- colnames(data)
        intr <- strsplit(vn[ idx ], ":")[[1]]
        data$intr <- as.numeric(unname(data[, intr[1] ] * data[, intr[2] ]))
        vn <- sub(":", "_", vn)
        intr <- paste(intr, collapse = "_")
        colnames(data) <- c(cnam, intr)

        if (length(idx) > 1)
            warnings("*** More than one interaction was requested.",
                    "Only the first interaction was considered")
    }

    pc <- prcomp(x = data[ vn ], retx = TRUE, center = center,
                 scale. = scale, tol = tol, rank. = max.pc)

    cn <- colnames(pc$x)
    if (ncol(pc$x) > n.pc) {
        ind.coord <- as.data.frame(pc$x[, seq_len(n.pc)])
        colnames(ind.coord) <- cn[seq_len(n.pc)]
    } else {
        ind.coord <- as.data.frame(pc$x)
        colnames(ind.coord) <- cn[seq_len(ncol(pc$x))]
    }

    resp <- as.character(formula)[2]
    res <- data[, resp]
    l <- levels(res)
    res <- as.character(res)
    res[res == l[1]] <- 0
    res[res == l[2]] <- 1
    res <- as.numeric(res)
    data <- data.frame(res, ind.coord)
    colnames(data) <- c(resp, colnames(ind.coord))
    predictors <- paste(colnames(ind.coord), collapse = " + ")
    formula <- paste0(resp, " ~ ", predictors)
    model <- suppressWarnings(glm(formula = formula,
                            family = binomial(link = "logit"), data = data))

    model <- structure(
                list(
                    logistic = model,
                    pca = pc,
                    reference.level = l[1],
                    positive.level = l[2],
                    interaction =   {if (interactions)
                                        intr
                                    else NULL}),
                class = "pcaLogisticR")
    return(model)
}
#'
#' @rdname pcaLogisticR
#' @name predict.pcaLogisticR
#' @title Logistic regression using Principal Component Analysis (PCA)
#' @description Logistic regression using Principal Components from PCA as
#'     predictor variables
#' @param object To use with function 'predict'. A 'pcaLogisticR' object
#'     containing a list of two objects: 1) an object of class inheriting from
#'     'glm' and 2) an object of class inheriting from 'prcomp'.
#' @param newdata To use with function 'predict'. New data for classification
#'     prediction
#' @param type To use with function 'predict'. The type of prediction required:
#' 'class', 'posterior', 'pca.ind.coord', or 'all'. If type = 'all', function
#' 'predict.pcaLogisticR' ('predict') returns a list with:
#'         1) 'class': individual classification.
#'         2) 'posterior': probabilities for the positive class.
#'         3) 'pca.ind.coord': PC individual coordinate.
#' Each element of this list can be requested independently using parameter
#' 'type'.
#' @param ... Not in use.
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @export
predict.pcaLogisticR <- function(object, ...) UseMethod("predict", object)
predict.pcaLogisticR <- function(
                                object,
                                newdata,
                                type = c("class", "posterior",
                                         "pca.ind.coord", "all"),
                                ...) {

    if (!is(object, "pcaLogisticR")) {
        stop("* Parameter 'object' must be a model from class 'pcaLogisticR'")
    }
    type <- match.arg(type)
    ## predictor names
    vn <- rownames(object$pca$rotation)

    if (!is.null(newdata) && inherits(newdata, c("pDMP", "InfDiv")))
        newdata <- unlist(newdata)
    if (inherits(newdata, "GRanges")) {
        if (is.element("pos", vn)) {
            position <- function(gr) {
                chrs <- split(gr, seqnames(gr))
                gr <- lapply(chrs, function(grc) {
                    x <- start(grc)
                    x.min <- min(x)
                    x.max <- max(x)
                    if (x.min == Inf)
                        x.min = 0
                    if (x.max == -Inf)
                        x.max = 1
                    delta <- max(c(x.max - x, 1))
                    return((x - x.min)/(delta))
                })
                return(unlist(gr))
            }
            newdata$pos <- position(newdata)
        }
        newdata$logP <- log10(newdata$wprob + 2.2e-308)
        newdata <- mcols(newdata)
    }

    if (!is.null(newdata)) {
        if (!is.null(object$interaction)) {
            idx <- match(object$interaction, vn)
            cnam <- colnames(newdata)
            intr <- strsplit(vn[ idx ], "_")[[1]]
            newdata$intr <- as.numeric(
                            unname(newdata[, intr[1] ] * newdata[, intr[2] ]))
            vn <- sub(":", "_", vn)
            colnames(newdata) <- c(cnam, paste(intr, collapse = "_"))
        }

        newdata <- newdata[vn]
        newdata <- as.matrix(newdata)

        ## Centering and scaling new individuals
        dt.scaled <- scale(newdata, center = object$pca$center,
                           scale = object$pca$scale)
        ## Coordinates of the individuals
        coord_func <- function(ind, loadings) {
            x <- loadings * ind
            return(apply(x, 2, sum))
        }
        nc <- ncol(object$pca$x)
        loadings <- object$pca$rotation[, seq_len(nc)]
        if (nc == 1)
            loadings <- as.matrix(loadings)

        ind.coord <- data.frame(t(apply(dt.scaled, 1, coord_func,
                                        loadings)))
        if (nc == 1) {
            ind.coord <- as.data.frame(t(ind.coord))
            colnames(ind.coord) <- "PC1"
            row.names(ind.coord) <- NULL
        }
        rm(dt.scaled, loadings); gc()
    }

    predictClass <- function(object, dt) {
        pred <- predict(object$logistic, newdata = dt,
                        type = "response")
        PredClass <- rep(object$reference.level, nrow(dt))
        PredClass[pred > 0.5] <- object$positive.level
        return(PredClass)
    }

    if (is.null(newdata))
        ind.coord <- object$pca$x
    else {
        rm(newdata)
        gc()
    }

    pred <- switch(type,
                posterior = predict(object$logistic,
                                    newdata = ind.coord,
                                    type = "response"),
                class = predictClass(object = object,
                                    dt = ind.coord),
                pca.ind.coord = ind.coord,
                all = list( class = predictClass(object = object,
                                                dt = ind.coord),
                            posterior = predict(object$logistic,
                                                newdata = ind.coord,
                                                type = "response"),
                            pca.ind.coord = ind.coord))
    return(pred)
}
