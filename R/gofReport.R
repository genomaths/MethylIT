#' @rdname gofReport
#' @title Report the Best Fitted Probability Distribution Model
#' @description In Methyl-IT,  the signal detection step requires for the
#'     knowledge of the probability distribution of the methylation noise, which
#'     is the null hypothesis used to discriminate the signal from the noise.
#'     Both, signal and noise, are expressed in terms of an information
#'     divergence of methylation levels, which (currently in Methyl-IT) could be
#'     the Hellinger divergence or the total variation distance.
#'     As shown in reference (1), on statistical physics basis, the probability
#'     distribution of the noise is a member of the generalized gamma
#'     distribution. In particular, if the methylation changes on the DNA
#'     molecule follow a Poisson process, then the noise follows a Weibull
#'     probability distribution.
#'
#' Function 'gofReport' search for the best fitted model between the set of
#' models requested by the user. Two goodness-of-fit (GoF) criteria are applied
#' to select the best fitted model: Akaike's information criterion (AIC) and the
#' correlation coefficient of cross-validations for the nonlinear regressions
#' (R.Cross.val) (2). These criteria evaluate different information inferred
#' from the models. AIC deals with the trade-off between the goodness of fit of
#' the model and the complexity of the model, while R.Cross.val provides
#' information on the prediction power/performance of the model when confronted
#' with external dataset.
#'
#' Although the numerical algorithms to accomplish the nonlinear fit are not
#' perfect, in general, the model with the lowest AIC must have the highest
#' R.Cross.val. If the model with the lowest AIC has not the highest
#' R.Cross.val, then further analyzes are required. These analyzes could include
#' the visualization of the graphics for the density distribution, evaluation of
#' whether the parameter values can be meaningful or not, etc. Nevertheless, the
#' best model will, in general, lead to the identification of the greater amount
#' of potential DMPs and DMPs, as well as, the highest classification accuracy
#' estimated with functions \code{\link{estimateCutPoint}} and
#' \code{\link{evaluateDIMPclass}}. In the worse scenario, these observations
#' can ultimately lead to a post-hoc decision on which the best model is.
#'
#' @param HD An 'InfDiv' object returned by function
#'     \code{\link{estimateDivergence}}.
#' @param model A character vector naming the models to fit. Default is
#'     model = c('Weibull2P', 'Weibull3P', 'Gamma2P', 'Gamma3P'). See
#'     \code{\link{nonlinearFitDist}} for more options.
#'@param column An integer number denoting the index of the GRanges column where
#'     the information divergence is given. Default column = 9, which is the
#'     column where the Hellinger divergence values are reported (by default) by
#'     function \code{\link{estimateDivergence}}.
#' @param absolute Logic (default, FALSE). Total variation (TV, the difference
#'     of methylation levels) is normally an output in the downstream MethylIT
#'     analysis. If 'absolute = TRUE', then TV is transformed into |TV|, which
#'     is an information divergence that can be fitted as well.
#' @param output If output == 'all', the table with the GoF statistics is
#'     returned in a list together with the best fitted model and the
#'     corresponding statistics. Default is 'best.model'
#' @param confl_model Logic. If TRUE, then the best model based on highest
#'     R.Cross.val is returned for those samples where the model(s) with lowest
#'     AIC has not the highest R.Cross.val.
#' @param num.cores The number of cores to use in the nonlinear fit step, i.e.
#'     at most how many child processes will be run simultaneously (see
#'     \code{\link[BiocParallel]{bplapply}} function from BiocParallel package).
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Further arguments to pass to \code{\link{nonlinearFitDist}}.
#' @importFrom matrixStats rowMins
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return If 'output = 'best.model'', a character vector with the name of the
#' best fitted model for each sample and model statistics, which can be used the
#' next step to get the potential DMPs with function
#' \code{\link{getPotentialDIMP}}. Otherwise, it will return a list with the
#' table carrying the GoF values and the previously mentioned data.
#'
#' The best model selection is based on the lowest AIC. However, if 'confl_model
#' = TRUE', then the returned list will contain a sublist named 'confl_model'
#' with the best model selected based on the highest R.Cross.val. This sublist
#' is returned only if at least there is one sample where the model with the
#' lowest AIC has not the highest R.Cross.val. This report is useful to analyze
#' any conflict between models. For example, some times the best model selected
#' based on AIC has a R.Cross.val = 0.99945, while the highest R.Cross.val is
#' 0.99950. In such a situation the model with lowest AIC is still fine.
#' However, some times the model with the best AIC has some meaningless
#' parameter value. For example,  scale = 0.0000000001 in 'GGamma3P' or in
#' 'Weibull2P' models. This last situation can result from the numerical
#' algorithm used in the parameter estimation. The numerical algorithms for
#' nonlinear fit estimation are not perfect!
#'
#' For distribution models that include a location parameter ('Weibull3P',
#' 'Gamma3P', and 'GGamma4P') there is an additional important constraint, which
#' cannot be evaluated with 'AIC' or 'R.Cross.val': \eqn{\mu > 0}. That is, the
#' information divergences are strictly positive magnitudes, therefore, any best
#' model in terms of AIC with values \eqn{\mu < 0} are meaningless and rejected.
#'
#' @export
#'
#' @author Robersy Sanchez 11/25/2019 <https://github.com/genomaths>
#' @examples
#' ## Loading information divergence dataset
#' data(HD)
#' ## Subsetting object HD (for the sake of runnig a faster example)
#' hd <- lapply(HD, function(x) x[seq_len(100)], keep.attr = TRUE)
#' ## The GoF report
#' dt <- gofReport(hd)
#'
#' ## To get the potential DMPs
#' ps_dmp <- getPotentialDIMP(LR = hd, nlms = dt$nlms, div.col = 9L,
#'                             dist.name = dt$bestModel)
#'
#' @references
#' \enumerate{
#'     \item R. Sanchez and S. A. Mackenzie, “Information Thermodynamics of
#'           Cytosine DNA Methylation,” PLoS One, vol. 11, no. 3, p. e0150427,
#'             Mar. 2016.
#'     \item Stevens JP. Applied Multivariate Statistics for the Social
#'             Sciences. Fifth Edit. Routledge Academic; 2009.
#' }
#' @seealso \code{\link{nonlinearFitDist}}

gofReport <- function(HD,
                    model = c("Weibull2P", "Weibull3P", "Gamma2P", "Gamma3P"),
                    column = 9, absolute = FALSE,
                    output = c("best.model", "all"),
                    confl_model = FALSE, num.cores = 1L,
                    verbose = FALSE, ...) {
    validateClass(HD)
    output <- match.arg(output)
    model <- unique(model)  # just in case

    idx <- match(model, c("Weibull2P", "Weibull3P",
        "Gamma2P", "Gamma3P", "GGamma3P", "GGamma4P"))

    if (any(is.na(idx))) {
        stop("*** At least one of the requested models is not valid \n",
            "The valid model names are: \n", "'Weibull2P', 'Weibull3P',
            'Gamma2P', 'Gamma3P', 'GGamma3P', and 'GGamma4P'")
    }
    nams <- c("w2p", "w3p", "g2p", "g3p", "gg3p", "gg4p")
    nams <- nams[idx]
    sn <- names(HD)

    stp <- seq_along(model)
    nlms <- list()
    pb <- txtProgressBar(max = max(stp), style = 3,
        char = "=")
    for (k in stp) {
        mdl <- suppressWarnings(nonlinearFitDist(LR = HD,
            column = column, absolute = absolute, num.cores = num.cores,
            dist.name = model[k], verbose = verbose))
        names(mdl) <- paste(names(mdl), nams[k], sep = "_")
        stat <- vapply(mdl, function(x) {
            AIC <- as.numeric(x$AIC[1])
            R.Cross.val <- as.numeric(x$R.Cross.val[1])
            if (is.element(x$model[1], c("Weibull3P",
                "Gamma3P", "GGamma4P"))) {
                ## An important constraint must hold: mu > 0
                x <- x$Estimate[match("mu", rownames(x))]
                if (x < 0) {
                    AIC <- Inf
                    R.Cross.val <- 0
                }
            }
            c(AIC = AIC, R.Cross.val = R.Cross.val)
        }, FUN.VALUE = numeric(2))

        colnames(stat) <- sn
        rownames(stat) <- paste(nams[k], rownames(stat),
            sep = "_")

        nlms[[k]] <- mdl
        if (k == 1)
            stats <- data.frame(t(stat))
        else stats <- cbind(stats, data.frame(t(stat)))
        setTxtProgressBar(pb, k)
    }
    close(pb)
    if (length(HD) > 1)
        nlms <- unlist(nlms, recursive = FALSE) else {
        nms <- sapply(nlms, function(x) names(x[1]))
        nlms <- lapply(nlms, function(x) x[[1]])
        names(nlms) <- nms
    }

    cat("\n *** Creating report ... \n")
    aic_col <- grep("AIC", colnames(stats))
    r_col <- grep("Cross", colnames(stats))

    mdl <- apply(stats, 1, function(x) {
        unique(c(nams[which.min(x[aic_col])], nams[which.max(x[r_col])]))
    })

    if (inherits(mdl, "matrix")) {
        ns <- colnames(mdl)
        ### Split the matrix into a list by columns
        mdl <- split(mdl, col(mdl))
        names(mdl) <- ns
    }

    conflict <- (inherits(mdl, "list"))
    if (conflict) {
        issue <- unlist(lapply(mdl, function(x) length(x) >
            1))
        if (confl_model) {
            mdl2 <- mdl
            mdl2[issue] <- vapply(mdl[issue], function(x) x[2],
                character(1))
            mdl2 <- unlist(mdl2)
            mdl2 <- mdl2[issue]
            issuekey <- paste(names(mdl2), mdl2, sep = "_")
            issue_nlms <- nlms[match(issuekey, names(nlms))]
        }
        mdl[issue] <- vapply(mdl[issue], function(x) x[1],
            character(1))

        bestAIC <- unlist(mdl)
        mdl[issue] <- "Needs revision"
        mdl <- unlist(mdl)
        warning("The best fitted model for sample(s) ",
            paste(sn[issue], collapse = ", "),
                " require(s) for further analysis. \n",
            "The model with the lowest AIC must have the highest R.Cross.val")
    } else bestAIC <- unlist(mdl)

    modelkey <- paste(names(bestAIC), bestAIC, sep = "_")
    nlms <- nlms[match(modelkey, names(nlms))]
    names(nlms) <- sn
    nlms <- structure(nlms, class = c("ProbDistrList",
        "list"))

    bestModel <- model[match(bestAIC, nams)]
    names(bestModel) <- names(HD)
    stats$bestModel <- mdl
    if (output == "best.model") {
        print(stats)
        cat("\n")
        res <- list(bestModel = bestModel, nlms = nlms)
    } else res <- list(stats = stats, bestModel = bestModel,
        nlms = nlms)
    if (confl_model && conflict)
        res$confl_model <- issue_nlms
    return(res)
}

