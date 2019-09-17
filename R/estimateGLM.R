#' @rdname estimateGLM
#'
#' @title Poisson and Negative Binomial regression analysis.
#' @description This function is called internally by countTest function. You
#'     would need to call it directly only in very special cases.
#' @description Perform Poisson and Negative Binomial regression analysis to
#'     compare the counts from different groups, treatment and control x:
#'     vector of counts groups: factor labeling the members from each group
#'     Evaluated models are "Poisson", "Quasipoisson", "Neg.Binomial.W", and
#'     "Neg.Binomial"
#'
#' @param x Matrix of counts.
#' @param groups Groups information derived from a
#'     \code{\link[MethylIT]{glmDataSet}} object.
#' @param baseMV Mean and variance of group counts. If
#'     baseMean >= baseVar*MVrate, then the nonlinear fit to "Poisson" and
#'     "QuasiPoisson" models are performed, otherwise only the nonlinear fit to
#'     "Neg.Binomial" and "Neg.Binomial with weights" models are performed
#' @param w group weights used in glm procedure
#' @param MVrate Minimum Mean/Variance rate.
#' @param test A character string matching one of "Wald" or "LRT". If test =
#'     "Wald", then the p-value of the Wald test for the coefficient of the
#'     independent variable (\emph{treatment group}) will be reported.
#'     If test = "LRT", then the p-value from a likelihood ratio test given by
#'     \code{\link[stats]{anova}} function from \emph{stats} packages will be
#'     the reported p-value for the group comparison when the best fitted model
#'     is the negative binomial. As suggested for \code{link[stats]{glm}}, if
#'     best fitted model is Poisson or quasi-Poisson, then the best test is
#'     'Chi-squared' or 'F-test', respectively. So, for the sake of simplicity,
#'     the corresponding suitable test will be applied when test = "LRT".
#'
#' @return GLM model of the group comparison for the given genomic region
#'
#' @importFrom stats glm poisson quasipoisson relevel glm.control
#' @importFrom MASS glm.nb negative.binomial
#' @keywords internal
.estimateGLM <- function(x, groups, baseMV, w, MVrate,
                         test = c("Wald", "LRT")) {

   neg_Bin <- function(data, weights, control) {

       res <- try(suppressWarnings(glm.nb(count ~ group,
                                           data=data, weights=weights,
                                           control=control)),
               silent=TRUE)
       if (inherits(res, "try-error")) {
           res <- try(suppressWarnings(glm(count ~ group,
                                       family=negative.binomial(theta = 1),
                                       data=data, weights=weights,
                                       control=control)),
                       silent=TRUE)
       }
       return(res)
   }

   model <- function(dt, model, weights=NULL) {
       if (!is.null(weights)) {
           lev <- levels(dt$group)
           weights <- c(rep(weights[1], sum(groups == lev[1])), rep(weights[2],
                               sum(groups == lev[2])))
       }
       controls <- glm.control(maxit=(10 ^ 3), epsilon=(1e-8), trace=FALSE)
       switch(model,
           Poisson=try(suppressWarnings(glm(count ~ group,
                       family=poisson(link="log"), data=dt, control=controls)),
                       silent=TRUE),
           QuasiPoisson=try(suppressWarnings(glm(count ~ group,
                       family=quasipoisson(link="log"), data=dt,
                       control=controls)), silent=TRUE),
           Neg.Binomial=try(suppressWarnings(glm.nb(count ~ group, data=dt,
                                                    control=controls)),
                            silent=TRUE),
           Neg.Binomial.W=neg_Bin(data=dt, weights = weights, control=controls))
   }

   levels(groups) <- c("CT", "TT") ## Control vs Treatment
   relevel(groups, ref="CT")
   # Add pseudocounts (1s) to all the individuals if at least one individual has
   # zero count.
   if(sum(x == 0) > 0) x <- x + 1

   dat <- data.frame(group=groups, count=x)
   mdl <- list(Eval=FALSE)
   mdl <- list()
   if (baseMV$baseMean >= baseMV$baseVar * MVrate) {
       mdls <- c("Poisson", "QuasiPoisson", "Neg.Binomial", "Neg.Binomial.W")
       mdl$P <- .evaluateModel(model(dat, "Poisson"), test = test[1])
       mdl$Q <- .evaluateModel(model(dat, "QuasiPoisson"), test = test[1])
       mdl$NB <- .evaluateModel(model(dat, "Neg.Binomial"), test = test[1])
       mdl$NBW <- .evaluateModel(model(dat, "Neg.Binomial.W", weights=w),
                                 test = test[1])
       aic <- c(mdl$P$AIC, mdl$Q$AIC, mdl$NB$AIC, mdl$NBW$AIC)

       aic <- aic[c(mdl$P$Eval, mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
       mdls <- mdls[c(mdl$P$Eval, mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
       mdl <- mdl[c(mdl$P$Eval, mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
       # print(mdls)
   } else {
       mdls <- c("QuasiPoisson", "Neg.Binomial", "Neg.Binomial.W")
       mdl$Q <- .evaluateModel(model(dat, "QuasiPoisson"), test = test[1])
       mdl$NB <- .evaluateModel(model(dat, "Neg.Binomial"), test = test[1])
       mdl$NBW <- .evaluateModel(model(dat, "Neg.Binomial.W", weights=w),
                                 test = test[1])

       aic <- c(mdl$Q$AIC, mdl$NB$AIC, mdl$NBW$AIC)

       aic <- aic[c(mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
       mdls <- mdls[c(mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
       mdl <- mdl[c(mdl$Q$Eval, mdl$NB$Eval, mdl$NBW$Eval)]
   }
   if (length(mdl) > 0) {
       Eval <- TRUE
       if (length(mdls) > 1) {
           ind <- which(aic == min(aic))
           mdl <- mdl[ind][[1]]
           mdls <- mdls[ind][1]
       } else mdl <- mdl[[1]]
   } else Eval <- FALSE
   if (Eval) {
       res <- data.frame(log2FC=mdl$log2FC, pvalue=mdl$coef.pval, model=mdls)
   } else {
       res <- data.frame(log2FC=NA, pvalue=NA, model=NA)
   }
   return(res)
}
