#' @rdname countTest2
#'
#' @title Regression Test for Count
#' @description Perform Poisson and Negative Binomial regression analysis to
#'     compare the counts from different groups, treatment and control. The
#'     difference between functions 'countTest2' and 'countTest' resides in the
#'     estimation of the prior weights used in Negative Binomial generalized
#'     linear model.
#' @details A pairwise group comparison, control versus treatment, is performed.
#'     The experimental design settings must be introduced using function
#'     \code{link{glmDataSet}} to provide dataset (DS) object.
#' @param DS A 'glmDataSet' object, which is created with function
#'     \code{link{glmDataSet}}.
#' @param countFilter whether or not to filter the counts according to the
#'     minimum count per region per each individual/sample, which is setting by
#'     "minCountPerIndv".
#' @param CountPerBp for each group the count per bp must be equal or greater
#'     than CountPerBp. The filter is applied if 'CountPerBp' is given and if
#'     'x' DESeqDataSet object has the rowRanges as a GRanges object on it
#' @param minCountPerIndv each gene or region must have more than
#'     'minCountPerIndv' counts (on average) per individual in at least one
#'     group.
#' @param maxGrpCV A numerical vector. Maximum coefficient of variance for each
#'     group. Defaul maxGrpCV = NULL. The numbers maxGrpCV[1] and maxGrpCV[2]
#'     will be taken as the maximun variances values permitted in control and
#'     in treatment groups, repectively. If only maxGrpCV[1] is provided, then
#'     maxGrpCV = c(maxGrpCV[1], maxGrpCV[1]). This parameter is addressed to
#'     prevent testing regions where intra-group variations are very large,
#'     e.g.: control = c(1,0,1,1) and traatment = c(1, 0, 1, 40). The
#'     coefficient of variance for the treatment group is 1.87, very high. The
#'     generalized linear regression analysis would yield statistical
#'     significant group differences, but evidently there is something wrong in
#'     one of the treatment samples. We would try the application of further
#'     statistical smoothing approach, but we prefer to leave the user decide
#'     which regions to test.
#' @param FilterLog2FC if TRUE, the results are filtered using the minimun
#'     absolute value of log2FoldChanges observed to accept that a gene in the
#'     treatment is differentially expressed in respect to the control
#' @param pAdjustMethod method used to adjust the results; default: BH
#' @param pvalCutOff cutoff used then a p-value adjustment is performed
#' @param MVrate Minimum Mean/Variance rate.
#' @param Minlog2FC minimum logarithm base 2 of fold changes.
#' @param test A character string matching one of "Wald" or "LRT". If test =
#'     "Wald", then the p-value of the Wald test for the coefficient of the
#'     independent variable (\emph{treatment group}) will be reported.
#'     If test = "LRT", then the p-value from a likelihood ratio test given by
#'     \code{\link[stats]{anova}} function from \emph{stats} packages will be
#'     the reported p-value for the group comparison when the best fitted model
#'     is the negative binomial. As suggested for \code{\link[stats]{glm}}, if
#'     best fitted model is Poisson or quasi-Poisson, then the best test is
#'     'Chi-squared' or 'F-test', respectively. So, for the sake of simplicity,
#'     the corresponding suitable test will be applied when test = "LRT".
#'
#' @param scaling integer (default 1). Scaling factor estimate the
#'     signal density as: scaling x "DIMP-Count-Per-Bp". For example,
#'     if scaling = 1000, then signal density denotes the number of DIMPs in
#'      1000 bp.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param saveAll if TRUE all the temporal results are returned
#' @param verbose if TRUE, prints the function log to stdout
#
#' @return a data frame or GRanges object (if the DS contain the GRanges
#'     information for each gene) with the test results and original count
#'     matrix, plus control and treatment signal densities and their variation.
#'
#' @examples
#' set.seed(133) # Set a seed
#' ## A GRanges ogbject with the count matrix in the metacolumns is created
#' countData <- matrix(sample.int(200, 500, replace = TRUE), ncol = 4)
#' colnames(countData) <- c("A1","A2","B1","B2")
#' start <- seq(1, 25e4, 2000)
#' end <- start + 1000
#' chr <- c(rep("chr1", 70), rep("chr2", 55))
#' GR <- GRanges(seqnames = chr, IRanges(start = start, end = end))
#' mcols(GR) <- countData
#' ## Gene IDs
#' names(GR) <- paste0("gene", 1:length(GR))
#'
#' ## An experiment design is set.
#' colData <- data.frame(condition = factor(c("A","A","B","B")),
#'                       c("A1","A2","B1","B2"),
#'                       row.names = 2)
#' ## A RangedGlmDataSet is created
#' ds <- glmDataSet(GR = GR, colData = colData)
#'
#' ## The gneralized linear model pairwise group comparison, group 'A'
#' ## ('control') versus 'B' (treatment) is performed.
#' countTest2(ds, num.cores = 1L, maxGrpCV = c(0.4, 0.4), verbose = FALSE)
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom stats p.adjust sd
#' @importFrom BiocGenerics counts
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom IRanges width
#' @importFrom stats var
#'
#' @export
countTest2 <- function(DS, num.cores=1, countFilter=TRUE, CountPerBp=NULL,
                      minCountPerIndv=3, maxGrpCV=NULL, FilterLog2FC=TRUE,
                      pAdjustMethod="BH", pvalCutOff=0.05, MVrate=0.98,
                      Minlog2FC=0.5, test = c("Wald", "LRT"), scaling =1L,
                      tasks=0L, saveAll=FALSE, verbose=TRUE ) {

   group <- DS$colData$condition
   lev <- DS$levels
   sample.names <- DS$sampleNames
   res <- GRanges()

   # ========================== Filtering Block ============================== #
   if (countFilter) {

       ## --------------------- Remove var == 0 positions -------------------- #
       dc <- DS$counts
       vars <- apply(dc, 1, var)
       DS <- DS[which(vars > 0)]

       # -------------------------- minCountPerIndv -------------------------- #
       ## Each gene must have more than 'minCountPerIndv' read-counts
       ## per individual in at least one group
       dc <- DS$counts
       g1 <- which(lev[1] ==  group)
       g2 <- which(lev[2] ==  group)
       m1 <- rowMeans(data.frame(dc[ ,g1]))
       m2 <- rowMeans(data.frame(dc[ ,g2]))
       idx <- which(m1 >= minCountPerIndv |
                   m2 >= minCountPerIndv)
       res <- NULL

       if (length(idx) == 0) {
           warning("* No genomic region passed the 'minCountPerIndv' \n",
               "filtering conditions")
       } else {
           DS <- DS[idx];
           dc <- DS$counts
           res <- GRanges()
       }

       # ---------------------------- CountPerBp ----------------------------- #
       if (!is.null(res) &&
           !is.null(CountPerBp) && class(DS$GR) == "GRanges") {
           ## For each group the count per bp must be equal or greater
           ## than CountPerBp
           size <- width(DS$GR)
           idx <- which((unname((rowMeans(dc[,g1])) / size) >= CountPerBp) |
                       (unname((rowMeans(dc[,g2])) / size) >= CountPerBp))
           if (length(idx) == 0) {
               warning("* No genomic region passed the 'CountPerBp' \n",
                   "filtering conditions")
               res <- NULL
           } else {
              DS <- DS[idx]
              dc <- DS$counts
           }
       }

       # ---------------------------- maxGrpCV ------------------------------- #
       if (!is.null(maxGrpCV) && !is.null(res)) {
           if (length(maxGrpCV) == 1 && is.numeric(maxGrpCV))
               maxGrpCV = c(maxGrpCV, maxGrpCV)
           dc <- DS$counts
           # Add Bayesian correction assuming uniform prior
           # distribution for count data.
           if (any(dc == 0)) dc <- dc + 1

           if (nrow(dc) > 1) {
               g1 <- which(lev[1] ==  group)
               g2 <- which(lev[2] ==  group)
               m1 <- rowMeans(dc[ ,g1])
               m2 <- rowMeans(dc[ ,g2])
               m1 <- sapply(m1, function(x) max(x,1)) # mean = 0 undefine the CV
               m2 <- sapply(m2, function(x) max(x,1))
               cv1 <- apply(dc[ ,g1], 1, sd)/m1
               cv2 <- apply(dc[ ,g2], 1, sd)/m2
           } else {
               m1 <- max(mean(dc[g1]), 1)
               m2 <- max(mean(dc[g2]), 1)
               cv1 <- sd(dc[ ,g1])/m1
               cv2 <- sd(dc[ ,g2])/m2
           }
           idx <- intersect(which(cv1 <= maxGrpCV[1]),
                           which(cv2 <= maxGrpCV[2]))

           if (length(idx) == 0) {
               warning("* No genomic region passed the 'maxGrpCV' \n",
                       "filtering conditions")
               res <- NULL
           } else {
               DS <- DS[idx]
               rm(m1,m2,cv1,cv2,dc); gc()
           }
       }
       if (verbose)
           message("*** Number of GR after filtering counts ", length(idx))

   }

   # ==================== Data preprocessing Block =========================== #
   if (!is.null(res)) {
       # == A rough estimation of prior weights ===
       X <- DS$counts
       # if only one range:
       if (is.null(nrow(X))) { # A trick
           baseMeanAndVar <- data.frame(baseMean = mean(X),
                                        baseVar = var(X))
           X <- data.frame(t(X))
       } else {
           baseMeanAndVar <- data.frame(baseMean=rowMeans(X),
                                       baseVar=rowVars(X))
       }

       m1 <- rowMeans(log(X[, group == lev[1]] + 1))
       m2 <- rowMeans(log(X[, group == lev[2]] + 1))
       v1 <- rowVars(as.matrix(log(X[, group == lev[1]] + 1)))
       v2 <- rowVars(as.matrix(log(X[, group == lev[2]] + 1)))
       v12 <- var(c(m1, m2), na.rm = TRUE, use = "everything")
       logBaseMean <- log(baseMeanAndVar[, 1] + 1)
       disp <- v1/(v12/(logBaseMean - m1)^2) + v2/(v12/(logBaseMean - m2)^2)

       ## Var = mean + dispersion * mean ^ 2
       ## weights
       ## w = 1 / (Mean + disp * Mean ^ 2)
       if(!is.null(disp)) {
          w1 <- 1 / (m1 + disp * m1 ^ 2) # weights
          w2 <- 1 / (m2 + disp * m2 ^ 2)
       } else {
          d1 <- apply(log(X[ ,group == lev[1]] + 1), 1, var)
          d2 <- apply(log(X[ ,group == lev[2]] + 1), 1, var)
          w1 <- 1 / (m1 + d1 * m1 ^ 2) # weights
          w2 <- 1 / (m2 + d2 * m2 ^ 2)
       }
       w1[is.infinite(w1)] <- 1
       w2[is.infinite(w2)] <- 1
       w <- cbind(w1, w2)

       # ============================= GLM Block ============================= #
       if (verbose) message("*** GLM...")
       if (num.cores == 1) {
           tests = c()
           for (k in 1:nrow(X)) {
               if (verbose) message("*** Processing sample", k, "\n")
               tests <- rbind(tests, .estimateGLM(x=X[k, ], groups=group,
                                               baseMV=baseMeanAndVar[k, ],
                                               w=w[k, ], MVrate=MVrate,
                                               test=test[1]))
         }
       }
       if (num.cores > 1) {
           if (Sys.info()['sysname'] == "Linux") {
               bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
           } else {
               bpparam <- SnowParam(workers = num.cores, type = "SOCK")
           }
           tests <- bplapply(1:nrow(X),
                           function(k) .estimateGLM(x=X[k, ], groups=group,
                                                   baseMV=baseMeanAndVar[k, ],
                                                   w=w[k, ], MVrate=MVrate,
                                                   test=test[1]),
                           BPPARAM=bpparam)
           tests=do.call(rbind, tests)
       }

       # ================= Post processing filtering Block =================== #
       DS$optionData <- DataFrame(tests)
       DS <- DS[!is.na(tests$log2FC)]
       if (FilterLog2FC && is.null(pvalCutOff) && !saveAll) {
           DS <- DS[abs(DS$optionData$log2FC) > Minlog2FC]
       }
       if (!is.null(pvalCutOff) && !saveAll) {
           if (FilterLog2FC)
               DS <- DS[ abs(DS$optionData$log2FC) > Minlog2FC ]
           DS$optionData$adj.pval <- p.adjust(DS$optionData$pvalue,
                                                   method=pAdjustMethod)
           DS <- DS[ DS$optionData$adj.pval < pvalCutOff ]
       } else {
           if (!is.null(pvalCutOff) && saveAll) {
               DS$optionData$adj.pval <- p.adjust(DS$optionData$pvalue
                                                   , method=pAdjustMethod)
               if (FilterLog2FC) {
                   idx <- which(abs(DS$optionData$log2FC) > Minlog2FC)
                   pval <- DS$optionData$pvalue[idx]
                   DS$optionData$adj.pval[idx] <- p.adjust(pval,
                                                           method=pAdjustMethod)
               }
           }
       }

       # ============================= Output Block ========================== #
       if (class(DS$GR) == "GRanges") {
           GR <- DS$GR
           dc <- DS$counts
           g1 <- which(lev[1] ==  group)
           g2 <- which(lev[2] ==  group)
           size <- width(GR)
           CT.CountPerBp <- unname((rowSums(dc[ ,g1]) / length(g1)) / size)
           TT.CountPerBp <- unname((rowSums(dc[ ,g2]) / length(g2)) / size)
           mcols(GR) <- data.frame(DS$counts, DS$optionData,
                                   CT.SignalDensity=(scaling * CT.CountPerBp),
                                   TT.SignalDensity=(scaling * TT.CountPerBp),
                                   SignalDensityVariation=scaling *
                                   (TT.CountPerBp - CT.CountPerBp))
           res <- GR[order(as.factor(seqnames(GR)), start(GR)), ]
       } else {
         res <- cbind(DS$counts, DS$optionData)
       }
   } else {
       # ------ If not individual sample passed the filtering conditions ----- #
       # An empty GRanges object will be returned
       if (class(DS$GR) == "GRanges") {
           res <- GRanges()
           x <- matrix(integer(), nrow = 1, ncol = length(sample.names))
           colnames(x) <- sample.names
           mcols(res) <- DataFrame(x[NULL,], log2FC=double(), pvalue=double(),
                                   model=character(), adj.pval=double(),
                                   CT.SignalDensity=double(),
                                   TT.SignalDensity=double(),
                                   SignalDensityVariation=double())
       } else {
           res <- GRanges()
           x <- matrix(integer(), nrow = 1, ncol = length(sample.names))
           colnames(x) <- sample.names
           mcols(res) <- DataFrame(x[NULL,], log2FC=double(), pvalue=double(),
                                   model=character(), adj.pval=double())
       }
   }
   return(res)
}
