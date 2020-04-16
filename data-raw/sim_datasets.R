# ============================================================================ #
#
# ======== Script used to generate the datasets used in the examples ========= #
#
# ============================================================================ #
library(MethylIT)
library(MethylIT.utils)

bmean <- function(alpha, beta) alpha/(alpha + beta)
alpha.ct <- 0.09
alpha.tt <- 0.2

# The number of cytosine sites to generate
sites = 10000
# Set a seed for pseudo-random number generation
set.seed(124)
control.nam <- c("C1", "C2", "C3")
treatment.nam <- c("T1", "T2", "T3")
strands <- sample(c("+", "-"), size = sites, replace = TRUE)


# Reference group
ref0 = simulateCounts(num.samples = 4, sites = sites, alpha = alpha.ct,
                       beta = 0.5, size = 50, theta = 4.5, strand = strands,
                       sample.ids = c("R1", "R2", "R3", "R4"))
# Control group
ctrl = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.ct,
                       beta = 0.5, size = 50, theta = 4.5, strand = strands,
                       sample.ids = control.nam)
# Treatment group
treat = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.tt,
                        beta = 0.5, size = 50, theta = 4.5, strand = strands,
                        sample.ids = treatment.nam)

# Reference sample
ref = poolFromGRlist(ref0, stat = "mean", num.cores = 4L, verbose = FALSE)

# Methylation level divergences
HD <- estimateDivergence(ref = ref, indiv = c(ctrl, treat), Bayesian = TRUE,
                       num.cores = 6L, percentile = 1, JD = TRUE,
                       verbose = FALSE)

critical.val <- do.call(rbind, lapply(HD, function(x) {
  hd.95 = quantile(x$hdiv, 0.95)
  tv.95 = quantile(abs(x$bay.TV), 0.95)
  return(c(tv = tv.95, hd = hd.95, num.sites.hd95 = sum(x$hdiv > hd.95),
           num.sites.tv95 = sum(x$bay.TV > tv.95)))}))
critical.val
#       tv.95%    hd.95% num.sites.hd95 num.sites.tv95
# C1 0.6768858  67.38872            340            336
# C2 0.6733885  66.11182            339            335
# C3 0.6678430  64.87552            341            338
# T1 0.9200526 130.33587            415            415
# T2 0.9319064 138.29089            413            413
# T3 0.9341989 138.70751            413            415

## The best nonlinear model
gof <- gofReport(HD, column = 9L, num.cores = 9L)
#      w2p_AIC w2p_R.Cross.val   w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val bestModel
# C1 -42646.25       0.9994076        NA              NA -41619.24       0.9993085       Inf       0.0000000       w2p
# C2 -43059.49       0.9994484        NA              NA -41989.47       0.9993685       Inf       0.0000000       w2p
# C3 -43622.49       0.9994884        NA              NA -43006.08       0.9994038       Inf       0.0000000       w2p
# T1 -52310.16       0.9993933 -52385.16       0.9993979 -45845.37       0.9986357 -46207.44       0.9986943       w3p
# T2 -50657.44       0.9992650 -50977.10       0.9992900 -44329.09       0.9983596 -44967.16       0.9984987       w3p
# T3 -54625.62       0.9995332 -54672.57       0.9995346 -46946.44       0.9987922 -47323.73       0.9988607       w3p

gof$bestModel
#         C1          C2          C3          T1          T2          T3
# "Weibull2P" "Weibull3P" "Weibull2P" "Weibull3P" "Weibull3P" "Weibull3P"

## Next, the potential signal can be estimated
PS <- getPotentialDIMP(LR = HD, nlms = gof$nlms, div.col = 9L, alpha = 0.05,
                       dist.name = gof$bestModel, tv.col = 7L, tv.cut = 0.68)

cutpoint = estimateCutPoint(LR = PS, simple = FALSE,
                          column = c(hdiv = TRUE, TV = TRUE,
                                     wprob = TRUE, pos = TRUE),
                          classifier1 = "qda",
                          control.names = control.nam,
                          treatment.names = treatment.nam,
                          clas.perf = TRUE, prop = 0.6,
                          div.col = 9L)


## DMPs are selected using the cupoints
dmps <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint,
                    tv.cut = 0.68)

## Classification of DMPs into two clases: DMPS from control and DMPs
## from treatment samples and evaluation of the classifier performance
## (for more details see ?evaluateDIMPclass).
lda_perf <- evaluateDIMPclass(LR = dmps,
                        column = c(hdiv = TRUE, TV = TRUE,
                                    wprob = TRUE, pos = TRUE),
                        classifier = 'lda', n.pc = 4L,
                        control.names =  c('C1', 'C2', 'C3'),
                        treatment.names = c('T1', 'T2', 'T3'),
                        center = TRUE, scale = TRUE, prop = 0.6)

qda_perf <- evaluateDIMPclass(LR = dmps,
                              column = c(hdiv = TRUE, TV = TRUE,
                                         wprob = TRUE, pos = TRUE),
                              classifier = 'qda', n.pc = 4L,
                              control.names =  c('C1', 'C2', 'C3'),
                              treatment.names = c('T1', 'T2', 'T3'),
                              center = TRUE, scale = TRUE, prop = 0.6)

pcaLda_perf <- evaluateDIMPclass(LR = dmps,
                              column = c(hdiv = TRUE, TV = TRUE,
                                         wprob = TRUE, pos = TRUE),
                              classifier = 'pca.lda', n.pc = 4L,
                              control.names =  c('C1', 'C2', 'C3'),
                              treatment.names = c('T1', 'T2', 'T3'),
                              center = TRUE, scale = TRUE, prop = 0.6)

pcaQda_perf <- evaluateDIMPclass(LR = dmps,
                              column = c(hdiv = TRUE, TV = TRUE,
                                         wprob = TRUE, pos = TRUE),
                              classifier = 'pca.qda', n.pc = 4L,
                              control.names =  c('C1', 'C2', 'C3'),
                              treatment.names = c('T1', 'T2', 'T3'),
                              center = TRUE, scale = TRUE, prop = 0.6)


logit_perf <- evaluateDIMPclass(LR = dmps,
                              column = c(hdiv = TRUE, TV = TRUE,
                                         wprob = TRUE, pos = TRUE),
                              classifier = 'logistic', n.pc = 4L,
                              control.names =  c('C1', 'C2', 'C3'),
                              treatment.names = c('T1', 'T2', 'T3'),
                              center = TRUE, scale = TRUE, prop = 0.6)

pcalogit_perf <- evaluateDIMPclass(LR = dmps,
                                column = c(hdiv = TRUE, TV = TRUE,
                                           wprob = TRUE, pos = TRUE),
                                classifier = 'pca.logistic', n.pc = 4L,
                                control.names =  c('C1', 'C2', 'C3'),
                                treatment.names = c('T1', 'T2', 'T3'),
                                center = TRUE, scale = TRUE, prop = 0.6)


# ========================= dataset for countTes2 =============================
set.seed(133) # Set a seed
## A GRanges object with the count matrix in the metacolumns is created
countData <- matrix(sample.int(200, 500, replace = TRUE), ncol = 4)
colnames(countData) <- c('A1','A2','B1','B2')

start <- seq(1, 25e4, 2000)
end <- start + 1000
chr <- c(rep('chr1', 70), rep('chr2', 55))
GR <- GRanges(seqnames = chr, IRanges(start = start, end = end))
mcols(GR) <- countData

## Gene IDs
names(GR) <- paste0('gene', 1:length(GR))

## An experiment design is set.
colData <- data.frame(condition = factor(c('A','A','B','B')),
                    c('A1','A2','B1','B2'), row.names = 2)

## A RangedGlmDataSet is created
ds <- glmDataSet(GR = GR, colData = colData)


usethis::use_data(PS, HD, cutpoint, gof, dmps, lda_perf, qda_perf,
                  logit_perf, pcalogit_perf, pcaLda_perf, pcaQda_perf, ds,
                  overwrite = TRUE )

data(PS, HD, cutpoint, gof, dmps, lda_perf, qda_perf,
     logit_perf, pcalogit_perf, pcaLda_perf, pcaQda_perf)












