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

# Reference group
ref0 = simulateCounts(num.samples = 4, sites = sites, alpha = alpha.ct,
                       beta = 0.5, size = 50, theta = 4.5,
                       sample.ids = c("R1", "R2", "R3", "R4"))
# Control group
ctrl = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.ct,
                       beta = 0.5, size = 50, theta = 4.5,
                       sample.ids = control.nam)
# Treatment group
treat = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.tt,
                        beta = 0.5, size = 50, theta = 4.5,
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
# C1 0.6770181  67.87938            340            333
# C2 0.6730444  67.43298            339            333
# C3 0.6639938  64.95108            341            337
# T1 0.9183686 129.36969            415            415
# T2 0.9307725 137.70360            415            415
# T3 0.9329020 137.59932            415            415

## The best nonlinear model
gof <- gofReport(HD, column = 9L, num.cores = 9L)
#      w2p_AIC w2p_R.Cross.val   w3p_AIC w3p_R.Cross.val   g2p_AIC g2p_R.Cross.val   g3p_AIC g3p_R.Cross.val bestModel
# C1 -47175.96       0.9996802        NA              NA -44083.56       0.9995605 -44083.47       0.9995638       w2p
# C2 -45656.31       0.9996065 -45656.54       0.9996072 -43515.39       0.9995038 -43564.18       0.9995185       w3p
# C3 -49292.07       0.9997625        NA              NA -47003.67       0.9997073 -47002.67       0.9997094       w2p
# T1 -52533.10       0.9994111 -52797.69       0.9994276 -45615.79       0.9986015 -46242.91       0.9987003       w3p
# T2 -51620.60       0.9993470 -52307.83       0.9993949 -44558.86       0.9984040 -45562.69       0.9986000       w3p
# T3 -55773.72       0.9995941 -55930.50       0.9995975 -47151.99       0.9988210 -47794.81       0.9989255       w3p

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

usethis::use_data(PS, HD, cutpoint, gof, dmps, lda_perf, qda_perf,
                  overwrite = TRUE )



















