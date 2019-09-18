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
ref0 = simulateCounts(num.samples = 4, sites = sites, alpha = alpha.ct, beta = 0.5,
                      size = 50, theta = 4.5, sample.ids = c("R1", "R2", "R3"))
# Control group
ctrl = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.ct, beta = 0.5,
                      size = 50, theta = 4.5, sample.ids = control.nam)
# Treatment group
treat = simulateCounts(num.samples = 3, sites = sites, alpha = alpha.tt, beta = 0.5,
                       size = 50, theta = 4.5, sample.ids = treatment.nam)

# Reference sample
ref = poolFromGRlist(ref0, stat = "mean", num.cores = 4L, verbose = FALSE)

# Methylation level divergences
HD <- estimateDivergence(ref = ref, indiv = c(ctrl, treat), Bayesian = TRUE,
                           num.cores = 6L, percentile = 1, verbose = FALSE)

## To perform the nonlinear regression analysisx
nlms <- nonlinearFitDist(HD, column = 4, verbose = FALSE)

## Next, the potential signal can be estimated
PS <- getPotentialDIMP(LR = HD, nlms = nlms, div.col = 9L, alpha = 0.05,
                       tv.col = 7L, tv.cut = 0.92)

cutpoint = estimateCutPoint(LR = PS, simple = FALSE,
                          column = c(hdiv = TRUE, TV = TRUE,
                                     wprob = TRUE, pos = TRUE),
                          classifier1 = "qda",
                          control.names = control.nam,
                          treatment.names = treatment.nam,
                          tv.cut = 0.92, clas.perf = TRUE, prop = 0.6,
                          div.col = 9L)

devtools::use_data(PS, HD, cutpoint, nlms, overwrite = TRUE )



















