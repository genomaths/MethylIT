library(testthat)
library(GenomicRanges)
library(MethylIT)

context("nonlinearFitDist tests")

test_that("nonlinearFitDist with Weibull distribution", {
   data(HD)
   nlms <- nonlinearFitDist(HD, npoints = 100, verbose = FALSE)
   expect_true(as.numeric(as.character(nlms$T1$R.Cross.val))[1] > 0.95)
})

test_that("nonlinearFitDist with Gamma2P distribution", {
   data(HD)
   nlms <- suppressWarnings(nonlinearFitDist(HD, npoints = 100,
                                           dist.name = "Gamma2P",
                                           verbose = FALSE))
   expect_true(as.numeric(as.character(nlms$T1$R.Cross.val))[1] > 0.95)
})
