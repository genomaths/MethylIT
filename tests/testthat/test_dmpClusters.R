library(testthat)
library(MethylIT)

context("dmpClusters tests")

test_that("dmpClusters test", {
   data(dmps)
   x1 = dmpClusters(GR = dmps, maxDist = 30, minNumDMPs = 10,
                    num.cores=2L, method = "fixed.int", verbose = FALSE)
   expect_true(length(x1) == 11)
})


